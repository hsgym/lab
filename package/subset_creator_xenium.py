# subset_creator_Xenium.py 
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Tuple
from PIL import Image
from tifffile import TiffFile
import dask.dataframe as dd

from package.subset_creator_abc import AbstractSubsetCreator
from codes.package.subset_creator_help import make_fov_info_json_xenium, image_size_xenium
from package.subset_creator_common import CommonSubsetCreator

# エラーハンドリング
def error_handler(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except FileNotFoundError as e:
            print(f"File not found error: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")
            raise
    return wrapper

class XeniumSubsetCreator(AbstractSubsetCreator, CommonSubsetCreator):
    def __init__(self, input_dir, output_dir, fov, width, height, resampling_factor):
        super().__init__(input_dir, output_dir, fov, width, height)
        self.resampling_factor = resampling_factor

    @error_handler
    def calculate_crop_coordinates(
            self, json_path:str, x_min_pixel:int, y_min_pixel:int, image_width:int, image_height:int
        ) -> Tuple[int, int, int, int]:
        x_init = int(max(0, min(x_min_pixel, image_width - self.width)))
        y_init = int(max(0, min(y_min_pixel, image_height - self.height)))
        x_last = x_init + self.width
        y_last = y_init + self.height
        
        data_to_write = { 
            f"subset{self.fov}": {
                "crop_x_pixel": [x_init, x_last],
                "crop_y_pixel": [y_init, y_last]
            }
        }
        self.write_crop_info_json(json_path, data_to_write)
        return x_init, y_init, x_last, y_last

    @error_handler
    def make_directory(self) -> str:
        subset_dir = os.path.join(self.output_dir, f"subset{self.fov}")
        os.makedirs(subset_dir, exist_ok=True)
        return subset_dir

    @error_handler       
    def decide_crop_area(self) -> Tuple[int, int, int, int]:
        crop_json_path = os.path.join(self.output_dir, f"crop_info.json")
        
        if os.path.exists(crop_json_path):
            x_init, y_init, x_last, y_last = self.read_crop_info_json(crop_json_path, self.fov)
            if (x_init == -1 & y_init == -1):
                fov_json_path = os.path.join(self.output_dir, "fov_info.json")
                x_init, y_init, _x_fov, _y_fov = self.read_crop_info_json(fov_json_path, self.fov)
                image_width, image_height = self.read_info_yaml(os.path.join(self.output_dir, "image_info.yaml"))
                x_init, y_init, x_last, y_last = self.calculate_crop_coordinates(self.fov, crop_json_path, x_init, y_init, self.width, self.height, image_width, image_height)
        else:
            fov_json_path = os.path.join(self.output_dir, "fov_info.json")
            x_init, y_init, _x_fov, _y_fov = self.read_crop_info_json(fov_json_path, self.fov)
            image_width, image_height = self.read_info_yaml(os.path.join(self.output_dir, "image_info.yaml"))
            x_init, y_init, x_last, y_last = self.calculate_crop_coordinates(self.fov, crop_json_path, x_init, y_init, self.width, self.height, image_width, image_height)
        
        return x_init, y_init, x_last, y_last

    @error_handler    
    def crop_table_data(self, subset_dir:str, x_init:int, y_init:int, x_last:int, y_last:int) -> None:
        parquet_path = os.path.join(self.input_dir, f"transcripts.parquet")
        try:
            data_ddf = dd.read_parquet(parquet_path) 
        except FileNotFoundError as e:
            raise FileNotFoundError(f'{e}: {parquet_path} not found')
        
        _x_row_name = "x_location"
        _y_row_name = "y_location"
        microns_per_pixel = self.resampling_factor
        bbox_micron_x = 0.0
        bbox_micron_y = 0.0
        
        x_min_micron, y_min_micron = self.pixel_to_microns(float(x_init), float(y_init), microns_per_pixel, bbox_micron_x, bbox_micron_y)
        x_width_micron, y_height_micron = self.pixel_to_microns(float(x_last), float(y_last), microns_per_pixel, bbox_micron_x, bbox_micron_y)
        
        selected_df = data_ddf.query(f'{x_min_micron} <= x_location & x_location <= {x_width_micron} & {y_min_micron} <= y_location & y_location <= {y_height_micron}').compute() 

        
        selected_df[["x_location_pixel", "y_location_pixel"]] = selected_df.apply(
            lambda row: self.microns_to_pixel(row["x_location"], row["y_location"], microns_per_pixel, bbox_micron_x, bbox_micron_y), 
            axis=1, result_type = "expand"
        )

        selected_df[["local_pixel_x", "local_pixel_y"]] = selected_df.apply(
            lambda row: self.global_to_local(row["x_location_pixel"], row["y_location_pixel"], x_init, y_init), 
            axis=1, result_type = "expand"
        )

        pkl_path = os.path.join(subset_dir, f"detected_transcripts_subset{self.fov}.pkl")
        selected_df.to_pickle(pkl_path)
        print(f"saved to {pkl_path}")
        return
    
    @error_handler    
    def crop_image(self, subset_dir:str, x_init:int, y_init:int, x_last:int, y_last:int, is_show:bool=True) -> None:
        image_output_dir = os.path.join(subset_dir, "images")
        os.makedirs(image_output_dir, exist_ok=True) 
        
        # 画像の読み込み
        ometif_path = os.path.join(self.input_dir, "morphology_focus.ome.tif") 
        with TiffFile(ometif_path) as tif:
            img = tif.series[0].levels[0].asarray().astype(np.uint8)
        
        cropped_img = img[y_init:y_last, x_init:x_last]
        
        if is_show:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
            
            rect = patches.Rectangle((x_init, y_init), self.width, self.height, linewidth=1, edgecolor='yellow', facecolor='none')
            ax1.add_patch(rect)
            ax1.set_title(f"fov{self.fov}: cropped area")
            ax1.imshow(img)
            
            ax2.set_title(f"fov{self.fov}: cropped image")
            ax2.imshow(cropped_img)
            
            plt.tight_layout()
            plt.show()
        
        # 出力ファイルパス
        file_name = os.path.basename(ometif_path).split('.', 1)[0]
        output_file = f"{file_name}_subset{self.fov}.png"
        output_path = os.path.join(image_output_dir, output_file)
     
        img_to_save = Image.fromarray(cropped_img)
        img_to_save.save(output_path)
        print(f"saved to {output_path}")
        return
    
    @error_handler  
    def create_subset(self) -> None:
        subset_dir = self.make_directory()
        x_init, y_init, x_last, y_last = self.decide_crop_area()
        
        self.crop_table_data(subset_dir, x_init, y_init, x_last, y_last)
        self.crop_image(subset_dir, x_init, y_init, x_last, y_last)

"""
# 使用例
if __name__ == "__main__":
    input_dir = "/work/datasets/Okada/output-XETG00130__0014491__Cont_3_1__20240403__094846"
    output_dir = "/work/output_Xenium"
    resampling_factor = 0.2125

    # xeniumの準備（一度だけ実行）
    input_json_path = os.path.join(input_dir, "aux_outputs/morphology_fov_locations.json") # micron座標
    output_json_path = os.path.join(output_dir, "fov_info.json")
    make_fov_info_json_xenium(input_json_path, output_json_path, resampling_factor)

    ometif_path = os.path.join(input_dir, "morphology_focus.ome.tif") 
    yaml_path = os.path.join(output_dir, "image_info.yaml")
    image_width, image_height = image_size_xenium(ometif_path, yaml_path)

    # subset分割
    fov = 31
    width = 4096
    height = 4096

    subset = XeniumSubsetCreator(input_dir, output_dir, fov, width, height, resampling_factor)
    subset.create_subset()
"""

