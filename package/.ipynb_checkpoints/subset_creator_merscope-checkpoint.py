"""
MERSCOPE用の子クラス
"""

import os
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Tuple, List
from PIL import Image
from tifffile import TiffFile
import dask.dataframe as dd

from package.subset_creator_abc import AbstractSubsetCreator
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

class MerscopeSubsetCreator(AbstractSubsetCreator, CommonSubsetCreator):
    def __init__(self, input_dir:int, output_dir:int, fov:int, width:int, height:int, image_keyword:List[str]=None, z:int = None):
        super().__init__(input_dir, output_dir, fov, width, height)
        self.image_keyword = image_keyword
        self.z = z

    @error_handler
    def read_manifest_json(self) -> Tuple[float, float, float, int, int]:
        json_path = os.path.join(self.input_dir, "images/manifest.json")
        with open(json_path, "r") as file:
            try:
                json_data = json.load(file)
            except FileNotFoundError as e:
                raise FileNotFoundError(f"{e}: {json_path} not found")
            except json.JSONDecodeError as e:
                raise json.JSONDecodeError(f"{e}: failed to decode JSON file")
        
        return json_data["microns_per_pixel"], json_data["bbox_microns"][0], json_data["bbox_microns"][1],\
            json_data["mosaic_width_pixels"], json_data["mosaic_height_pixels"]

    @error_handler
    def calculate_crop_coordinates(self, json_path:str, x_min_micron:float, y_min_micron:float)-> Tuple[int, int, int, int]:
        microns_per_pixel, bbox_micron_x, bbox_micron_y, image_width, image_height = self.read_manifest_json(self.input_dir)
        x_min_pixel, y_min_pixel = self.microns_to_pixel(x_min_micron, y_min_micron, microns_per_pixel, bbox_micron_x, bbox_micron_y)

        x_init = int(max(0, min(x_min_pixel, image_width - self.width)))
        y_init = int(max(0, min(y_min_pixel, image_height - self.height)))
        x_last = x_init + self.width
        y_last = y_init + self.height
        
        data_to_write = { 
            f"subset{self.fov}": {
            "crop_x_pixel": [x_init, x_last],
            "crop_y_pixel": [y_init, y_last] }
        }
        self.write_crop_info_json(json_path, data_to_write)   
        
        return x_init, y_init, x_last, y_last

    @error_handler
    def get_min_coordinates(self) -> Tuple[float, float]:
        csv_path = os.path.join(self.input_dir, "detected_transcripts.csv")

        try:
            cols = ['global_x', 'global_y', 'fov']
            data_ddf = dd.read_csv(csv_path, usecols=cols) 
        except FileNotFoundError as e:
            raise FileNotFoundError(f'{e}: {csv_path} not found')

        selected_ddf = data_ddf.query(f'fov == {self.fov}').compute() 

        if selected_ddf.empty:
            raise Exception(f"no match data for fov={self.fov}")

        x_min_micron = selected_ddf["global_x"].min()
        y_min_micron = selected_ddf["global_y"].min()
        return x_min_micron, y_min_micron

    @error_handler
    def make_directory(self)->str:
        subset_dir = os.path.join(self.output_dir, f"subset{self.fov}")
        os.makedirs(subset_dir, exist_ok=True)
        return subset_dir

    @error_handler
    def decide_crop_area(self)->Tuple[int, int, int, int]:
        crop_json_path = os.path.join(self.output_dir, f"crop_info.json")
        if os.path.exists(crop_json_path):
            x_init, y_init, x_last, y_last = self.read_crop_info_json(crop_json_path, self.fov)
            if x_init == -1 & y_init == -1:
                x_min_micron, y_min_micron = self.get_min_coordinates(self.input_dir, self.fov)
                x_init, y_init, x_last, y_last = self.calculate_crop_coordinates(self.input_dir, self.fov, crop_json_path, x_min_micron, y_min_micron, self.width, self.height)
        else:
            x_min_micron, y_min_micron = self.get_min_coordinates(self.input_dir, self.fov)
            x_init, y_init, x_last, y_last = self.calculate_crop_coordinates(self.input_dir, self.fov, crop_json_path, x_min_micron, y_min_micron, self.width, self.height)
        
        return x_init, y_init, x_last, y_last

    @error_handler
    def crop_table_data(self, subset_dir:str, x_init:int, y_init:int, x_last:int, y_last:int) -> None:
        csv_path = os.path.join(self.input_dir, "detected_transcripts.csv")

        try:
            data_ddf = dd.read_csv(csv_path)#, usecols=cols) 
        except FileNotFoundError as e:
            raise FileNotFoundError(f'{e}: {csv_path} not found')
        
        # ddf から x_min_pixel ~ x_min_pixel+widht かつ y_min_pixel ~ y_min_pixel + height を抽出する 
        microns_per_pixel, bbox_micron_x, bbox_micron_y, _image_width, _image_height = self.read_manifest_json(self.input_dir)
        x_width_micron, y_height_micron = self.pixel_to_microns(x_last, y_last, microns_per_pixel, bbox_micron_x, bbox_micron_y)
        x_min_micron, y_min_micron = self.get_min_coordinates(self.input_dir, self.fov)
        selected_df = data_ddf.query(f'{x_min_micron} <= global_x & global_x <= {x_width_micron} & {y_min_micron} <= global_y & global_y <= {y_height_micron}').compute() 

        microns_per_pixel, bbox_micron_x, bbox_micron_y, _image_width, _image_height = self.read_manifest_json(self.input_dir)
        selected_df[["global_pixel_x", "global_pixel_y"]] = selected_df.apply(
            lambda row: self.microns_to_pixel(row["global_x"], row["global_y"], microns_per_pixel, bbox_micron_x, bbox_micron_y), 
            axis=1, result_type = "expand"
        )

        selected_df[["local_pixel_x", "local_pixel_y"]] = selected_df.apply(
            lambda row: self.global_to_local(row["global_pixel_x"], row["global_pixel_y"], x_init, y_init), 
            axis=1, result_type = "expand"
        )

        # 出力ファイルパス
        pkl_path = os.path.join(subset_dir, f"detected_transcripts_subset{self.fov}.pkl")
        selected_df.to_pickle(pkl_path)
        print(f"saved to {pkl_path}")
        return

    @error_handler
    def get_imagepath(self, image_keywords) -> List[str]:
        path_list = []    
        image_input_dir = os.path.join(self.input_dir, "images")

        if not os.path.exists(image_input_dir):
            raise Exception(f"{image_input_dir} not found")

        if not image_keywords:
            tif_files = [f for f in os.listdir(image_input_dir) if f.endswith(".tif")]
        else:
            tif_files = [f for f in os.listdir(image_input_dir) if f.endswith(".tif") and all(word in f for word in image_keywords)]

        for tif_file in tif_files:
            image_path = os.path.join(image_input_dir, tif_file)
            path_list.append(image_path)

        if not path_list:
            raise Exception(f"no .tif files {'containing' if image_keywords else 'found'} {' and '.join(image_keywords)} in {image_input_dir}")

        return path_list
        
    @error_handler
    def crop_image(
            self, subset_dir:str, x_init:int, y_init:int, x_last:int, y_last:int, is_show:bool=True) -> None:
        
        image_output_dir = os.path.join(subset_dir, "images")
        os.makedirs(image_output_dir, exist_ok=True) 

        path_list = self.get_imagepath(self.input_dir, self.image_keyword)
        
        for index, tif_path in enumerate(path_list):
            with TiffFile(tif_path) as tif:
                img = tif.asarray(out="memmap") # .astype(np.uint8)
                # print(img.dtype)

            # 列->行の順番なのでyが先になることに注意
            cropped_img = img[y_init:y_last, x_init:x_last]
            scaled_img = (cropped_img - cropped_img.min()) / (cropped_img.max() - cropped_img.min()) * 255
            scaled_img = scaled_img.astype(np.uint8)

            if is_show and index == 0:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
                
                rect = patches.Rectangle((x_init, y_init), self.width, self.height, linewidth=1, edgecolor='yellow', facecolor='none')
                ax1.add_patch(rect)
                ax1.set_title(f"fov{self.fov}: cropped area")
                ax1.imshow(img)
                
                ax2.set_title(f"fov{self.fov}: cropped image")
                ax2.imshow(scaled_img)
                
                plt.tight_layout()
                plt.show()
            
            # 出力ファイルパス
            file_name = os.path.basename(tif_path)
            file_name_without_extension, _file_extension = os.path.splitext(file_name)
            output_file = file_name_without_extension + f"_subset{self.fov}.png"
            output_path = os.path.join(image_output_dir, output_file)

            img_to_save = Image.fromarray(scaled_img)
            img_to_save.save(output_path)
            print(f"saved to {output_path}")
        return 
        
    @error_handler
    def temporary_test(self)->None:
        pass

    @error_handler
    def create_subset(self) -> None:
        subset_dir = self.make_directory()
        x_init, y_init, x_last, y_last = self.decide_crop_area()
        
        self.crop_table_data(subset_dir, x_init, y_init, x_last, y_last)
        self.crop_image(subset_dir, x_init, y_init, x_last, y_last)


"""
# 使用例
if __name__ == "__main__":
    input_dir = "/work/datasets"
    output_dir = "/work/output_MERSCOPE"

    # subset分割
    fov = 1
    width = 2048
    height = 2048
    image_keyword = ["DAPI", "z0"]
    
    subset = MerscopeSubsetCreator(input_dir, output_dir, fov, width, height, image_keyword)
    subset.create_subset()
"""