# 新規subset作成用
import os
import json
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from typing import Tuple, List
from pandas import DataFrame
from PIL import Image
from tifffile import TiffFile
import dask.dataframe as dd
from tqdm import tqdm


# 画像の切り抜く始点の情報（オフセットに相当）
def read_crop_info_json(output_dir:str, fov:int) -> Tuple[int, int]:
    json_path = os.path.join(output_dir, f"crop_info.json")
    with open(json_path, "r") as file:
        try:
            json_data = json.load(file)
        except FileNotFoundError as e:
            raise FileNotFoundError(f"{e}: {json_path} not found")
        except json.JSONDecodeError as e:
            raise json.JSONDecodeError(f"{e}: failed to decode JSON file")
    
    return json_data[f"subset{fov}"]["crop_x_pixel"][0], json_data[f"subset{fov}"]["crop_y_pixel"][0]



# 該当fovの中で最小のxyを返す 
def get_min_coordinates(input_dir:str, fov:int) -> Tuple[float, float]:
    csv_path = os.path.join(input_dir, "detected_transcripts.csv")

    try:
        cols = ['global_x', 'global_y', 'fov']
        data_ddf = dd.read_csv(csv_path, usecols=cols) 
    except FileNotFoundError as e:
        raise FileNotFoundError(f'{e}: {csv_path} not found')

    selected_ddf = data_ddf.query(f'fov == {fov}').compute() 

    if selected_ddf.empty:
        raise Exception(f"no match data for fov={fov}")

    x_min = selected_ddf["global_x"].min()
    y_min = selected_ddf["global_y"].min()
    return x_min, y_min

# manifest.jsonより必要情報の取得
def read_manifest_json(input_dir:str) -> Tuple[float, float, float, int, int]:
    json_path = os.path.join(input_dir, "images/manifest.json")
    with open(json_path, "r") as file:
        try:
            json_data = json.load(file)
        except FileNotFoundError as e:
            raise FileNotFoundError(f"{e}: {json_path} not found")
        except json.JSONDecodeError as e:
            raise json.JSONDecodeError(f"{e}: failed to decode JSON file")
    
    return json_data["microns_per_pixel"], json_data["bbox_microns"][0], json_data["bbox_microns"][1],\
           json_data["mosaic_width_pixels"], json_data["mosaic_height_pixels"]


def microns_to_pixel(
        x_micron:float, y_micron:float, microns_per_pixel:float, bbox_x:float, bbox_y:float
        ) -> Tuple[float, float]:
    x_pixel = (x_micron - bbox_x) / microns_per_pixel
    y_pixel = (y_micron - bbox_y) / microns_per_pixel
    return x_pixel, y_pixel


# 変換式： 0.108 * xp + bbox[0] = xm / 0.108 * yp + bbox[1] = ym
def pixel_to_microns(
        x_pixel:float, y_pixel:float, microns_per_pixel:float, bbox_x:float, bbox_y:float
        ) -> Tuple[float, float]:
    x_micron = microns_per_pixel * x_pixel + bbox_x
    y_micron = microns_per_pixel * y_pixel + bbox_y
    return x_micron, y_micron


def global_to_local(global_x:float, global_y:float, x_offset, y_offset) -> Tuple[float, float]:
    local_x = global_x - x_offset
    local_y = global_y - y_offset
    return local_x, local_y


# dfにpixelの列を追加する
def add_pixel_columns(input_dir:str, df:DataFrame, x_init:int, y_init:int) -> DataFrame:
    microns_per_pixel, bbox_micron_x, bbox_micron_y, _image_width, _image_height = read_manifest_json(input_dir)
    df[["global_pixel_x", "global_pixel_y"]] = df.apply(
        lambda row: microns_to_pixel(row["global_x"], row["global_y"], microns_per_pixel, bbox_micron_x, bbox_micron_y), 
        axis=1, result_type = "expand"
    )

    df[["local_pixel_x", "local_pixel_y"]] = df.apply(
        lambda row: global_to_local(row["global_pixel_x"], row["global_pixel_y"], x_init, y_init), 
        axis=1, result_type = "expand"
    )
    print("Add local/global pixel columns. Dataframe shape is ", df.shape) # 14列になる
    return df

# 2048*2048でcsvを切り抜き、保存
def crop_csv(input_dir:str, fov:int, subset_dir:str, x_min_micron:float, y_min_micron:float, 
             x_init:int, y_init:int, width:int = 2048, height:int = 2048) -> None:
    csv_path = os.path.join(input_dir, "detected_transcripts.csv")

    try:
        # cols = ['global_x', 'global_y', 'global_z', 'x', 'y', 'fov', 'gene', 'transcript_id']
        data_ddf = dd.read_csv(csv_path)#, usecols=cols) 
    except FileNotFoundError as e:
        raise FileNotFoundError(f'{e}: {csv_path} not found')
    
    # ddf から x_min_pixel ~ x_min_pixel+widht かつ y_min_pixel ~ y_min_pixel + height を抽出する -> 120sec
    microns_per_pixel, bbox_micron_x, bbox_micron_y, _image_width, _image_height = read_manifest_json(input_dir)
    x_width_micron, y_height_micron = pixel_to_microns(x_init+width, y_init+height, microns_per_pixel, bbox_micron_x, bbox_micron_y)
    selected_ddf = data_ddf.query(f'{x_min_micron} <= global_x & global_x <= {x_width_micron} & \
                                  {y_min_micron} <= global_y & global_y <= {y_height_micron}').compute() 

    print(f"global_x {x_min_micron}:{x_width_micron} \nglobal_y {y_min_micron}:{y_height_micron}")

    selected_ddf = add_pixel_columns(input_dir, selected_ddf, x_init, y_init)

    # 出力ファイルパス
    pkl_path = os.path.join(subset_dir, f"detected_transcripts_subset{fov}.pkl")
    selected_ddf.compute().to_pickle(pkl_path)
    # print(f"Selected data saved to '{pkl_path}'")

    return


# .tifのパスを取得
def get_imagepath(input_dir:str) -> List[str]:
    path_list = []    
    image_input_dir = os.path.join(input_dir, "images")

    if not os.path.exists(image_input_dir):
        raise Exception(f"{image_input_dir} not found")

    tif_files = [f for f in os.listdir(image_input_dir) if f.endswith(".tif")]

    for tif_file in tif_files:
        image_path = os.path.join(image_input_dir, tif_file)
        path_list.append(image_path)

    if not path_list:
        raise Exception(f"no .tif files in {image_input_dir}")
    
    return path_list


# 画像を切り抜く
def crop_images(
        fov:int, subset_dir:str, path_list:List[str], x_init:int, y_init:int, is_show:bool, w:int=2048, h:int=2048
        ) -> None:
    
    image_output_dir = os.path.join(subset_dir, "images")
    os.makedirs(image_output_dir, exist_ok=True) 

    for index, tif_path in enumerate(path_list):
        with TiffFile(tif_path) as tif:
            img = tif.asarray(out="memmap")

        if(is_show):
            if index == 0:
                num_columns = 5
                num_images = len(path_list)
                num_rows = math.ceil(num_images / num_columns)

                img0 = np.asarray(img)
                fig, ax = plt.subplots(figsize=(6, 6))

                rect = patches.Rectangle((x_init, y_init), w, h, linewidth=1, edgecolor='yellow', facecolor='none')
                ax.add_patch(rect)

                ax.set_title(f"{tif_path}: cropped area is shown.")
                ax.imshow(img0)

                fig, ax = plt.subplots(nrows=num_rows, ncols=num_columns, figsize=(20, 5*num_rows), tight_layout=True)         

        # 列->行の順番なのでyが先になることに注意
        cropped_img = img[y_init:y_init+h, x_init:x_init+w]

        # uint16 -> uint8へ正規化して型変換
        scaled_img = (cropped_img - cropped_img.min()) / (cropped_img.max() - cropped_img.min()) * 255
        scaled_img = scaled_img.astype(np.uint8)
        
        # 出力ファイルパス
        file_name = os.path.basename(tif_path)
        file_name_without_extension, _file_extension = os.path.splitext(file_name)
        output_file = file_name_without_extension + f"_subset{fov}.png"
        output_path = os.path.join(image_output_dir, output_file)
        
        if(is_show):
            row, col = divmod(index, num_columns)
            ax[row, col].imshow(scaled_img)
            ax[row, col].set_title(file_name_without_extension)

        img_to_save = Image.fromarray(scaled_img)
        img_to_save.save(output_path)
        # print(f"Saved {output_path}")

    if(is_show):
        for index in range(num_images, num_rows * num_columns):
            row, col = divmod(index, num_columns)
            ax[row, col].axis('off')
        plt.show()
    
    # print("all cropped images saved.")
    
    return 

def create_subset(input_dir:str, output_dir:str, fov:int, is_show:bool = False, width:int = 2048, height:int = 2048):
    if not os.path.isdir(input_dir):
        print("no input directory")
        return
    try:
        # 出力ディレクトリ作成
        subset_dir = os.path.join(output_dir, f"subset{fov}")
        if not os.path.isdir(subset_dir):
            os.makedirs(subset_dir)
            # print(f"create {subset_dir}")

        x_min_micron, y_min_micron = get_min_coordinates(input_dir, fov)
        x_init, y_init = read_crop_info_json(output_dir, fov)
        #calculate_crop_coordinates(input_dir, subset_dir, x_min_micron, y_min_micron)
        crop_csv(input_dir, fov, subset_dir, x_min_micron, y_min_micron, x_init, y_init)
        
        tif_path_list = get_imagepath(input_dir)
        crop_images(fov, subset_dir, tif_path_list, x_init, y_init, is_show)
        # print(f"Create subset{fov} and saved to {subset_dir}.")
        
    except Exception as e:
        print(f"Error: {e}")
        # エラーが発生した場合、subset_dir が空かどうかを確認して削除
        if os.path.isdir(subset_dir) and not os.listdir(subset_dir):
            try:
                os.rmdir(subset_dir)  # ディレクトリが空の場合削除
                print(f"{subset_dir} is empty and remove it")
            except OSError as e:
                print(f"Error while deleting directory: {e}")
        return


if __name__ == "__main__":
    input_dir = "/work/datasets/"
    output_dir = "/work/outputs/"
    # target_str = "" # "z0"だけを対象とするなどしてテストが可能
    target_fov = 1
    for fov in tqdm(target_fov):
        create_subset(input_dir, output_dir, target_fov) # 画像表示が必要な場合はis_show=True