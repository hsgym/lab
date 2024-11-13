"""
子クラスにおいて共通して利用する関数をまとめたクラス：
"""

import json
import os
import yaml

from abc import ABC
from typing import Tuple


# ファクトリークラス（いつか作りたい） --------------------------------------------------------------------------
"""
class SubsetCreatorFactory:
@staticmethod
def create_subset_creator(data_format: str) -> AbstractSubsetCreator:
    if data_format == 'xenium':
        return XeniumSubsetCreator()
    elif data_format == 'merscope':
        return MERSCOPESubsetCreator()
    else:
        raise ValueError(f"Unsupported data format: {data_format}")
"""        
  
class CommonSubsetCreator():
    """
    ファイルの読み書きや座標変換用の関数を持つ
    """
    def read_crop_info_json(self, json_path:str, fov:int) -> Tuple[int, int, int, int]:
        with open(json_path, "r") as file:
            try:
                json_data = json.load(file)
            except FileNotFoundError as e:
                raise FileNotFoundError(f"{e}: {json_path} not found")
            except json.JSONDecodeError as e:
                raise json.JSONDecodeError(f"{e}: failed to decode JSON file")
            
        if f"subset{fov}" in json_data:
            return json_data[f"subset{fov}"]["crop_x_pixel"][0], json_data[f"subset{fov}"]["crop_y_pixel"][0],\
                json_data[f"subset{fov}"]["crop_x_pixel"][1], json_data[f"subset{fov}"]["crop_y_pixel"][1]
        else:
            return -1, -1, -1, -1


    def write_crop_info_json(self, json_path:str, data_to_write) -> None:
        if not os.path.exists(json_path):
            with open(json_path, "w") as f:
                json.dump(data_to_write, f)
            return
        
        with open(json_path, "r") as file:
            try:
                json_data = json.load(file) 
            except json.JSONDecodeError as e:
                raise json.JSONDecodeError(f"{e}: failed to decode JSON file")

        json_data.update(data_to_write)

        with open(json_path, 'w') as f:
            json.dump(json_data, f)
        return


    def read_info_yaml(self, yaml_path:str) -> Tuple[int, int]:
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        return data['imgae_data']['width'],  data['imgae_data']['height']


    def write_info_yaml(self, yaml_path, data_to_write) -> None:
        with open(yaml_path,'w')as f:
            yaml.dump(data_to_write, f, default_flow_style=False, allow_unicode=True)
        return


    def global_to_local(self, global_x:float, global_y:float, x_offset, y_offset) -> Tuple[float, float]:
        return (global_x - x_offset), (global_y - y_offset)

    def microns_to_pixel(self, x_micron:float, y_micron:float, microns_per_pixel:float, bbox_x:float, bbox_y:float) -> Tuple[float, float]:
        x_pixel = (x_micron - bbox_x) / microns_per_pixel
        y_pixel = (y_micron - bbox_y) / microns_per_pixel
        return x_pixel, y_pixel

    def pixel_to_microns(self, x_pixel:float, y_pixel:float, microns_per_pixel:float, bbox_x:float, bbox_y:float) -> Tuple[float, float]:
        x_micron = microns_per_pixel * x_pixel + bbox_x
        y_micron = microns_per_pixel * y_pixel + bbox_y
        return x_micron, y_micron
