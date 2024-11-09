import json
import numpy as np 
import os
from tifffile import TiffFile
from typing import Optional, Tuple
from package.subset_creator_common import CommonSubsetCreator

"""
補助的に利用する可能性のある関数の置き場
"""


# 前準備（ユーザ自身が決める値のヒント）
def make_fov_info_json_xenium(input_json_path: str, output_json_path: str, resampling_factor: float) -> None:
    common = CommonSubsetCreator()
    with open(input_json_path, 'r') as f:
        fov_conf = json.load(f)

    data_to_write = {}
    fov_index = 1  # FOVのインデックス

    if fov_conf['fov_locations']:
        first_fov_info = next(iter(fov_conf['fov_locations'].values()))
        width = int(first_fov_info['width'] / resampling_factor)
        height = int(first_fov_info['height'] / resampling_factor)

    for fov_name, fov_info in fov_conf['fov_locations'].items():
        x_init = int(fov_info['x'] / resampling_factor)
        y_init = int(fov_info['y'] / resampling_factor)

        data_to_write[f"subset{fov_index}"] = {
            "crop_x_pixel": [x_init, x_init + width],
            "crop_y_pixel": [y_init, y_init + height],
            "original_fov_name": fov_name
        }
        fov_index += 1
    common.write_crop_info_json(output_json_path, data_to_write)
    return

def image_size_xenium(tif_path:str, yaml_path:str)->Tuple[int, int]:
    common = CommonSubsetCreator()
    with TiffFile(tif_path) as tif:
        img = tif.series[0].levels[0].asarray().astype(np.uint8)
        image_height, image_width = img.shape
    
    data_to_write = {
        "imgae_data": {
            "width": image_width,
            "height": image_height
        }
    }
    common.write_info_yaml(yaml_path, data_to_write)
    return image_width, image_height


def decide_fov_hint(data_format:str, input_dir:str, output_dir:str, resampling_factor:Optional[float]=None)->None:
    common = CommonSubsetCreator()
    if data_format == 'xenium':
        input_json_path = os.path.join(input_dir, "aux_outputs/morphology_fov_locations.json") # micron座標
        output_json_path = os.path.join(output_dir, "fov_info.json")
        make_fov_info_json_xenium(input_json_path, output_json_path, resampling_factor)
        
        ometif_path = os.path.join(input_dir, "morphology_focus.ome.tif") 
        yaml_path = os.path.join(output_dir, "image_info.yaml")
        image_width, image_height = image_size_xenium(ometif_path, yaml_path)
        print("image size (width, height): ", image_width, image_height)
    elif data_format == 'merscope':
        _microns_per_pixel, _bbox_micron_x, _bbox_micron_y, image_width, image_height = common.read_manifest_json(input_dir)
        print("image size (width, height): ", image_width, image_height)
    else:
        raise ValueError(f"Unsupported data format: {data_format}")