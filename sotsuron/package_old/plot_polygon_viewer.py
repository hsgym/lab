# 新規plot_polygon_viewer
import os
from typing import Tuple, List, Dict, Optional

import dask.dataframe as dd
import h5py
import json
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from pandas import DataFrame
from PIL import Image
from tifffile import TiffFile
from matplotlib.axes._axes import Axes

from package_old.subset_creater import (
    read_manifest_json, 
    microns_to_pixel, 
    global_to_local,
    read_crop_info_json
    )


def convert_to_local_pixel(
    input_dir:str, output_dir:str, fov, df:DataFrame, original:str, result:str
    ) -> pd.DataFrame:

    microns_per_pixel, bbox_micron_x, bbox_micron_y, _image_width, _image_height = read_manifest_json(input_dir)
    x_init, y_init = read_crop_info_json(output_dir, fov)

    df[[f"global_pixel_{original}_x", f"global_pixel_{original}_y"]] = df.apply(
        lambda row: microns_to_pixel(row[f"{original}_x"], row[f"{original}_y"], \
            microns_per_pixel, bbox_micron_x, bbox_micron_y),
        axis=1, result_type="expand"
    )

    df[[f"{result}_x", f"{result}_y"]] = df.apply(
        lambda row: global_to_local(row[f"global_pixel_{original}_x"], row[f"global_pixel_{original}_y"], x_init, y_init),
        axis=1, result_type="expand"
    )
    return df


# メタデータの読み込み
def load_meta_csv(input_dir:str, output_dir:str, fov:int) -> DataFrame:
    csv_path = os.path.join(input_dir, "cell_metadata.csv")

    try:
        data_df = pd.read_csv(csv_path) 
    except FileNotFoundError as e:
        raise FileNotFoundError(f'{e}: {csv_path} not found')

    selected_df = data_df.query(f"fov=={fov}").copy()
    
    if selected_df.empty:
        raise Exception(f"no match data for fov={fov}")

    selected_df = convert_to_local_pixel(input_dir, output_dir, fov, selected_df, "center", "local_pixel_center")
    selected_df = convert_to_local_pixel(input_dir, output_dir, fov, selected_df, "max", "local_pixel_max")
    selected_df = convert_to_local_pixel(input_dir, output_dir, fov, selected_df, "min", "local_pixel_min")

    return selected_df


# HDF5ファイルからデータを読み込む関数
def load_cell_data_from_hdf5(cell_id:str, hdf5_file:h5py.File) -> Dict[str, pd.DataFrame]:
    cell_data = {}
    cell_path = f"featuredata/{cell_id}/zIndex_0/p_0/coordinates"
    if cell_path in hdf5_file:
        cell_df = hdf5_file[cell_path][0]
        cell_data = pd.DataFrame(data=cell_df, columns=["global_x", "global_y"])
    return cell_data

# データ読み込みと変換  <-- mainプログラム参照
def load_and_convert_data(input_dir:str, output_dir:str, fov:int) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    selected_cell_data = load_meta_csv(input_dir, output_dir, fov)
    selected_cell_data.rename(columns={"Unnamed: 0": "cell_id"}, inplace = True)
    selected_cell_data.set_index("cell_id", inplace=True)

    cell_data_dict = {}

    hdf5_path = os.path.join(input_dir, f"cell_boundaries/feature_data_{fov}.hdf5")
    if not os.path.isfile(hdf5_path):
            raise FileNotFoundError(f'HDF5 file not found: {hdf5_path}')

    with h5py.File(hdf5_path, 'r') as f:
        cell_ids = f["featuredata"].keys()
        for cell_id in cell_ids:
            cell_data_dict[cell_id] = load_cell_data_from_hdf5(cell_id, f)
            cell_data_dict[cell_id] = convert_to_local_pixel(input_dir, output_dir, fov, cell_data_dict[cell_id], "global", "local_pixel")
    return cell_data_dict, selected_cell_data


def read_select_transcripts_csv(subset_dir:str,fov:int, z:int = None) -> Tuple[DataFrame, List[str]]:
    pkl_path = os.path.join(subset_dir, f"detected_transcripts_subset{fov}.pkl")
    try:
        df = pd.read_pickle(pkl_path)
    except FileNotFoundError as e:
        raise FileNotFoundError(f'{e}: {pkl_path} not found')
    
    if z is None: # zが無い時
        selected_df = df.copy()
    else:
        selected_df = df.query(f"global_z == {z}")

    if selected_df.empty:
        raise Exception(f"DataFrame is empty")

    gene_name_list = selected_df['gene'].unique().tolist()
    return selected_df, gene_name_list



# 画像を読み込む
def load_image(img_dir:str, fov:int, z:int) -> Image:
    png_path = os.path.join(img_dir, f"mosaic_DAPI_z{z}_subset{fov}.png")
    if not os.path.isfile(png_path):
        raise FileNotFoundError(f"File not found: {png_path}")
    return Image.open(png_path)


def draw_gene_plot(
    ax:Axes, df:pd.DataFrame, gene_list: List[str], total_gene_num:int, size:int=2
    ) -> None:
    
    if len(gene_list) != total_gene_num:
        selected_df = df[df["gene"].isin(gene_list)]
        sns.scatterplot(data=selected_df, x="local_pixel_x", y="local_pixel_y", hue="gene", palette="hls", s=size, legend=True, ax=ax)
    else:
        sns.scatterplot(data=df, x="local_pixel_x", y="local_pixel_y", hue="gene", palette="hls", s=size, legend=False, ax=ax)
    return


# ポリゴンを1つ表示
def draw_polygon(ax:plt.Axes, x_data:float, y_data:float) -> None:
    ax.plot(x_data, y_data, linestyle="-", alpha=0.9, c = "#d62728")
    ax.fill(x_data, y_data, alpha=0.1, c="#d62728")
    return 


# ポリゴンを表示（ループを回す）
def draw_all_polygon(ax: plt.Axes, cell_data_dict: Dict[str, pd.DataFrame], metadata_df: DataFrame) -> None:
    for cell_id in cell_data_dict:
        if cell_id in metadata_df.index:
            cell_data = cell_data_dict[cell_id]
            draw_polygon(ax, cell_data["local_pixel_x"], cell_data["local_pixel_y"])
    return


# ポリゴンを表示(位置合わせ付き)
def draw_id_polygon(
    cell_id: str, ax: plt.Axes, cell_data_dict: Dict[str, pd.DataFrame], metadata_df: pd.DataFrame, w: int = 512, h: int = 512
    ) -> Tuple[float, float]:

    if not cell_id:
        raise ValueError("Cell ID not provided.")
    if cell_id not in cell_data_dict:
        raise ValueError(f"Cell with ID {cell_id} not found.")
    
    cell_data = cell_data_dict[cell_id]

    selected_metadata_df = metadata_df.loc[cell_id]
    if selected_metadata_df.empty:
        raise ValueError(f"Cell metadata with ID {cell_id} not found.")

    x = cell_data["local_pixel_x"]
    y = cell_data["local_pixel_x"]
    x_center = selected_metadata_df["local_pixel_center_x"]
    y_center = selected_metadata_df["local_pixel_center_y"]

    # crop位置に合わせて座標を計算し直す
    x_init = x_center - w/2
    y_init = y_center - h/2
    x_conv = x - x_init
    y_conv = y - y_init

    draw_polygon(ax, x_conv, y_conv)
    return x_init, y_init


# fov全体を表示する
def show_whole_fov(
    ax:plt.Axes, fov: int, z: int, cell_data_dict: Dict[str, pd.DataFrame], metadata_df: pd.DataFrame, selected_df:DataFrame,
    gene_list:List[str], total_gene_num:int, show_plot: bool = True, show_polygon: bool = True
    ) -> None:
    

    if show_plot: # plot用の関数を呼ぶ
        selected_z_df = selected_df.query(f"global_z=={z}")
        if selected_z_df.empty:
            raise ValueError("match z data none")
        if len(gene_list) != total_gene_num:
            draw_gene_plot(ax, selected_z_df, gene_list, total_gene_num, size=4) # 一部の遺伝子を表示
        else:
            draw_gene_plot(ax, selected_z_df, gene_list, total_gene_num) # 全遺伝子表示
    
    if show_polygon:
        # polygon用の関数を呼ぶ
        draw_all_polygon(ax, cell_data_dict, metadata_df)
    return 


def calc_df(
    df:DataFrame, metadata_df, cell_id, w:int = 512, h:int=512
    ) -> Tuple[DataFrame, float, float]:

    selected_metadata_df = metadata_df.loc[cell_id]
    
    x_init = selected_metadata_df["local_pixel_center_x"] - w/2
    y_init = selected_metadata_df["local_pixel_center_y"] - h/2

    new_df = df.copy() 
    new_df["local_pixel_x"] = df["local_pixel_x"] - x_init
    new_df["local_pixel_y"] = df["local_pixel_y"] - y_init

    return new_df, x_init, y_init


# cell_idの周辺領域を表示する
def show_id_surroundings(
    ax:plt.Axes, 
    cell_data_dict: Dict[str, pd.DataFrame], 
    metadata_df: pd.DataFrame, 
    selected_df:DataFrame, 
    cell_id: str,
    gene_list:List[str], 
    total_gene_num: int,
    show_plot: bool = True, 
    show_polygon: bool = True, 
    w: int = 512, 
    h: int = 512
    ) -> Tuple[float, float]:

    new_df, x_offset, y_offset = calc_df(selected_df, metadata_df, cell_id)

    if show_plot: # plot用の関数を呼ぶ
        draw_gene_plot(ax, new_df, gene_list, total_gene_num, size=4)
    if show_polygon: # polygon用の関数を呼ぶ
        x_offset, y_offset = draw_id_polygon(cell_id, ax, cell_data_dict, metadata_df)

    return x_offset, y_offset



def display_area(
    img_dir: str, 
    fov: int, 
    z: int, 
    cell_data_dict: Dict[str, pd.DataFrame], 
    metadata_df: pd.DataFrame,  
    selected_df:DataFrame, 
    gene_list:List[str],
    total_gene_num:int,
    cell_id: Optional[str] = None, 
    w: int = 512, 
    h: int = 512,
    show_plot: bool = True, 
    show_polygon: bool = True, 
    is_save: bool = False,
    ax: Optional[plt.Axes] = None  # ax パラメータを追加
    ) -> None:
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    img = load_image(img_dir, fov, z)

    if cell_id: # 一部を表示
        x_offset, y_offset = show_id_surroundings(ax, cell_data_dict, metadata_df, selected_df, gene_list, total_gene_num, cell_id, show_plot=show_plot, show_polygon=show_polygon)
        if x_offset == -1 or y_offset == -1:
            raise ValueError(f"(x_offset, y_offset) = ({x_offset}, {y_offset})")
        img_crop = img.crop((x_offset, y_offset, x_offset + w, y_offset + h))
        ax.set_title(f"fov={fov}, z={z}, id={cell_id}")
        ax.set_xlim([0, w]) 
        ax.set_ylim([h, 0]) 
    else: # fov全体を表示
        show_whole_fov(ax, fov, z, cell_data_dict, metadata_df, selected_df, gene_list, total_gene_num, show_plot=show_plot, show_polygon=show_polygon)
        img_crop = img
        ax.set_title(f"fov={fov}, z={z}")
        ax.set_xlim([0, 2048]) 
        ax.set_ylim([2048, 0]) 

    ax.imshow(img_crop, cmap="gray")
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    if is_save:
        tmp_dir = "/work/tmp_output/"
        filename = os.path.join(tmp_dir, f"fov{fov}_z{z}_display_area.png")
        plt.savefig(filename, format="png")
        print(f"saved {filename}")

    if ax is None:
        plt.show()


def display_focus(
    img_dir:str, fov:int, z:int, cell_data_dict:Dict[str, pd.DataFrame], metadata_df:pd.DataFrame,  selected_df:DataFrame, 
    cell_id:str, gene_list:Optional[List[str]]=None, show_plot:bool = False, w:int=512, h:int=512
    ) -> None: 
    fig, ax = plt.subplots(figsize=(8, 8))
    img = load_image(img_dir, fov, z)

    if not cell_id:
        raise ValueError("Cell ID not provided.")
    if cell_id not in cell_data_dict:
        raise ValueError(f"Cell with ID {cell_id} not found.")
    
    cell_data = cell_data_dict[cell_id]

    selected_metadata_df = metadata_df.loc[cell_id]
    if selected_metadata_df.empty:
        raise ValueError(f"Cell metadata with ID {cell_id} not found.")

    if show_plot: # plot用の関数を呼ぶ
        draw_gene_plot(ax, selected_df, gene_list, total_gene_num)

    draw_polygon(ax, cell_data["local_pixel_x"], cell_data["local_pixel_y"])

    x_offset = selected_metadata_df["local_pixel_center_x"] - w/2
    y_offset = selected_metadata_df["local_pixel_center_y"] - h/2
    rect = patches.Rectangle((x_offset, y_offset), w, h, linewidth=1, edgecolor='red', facecolor='none')
    ax.add_patch(rect)

    ax.set_title(f"fov={fov}, id={cell_id}")
    ax.imshow(img, cmap="gray")
    ax.set_xlim([0, 2048]) 
    ax.set_ylim([2048, 0]) 
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    plt.show()
    
    return


if __name__ == "__main__":

    input_dir = "/work/datasets/"
    output_dir = "/work/outputs/"
    fov = 933
    z = 3 # 仮
    try:
        subset_dir = os.path.join(output_dir, f"subset{fov}")
        img_dir = os.path.join(subset_dir, "images")
        
        # selected_df:表示する転写物のdf
        selected_df, gene_list = read_select_transcripts_csv(subset_dir, fov, z)
        # cell_data_dict:細胞のポリゴンデータ、metadata_df:細胞のメタ情報 
        cell_data_dict, metadata_df = load_and_convert_data(input_dir, output_dir, fov)
        
        total_gene_num = len(gene_list)

        # fov全体
        display_area(img_dir, fov, z, cell_data_dict, metadata_df, selected_df, gene_list, total_gene_num)
        display_area(img_dir, fov, z, cell_data_dict, metadata_df, selected_df, gene_list, total_gene_num, show_plot=False)
        display_area(img_dir, fov, z, cell_data_dict, metadata_df, selected_df, gene_list, total_gene_num, show_polygon=False)
    except Exception as e:
        print(f"{e}")
    