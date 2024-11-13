#!/usr/bin/env python
# coding: utf-8
import os
from typing import Tuple, List, Dict, Optional, Union, Callable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from shapely.geometry import Point, Polygon
from package.plot_polygon_viewer_v2 import (
    load_and_convert_data,
    read_select_transcripts_csv,
    load_image,
    show_whole_fov,
)

from rtree import index 


def process_cells(metadata_df: pd.DataFrame, trans_df: pd.DataFrame) -> pd.DataFrame:
    result_df = pd.concat(process_cell_id(metadata_df, trans_df, target_id) for target_id in metadata_df["cell_id"].unique() if target_id != "0")
    return result_df


def process_cell_id(meta_df: pd.DataFrame, trans_df: pd.DataFrame, target_id:str) -> Union[pd.DataFrame, None]:
    if target_id == "0":
        return None

    cell_info = meta_df[meta_df["cell_id"] == target_id].iloc[0]
    bbox_xmin, bbox_xmax = cell_info["local_pixel_min_x"], cell_info["local_pixel_max_x"]
    bbox_ymin, bbox_ymax = cell_info["local_pixel_min_y"], cell_info["local_pixel_max_y"]
    bbox_xcenter, bbox_ycenter = cell_info["local_pixel_centroid_x"], cell_info["local_pixel_centroid_y"]

    bbox_df = trans_df[(trans_df["local_pixel_x"] >= bbox_xmin) & (trans_df["local_pixel_x"] <= bbox_xmax) &
                       (trans_df["local_pixel_y"] >= bbox_ymin) & (trans_df["local_pixel_y"] <= bbox_ymax)]

    return calculate_distances(bbox_df, bbox_xcenter, bbox_ycenter)


def calculate_distances(bbox_df:pd.DataFrame, bbox_xcenter:float, bbox_ycenter:float) -> pd.DataFrame:
    bbox_df_copy = bbox_df.copy()
    bbox_df_copy["x_diff"] = np.abs(bbox_df_copy["local_pixel_x"] - bbox_xcenter)
    bbox_df_copy["y_diff"] = np.abs(bbox_df_copy["local_pixel_y"] - bbox_ycenter)
    return bbox_df_copy


def calc_norm(x:float, y:float) -> float:
    return np.hypot(x, y)


def calculate_plot_data(result_df: pd.DataFrame, trans_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
    grouped_df = result_df.groupby(["gene"], dropna=False)
    centroids_df = grouped_df[["x_diff", "y_diff"]].mean()

    ratio, total_count = calculate_gene_ratio(trans_df)
    centroids_df["center_diff"] = centroids_df.apply(lambda row: calc_norm(row["x_diff"], row["y_diff"]), axis=1)
   
    plot_df = centroids_df.copy()
    plot_df['ratio'] = ratio
    plot_df['total_count'] = total_count
    return plot_df, ratio


def calculate_gene_ratio(trans_df: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    zero_count = trans_df[trans_df["cell_id"] == "0"].groupby("gene")["cell_id"].count()
    total_count = trans_df.groupby("gene")["cell_id"].count()

    # 0除算防止
    ratio = zero_count / total_count if not total_count.empty else pd.Series(dtype=float)
    return ratio, total_count


def show_bubble_chart(
    ax: plt.Axes,  # 外部からAxesを渡す
    df: pd.DataFrame,
    x_metric: str,
    y_metric: str,
    size_metric: str,
    title: str,
    fov: int,
    filter_condition: Union[None, Callable[[pd.DataFrame], pd.Index]] = None,
    is_save: bool = False
 ) -> None:

    # フィルタリングが指定されていれば適用
    if filter_condition is not None:
        if callable(filter_condition): # フィルタリング条件が関数の場合
            selected_indices = filter_condition(df)
        else: # filter_conditionが直接Indexまたはリストを指定する場合
            selected_indices = filter_condition

        df = df.loc[selected_indices]

    x_values = df[x_metric]
    y_values = df[y_metric]
    bubble_sizes = df[size_metric]

    _scatter = ax.scatter(x_values, y_values, s=bubble_sizes, alpha=0.2, c="#0C7CC1")

    # バブルにラベルを追加する
    label_list = df .nlargest(50, size_metric).index.tolist()
    for i, label in enumerate(df.index):
        if label in label_list:
            ax.text(x_values[i], y_values[i], label, fontsize=5, ha='center', va='center')

    ax.set_xlabel(f"{x_metric}")
    ax.set_ylabel(f"{y_metric}")
    ax.set_title(title)

    # x軸およびy軸の範囲を手動で設定して左下を(0, 0)に
    y_min, y_max = ax.get_ylim()  # 現在のy軸の最小値と最大値を取得
    x_min, x_max = ax.get_xlim()  # 現在のy軸の最小値と最大値を取得
    # ax.set_xlim(0, x_max) 
    # ax.set_ylim(0, y_max) 

    if is_save:
        tmp_dir = "/work/tmp_output/"
        filename = os.path.join(tmp_dir, f"fov{fov}_show_bubble_chart.png")
        plt.savefig(filename, format="png")
        print(f"saved {filename}")

    return 


"""
def show_selected_genes(
    img_dir: str, 
    fov: int, 
    z: int, 
    cell_data_dict: Dict[str, Dict[int, pd.DataFrame]], 
    metadata_df: pd.DataFrame, 
    selected_df: pd.DataFrame, 
    gene_list: List[str]) -> None:
    fig, ax = plt.subplots(figsize=(8, 8))
    img = load_image(img_dir, fov, z, "DAPI")

    show_whole_fov(ax, fov, z, cell_data_dict, metadata_df, selected_df, gene_list, show_plot=True)
    img_crop = img

    ax.set_title("Filtered & Top-5 Center Diff")
    ax.imshow(img_crop, cmap="gray")
    ax.set_xlim([0, 2048]) 
    ax.set_ylim([2048, 0]) 
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.show()

"""


def calc_cell_centroid(
        cell_data_dict: Dict[str, pd.DataFrame], metadata_df: pd.DataFrame
        ) -> pd.DataFrame:    
    # metadataに対して、各細胞の重心を計算→保存
    centroids = []
    
    # ポリゴンのデータを使って重心を計算
    for cell_id, polygon_data in cell_data_dict.items():
        polygon_coords = polygon_data[["local_pixel_x", "local_pixel_y"]].to_numpy()
        polygon = Polygon(polygon_coords)
        centroid = polygon.centroid  # ポリゴンの重心を計算
        centroids.append((cell_id, centroid))
    
    # ポリゴンID、x座標、y座標を抽出
    cell_ids = [cell_id for cell_id, _ in centroids]
    x_coords = [centroid.x for _, centroid in centroids]
    y_coords = [centroid.y for _, centroid in centroids]
    
    centroid_coords_df = pd.DataFrame({'cell_id': cell_ids, 'local_pixel_centroid_x': x_coords, 'local_pixel_centroid_y': y_coords})
    if centroid_coords_df.empty:
        print("Warning: No centroids calculated.")
        return
        
    meta = metadata_df.merge(centroid_coords_df, on="cell_id", how="left")    
    meta.set_index("cell_id", inplace=True)
    
    meta_df = meta.copy()
    meta_df.reset_index(inplace=True)

    return meta_df


#  遺伝子に細胞IDを割り当てる関数
def assign_cell_id_to_genes(
        selected_df: pd.DataFrame, cell_data_dict: Dict[str, pd.DataFrame]
        ) -> pd.DataFrame:
    points = np.array(selected_df[["local_pixel_x", "local_pixel_y"]])
    point_objects = [Point(p) for p in points]
    
    _id_mapping, reverse_id_mapping, idx = build_polygon_index(cell_data_dict)
    results = query_polygon_membership(cell_data_dict, point_objects, idx, reverse_id_mapping)

    # 結果のリストをselected_dfに追加
    selected_df["cell_id"] = results
    return selected_df


# ポリゴンの境界ボックスを構築し、R-treeインデックスを作成する関数
def build_polygon_index(cell_data_dict: Dict[str, pd.DataFrame]) -> Tuple[Dict[str, int], Dict[int, str], index.Index]:
    id_mapping = {}
    id_counter = 1
    idx = index.Index()

    for polygon_id, polygon_data in cell_data_dict.items():
        polygon_coords = polygon_data[["local_pixel_x", "local_pixel_y"]].to_numpy()
        polygon = Polygon(polygon_coords)
        minx, miny, maxx, maxy = polygon.bounds
        id_mapping[polygon_id] = id_counter
        idx.insert(id_counter, (minx, miny, maxx, maxy))
        id_counter += 1

    reverse_id_mapping = {v: k for k, v in id_mapping.items()}
    return id_mapping, reverse_id_mapping, idx


# ポリゴンの境界ボックスをクエリし、遺伝子に細胞IDを割り当てる関数
def query_polygon_membership(
        cell_data_dict: Dict[str, pd.DataFrame], 
        point_objects: List[Point], 
        idx: index.Index, 
        reverse_id_mapping: Dict[int, str]
        ) -> List[Union[str, int]]:
    results = []

    for point_obj in point_objects:
        found = False
        for polygon_id_int in idx.intersection((point_obj.x, point_obj.y)):
            polygon_id = reverse_id_mapping[polygon_id_int]
            polygon_coords = cell_data_dict[polygon_id][["local_pixel_x", "local_pixel_y"]].to_numpy()
            polygon = Polygon(polygon_coords)
            if polygon.contains(point_obj):
                results.append(polygon_id)
                found = True
                break
        if not found:
            results.append("0")  # 0を文字列として格納

    return results


def display_bubble_chart(
        selected_df: pd.DataFrame,
        cell_data_dict: Dict[str, pd.DataFrame], 
        metadata_df: pd.DataFrame,
        fov: int,
        is_save: bool = False,
        ax: Optional[plt.Axes] = None  # ax パラメータを追加
        ) -> None:
    trans_df = assign_cell_id_to_genes(selected_df, cell_data_dict)
    meta_df = calc_cell_centroid(cell_data_dict, metadata_df)

    # 細胞ごとに重心からのずれを計算
    result_df = process_cells(meta_df, trans_df)
    
    # total count列とoutter ratio列の追加
    plot_df, _ratio = calculate_plot_data(result_df, trans_df)
    
    # バブルチャートの作成
    custom_filter_condition = lambda df: df['total_count'].nlargest(50).index
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))    
    show_bubble_chart(ax, plot_df, 'center_diff', 'ratio', 'total_count', f'ubble Chart - fov{fov}', fov, filter_condition=custom_filter_condition, is_save=is_save)
    if ax is None:
        plt.show()


# ====================================================================
if __name__ == "__main__":
    input_dir = "/work/datasets/"
    output_dir = "/work/outputs/"

    target_fov = 801
    target_z = 2  # ターゲットの z 座標を設定
    
    # データの用意。
    sub_dir = os.path.join(output_dir, f"subset{target_fov}")
    img_dir = os.path.join(sub_dir, "images")

    # cell_data_dict -> polygonデータ（辞書型）
    # selected_df -> pointデータ（DataFrame）
    cell_data_dict, metadata_df = load_and_convert_data(input_dir, sub_dir, target_fov)
    selected_df, gene_name_list = read_select_transcripts_csv(sub_dir, target_fov, target_z)
    
    display_bubble_chart(selected_df, cell_data_dict, metadata_df, target_fov)

