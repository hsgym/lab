import os

import cv2
import numpy as np
import matplotlib.pyplot as plt

from itertools import product
from joblib import Parallel, delayed
from typing import List, Tuple

from skimage.metrics import (adapted_rand_error,
                              variation_of_information)

from package_old.save_read_file import read_npy, read_pkl
from package_old.graphcut_masks_utils import load_image, visualize_polygons
from package_old.execute_graphcut import filtering_polygon

# iou計算
def calculate_iou(args: Tuple[np.ndarray, np.ndarray]) -> float:
    mask1, mask2 = args
    intersection_area = np.logical_and(mask1, mask2).sum()
    union_area = np.logical_or(mask1, mask2).sum()
    iou = intersection_area / union_area if union_area > 0 else 0
    return iou

# マスクを作成する関数
def generate_masks_from_polygon(polygons: List[np.ndarray], image_shape: Tuple[int, int]) -> Tuple[np.ndarray, np.ndarray]:
    length = len(polygons)
    masks = np.zeros((len(polygons), image_shape[0], image_shape[1]), dtype=np.uint8)
    all_masks = np.zeros((image_shape[0], image_shape[1]), dtype=np.uint8)

    for i, poly in enumerate(polygons):
        try:
            # ポリゴンが閉じた形状になるように最初の座標を追加
            poly = np.vstack([poly, poly[0]])
            poly = poly.astype(np.int32)
            cv2.fillPoly(masks[i], [poly], i+length+1)
            all_masks += masks[i]
        except Exception as e:
            continue

    return masks, all_masks


def calculate_average_iou(masks_ans: np.ndarray, masks_method: np.ndarray) -> float:
    iou_values = Parallel(n_jobs=-1)(delayed(calculate_iou)(args) for args in product(masks_ans, masks_method))
    average_iou = np.mean([iou for iou in iou_values if iou > 0])
    return average_iou


def print_result(method: str, polygons: List[np.ndarray], all_masks: np.ndarray, gene_count: np.ndarray, iou: float = 0) -> None:
    total_gene_count = np.sum(gene_count)
    outside_count = np.sum(gene_count[all_masks == 0])
    inside_count = np.sum(gene_count[all_masks >= 1])
    cell_count = len(polygons)
    inside_rate = inside_count / total_gene_count
    gene_per_cell = inside_count / cell_count

    if method == "answer":
        print("total_gene_count: ", total_gene_count)
        print("# ---------------------- \n  answer:")
        print("    cell_count: ", cell_count)
        print("    gene_count: ")
        print("      inside_count:  ", inside_count)
        print("      outside_count: ", outside_count)
    else:
        print(f"# --------------------- \n  {method}:")
        print("    cell_count: ", cell_count)
        print("    gene_count: ")
        print("      inside_count:  ", inside_count)
        print("      outside_count: ", outside_count)
        print("    average_IoU:      ", iou)
    
    print("    inside_gene_rate: ", inside_rate)
    print("    gene_per_cell:    ", gene_per_cell)

    return


def evaluate_skimage(all_masks_ans: np.ndarray, all_masks_method: np.ndarray) -> None:
    error, precision, recall = adapted_rand_error(all_masks_ans, all_masks_method)
    splits, merges = variation_of_information(all_masks_ans, all_masks_method)
    
    print(f'    skimage:')
    print(f'      Adapted_Rand_error: {error}')
    print(f'      Adapted_Rand_precision: {precision}')
    print(f'      Adapted_Rand_recall: {recall}')
    print(f'      False_Splits: {splits}')
    print(f'      False_Merges: {merges}')
    return 

def show_polygons_imgae(img, title, width, height, polygons, save_fig, save_path):
    fig, ax = plt.subplots(figsize=(8,8))
    ax.imshow(img, "gray")
    visualize_polygons(polygons, title, width, height, ax)

    if save_fig:
        ax.set_title("")  
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        plt.savefig(save_path, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()
    return

def evaluate_all(polygons, height, width, method, gene_count):
    masks, all_masks = generate_masks_from_polygon(polygons, (height, width))
    print_result(method, polygons, all_masks, gene_count)
    evaluate_skimage(all_masks, all_masks) # test用(1回実行したら消していいかも)
    return 


if __name__ == "__main__":
    fov = 100  #100, 237, 845, 63
    output_dir = "/work/outputs/a/"
    img_dir = os.path.join(output_dir, f"subset{fov}/images")
    input_dir = "/work/datasets"
    experiment_dir = f"/work/experiment/fov{fov}"
    ans_dir = "/work/answers"
    width, height = 2048, 2048
    save_dir = os.path.join(experiment_dir, "image")

    save_fig = False ### 
    base = "DAPI"  # or "Cellbound3"

    # 画像の確認
    img = load_image(img_dir, fov, z=3, base=base)
    fig, ax = plt.subplots(figsize=(8,8))
    ax.imshow(img, "gray")
    ax.set_xlim([0, width])
    ax.set_ylim([height, 0])
    ax.set_title(f"fov{fov}")
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.show()


    gene_count_path = os.path.join(experiment_dir, f"gene_count_matrix_fov{fov}.npy")
    gene_count = read_npy(gene_count_path)
    gene_count = gene_count.T


    # 正解データについて
    method = "answer"
    polygon_ans_path = os.path.join(experiment_dir, f"polygon_list_answer_{base}_fov{fov}.pkl")
    polygons_ans = read_pkl(polygon_ans_path)
    path = os.path.join(save_dir, f"fov{fov}_{base}_{method}.png")
    show_polygons_imgae(img, f"fov{fov} - {method}", width, height, polygons_ans, save_fig, path)
    print(f"== fov{fov}, {base} ==")
    evaluate_all(polygons_ans, height, width, method, gene_count)

    # 提案手法について
    method =  "propose" 
    polygon_method_path = os.path.join(experiment_dir, f"polygon_list_{method}_{base}_fov{fov}_v3.pkl")
    polygons_propose_0 = read_pkl(polygon_method_path)
    polygons_propose = filtering_polygon(polygons_propose_0, height, width, 1000, 100)
    path = os.path.join(save_dir, f"fov{fov}_{base}_{method}.png")
    show_polygons_imgae(img, f"fov{fov} - {method}", width, height, polygons_propose, save_fig, path)
    print(f"== fov{fov}, {base} ==")
    evaluate_all(polygons_propose, height, width, method, gene_count)


    # cellposeについて
    method = "cellpose"
    polygon_method_path = os.path.join(experiment_dir, f"polygon_list_{method}_{base}_fov{fov}.pkl")
    polygons_cellpose = read_pkl(polygon_method_path)
    path = os.path.join(save_dir, f"fov{fov}_{base}_{method}.png")
    show_polygons_imgae(img, f"fov{fov} - {method}", width, height, polygons_cellpose, save_fig, path)
    print(f"== fov{fov}, {base} ==")
    evaluate_all(polygons_cellpose, height, width, method, gene_count)
    
    
    # merscope について
    method = "merscope"
    polygon_method_path = os.path.join(experiment_dir, f"polygon_list_{method}_{base}_fov{fov}.pkl")
    polygons_merscope = read_pkl(polygon_method_path)
    path = os.path.join(save_dir, f"fov{fov}_{base}_{method}.png")
    show_polygons_imgae(img, f"fov{fov} - {method}", width, height, polygons_merscope, save_fig, path)
    print(f"== fov{fov}, {base} ==")
    evaluate_all(polygons_merscope, height, width, method, gene_count)

