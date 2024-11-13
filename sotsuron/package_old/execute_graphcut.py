import os
import time
import cv2
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from shapely.geometry import Polygon, LineString

import math
from PIL import Image
from typing import Optional, List

from package_old.save_read_file import read_npy, save_list_as_pkl
from package_old.graphcut_masks_utils import visualize_polygons, load_image, show_ndarray

def read_cellpose_npy(cellpose_path):
    cellpose_data = np.load(cellpose_path, allow_pickle=True).item()
    flows = cellpose_data.get("flows")
    cell_prob = flows[1][0]
    return cell_prob

def exec_laplacian(image, ksize=7):
    image_filtered = cv2.Laplacian(image, cv2.CV_8U, ksize=ksize)
    show_ndarray(image_filtered, f"laplacian ksize={ksize}")
    return image_filtered

def exec_median(image, ksize):
    image_filtered = cv2.medianBlur(image, ksize=ksize)
    show_ndarray(image_filtered, f"median ksize={ksize}")
    return image_filtered

def show_graphcut(partition, height, width):
    binary_image = np.zeros((height, width), dtype=np.uint8)

    source_nodes = set(partition[0])

    # ソースセットの領域を白に、ターゲットセットの領域を黒に設定
    for i in range(height):
        for j in range(width):
            if (i, j) in source_nodes:
                binary_image[i, j] = 255  # 白
    
    show_ndarray(binary_image, "graphcut result")
    return


def make_polygon(G, partition, height, width) -> List[np.ndarray]:
    ignored_nodes = ["obj_terminal", "bkg_terminal"]
    source_graph_ig = G.subgraph(partition[0]).copy()
    source_graph_ig.remove_nodes_from(ignored_nodes)

    polygon_data = []

    # 各連結成分ごとに多角形を生成し描画
    for component in nx.connected_components(source_graph_ig):
        source_nodes = list(component)
        source_positions = [(j, i) for i, j in source_nodes]

        if len(source_positions) < 4:
            continue

        # 凸包を求める
        convex_hull = Polygon(source_positions).convex_hull

        if isinstance(convex_hull, LineString):
            continue
        
        # 凸包の頂点座標をリストに追加
        polygon_data.append(np.array(convex_hull.exterior.xy).T.astype(int))
        
    visualize_polygons(polygon_data, "polygons", height, width, color="blue")
    return polygon_data


def filtering_polygon(polygon_data, height, width, area_threshold1=500, area_threshold2 = 500) -> List[np.ndarray]:
    filtered_polygon_data = [polygon for polygon in polygon_data if Polygon(polygon).area > area_threshold1]
    
    result = []

    for i, poly_A in enumerate(filtered_polygon_data):
        # 内外判定
        is_inside = any(Polygon(poly_A).within(Polygon(poly_B)) for j, poly_B in enumerate(filtered_polygon_data) if i != j)

        # 重なり判定と面積の比較: 
        overlaps_and_larger = any(
            Polygon(poly_A).intersects(Polygon(poly_C)) and Polygon(poly_A).area < Polygon(poly_C).area
            and Polygon(poly_A).intersection(Polygon(poly_C)).area >= area_threshold2
            for k, poly_C in enumerate(filtered_polygon_data) if i != k
        )

        # 条件を満たさない場合、resultに追加
        if not is_inside and not overlaps_and_larger:
            result.append(np.array(poly_A))

    visualize_polygons(result, "filtering", height, width, color="blue")    

    return filtered_polygon_data


def smooth_term(p_att:float, q_att:float, lam=0.001, kappa=0.1) -> float:
    if p_att == q_att:
        return 0.0

    p = max(0,1.0 / (1.0 + math.exp(-10*p_att/255.0  + 5)))
    q = max(0,1.0 / (1.0 + math.exp(-10*q_att/255.0  + 5)))
    
    cost = lam * math.exp(-kappa * (p - q)**2) # 本当は２点の距離で割る -> 今は必ず１なので割愛
    return cost


def calc_intensity(cell_prob, edge_fil, density, alpha=1.0, beta=1.0, gamma=40) -> np.ndarray:
    X = alpha*cell_prob
    Y = beta*edge_fil 
    Z = gamma* density

    tmp = np.zeros((500, 500), dtype=np.float64)
    tmp = np.minimum(X+Z, 255)
    tmp = np.maximum(tmp-Y, 0).astype(int)

    show_ndarray(tmp, "calculation_intensity")
    return tmp


# T -> p：コストが小さいほどobj
def calc_data_obj_matrix(prob:np.ndarray) -> np.ndarray: 
    epsilon = 0.001
    data_obj = np.maximum(0.0, -np.log(prob + epsilon)) 
    data_obj = np.clip(data_obj, 0.0, np.finfo(float).max) 
    data_obj = data_obj.astype(float)
    show_ndarray(data_obj, "data_obj")
    return data_obj

#  S -> p：コストが小さいほどbkg
def calc_data_bkg_matrix(prob:np.ndarray) -> np.ndarray:
    epsilon = 0.001
    data_bkg = np.maximum(0, -np.log(1 - prob + epsilon)) 
    data_bkg = np.clip(data_bkg, 0.0, np.finfo(float).max) 
    data_bkg = data_bkg.astype(float)
    show_ndarray(data_bkg, "data_bkg")
    return data_bkg

def make_graph(origin_map:np.ndarray, intensity:np.ndarray):
    height, width = origin_map.shape
    # グラフの作成
    G = nx.Graph()

    # ピクセルのノードを追加
    for i in range(height):
        for j in range(width):
            node_id = (i, j)         
            G.add_node(node_id, prob=origin_map[i][j], intensity=intensity[i][j])
    # 特別な頂点の作成：始点(source)、終点(target)の作成
    source = "obj_terminal"
    target = "bkg_terminal"
    G.add_node(source, prob=255) # 前景側
    G.add_node(target, prob=0)   # 背景側


    # n-linkを張る：隣接ピクセル間にエッジを張る
    for i in range(height):
        for j in range(width):
            node_id = (i, j)
            neighbors = [(i-1, j), (i+1, j), (i, j-1), (i, j+1), (i-1, j-1), (i-1, j+1), (i+1, j-1), (i+1, j+1)]
            for neighbor in neighbors:
                if 0 <= neighbor[0] < height and 0 <= neighbor[1] < width:
                    if neighbor in G[node_id]:
                        continue
                    G.add_edge(node_id, neighbor, weight=smooth_term(G.nodes[node_id]["intensity"], G.nodes[neighbor]["intensity"]))


    # t-linkを張る：終点とのエッジ
    # sigmoid
    sigmoid_i = np.maximum(0,1.0 / (1.0 + np.exp(-10*intensity/255.0  + 5))).astype(np.float64)

    data_bkg = calc_data_bkg_matrix(sigmoid_i)
    data_obj = calc_data_obj_matrix(sigmoid_i)

    inf = float(2.0*height*width)

    for i in range(height):
        for j in range(width):
            node_id = (i, j)
            node_prob = int(G.nodes[node_id]['intensity'])

            if node_prob <= 0.01:   # bkg
                G.add_edge(node_id, "obj_terminal", weight=0.0)
                G.add_edge(node_id, "bkg_terminal", weight=inf)            
            elif node_prob >= 0.99: # obj
                G.add_edge(node_id, "obj_terminal", weight=inf)
                G.add_edge(node_id, "bkg_terminal", weight=0.0)    
            else:                  # other
                G.add_edge(node_id, "obj_terminal", weight=data_bkg[i][j])
                G.add_edge(node_id, "bkg_terminal", weight=data_obj[i][j])
    print(G.number_of_nodes(), G.number_of_edges()) 
    return G

def make_filter_polygons(G, partition, height, width, thershold1, thershold2) -> List[np.ndarray]:
    polygon_data = make_polygon(G, partition, height, width)
    filtered_polygon_data = filtering_polygon(polygon_data, height, width, thershold1, thershold2)
    return filtered_polygon_data

def show_polygons_with_image(fov, polygons, left, right, upper, lower):
    output_dir = "/work/outputs/a/"
    img_dir = os.path.join(output_dir, f"subset{fov}/images")
    img = load_image(img_dir, fov, z=3, base="DAPI")
    fig, ax = plt.subplots(figsize=(4,4))
    img = img.crop((left, upper, right, lower))

    ax.imshow(img, "gray")
    visualize_polygons(polygons, f"fov{fov}", right-left, lower-upper, ax)
    plt.show()
    return


if __name__ == "__main__":
    # 定数を宣言する
    fov = 237 # 237, 845, 63
    output_dir = "/work/outputs/a/"
    _input_dir = "/work/datasets"
    experiment_dir = f"/work/experiment/fov{fov}"

    base = "DAPI"
    cellpose_path = os.path.join(output_dir, f"subset{fov}/cellpose_mask_{base}_cyto/mosaic_{base}_z3_subset{fov}_seg.npy")
    cell_prob_map = read_cellpose_npy(cellpose_path)
    show_ndarray(cell_prob_map, f"cell probabiliry map({base})")

    gene_count_path = os.path.join(experiment_dir, f"gene_count_matrix_fov{fov}.npy")
    density_array = read_npy(gene_count_path)
    density = density_array.T

    left, right  = 0, 2048
    upper, lower = 0, 2048
    width, height = right-left, lower-upper


    origin_map = cell_prob_map[left:right, upper:lower]
    density_crop = density[left:right, upper:lower]

    # フィルタ処理
    med1 = exec_median(origin_map, ksize = 1)
    lap = exec_laplacian(med1, ksize = 7)
    med2 = exec_median(lap, ksize = 11)


    # I_i,j：ピクセル(i,j)での細胞領域らしさ
    intensity = calc_intensity(origin_map, med2, density_crop, alpha=1.0, beta=1.2, gamma=40)

    show_ndarray(intensity, "test")

    # グラフ構築 --------------------------------
    G = make_graph(origin_map, intensity)

    # グラフカットの実行 -------------------------------
    start = time.time()
    _cut_value, partition = nx.minimum_cut(G, "obj_terminal", "bkg_terminal", capacity='weight')
    end = time.time()
    print("time: ", end-start)
    # source_graph = G.subgraph(partition[0])
    # target_graph = G.subgraph(partition[1])
    show_graphcut(partition, height, width)

    result_polygons = make_filter_polygons(G, partition, height, width, 0, 0)

    show_polygons_with_image(fov, result_polygons, left, right, upper, lower)

    # リストの保存
    save_list_path = os.path.join(experiment_dir, f"polygon_list_propose_{base}_fov{fov}.pkl")
    # save_list_as_pkl(save_list_path, result_polygons)