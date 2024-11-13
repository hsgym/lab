import os
import pickle
from typing import Optional, List
import networkx as nx
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from package.myLP import my_LP_test
from package_v2.process_data_for_LP import main_merscope
from package.viewer_common import CommonPlotPolygonViewer
from package_v2.visualize_celltype_with_graph import visualize_predicted_type

import csv

# グラフ構築
# diff[A]がすでにあるならばそれを返す関数
def check_diff_area(node_name: str, pd_df: pd.DataFrame, md_df: pd.DataFrame) -> Optional[float]:
    if node_name not in pd_df.columns:
        return None
    return pd_df[node_name].sum() + md_df[node_name].sum()


# グラフ構築
def create_graph(num_x_bins:int, num_y_bins:int, celltype_df:pd.DataFrame, pd_df, md_df) -> nx.Graph:
    G = nx.Graph()

    # ノードの追加
    for x in range(num_x_bins):
        for y in range(num_y_bins):
            cell_name = f"{x}_{y}"
            G.add_node(cell_name, pos=[(x, y)], cell_type=celltype_df.at[cell_name, "predicted_type"], diff=check_diff_area(cell_name, pd_df, md_df))

    # 隣接エッジの追加
    neighbors_offsets = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]  # 隣接する領域のオフセット
    for x in range(num_x_bins):
        for y in range(num_y_bins):
            for dx, dy in neighbors_offsets:
                xx, yy = x + dx, y + dy
                if (0 <= xx < num_x_bins) and (0 <= yy < num_y_bins) and ((x < xx) or (x == xx and y < yy)):
                        G.add_edge(f"{x}_{y}", f"{xx}_{yy}", att="neighbor")

    print("Node: ", G.number_of_nodes(), " / Edge: ", G.number_of_edges())
    return G


# グラフ集約の候補を挙げるための処理
# 領域Xに対する空間データの取得（領域X=A+Bとかになっている）
def get_merge_spe(G, node:str, neighbor:str, spe_df_com:pd.DataFrame) -> pd.DataFrame:
    # 空間データの取得
    pos_list1 = G.nodes[node]["pos"]
    pos_list2 = G.nodes[neighbor]["pos"]

    # data1とdata2をリストに収集して一度に合計
    data1_list = [spe_df_com[f"{x}_{y}"] for (x, y) in pos_list1]
    data2_list = [spe_df_com[f"{x}_{y}"] for (x, y) in pos_list2]

    data1 = sum(data1_list)
    data2 = sum(data2_list)
    
    spe_merge = pd.DataFrame(data1 + data2, columns=[f"{node}_{neighbor}"])
    return spe_merge


def merge_dataframes(x_merge_df, pd_merge_df, md_merge_df, x_tmp, pd_tmp, md_tmp):
    if x_merge_df is None:
        return x_tmp.copy(), pd_tmp.copy(), md_tmp.copy()
    return (pd.concat([x_merge_df, x_tmp]),
            pd.concat([pd_merge_df, pd_tmp], axis=1),
            pd.concat([md_merge_df, md_tmp], axis=1))


def get_merge_candidates_list(G:nx.Graph, pd_df, md_df, ref_df_com , spe_df_com, save_dir):
    checked_pair = set()
    merge_candidates = []
    x_merge_df, pd_merge_df, md_merge_df = None, None, None

    for node, data in tqdm(G.nodes(data=True)):
        if data["cell_type"] != "Background":
            continue
            
        neighbors = list(nx.all_neighbors(G, node))
        if any(G.nodes[neighbor]["cell_type"] != "Background" for neighbor in neighbors):
            continue

        diff1 = G.nodes[node]["diff"]

        for neighbor in neighbors:
            if (node, neighbor) in checked_pair or (neighbor, node) in checked_pair:
                continue
            checked_pair.add((node, neighbor))
            
            diff2 = G.nodes[neighbor]["diff"]
            diff_sum = diff1 + diff2
            
            name = f"{node}_{neighbor}"
            diff_merge = check_diff_area(name, pd_df, md_df)
            if diff_merge is None:
                spe_merge = get_merge_spe(G, node, neighbor, spe_df_com)
                x_tmp, pd_tmp, md_tmp = my_LP_test(ref_df_com, spe_merge, save_dir)
                x_merge_df, pd_merge_df, md_merge_df = merge_dataframes(x_merge_df, pd_merge_df, md_merge_df, x_tmp, pd_tmp, md_tmp)
                diff_merge = check_diff_area(name, pd_tmp, md_tmp)
        
            if diff_merge < diff_sum:
                merge_candidates.append((node, neighbor, diff_sum - diff_merge))
    return x_merge_df, pd_merge_df, md_merge_df, merge_candidates

def save_results(save_dir, x_merge_df, pd_merge_df, md_merge_df, merge_candidates):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    x_merge_path = os.path.join(save_dir, "x_merge_df.pickle")
    pd_merge_path = os.path.join(save_dir, "pd_merge_df.pickle")
    md_merge_path = os.path.join(save_dir, "md_merge_df.pickle")
    list_merge_path = os.path.join(save_dir, "merge_candidates.pickle")

    with open(x_merge_path, mode='wb') as fo:
        pickle.dump(x_merge_df, fo)

    with open(pd_merge_path, mode='wb') as fo:
        pickle.dump(pd_merge_df, fo)

    with open(md_merge_path, mode='wb') as fo:
        pickle.dump(md_merge_df, fo)

    with open(list_merge_path, mode='wb') as fo:
        pickle.dump(merge_candidates, fo)

    return 


def transfer_to_celltype_df(lp_df: pd.DataFrame) -> pd.DataFrame:
    celltype_df = pd.DataFrame(index=lp_df.index, columns=["predicted_type"])
    for i in range(lp_df.shape[0]):
        row = lp_df.iloc[i]
        max_type = row.idxmax()
        if row[max_type] == 0:
            celltype_df.iloc[i] = "Background"
        else:
            celltype_df.iloc[i] = max_type
    return celltype_df


# グラフ集約操作
def merge_node(G, node1:str, node2:str, diff, celltype_df_merge):
    new_node = f"{node1}_{node2}"
    new_pos = G.nodes[node1]['pos'] + G.nodes[node2]['pos']
    try:
        new_type = celltype_df_merge.at[f"{node1}_{node2}", "predicted_type"]
    except KeyError:
        new_type = celltype_df_merge.at[f"{node2}_{node1}", "predicted_type"]

    # 新しい頂点の追加
    G.add_node(new_node, pos=new_pos, cell_type=new_type, diff=diff)

    # 隣接エッジの張り直し
    for neighbor in set(G.neighbors(node1)).union(set(G.neighbors(node2))):
        if neighbor != new_node:  # 自分自身にはエッジを張らない
            G.add_edge(new_node, neighbor, att="neighbor")
    
    # 元々の頂点とエッジの削除
    G.remove_node(node1)
    G.remove_node(node2)
    return new_node


def process_merge_nodes(G:nx.Graph, merge_candidates:List, celltype_df_merge:pd.DataFrame):
    merge_candidates.sort(key=lambda x: x[2], reverse=True)
    
    merged_nodes = set()
    
    for node1, node2, diff in merge_candidates:
        if node1 in merged_nodes or node2 in merged_nodes:
            continue

        new_node_name = merge_node(G, node1, node2, diff, celltype_df_merge)
        merged_nodes.add(new_node_name)
        merged_nodes.add(node1)
        merged_nodes.add(node2)

    return G

def generate_seaborn_colors(ref_df_com):
    # Seabornのカラーパレットを取得
    palette = sns.color_palette("hls", 21)  # "hls"は色相環に基づいたパレット

    # 特別なカテゴリに対する色を指定
    celltype_list = ref_df_com.columns.tolist()
    osteo_celltype_list = ["Osteoblast", "Osteoclast", "Osteocyte"]
    remain_celltype_list = [i for i in celltype_list if i not in osteo_celltype_list]
    special_categories = {
        osteo_celltype_list[0]: palette[0],  # 特別なカテゴリ1
        osteo_celltype_list[1]: palette[7],  # 特別なカテゴリ2
        osteo_celltype_list[2]: palette[14], # 特別なカテゴリ3
    }

    # 残りのカテゴリの色を設定
    colors = {}
    for index, item in enumerate(remain_celltype_list):
        colors[item] = palette[index]

    # 特別なカテゴリを追加
    colors.update(special_categories)

    return colors

def save_colors_to_csv(colors, csv_path):
    # カテゴリ名をアルファベット順にソート
    sorted_colors = dict(sorted(colors.items()))
    sorted_colors["Background"] = "#AAAAAA"
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['cell_type', 'color'])
        for cell_type, color in sorted_colors.items():
            if cell_type == "Background":
                color_hex = color
            else:
                color_rgb = tuple(int(c * 255) for c in color)
                color_hex = '#{:02x}{:02x}{:02x}'.format(*color_rgb)
            writer.writerow([cell_type, color_hex])

def load_colors_from_csv(csv_path):
    colors = {}
    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # ヘッダーをスキップ
        for row in reader:
            cell_type, color_hex = row
            # 16進数カラーコードをRGB値に変換
            color_rgb = tuple(int(color_hex[i:i+2], 16) / 255.0 for i in (1, 3, 5))
            colors[cell_type] = color_rgb
    return colors


if __name__ == "__main__":
    # グローバル変数の宣言
    input_dir = "/work/datasets/Yahara/202304161129_MsFetusHumerus-23003-2-VS71-YS_VMSC07201/region_0"
    output_dir="/work/output_MERSCOPE"
    ref_dir = "/work/my_data/reference"
    fov=236
    width=2048
    height=2048
    z=3
    reference_data_path="/work/datasets/reference/adata.h5ad"
    num_x_bins, num_y_bins = 64, 64

    lp_dir = f"/work/tmpout/merscope_{num_x_bins}x{num_y_bins}/LP_outputs"
    md_path = os.path.join(lp_dir, "md.csv")
    pd_path = os.path.join(lp_dir, "pd.csv")
    lp_path = os.path.join(lp_dir, "x.csv")


    # 既にあるデータの読み込み
    md_df = pd.read_csv(md_path, index_col=["gene"])
    pd_df = pd.read_csv(pd_path, index_col=["gene"])
    lp_df = pd.read_csv(lp_path, index_col=["bin_x_y"])

    ref_df_com, spe_df_com = main_merscope(ref_dir, input_dir, output_dir, fov, width, height, z, num_x_bins, num_y_bins, reference_data_path)

    celltype_df = transfer_to_celltype_df(lp_df)

    G = create_graph(num_x_bins, num_y_bins, celltype_df, pd_df, md_df)
    G_init = G.copy()

    # 色を生成
    color_dict_csv_path = "/work/category_to_color.csv"
    if os.path.exists(color_dict_csv_path):
        category_to_color = load_colors_from_csv(color_dict_csv_path)
    else:
        category_to_color = generate_seaborn_colors(ref_df_com)
        save_colors_to_csv(category_to_color, color_dict_csv_path)
    

    img_dir = os.path.join(output_dir, f"subset{fov}/images")
    plotter1 = CommonPlotPolygonViewer(img_dir, fov, z, f"bin{num_x_bins}x{num_y_bins}", width, height)
    visualize_predicted_type(G_init, num_x_bins, num_y_bins, category_to_color, False, plotter1)
    visualize_predicted_type(G_init, num_x_bins, num_y_bins, category_to_color, True, plotter1)

    x_merge_df, pd_merge_df, md_merge_df, merge_candidates = get_merge_candidates_list(G, pd_df, md_df, ref_df_com , spe_df_com)
    save_dir = "/work/tmpout"
    save_results(save_dir, x_merge_df, pd_merge_df, md_merge_df, merge_candidates)
    celltype_df_merge = transfer_to_celltype_df(x_merge_df)

    process_merge_nodes(G, merge_candidates, celltype_df_merge)

    plotter2 = CommonPlotPolygonViewer(img_dir, fov, z, f"bin{num_x_bins}x{num_y_bins} merged", width, height)
    visualize_predicted_type(G, num_x_bins, num_y_bins, category_to_color, False, plotter2)
    visualize_predicted_type(G, num_x_bins, num_y_bins, category_to_color, True, plotter2)
    print("Node: ", G.number_of_nodes(), " / Edge: ", G.number_of_edges())

