import pandas as pd
import numpy as np
import os
import scanpy as sc
from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
import os
import os
from typing import Optional, List
import seaborn as sns


from package.viewer_common import CommonPlotPolygonViewer

from package.myLP import my_LP_test
from package.reader_merscope import MerscopeDataReader
from package.reader_xenium import XeniumDataReader

# 参照データの平均化
def process_reference_data(path:Optional[str]=None) -> pd.DataFrame:
    if path is None:
        path = "/work/datasets/reference/adata.h5ad"

    adata = sc.read_h5ad(path)
    adata_df = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)
    adata_df["cell"] = adata.obs["CellAssign_prediction"]
    group_counts = adata.obs.groupby("CellAssign_prediction", as_index=False).size()
    group_counts = group_counts.set_index('CellAssign_prediction')
    
    adata_df = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)
    adata_df["CellAssign_prediction"] = adata.obs["CellAssign_prediction"]
    celltype_means = adata_df.groupby("CellAssign_prediction").mean()
    return celltype_means.T


def bin_data(gene_df:pd.DataFrame, width, height, num_x_bins:int, num_y_bins:int) -> pd.DataFrame:
    x_bins = np.linspace(0, width, num_x_bins + 1)
    y_bins = np.linspace(0, height, num_y_bins + 1)
    # 以下を警告を回避するために変更
    gene_df['x_bin'] = np.digitize(gene_df['local_pixel_x'], bins=x_bins) - 1
    gene_df['y_bin'] = np.digitize(gene_df['local_pixel_y'], bins=y_bins) - 1
    return gene_df

# binごとに遺伝子発現量を集計
def aggregate_expression(gene_df_bin:pd.DataFrame, num_x_bins, num_y_bins) -> pd.DataFrame:
    """
    region_expression = gene_df_bin.groupby(['x_bin', 'y_bin', 'gene']).size().unstack(fill_value=0)
    region_expression_reset = region_expression.reset_index()
    region_expression_pivot = region_expression_reset.melt(id_vars=['x_bin', 'y_bin'], var_name='gene', value_name='expression')
    region_expression_final = region_expression_pivot.pivot_table(index=['x_bin', 'y_bin'], columns='gene', values='expression', fill_value=0)
    region_expression_final.index = region_expression_final.index.map(lambda x: f"{x[0]}_{x[1]}")
    new_region_expression_common = region_expression_final.set_index(region_expression_final.index)
    new_region_expression_common.index.name = "bin_x_y"
    """
    region_expression = gene_df_bin.groupby(['x_bin', 'y_bin', 'gene']).size().unstack(fill_value=0)
    region_expression_reset = region_expression.reset_index()
    region_expression_pivot = region_expression_reset.melt(id_vars=['x_bin', 'y_bin'], var_name='gene', value_name='expression')
    region_expression_final = region_expression_pivot.pivot_table(index=['x_bin', 'y_bin'], columns='gene', values='expression', fill_value=0)

    # 全ての可能なbinの組み合わせを作成
    all_bins = pd.MultiIndex.from_product([range(num_x_bins), range(num_y_bins)], names=['x_bin', 'y_bin'])
    region_expression_final = region_expression_final.reindex(all_bins, fill_value=0)

    region_expression_final.index = region_expression_final.index.map(lambda x: f"{x[0]}_{x[1]}")
    new_region_expression_common = region_expression_final.set_index(region_expression_final.index)
    new_region_expression_common.index.name = "bin_x_y"

    return new_region_expression_common


def process_spatial_data(gene_df:pd.DataFrame, num_x_bins:int, num_y_bins:int, width:int, height:int) -> pd.DataFrame:
    gene_df_bin = bin_data(gene_df, width, height, num_x_bins, num_y_bins)
    aggregate_bins_data = aggregate_expression(gene_df_bin,  num_x_bins, num_y_bins)
    return aggregate_bins_data

def common_rows(celltype_df:pd.DataFrame, bin_df:pd.DataFrame):
    spe_genes = bin_df.columns.values
    ref_genes = celltype_df.columns.values
    ref_genes = np.delete(ref_genes, np.where(ref_genes == "groupby_size")[0], axis=0)
    common_genes = list(set(spe_genes).intersection(set(ref_genes)))
    bin_df_com = bin_df[common_genes]
    celltype_df_com = celltype_df[common_genes]
    return bin_df_com, celltype_df_com

def read_or_process_reference_data(ref_path, reference_data_path:Optional[str]=None):
    
    if os.path.exists(ref_path):
        return pd.read_pickle(ref_path).T
    if reference_data_path is None:
        raise ValueError("reference_data_path is required.")
    return process_reference_data(reference_data_path)

def read_or_process_merscope_data(spe_path, input_dir, output_dir, fov, width, height, z, num_x_bins, num_y_bins):
    # if os.path.exists(spe_path):
    #    return pd.read_pickle(spe_path).T
    reader = MerscopeDataReader(input_dir, output_dir, fov, width, height)
    gene_df, _gene_name_list = reader.read_gene_data()
    return process_spatial_data(gene_df, num_x_bins, num_y_bins, width, height)

def read_or_process_xenium_data(spe_path, input_dir, output_dir, fov, width, height, resampling_factor, num_x_bins, num_y_bins):
    # if os.path.exists(spe_path):
    #     return pd.read_pickle(spe_path).T
    reader = XeniumDataReader(input_dir, output_dir, fov, width, height, resampling_factor)
    gene_df, _gene_name_list = reader.read_transcripts_data()
    return process_spatial_data(gene_df, num_x_bins, num_y_bins, width, height)

def main_merscope(ref_dir, input_dir, output_dir, fov, width, height, z, num_x_bins, num_y_bins, reference_data_path:Optional[str]=None):
    ref_path = os.path.join(ref_dir, f"ref_input_merscope_subset{fov}.pkl")
    spe_path = os.path.join(output_dir, f"subset{fov}/detected_transcripts_subset{fov}_bin{num_x_bins}x{num_y_bins}.pkl")

    if os.path.exists(ref_path) and os.path.exists(spe_path):
        ref_df = pd.read_pickle(ref_path)
        spe_df = pd.read_pickle(spe_path)
        return ref_df, spe_df

    ref_df = read_or_process_reference_data(ref_path, reference_data_path)
    spe_df = read_or_process_merscope_data(spe_path, input_dir, output_dir, fov, width, height, z, num_x_bins, num_y_bins)

    bin_df_com, celltype_df_com = common_rows(ref_df, spe_df)
    ref_df_com = celltype_df_com.copy().T
    spe_df_com = bin_df_com.copy().T

    # if not os.path.exists(ref_path):
    ref_df_com.to_pickle(ref_path)
    # if not os.path.exists(spe_path):
    spe_df_com.to_pickle(spe_path)
    return ref_df_com, spe_df_com


def main_xenium(ref_dir, input_dir, output_dir, fov, width, height, resampling_factor,  num_x_bins, num_y_bins, reference_data_path:Optional[str]=None):
    ref_path = os.path.join(ref_dir, f"ref_input_xenium_subset{fov}.pkl")
    spe_path = os.path.join(output_dir, f"subset{fov}/detected_transcripts_subset{fov}_bin{num_x_bins}x{num_y_bins}.pkl")
    
    if os.path.exists(ref_path) and os.path.exists(spe_path):
        ref_df = pd.read_pickle(ref_path)
        spe_df = pd.read_pickle(spe_path)
        return ref_df, spe_df

    ref_df = read_or_process_reference_data(ref_path, reference_data_path)
    spe_df = read_or_process_xenium_data(spe_path, input_dir, output_dir, fov, width, height, resampling_factor, num_x_bins, num_y_bins)

    bin_df_com, celltype_df_com = common_rows(ref_df, spe_df)
    ref_df_com = celltype_df_com.copy().T
    spe_df_com = bin_df_com.copy().T

  #  if not os.path.exists(ref_path):
    ref_df_com.to_pickle(ref_path)
    
  #  if not os.path.exists(spe_path):
    spe_df_com.to_pickle(spe_path)
    return ref_df_com, spe_df_com


# LP の結果を可視化する
def visualize_heatmap(LP_df:pd.DataFrame):
    fig, ax = plt.subplots(figsize=(10,10))
    sns.heatmap(LP_df,  cmap="RdBu", cbar=True)
    ax.set_title('Heatmap of Proportion')
    ax.set_xlabel('cell type')
    ax.set_ylabel('grid')
    plt.tight_layout()
    plt.show()
    return

def predict_celltype(LP_df:pd.DataFrame):
    LP_type_df = pd.DataFrame(index=LP_df.index, columns=["predicted_type"])
    for i in range(LP_df.shape[0]):
        row = LP_df.iloc[i]
        max_type = row.idxmax()
        if row[max_type] == 0:
            LP_type_df.iloc[i] = "Background"
        else:
            LP_type_df.iloc[i] = max_type
    return LP_type_df

def visualize_tile(LP_type_df:pd.DataFrame, num_x_bins, num_y_bins, plotter:Optional[object]=None):
    title = "Predicted Cell Types"
    unique_cell_types = LP_type_df['predicted_type'].unique()
    category_colors = dict(zip(unique_cell_types, sns.color_palette('hls', len(unique_cell_types))))

    fig, ax = plt.subplots()# figsize=(8, 6))

    for id, row in LP_type_df.iterrows():
        i, j = id.split("_")
        x_corner = int(i) * (plotter.width / num_x_bins)
        y_corner = int(j) * (plotter.height / num_y_bins)
        
        color = category_colors[row['predicted_type']]
        
        ax.fill([x_corner, x_corner + plotter.width / num_x_bins, x_corner + plotter.width / num_x_bins, x_corner],
                [y_corner, y_corner, y_corner + plotter.height / num_y_bins, y_corner + plotter.height / num_y_bins],
                color=color, edgecolor='black')

    handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in category_colors.values()]
    labels = list(category_colors.keys())
    ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1), title='Predicted Cell Type')

    ax.set_xlim(0, plotter.width)
    ax.set_ylim(plotter.height, 0)
    ax.set_title(title)
    fig.tight_layout()
    fig.subplots_adjust(right=0.75)
    plt.show()
    return

def visualize_piechart(LP_df:pd.DataFrame, num_x_bins, num_y_bins, is_background:bool=False, plotter:Optional[object]=None, png_path:Optional[str]=None):
    
    # 列名と色の対応付け
    num_columns = len(LP_df.columns)
    colors = sns.color_palette("hls", num_columns)
    column_color_map = dict(zip(LP_df.columns, colors))

    # パイチャートの可視化
    fig, axes = plt.subplots(num_x_bins, num_y_bins, figsize=(14, 14), gridspec_kw={'wspace': 0, 'hspace': 0})
    fig.suptitle("pie chart of cell type")
    
    if is_background:
        img_fov = plotter.load_image(png_path)
        fig_background = fig.add_subplot(111, zorder=-1)
        fig_background.imshow(img_fov, cmap="gray")
        fig_background.axis('off')  # 背景の軸を非表示に

    for id, row in LP_df.iterrows():
        i, j = id.split("_")
        ax = axes[int(j), int(i)]
        values = row.values

        # フィルタリング
        filtered_values = [value for value in values if value != 0]
        filtered_labels = [label for value, label in zip(values, LP_df.columns) if value != 0]
        filtered_colors = [column_color_map[label] for label in filtered_labels]

        if filtered_values:
            ax.pie(filtered_values, colors=filtered_colors, startangle=90, counterclock=False)
        else:
            ax.axis("off")

    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=column_color_map[col], markersize=10) for col in LP_df.columns]
    labels = LP_df.columns

    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), title="cell type")

    plt.subplots_adjust(right=0.85) 
    plt.tight_layout()
    plt.show()
    return

"""
if __name__ == "__main__":
    # MERSCOPEの場合 -----
    input_dir = "/work/datasets/Yahara/202304161129_MsFetusHumerus-23003-2-VS71-YS_VMSC07201/region_0"
    output_dir="/work/output_MERSCOPE"
    ref_dir = "/work/my_data/reference"
    fov=236
    width=2048
    height=2048
    z=3
    num_x_bins=8
    num_y_bins=8
    reference_data_path="/work/datasets/reference/adata.h5ad"
    ref_df_com, spe_df_com = main_merscope(ref_dir, input_dir, output_dir, fov, width, height, z, num_x_bins, num_y_bins, reference_data_path)

    # Xeniumの場合 -----
    input_dir = "/work/datasets/Okada/output-XETG00130__0014491__Cont_3_1__20240403__094846"
    output_dir = "/work/output_Xenium"
    ref_dir = "/work/my_data/reference"
    resampling_factor = 0.2125
    fov=30
    width=4096
    height=4096
    num_x_bins=8
    num_y_bins=8
    reference_data_path="/work/datasets/reference/adata.h5ad"
    
    # LP の実行 ------
    ref_df_com, spe_df_com = main_xenium(ref_dir, input_dir, output_dir, fov, width, height, resampling_factor, num_x_bins, num_y_bins)
    result = my_LP_test(ref_df_com, spe_df_com , "/work/tmpout")
    LP_df = result[0]
    print(LP_df.head(5))

    # 可視化 -----
    img_dir = os.path.join(output_dir, f"subset{fov}/images")
    plotter = CommonPlotPolygonViewer(img_dir, fov, z, "default", width, height)
    visualize_heatmap(LP_df)
    LP_type_df = predict_celltype(LP_df)
    visualize_tile(LP_type_df)
    visualize_piechart(LP_df, is_background=True, plotter=plotter)
"""