import os
import scanpy as sc
import pandas as pd
import numpy as np


from package_old.plot_polygon_viewer import read_select_transcripts_csv
from package_old.evaluate_segmentation import generate_masks_from_polygon
from package_old.save_read_file import read_pkl


def add_row_cell_id(df: pd.DataFrame, mask: np.ndarray) -> pd.DataFrame:
    # 新しい列 'object_id' を初期化
    df['cell_id'] = 0

    # 各プロット点について object_id を設定
    df['x_round'] = df['local_pixel_x'].round().astype(int).clip(0, mask.shape[1] - 1)
    df['y_round'] = df['local_pixel_y'].round().astype(int).clip(0, mask.shape[0] - 1)

    # マスクの値を取得して新しい列に格納
    df['cell_id'] = mask[df['y_round'], df['x_round']]
    return df


def transfer_anndata(df:pd.DataFrame):
    obs = pd.DataFrame()
    obs["cell_id"] = df.index
    var_names = df.columns
    var = pd.DataFrame(index=var_names)
    X = df.values
    adata = sc.AnnData(X=X, obs=obs, var=var, dtype="int32")
    print(adata)
    return adata

def process_adata(adata, is_show:bool):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=50)
    sc.tl.leiden(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color="leiden")
    sc.pp.log1p(adata)

    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', key_added='ttest_ranking')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key="ttest_ranking")

    if is_show:
        result = adata.uns['ttest_ranking']
        groups = result['names'].dtype.names
        d = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names']}).head(5)
        print(d)

    sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True, cmap = "viridis", swap_axes=False,  key="ttest_ranking", dendrogram=True)
    return adata

if __name__ == "__main__":
    fov = 63
    output_dir = "/work/outputs/a/"
    sub_dir = os.path.join(output_dir, f"subset{fov}")
    img_dir = os.path.join(output_dir, f"subset{fov}/images")
    input_dir = "/work/datasets"
    experiment_dir = f"/work/experiment/fov{fov}"
    ans_dir = "/work/answers"
    width, height = 2048, 2048

    base = "DAPI" # "DAPI" or "Cellbound3"
    method = "propose"

    selected_df, gene_name_list = read_select_transcripts_csv(sub_dir, fov)
    raw_df = selected_df[["gene", "local_pixel_x", "local_pixel_y", "transcript_id"]]


    polygon_path = os.path.join(experiment_dir, f"polygon_list_{method}_{base}_fov{fov}.pkl")
    polygon_list = read_pkl(polygon_path)
    masks, all_masks = generate_masks_from_polygon(polygon_list, (height, width))

    masks_df = add_row_cell_id(raw_df, all_masks)

    # 行が細胞，列が遺伝子になるはず
    cell_df = masks_df.groupby(["cell_id", "gene"]).size().unstack(fill_value=0)
    print(cell_df)


    # 行indexの値で抽出
    length = len(polygon_list)
    filtered_df = cell_df[(cell_df.index.astype(int) >= 1) & (cell_df.index.astype(int) <= 2 * length + 1)]
    # 結果を表示
    print(filtered_df)

    # anndataへの変換
    adata = transfer_anndata(filtered_df)
    adata = process_adata(adata, is_show=True)