# 必要なライブラリのimport 
import os
import numpy as np
import pandas as pd

from package_old.plot_polygon_viewer import (
    read_select_transcripts_csv
)

from package_old.save_read_file import save_matrix

# ----------
def calc_gene_count(df: pd.DataFrame) -> np.ndarray:
    tmp = df.copy()
    tmp = tmp[["local_pixel_x", "local_pixel_y"]]
    # float型→四捨五入してint型へ（0-2048の範囲に無理矢理clipで入れる）
    tmp['x_pixel_int'] = tmp['local_pixel_x'].round().astype(int).clip(0, 2048 - 1)
    tmp['y_pixel_int'] = tmp['local_pixel_y'].round().astype(int).clip(0, 2048 - 1)
    grouped = tmp.groupby(["x_pixel_int", "y_pixel_int"], as_index=False).size()
    grouped = grouped.rename({"size":"density"}, axis=1)

    x_range = range(0, 2048)
    y_range = range(0, 2048)

    gene_count_matrix = np.zeros((len(x_range), len(y_range)))

    for index, row in grouped.iterrows():
        x_index = x_range.index(row['x_pixel_int'])
        y_index = y_range.index(row['y_pixel_int'])
        gene_count_matrix[x_index, y_index] = row['density']

    gene_count_matrix = gene_count_matrix.astype(int)
    return gene_count_matrix



if __name__ == "main":
    fov_list = [800, 801, 1101, 1201, 1301, 1401, 1501, 1701] 

    input_dir = "/work/datasets"
    output_dir = "/work/outputs/a"

    for fov in fov_list:
        subset_dir = os.path.join(output_dir, f"subset{fov}")
        img_dir = os.path.join(subset_dir, "images")

        # selected_df:表示する転写物のdf
        selected_df, _gene_list = read_select_transcripts_csv(subset_dir, fov)
        
        row, _col = selected_df.shape
        total_gene_count = row
        print(total_gene_count) # logに保存

        gene_count_matrix = calc_gene_count(selected_df)
        print(gene_count_matrix.max().max()) # 検証用
        print(gene_count_matrix.shape) # (2048, 2048)であることを確認

        experiment_dir = f"/work/experiment/fov{fov}"
        save_path = os.path.join(experiment_dir, f"gene_count_matrix_fov{fov}.npy")
        save_matrix(save_path, gene_count_matrix)
        print(f"save to {save_path}")