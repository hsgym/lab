"""
MERSCOPE用の子クラス
"""

import os
import geopandas as gpd
import pandas as pd
from pandas import DataFrame
from shapely.geometry import box
from shapely.geometry import Point
from shapely.ops import transform
from typing import List, Tuple, Optional
import pandas as pd
from typing import Optional, Dict
import json
import os
import h5py
from package.reader_abc import AbstractDataReader
from package.subset_creator_common import CommonSubsetCreator

class MerscopeDataReader(AbstractDataReader):
    def __init__(self, input_dir, output_dir, fov, width, height, z = None):
        super().__init__(input_dir, output_dir, fov, width, height)
        self.z = z
        self.common = CommonSubsetCreator()

    def read_gene_data(self) -> Tuple[DataFrame, List[str]]:
        subset_dir, img_dir = self.get_path()
        pkl_path = os.path.join(subset_dir, f"detected_transcripts_subset{self.fov}.pkl")
        
        df = pd.read_pickle(pkl_path)
    
        if self.z is None: # zが無い時
            selected_df = df.copy()
        else:
            selected_df = df.query(f"global_z == {self.z}")

        if selected_df.empty:
            raise Exception(f"DataFrame is empty")
        
        if 'feature_name' in selected_df.columns:
            gene_df = selected_df.rename(columns={'feature_name': 'gene'})  
        else:
            gene_df = selected_df.copy()

        gene_name_list = gene_df['gene'].unique().tolist()
        return gene_df, gene_name_list
            
    def transfer_position(self, x, y, x_init, y_init):
        return x-x_init, y-y_init

    def transform_multipolygon(self, multipolygon, x_init, y_init):
        return multipolygon.apply(lambda geom: transform(lambda x, y, z=None: self.transfer_position(x, y, x_init, y_init), geom))
        
    def read_cell_data(self, mosaic_path:Optional[str]=None) -> pd.DataFrame:
        subset_dir, img_dir = self.get_path()
        
        if mosaic_path is None:
            mosaic_path = os.path.join(self.input_dir, "Cellpose_DAPI_CB3/cellpose_mosaic_space.parquet")

        cell_df = gpd.read_parquet(mosaic_path)  

        if self.z is None:
            z_df = cell_df.copy()
        else:
            z_df = cell_df.query(f"ZIndex == {self.z}")

        crop_json_path = os.path.join(self.output_dir, f"crop_info.json")
        x_init, y_init, x_last, y_last = self.common.read_crop_info_json(crop_json_path, self.fov)
        search_area = box(x_init, y_init, x_last, y_last)
        
        # 範囲内の行のみを抽出
        fov_cell_df = z_df[z_df['Geometry'].apply(lambda geom: geom.intersects(search_area))].copy()

        print("cell_count: ", fov_cell_df.shape)
        
        # fov単位での座標に変換（ピクセル単位系）
        fov_cell_df["Geometry_local"] = self.transform_multipolygon(fov_cell_df["Geometry"], x_init, y_init)
        
        selected_fov_cell_df = fov_cell_df.filter(items=['ID', 'ZIndex', "Geometry", "Geometry_local"])
        return selected_fov_cell_df # , selected_cell_data

"""
# 使用例
if __name__ == "__main__":
    input_dir = "/work/datasets/"
    output_dir = "/work/output_MERSCOPE"
    
    fov = 236
    width=2048
    height=2048
    z=3
    reader = MerscopeDataReader(input_dir, output_dir, fov, width, height, z)
    gene_df, gene_name_list = reader.read_gene_data()
    cell_df = reader.read_cell_data()
"""
