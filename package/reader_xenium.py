"""
Xenium用の子クラス
"""

import os
import dask.dataframe as dd
import geopandas as gpd
import pandas as pd
from pandas import DataFrame
from shapely.geometry import Polygon
from typing import Tuple, List, Optional

from package.reader_abc import AbstractDataReader
from package.subset_creator_common import (
    read_crop_info_json, global_to_local, microns_to_pixel, pixel_to_microns
)


class XeniumDataReader(AbstractDataReader):
    def __init__(self, input_dir, output_dir, fov, width, height, resampling_factor:float =0.2125, z:int = None):
        super().__init__(input_dir, output_dir, fov, width, height)
        self.resampling_factor = resampling_factor
        self.z = z

    def read_gene_data(self) -> Tuple[DataFrame, List[str]]:
        subset_dir, img_dir = self.get_path()
        pkl_path = os.path.join(subset_dir, f"detected_transcripts_subset{self.fov}.pkl")
        
        df = pd.read_pickle(pkl_path)
    
        if self.z is None: # zが無い時
            selected_df = df.copy()
        else:
            selected_df = df.query(f"z_location == {self.z}")

        if selected_df.empty:
            raise Exception(f"DataFrame is empty")
        
        if 'feature_name' in selected_df.columns:
            gene_df = selected_df.rename(columns={'feature_name': 'gene'})  
        else:
            gene_df = selected_df.copy()

        gene_name_list = gene_df['gene'].unique().tolist()
        return gene_df, gene_name_list

    def filter_and_create_polygons(self, filtered_df: pd.DataFrame) -> gpd.GeoDataFrame:
        polygons = []
        grouped = filtered_df.groupby('cell_id')

        for cell_id, group in grouped:
            vertices = list(zip(group['local_pixel_x'], group['local_pixel_y']))
            if len(vertices) > 2:  # ポリゴンを構成するのに最低3点必要
                polygon = Polygon(vertices)
                polygons.append({'cell_id': cell_id, 'Geometry_local': polygon})
        return gpd.GeoDataFrame(polygons, geometry='Geometry_local', crs="EPSG:4326")
    
    def read_cell_data(self, parquet_path:Optional[str]=None) -> pd.DataFrame:
        if parquet_path is None:
            output_parquet = os.path.join(self.output_dir, 'cell_boundaries_pixel.parquet')

        # fov_cell_dfを作る関数
        cell_ddf = dd.read_parquet(parquet_path)
        crop_json_path = os.path.join(self.output_dir, f"crop_info.json")
        x_init, y_init, x_last, y_last = read_crop_info_json(crop_json_path, self.fov)
        # pixel -> micron
        microns_per_pixel = self.resampling_factor
        bbox_micron_x = 0.0
        bbox_micron_y = 0.0

        print(self.resampling_factor)
        print(float(x_init), float(y_init), microns_per_pixel, bbox_micron_x, bbox_micron_y)
        x_min_micron, y_min_micron = pixel_to_microns(float(x_init), float(y_init), microns_per_pixel, bbox_micron_x, bbox_micron_y)
        x_width_micron, y_height_micron = pixel_to_microns(float(x_last), float(y_last), microns_per_pixel, bbox_micron_x, bbox_micron_y)

        filltered_df = cell_ddf.query(f'{x_min_micron} <= vertex_x & vertex_x <= {x_width_micron} & {y_min_micron} <= vertex_y & vertex_y <= {y_height_micron}').compute() 

        filltered_df[["vertex_x_pixel", "vertex_y_pixel"]] = filltered_df.apply(
            lambda row: microns_to_pixel(row["vertex_x"], row["vertex_y"], microns_per_pixel, bbox_micron_x, bbox_micron_y), 
            axis=1, result_type = "expand"
        )

        filltered_df[["local_pixel_x", "local_pixel_y"]] = filltered_df.apply(
            lambda row: global_to_local(row["vertex_x_pixel"], row["vertex_y_pixel"], x_init, y_init), 
            axis=1, result_type = "expand"
        )

        cell_df = self.filter_and_create_polygons(filltered_df)
        return cell_df # , selected_cell_data
    

""" 
# 使用例
if __name__ == "__main__":
    input_dir = "/work/datasets/Okada/output-XETG00130__0014491__Cont_3_1__20240403__094846"
    output_dir = "/work/output_Xenium"
    resampling_factor = 0.2125

    # xenium data の分割
    fov = 31
    width = 4096
    height = 4096

    reader = XeniumDataReader(input_dir, output_dir, fov, width, height, resampling_factor)
    gene_df, gene_name_list = reader.read_gene_data()
    cell_df = reader.read_cell_data()
"""