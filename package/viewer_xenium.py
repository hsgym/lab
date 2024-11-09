"""
Xenium用の子クラス
"""

import os
from typing import List

from package.reader_xenium import XeniumDataReader
from package.viewer_abc import AbstractPlotPolygonViewer
from package.viewer_common import CommonPlotPolygonViewer

class XeniumViewer(AbstractPlotPolygonViewer):
    def __init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int, resampling_factor:float = 0.2125, z = None):
        super().__init__(input_dir, output_dir, fov, width, height)
        self.resampling_factor = resampling_factor
        self.z = z
       
    def view_area(self, title:str="", show_plot:bool=True, show_polygon:bool=True, image_keywords:List[str]=None) -> None:
        reader = XeniumDataReader(self.input_dir, self.output_dir, self.fov, self.width, self.height, self.resampling_factor, self.z)
        gene_df, gene_list = reader.read_transcripts_data()
        
        output_parquet = os.path.join(self.output_dir, 'cell_boundaries_pixel.parquet')
        cell_df = reader.read_cell_data(output_parquet)
        
        subset_dir, img_dir = reader.get_path()
        total_gene_num = len(gene_list)
        plotter = CommonPlotPolygonViewer(img_dir, self.fov, self.z, title, self.width, self.height, gene_list=gene_list)
        plotter.display_area(cell_df, gene_df, show_plot, show_polygon, image_keywords)
        # plotter.display_area(fov_cell_df, selected_df, show_plot=False)
        # plotter.display_area(fov_cell_df, selected_df, show_polygon=False)

"""
if __name__ == "__main__":
    input_dir = "/work/datasets/Okada/output-XETG00130__0014491__Cont_3_1__20240403__094846"
    output_dir = "/work/output_Xenium"
    resampling_factor = 0.2125

    # xenium data の分割
    fov = 31
    width = 4096
    height = 4096
    
    viewer = XeniumViewer(input_dir, output_dir, fov, width, height, resampling_factor)
    title = "default"
    viewer.view_area(title, True, True)
"""