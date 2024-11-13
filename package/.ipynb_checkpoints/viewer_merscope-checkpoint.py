"""
MERSCOPE用の子クラス
"""


import os
from typing import List

from package.reader_merscope import MerscopeDataReader
from package.viewer_abc import AbstractPlotPolygonViewer
from package.viewer_common import CommonPlotPolygonViewer


class MerscopeViewer(AbstractPlotPolygonViewer):
    def __init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int, z = None):
        super().__init__(input_dir, output_dir, fov, width, height)
        self.z = z
            
    def view_area(self, title:str = "", show_plot:bool=True, show_polygon:bool=True, image_keywords:List[str]=None) -> None:
        reader = MerscopeDataReader(self.input_dir, self.output_dir, self.fov, self.width, self.height, self.z)
        gene_df, gene_list = reader.read_gene_data()
        mosaic_path = os.path.join(self.input_dir, "Cellpose_DAPI_CB3/cellpose_mosaic_space.parquet")
        cell_df = reader.read_cell_data(mosaic_path)
        
        subset_dir, img_dir = reader.get_path()

        total_gene_num = len(gene_list)
        plotter = CommonPlotPolygonViewer(img_dir, self.fov, self.z, title, self.width, self.height, gene_list=gene_list)
        plotter.display_area(cell_df, gene_df, show_plot, show_polygon, image_keywords)


"""
if __name__ == "__main__":
    input_dir = "/work/datasets/"
    output_dir = "/work/output_MERSCOPE"
    
    fov = 236
    width=2048
    height=2048
    z=3
    viewer = MerscopeViewer(input_dir, output_dir, fov, width, height, z)
    title = "default"
    viewer.view_area(title, True, True, ["DAPI","z3"])
"""