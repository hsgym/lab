"""
共通して利用する関数のクラス
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns

from pandas import DataFrame
from PIL import Image
from typing import List, Optional

class CommonPlotPolygonViewer:
    def __init__(self, img_dir: str, fov: int, z: int, title: str, width: int = 2048, height: int = 2048, gene_list: List[str] = []):
        self.img_dir = img_dir
        self.fov = fov
        self.z = z
        self.title = title
        self.width = width
        self.height = height
        self.gene_list = gene_list

    def load_image(self, png_path:str) -> Image:
        # png_path = os.path.join(self.img_dir, f"mosaic_{base}_z{self.z}_subset{self.fov}.png")
        if not os.path.isfile(png_path):
            raise FileNotFoundError(f"File not found: {png_path}")
        return Image.open(png_path)

    def draw_gene_plot(self, ax: plt.Axes, df: DataFrame, show_legend:bool, size: int = 2) -> None:
        selected_df = df[df["gene"].isin(self.gene_list)]
        sns.scatterplot(data=selected_df, x="local_pixel_x", y="local_pixel_y", hue="gene", palette="hls", s=size, ax=ax, alpha=0.8, legend = show_legend)
        return

    """
    def draw_polygon(self, ax: plt.Axes, x_data: float, y_data: float) -> None:
        ax.plot(x_data, y_data, linestyle="-", alpha=0.7, c="#d62728")
    """
    def draw_poly(self, ax: plt.Axes, x_data: float, y_data: float) -> None:
        ax.plot(x_data, y_data, linestyle="-", alpha=0.7, c="#d62728")

    def draw_multipolygon(self, ax: plt.Axes, fov_cell_df: DataFrame) -> None:
        for poly in fov_cell_df["Geometry_local"]:
            if poly.geom_type == 'Polygon':
                x, y = poly.exterior.xy
                self.draw_poly(ax, x, y)
            elif poly.geom_type == 'MultiPolygon':
                x, y = poly.geoms[0].exterior.xy
                self.draw_poly(ax, x, y)


    def show_whole_fov(self, ax: plt.Axes, fov_cell_df: DataFrame, selected_df: DataFrame, show_plot: bool = True, show_polygon: bool = True) -> None:
        if show_plot:
            if self.z is None:
                selected_z_df = selected_df.copy()
            else:
                selected_z_df = selected_df.query(f"global_z=={self.z}")
            if selected_z_df.empty:
                raise ValueError("No data matching the z value")

            size = 4
            show_legend = (len(self.gene_list) <= 20)
            self.draw_gene_plot(ax, selected_z_df, show_legend, size=size)

        if show_polygon:
            self.draw_multipolygon(ax, fov_cell_df)

    def get_imagepath(self, input_dir: str, image_keywords: List[str] = []) -> List[str]:
        # image_input_dir = os.path.join(input_dir, "images")

        if not os.path.exists(input_dir):
            raise Exception(f"{input_dir} not found")
        
        # 以下、未テストのコード
        all_png_files = [f for f in os.listdir(input_dir) if f.endswith(".png")]

        if image_keywords:
            filtered_files = [f for f in all_png_files if all(word in f for word in image_keywords)]
        else:
            filtered_files = all_png_files

        path_list = [os.path.join(input_dir, f) for f in filtered_files]

        if not path_list:
            raise FileNotFoundError(f"No .png files {'containing' if image_keywords else 'found'} {' and '.join(image_keywords) if image_keywords else ''} in {input_dir}")

        return path_list
        
        """
        # 以下はテスト済みのコード（動かなかったらこれに戻す）
        path_list = []
        if not image_keywords:
            tif_files = [f for f in os.listdir(input_dir) if f.endswith(".png")]
        else:
            tif_files = [f for f in os.listdir(input_dir) if f.endswith(".png") and all(word in f for word in image_keywords)]

        for tif_file in tif_files:
            image_path = os.path.join(input_dir, tif_file)
            path_list.append(image_path)

        if not path_list:
            raise Exception(f"no .tif files {'containing' if image_keywords else 'found'} {' and '.join(image_keywords)} in {input_dir}")

        return path_list        
        """

        

    def display_area(self, fov_cell_df: DataFrame, selected_df: DataFrame ,show_plot: bool = True, show_polygon: bool = True, image_keywords:List[str]=[], is_save: bool = False, ax: Optional[plt.Axes] = None) -> None:
        if ax is None:
            fig, ax = plt.subplots(figsize=(16, 16))
        
        png_path = self.get_imagepath(self.img_dir, image_keywords)
        img = self.load_image(png_path[0])
        self.show_whole_fov(ax, fov_cell_df, selected_df, show_plot=show_plot, show_polygon=show_polygon)
        
        ax.set_title(f"fov={self.fov}, z={self.z}, {self.title}")
        ax.set_xlim([0, self.width])
        ax.set_ylim([self.height, 0])
        ax.imshow(img, cmap="gray")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        if is_save:
            tmp_dir = "/work/tmp_output/"
            os.makedirs(tmp_dir, exist_ok=True)
            filename = os.path.join(tmp_dir, f"fov{self.fov}_z{self.z}_display_area.png")
            plt.savefig(filename, format="png")
            print(f"saved {filename}")

        if ax is None:
            fig.tight_layout()
            plt.show()
    