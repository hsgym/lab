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

def show_ndarray(matrix:np.ndarray, title:str) -> None:
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.imshow(matrix, cmap="binary_r")
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_title(title)
    plt.show()

def load_image(img_dir:str, fov:int, z:int, base:str="DAPI") -> Image:
    png_path = os.path.join(img_dir, f"mosaic_{base}_z{z}_subset{fov}.png")
    if not os.path.isfile(png_path):
        raise FileNotFoundError(f"File not found: {png_path}")
    return Image.open(png_path)

def visualize_polygons(polygons, title, height,width,  ax: Optional[plt.Axes] = None, color: Optional[str] = "yellow"):
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 4))

    for poly1 in polygons:
        # 4未満の座標がある場合は無視
        if len(poly1) < 4:
            continue
        ax.plot(*Polygon(poly1).exterior.xy, color=color)

    ax.set_xlim([0, width])
    ax.set_ylim([height, 0])
    ax.set_title(title)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.tight_layout()  