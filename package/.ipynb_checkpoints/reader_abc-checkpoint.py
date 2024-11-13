"""
抽象クラス
＋共通して利用する関数
"""

import os
from abc import ABC, abstractmethod

class AbstractDataReader(ABC):
    def __init__(self, input_dir, output_dir, fov, width, height):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.fov = fov
        self.width = width
        self.height = height    
        
    @abstractmethod
    def read_gene_data(self):
        pass

    @abstractmethod
    def read_cell_data(self):
        pass

    def get_path(self):
        subset_dir = os.path.join(self.output_dir, f"subset{self.fov}")
        img_dir = os.path.join(subset_dir, "images")
        return subset_dir, img_dir