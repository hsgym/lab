"""
抽象クラス
"""

from abc import ABC, abstractmethod

class AbstractPlotPolygonViewer(ABC):
    def __init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.fov = fov
        self.width = width
        self.height = height

    @abstractmethod
    def view_area(self):
        pass