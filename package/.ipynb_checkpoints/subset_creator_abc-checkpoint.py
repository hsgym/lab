"""
抽象クラス：
"""

from abc import ABC, abstractmethod

class AbstractSubsetCreator(ABC):
    def __init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.fov = fov
        self.width = width
        self.height = height    

    @abstractmethod
    def make_directory(self):
        pass
    
    @abstractmethod
    def decide_crop_area(self):
        pass
    
    @abstractmethod
    def crop_table_data(self):
        pass
    
    @abstractmethod
    def crop_image(self):
        pass
    
    @abstractmethod
    def create_subset(self):
        pass


