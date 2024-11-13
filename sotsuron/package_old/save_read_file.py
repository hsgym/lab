import pickle
import numpy as np
from typing import List



# listを.pklで保存する
def save_list_as_pkl(save_path: str, data_list: List[np.ndarray]) -> None:
    try:
        with open(save_path, mode='wb') as f:
            pickle.dump(data_list, f)
        print(f"Saved to {save_path}")
    except Exception as e:
        print(f"Error occurred while `save_list_as_pkl`: {e}")


def read_pkl(pkl_path) -> List[np.ndarray]:
    try:
        with open(pkl_path, 'rb') as file:
            data = pickle.load(file)
        return data
    except Exception as e:
        print(f"Error occurred while `read_pkl()`: {e}")
        return []



# numpy行列を.npyで保存する
def save_matrix(npy_path: str, matrix: np.ndarray) -> None:
    try:
        np.save(npy_path, matrix)
        print(f"Saved to {npy_path}")
    except Exception as e:
        print(f"Error occurred while `save_matrix()`: {e}")

def read_npy(npy_path:str) -> np.ndarray:
    try:
        data = np.load(npy_path)
        return data
    except Exception as e:
        print(f"Error occurred `read_npy()`: {e}")
        return None


