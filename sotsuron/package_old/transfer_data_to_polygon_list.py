import os
import numpy as np
import cv2
import h5py
import geojson
from typing import List
from package_old.save_read_file import save_list_as_pkl
from package_old.plot_polygon_viewer import convert_to_local_pixel, load_cell_data_from_hdf5
from package_old.execute_graphcut import visualize_polygons



def read_geojson(file_path:str) -> List[np.ndarray]:
    with open(file_path, 'r') as f:
        geojson_data = geojson.load(f)

    polygon_list = []

    for geometry in geojson_data['geometries']:
        if geometry['type'] == 'Polygon':
            coordinates = geometry['coordinates'][0]  # 最初のリストが座標のリストです
            polygon = np.array(coordinates, dtype=int)
            polygon_list.append(polygon)

    return polygon_list


def read_cellpose_npy_outlines(cellpose_path):
    cellpose_data = np.load(cellpose_path, allow_pickle=True).item()
    outlines = cellpose_data.get("outlines")
    return outlines
    

def extract_polygons(outlines) -> List[np.ndarray]:
    polygons_data = []

    for object_id in np.unique(outlines):
        if object_id == 0:
            continue
            
        mask = outlines == object_id
        contours, _ = cv2.findContours(mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        # contoursには複数の輪郭が含まれている可能性があるため、各輪郭ごとに処理
        for contour in contours:
            coordinates = contour.squeeze()
            polygons_data.append(coordinates)

    return polygons_data


def load_and_convert_polygon_data(hdf5_path:str) -> List[np.ndarray]:
    cell_data_dict = {}
    polygon_list = []

    with h5py.File(hdf5_path, 'r') as f:
        cell_ids = f["featuredata"].keys()
        for cell_id in cell_ids:
            cell_data_dict[cell_id] = load_cell_data_from_hdf5(cell_id, f)        
            cell_data_dict[cell_id] = convert_to_local_pixel(input_dir, output_dir, fov, cell_data_dict[cell_id], "global", "local_pixel")
            # Extract 'local_pixel_x' and 'local_pixel_y' columns from DataFrame
            local_pixel_x = cell_data_dict[cell_id]['local_pixel_x'].values.astype(int)
            local_pixel_y = cell_data_dict[cell_id]['local_pixel_y'].values.astype(int)

            # Combine 'local_pixel_x' and 'local_pixel_y' columns into a NumPy array of tuples
            local_polygon_data = np.column_stack((local_pixel_x, local_pixel_y))

            # Append the NumPy array to the polygon_list
            polygon_list.append(local_polygon_data)
    return polygon_list


def transfer_to_polygon_list(path:str, method:str) -> List[np.ndarray]:
    try:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"file not found: {path}")
        if method == "baysor":    
            polygon_list = read_geojson(path)
        elif method == "answer" or "cellpose":
            outlines = read_cellpose_npy_outlines(path)
            polygon_list = extract_polygons(outlines)
        elif method == "merscope":
            polygon_list = load_and_convert_polygon_data(path)
        else:
            raise ValueError(f"invalid method: {method}")

        cell_count = len(polygon_list)
        print(f"{method} cell_count : ", cell_count) # 記録を取ること
        return polygon_list

    except Exception as e:
        print(f"exception ocuured: {e}")
        raise


if __name__ == "__main__":
    input_dir = "/work/datasets"
    output_dir = "/work/output/a"
    ans_dir = "/work/answers"
    fov = 800
    base = "DAPI" # or "Cellbound3"
    width, height = 2048, 2048
    fov = 800
    z = 3 # 仮
    experiment_dir = f"/work/experiment/fov{fov}"

    method = "baysor" # "cellpose", "answer", "merscope"
    
    if method == "baysor":
        polygon_data_path = f"/work/baysor/segmentation_polygons_fov{fov}.json"
        save_list_path = os.path.join(experiment_dir, f"polygon_list_{method}_{base}_fov{fov}.pkl")
    elif method == "answer":
        papolygon_data_pathth = os.path.join(ans_dir, f"{base}_z3_subset{fov}_seg.npy")
        save_list_path = os.path.join(experiment_dir, f"polygon_list_{method}_{base}_fov{fov}.pkl")
    elif method == "cellpose":
        polygon_data_path = os.path.join(output_dir, f"subset{fov}/cellpose_mask_{base}_cyto/mosaic_{base}_z3_subset{fov}_seg.npy")
        save_list_path = os.path.join(experiment_dir, f"polygon_list_{method}_{base}_fov{fov}.pkl")
    else:
        polygon_data_path = os.path.join(input_dir, f"cell_boundaries/feature_data_{fov}.hdf5")
        save_list_path = os.path.join(experiment_dir, f"polygon_list_merscope_{base}_fov{fov}.pkl")
    
    polygon_list = transfer_to_polygon_list(polygon_data_path, method)
    visualize_polygons(polygon_list, f"fov{fov} - {method}", width, height)
    save_list_as_pkl(save_list_path, polygon_list)

    

