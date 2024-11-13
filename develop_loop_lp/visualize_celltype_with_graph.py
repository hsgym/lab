from package.viewer_common import CommonPlotPolygonViewer
from typing import Optional
import networkx as nx
import matplotlib.pyplot as plt

# グラフ描画
def visualize_predicted_type(G: nx.Graph, num_x_bins:int, num_y_bins:int, category_to_color:dict, show_image:bool, png_path: Optional[str]=None, plotter: Optional[CommonPlotPolygonViewer]=None):
    # 描画に必要な情報の設定
    radius = int(plotter.width / (2 * num_x_bins)) - 2
    w2 = (plotter.width / num_x_bins) / 2
    h2 = (plotter.height / num_y_bins) / 2 
    
    fig, ax = plt.subplots()
    
    cell_type_set = set()
    for node, data in G.nodes(data=True):
        cell_type = data['cell_type']
        # if cell_type == "Background":
        #    continue

        color = category_to_color[cell_type]
        cell_type_set.add(cell_type)
        pos_prev = None
        
        for (x, y) in data["pos"]:
            x_center = x * (plotter.width / num_x_bins) + w2 + 1
            y_center = y * (plotter.height / num_y_bins) + h2 + 1
            # ax.scatter(x_center, y_center, s=radius, c=[color])
            if cell_type == "Background":
                circle = plt.Circle((x_center, y_center), radius, facecolor=None, edgecolor = color, linewidth=1)
            else:
                circle = plt.Circle((x_center, y_center), radius, facecolor=color)
            ax.add_patch(circle)
            if pos_prev is not None:
                ax.plot([x_center, pos_prev[0]], [y_center, pos_prev[1]], color='gray')
            pos_prev = (x_center, y_center)
    
    if show_image:
        if png_path is None:
            png_path = plotter.get_imagepath(plotter.img_dir, image_keywords=["DAPI","z3"])[0]
        img = plotter.load_image(png_path)
        ax.imshow(img, cmap="gray")

    # 凡例描画
    sorted_cell_types = sorted(cell_type_set)
    handles = [plt.Rectangle((0, 0), 1, 1, color=category_to_color[cat]) for cat in sorted_cell_types]
    labels = list(sorted_cell_types)
    ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1), title='Cell Type')

    ax.set_xlim(0, plotter.width)
    ax.set_ylim(plotter.height, 0)
    ax.set_title(f"cell_type{' + image' if show_image else ''} : {plotter.title}")
    fig.tight_layout()
    fig.subplots_adjust(right=0.75)
    plt.show()