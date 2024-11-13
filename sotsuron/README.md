# 卒業研究用のプログラム置き場

-> ファイルの読み書き部分などは正式なパッケージへ移行する予定

現在は使用しない予定＆誤ってこのパッケージを参照しないように名称を変更している

-----------

# プログラム構造のメモ

## subset.py
### 目的
MERSCOPEデータからsubsetを制作する

### 留意事項
変数の命名について：
ただの`x/y`や`micron_x/y`はマイクロメートル単位系、`pixel_x/y`はピクセル単位系。
`global`は全体での座標、`local`はfov内の相対的な座標。
`x/y_init`はピクセル単位系であり、特に**切り抜くとき**の始点となる座標。(int)
切り抜きに使用する座標はintだが、他はfloat型

パスの渡し方について：
モジュールとして利用できるように、パスは極力関数内で作成し、引数として渡さずに済むようにしたい

例外処理について：
できるだけraiseさせ、呼び出しもとでtry-except


### 入力
- input_dir = "/work/datasets/"
- output_dir = "/work/outputs/"
- target_str = "" # "z0"だけを対象とするなどしてテストが可能
- target_fov = 801

### 出力
- output_dir/subsetXX内にファイルを出力する
- 画像はさらにoutput_dir/subsetXX/images内に出力する
- 切り抜くときに使用した情報はjsonファイルに書き込む

### 関数一覧
`get_min_coordinates`
- 引数：なし
- 戻り値：Tuple[float, float]
- 用途：detected_transcripts.csvを開き、global_x/yの最小値を取得する

`read_manifest_json`
- 引数：なし
- 戻り値：Tuple[float, float, float, int, int]:
- 用途：images/manifest.jsonを開き、"microns_per_pixel"、"bbox_microns"のxy、"mosaic_width/height_pixels"の5つの値を取得する
それぞれ変換比率、x/yのオフセット、画像の最大サイズに相当する

`microns_to_pixel`
- 引数：x_micron:float, y_micron:float, microns_per_pixel:float, bbox_x:float, bbox_y:float
- 戻り値：Tuple[float, float]
- 用途：マイクロメートルからピクセル単位系へ変換し、x_pixel、y_pixelを返す

`pixel_to_microns`
- 引数：x_pixel:float, y_pixel:float, microns_per_pixel:float, bbox_x:float, bbox_y:float
- 戻り値：Tuple[float, float]:
- 用途：ピクセルからマイクロメートル単位系へ変換し、x_micron、y_micronを返す

`global_to_local_pixel`
- 引数：global_x:float, global_y:float, x_offset, y_offset
- 戻り値：Tuple[float, float]:
- 用途：globalの座標を、fov内での相対的な座標に変換し、local_x/yを返す。
  
`write_coordinates_json`
- 引数：subset_dir: str, data: dict
- 戻り値：なし
- 用途：subsetXX/crop_info.jsonにdataを保存する。
`data = {"_comment": "適当なコメント", "key1": value1, "key2": value2 }`などとする

`calculate_crop_coordinates`
- 引数：subset_dir:str, x_min_micron:float, y_min_micron:float, width:int=2048, height:int=2048
- 戻り値：Tuple[int, int]
- 用途：画像を切り抜く際の座標を計算し、write_coordinates_jsonを呼び出し書き込む。
x_init, y_initを返す。
    

`add_pixel_columns`
- 引数：df:DataFrame, x_init:int, y_init:int
- 戻り値：DataFrame
- 用途：引数で渡したdfに、global_pixel_x/y、local_pixel_x/yの列を追加する。
  

`crop_csv`
- 引数：subset_dir:str, x_min_micron:float, y_min_micron:float,  x_init:int, y_init:int, width:int = 2048, height:int = 2048
- 戻り値：なし
- 用途：detected_transcripts.csvを開く。
x_init:x_init+w/ y_init:y_init+h をマイクロメートル単位に変換し、その範囲で抽出する。
抽出したdfをdetected_transcripts_subsetXX.csvとして保存する。


`get_imagepath`
- 引数：なし
- 戻り値：List[str]
- 用途：input_dir/images内.tifファイルのパスのリストを取得する。

`crop_images`
- 引数：subset_dir:str, path_list:List[str], x_init:int, y_init:int, is_show:bool, w:int=2048, h:int=2048
- 戻り値：なし
- 用途：subsetXX/images内に切り抜いた後の画像を_subsetXX.pngとして保存する。
is_show=Trueの時、サムネイルとして切り抜いた後の画像や全体像を表示する。


`create_subset`
- 引数：fov:int, is_show:bool = False, width:int = 2048, height:int = 2048
- 戻り値：なし
- 用途：上記の関数を利用し、fovごとにデータをサブセット化する

--------------------------------------------------------------------------------


## gene_plot.py
### 目的
subsetを読み込み、画像と遺伝子発現のプロットを重ねて表示する

### 関数一覧
`read_select_transcripts_csv`
- 引数：fov:int, z:int
- 戻り値：Tuple[DataFrame, List[str]]
- 用途：detected_transcripts_subsetXX.csvから指定されたz座標のデータを抽出する。
抽出したdfと含まれる遺伝子名のリストを返す。

`load_image`
- 引数：img_dir:str, fov:int, z:int, img_type:str ="DAPI"
- 戻り値：Image
- 用途：img_dir/mosaic_{img_type}_z{z}_subset{fov}.pngの画像を読み込む。


`draw_gene_plot`
- 引数：ax:Axes, df:pd.DataFrame, gene_list: List[str] = None
- 戻り値：なし
- 用途：ax上にdfからgene_listが指定されているならばリストに含まれている行を抽出する。（そうでないならそのまま）
その後local_pixel_x/y列の情報からプロットを行う。

`show_gene_plot`
- 引数：img_dir:str, fov:int, z:int, df:pd.DataFrame, gene_list:List[str]
- 戻り値：なし
- 用途：load_imageで読み込んだ画像上に、draw_gene_plotでのプロット結果を重ねて表示する


## cell_polygon.py
### 目的
hdf5やmetadataを読み込みポリゴンを表示する
subsetから画像を読み込む。プロットもする。


### 関数一覧
`read_crop_info_json`
- 引数：なし
- 戻り値：Tuple[int, int]
- 用途：subsetXX/crop_info.jsonを読み込み、画像を切り抜いた時の始点を獲得する

`load_cell_data_from_hdf5`
- 引数：cell_id:str, hdf5_file:h5py.File
- 戻り値：Dict[int, DataFrame]
- 用途：hdf5ファイルの内容を辞書型に直していく関数

`convert_to_local_pixel`
- 引数：df:DataFrame, original:str, result:str
- 戻り値：DataFrame
- 用途：original_x/yという名称のglobalでマイクロメートル単位の値を、result_x/yという名称のlocalピクセル単位に変換する。dfに列を追加する。

`load_meta_csv`
- 引数：fov:int
- 戻り値：DataFrame
- 用途：cell_metadata.csvを読み込み、convert_to_local_pixelでlocalなピクセル値に変換する
  
`load_and_convert_data`
- 引数：fov:int
- 戻り値：Tuple[Dict[str, Dict[int, DataFrame]], DataFrame]
- 用途：上記の関数を利用して、hdf5ファイルの内容を辞書型に、metadaa.csvの内容をdfに直して返す。
座標計算済み。

`draw_all_polygon`
- 引数：ax:plt.Axes, cell_data_dict:Dict[str, pd.DataFrame], metadata_df:DataFrame, z:int
- 戻り値：なし
- 用途：fov全体のポリゴンを表示する。コードを書き直すことで、中央値やboxを表示可能。


`draw_selected_polygon`
- 引数：cell_id: str, ax: plt.Axes, cell_data_dict: Dict[str, pd.DataFrame],metadata_df: pd.DataFrame, z: int, w: int = 512, h: int = 512
- 戻り値：Tuple[float, float]
- 用途：指定したcell_idのポリゴンとその周辺領域を表示する。コードを書き直すことで、中央値やboxを表示可能。表示領域の座標を変換できるように始点の位置を返す。


`show_polygon`
- 引数：img_dir:str, fov:int, z:int, cell_data_dict: Dict[str, pd.DataFrame],metadata_df: pd.DataFrame, cell_id:Optional[str] = None, w:int=512, h:int=512
- 戻り値：なし
- 用途：画像上にポリゴンを表示する。
cell_idがない時はfov全体を、ある場合はその周辺のみを切り抜いて表示する。


`show_polygon_and_plot`
- 引数：img_dir:str, fov:int, z:int, cell_data_dict:Dict, metadata_df:pd.DataFrame, cell_id:Optional[str] = None, w:int=512, h:int=512
- 戻り値：なし
- 用途：画像上にポリゴンと遺伝子発現プロットを表示する。
余裕があったら、cell_idの周辺のみを表示できるように拡張したいと思っている。


`show_area_of_focus`
- 引数：cell_id, cell_data_dict, selected_data, target_z, img_dir, target_fov
- 戻り値：なし
- 用途：全体像のうち切り抜いた位置を示す画像が表示できると思っている。（切り抜かない予定）


## test.py
fov全体の遺伝子発現プロット、ポリゴン、両方表示
cell_id周辺の遺伝子発現プロット、ポリゴン、両方表示　を可能にする

### 関数一覧
`read_crop_info_json`
- sub_dir:str
- Tuple[int, int]
- 用途：subsetXX/crop_info.jsonを読み込み、画像を切り抜いた時の始点を獲得する

`load_cell_data_from_hdf5`
- cell_id:str, hdf5_file:h5py.File
- Dict[int, DataFrame]
- 用途：hdf5ファイルの内容を辞書型に直していく関数

`convert_to_local_pixel`
- 引数：sub_dir:str, df:DataFrame, original:str, result:str
- 戻り値：DataFrame
- 用途：original_x/yという名称のglobalでマイクロメートル単位の値を、result_x/yという名称のlocalピクセル単位に変換する。dfに列を追加する。

`load_meta_csv`
- 引数：input_dir:str, sub_dir:str, fov:int
- 戻り値：DataFrame
- 用途：cell_metadata.csvを読み込み、convert_to_local_pixelでlocalなピクセル値に変換する
  
`load_and_convert_data`
- 引数：input_dir:str, fov:int
- 戻り値：Tuple[Dict[str, Dict[int, DataFrame]], DataFrame]
- 用途：上記の関数を利用して、hdf5ファイルの内容を辞書型に、metadaa.csvの内容をdfに直して返す。
座標計算済み。

`read_select_transcripts_csv`
- 引数：sub_dir:str, fov:int, z:int
- 戻り値：Tuple[DataFrame, List[str]]
- 用途：detected_transcripts_subsetXX.csvから指定されたz座標のデータを抽出する。
抽出したdfと含まれる遺伝子名のリストを返す。

`load_image`
- 引数：img_dir:str, fov:int, z:int, img_type:str ="DAPI"
- 戻り値：Image
- 用途：img_dir/mosaic_{img_type}_z{z}_subset{fov}.pngの画像を読み込む。


`draw_gene_plot`
- 引数：ax:Axes, df:pd.DataFrame, gene_list: List[str] = None
- 戻り値：なし
- 用途：gene_listが指定されているならばリストに含まれている行を抽出する。（そうでないならそのまま）
その後local_pixel_x/y列の情報からプロットを行う。

`draw_polygon`
- ax:plt.Axes, x_data:float, y_data:float
- None
- ポリゴンを1つ表示する

`draw_all_polygon`
- ax:plt.Axes, cell_data_dict:Dict[str, pd.DataFrame], metadata_df:DataFrame, z:int
- None
- ループを回してポリゴンをすべて表示する

`draw_id_polygon`
- cell_id: str, ax: plt.Axes, cell_data_dict: Dict[str, pd.DataFrame], metadata_df: pd.DataFrame, z: int, w: int = 512, h: int = 512
- Tuple[float, float]:
- cell_id周辺に位置合わせを行った上でポリゴンを表示する
 
`show_whole_fov`
- ax:plt.Axes, fov: int, z: int, cell_data_dict: Dict, metadata_df: pd.DataFrame, gene_list: Optional[List[str]] = None, show_plot: bool = True, show_polygon: bool = True
- None
- fov全体を表示する
show_plotがTrue⇒遺伝子のプロットを表示する。
さらに、gene_listが与えられている⇒その遺伝子のみを抽出し、そうでない⇒全遺伝子を表示
show_polygonがTrue⇒ポリゴンを表示する。

`calc_df`
- df:DataFrame, metadata_df, cell_id, w:int = 512, h:int=512
- Tuple[DataFrame, float, float]
- cell_id周辺を表示する際の遺伝子発現データの位置合わせ関数。
位置合わせ後の値でdfのlocal_pixel_x/yを更新して返す。
このときのoffsetを返す。


`show_id_surroundings`
-  ax:plt.Axes, fov: int, z: int, cell_data_dict: Dict, metadata_df: pd.DataFrame, cell_id: str, gene_list: Optional[List[str]] = None, show_plot: bool = True, show_polygon: bool = True, w: int = 512, h: int = 512
- Tuple[float, float]
- id周辺領域を表示する。
show_plotがTrue⇒遺伝子のプロットを表示する。
さらに、gene_listが与えられている⇒その遺伝子のみを抽出し、そうでない⇒全遺伝子を表示
show_polygonがTrue⇒ポリゴンを表示する。
切り抜くときの始点をoffsetとして返す。


`display_area`
- img_dir: str, fov: int, z: int, cell_data_dict: Dict[str, pd.DataFrame], metadata_df: pd.DataFrame, cell_id: Optional[str] = None, w: int = 512, h: int = 512, gene_num:Optional[int] = None , show_plot: bool = True, show_polygon: bool = True
- None:
- cell_idが指定されているならばその周辺のみを表示する
そうでないならばfov全体を表示する
plt.showで最終的な画像に重ね合わせて全てを表示する。