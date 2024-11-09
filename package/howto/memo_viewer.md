# viewer系のプログラム用の設定書
## 要望
基本的な空間トランスクリプトームデータを可視化したい
（画像、遺伝子発現点、細胞領域）

### 現状
技術によって測定の仕方が異なるため、統一的な方法では読み込めない

### 理想
subsetを指定することで、画像、遺伝子発現点、細胞領域を可視化できる
（これは３つのうち任意のものを重ねて可視化できる：画像only、画像+細胞領域とか）
後で他のデータを読み込みたくなる可能性があるため拡張の余地がある
論文や発表資料用に画像をきれいな形で保存できるようにしたい

## 要求
- 指定されたfovに対して必要なデータの読み込み
- 描画
- 要求があればきれいな形で保存する

### input
must
- 入力データのあるディレクトリのパス
- 出力先ディレクトリのパス
- 測定方法（Xenium or MERSCOPE）が分かるような何か
- 読み込みたい視野（fov）
- 描画したいもの（画像、発現点、細胞領域）
- 保存の有無

### output
- 指定された画像

## 要件
1. 引数から描画したいものを解析
2. 必要なデータの読み込み --> Readerクラスとして実装
3. 可視化

=================================

# reader系のプログラム用設定書
## 要望
空間トランスクリプトームデータを各小領域について読み込みたい
※最終的に可視化することを想定するため、可視化に必要なデータの読み込みというイメージ

### 現状
技術によって測定形式が異なるため統一的な方法では読み込めない

### 理想
subsetを指定した時、その領域の遺伝子発現のテーブルデータ、細胞領域のテーブルデータを読み込む
画像も読み込む？
後で他のデータを読み込みたくなる可能性があるため、できるだけ汎用性が高い形にしてほしい

## 要求
- 指定されたfovの空間トランスクリプトームデータ（MERSCOPE、Xenium）の読み込み
- 〃細胞領域データの読み込み（デフォルトの形を読みたい）

### input
must
- 入力データのあるディレクトリのパス
- 出力先ディレクトリのパス
- 測定方法（Xenium or MERSCOPE）が分かるような何か
- 読み込みたい視野（fov）

option
- 分割領域のサイズ


### output
- 読み込んだ後のデータ（すべてローカルピクセル単位系に変換する）


## 要件
1. 出力先ディレクトリの定義
2. 分割座標の決定
3. 遺伝子発現データの分割
4. 画像データの分割


### 1. 出力先ディレクトリの定義
- 出力先パスの存在確認
- subset用パスの存在確認+作成

### 2. 分割座標の決定
- 何らかのデータの読み込み
- 分割座標の決定
- 分割座標を記録

### 3. 遺伝子発現データの分割
- テーブルデータ全体を必要な列だけ読み込み
- 分割座標を基に抽出
- 座標変換（最終的にローカルのピクセル単位系が欲しい）
- 保存

単位系の変換
ローカル←→グローバルの変換
- オフセット(分割の始点となる座標の値)を足し引きする
- ピクセル単位系かマイクロメートル単位系かを統一すること

ピクセル←→マイクロメートル
- 何らかの倍率をかける

※Xeniumの場合
- 倍率が判明している

※MERSCOPEの場合
- `images/manifest.json`内の"microns_per_pixel", "bbox_microns"を用いて変換可能


### 4. 画像データの分割
- 画像データ全体の読み込み
- 分割座標を元にクリップ
- 保存

---

# viewer系統

## viewer_abc.py

### 抽象クラス: `AbstractDisplayArea`
`__init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int)`
param:
- `self`
- `input_dir:str` : 入力データのディレクトリパス
- `output_dir:str` : 出力データのディレクトリパス
- `fov:int` : 対象視野の番号
- `width:int, height:int` : 分割視野の幅と高さ

---

`view_area()`

---------------------------------------------------------------------------

## viewer_merscope.py

### 子クラス: `MerscopeDisplayArea(AbstractDisplayArea)`

`__init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int, z:int=None)`
param:
- `self` ~ `height` : 継承元のクラスと同じ
- `z:int` : z座標指定用（デフォルト値は `None`）

---

`view_area(self, title:str = "default", show_plot:bool = True, show_polygon:bool = True, image_keywords:List[str] = None)`
param:
- `self` : MerscopeDisplayAreaクラスのインスタンス
- `title:str`：描画画像のタイトルを指定する（デフォルトタイトルは"default"）
- `show_plot:bool`：遺伝子発現点のプロットを描画するか否かを指示するフラグ（デフォルト値はTrue）
- `show_polygon:bool`：細胞涼気の多角形を〃
- `image_keywords:List[str]`：謎のリスト（デフォルト値は `None`）

------------------------------------------------------------------------------

## viewer_xenium.py

### 子クラス: `XeniumDisplayArea(AbstractDisplayArea)`

`__init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int, resampling_factor:float=0.2125, z:int=None)`
param:
- `self` ~ `height` : 継承元のクラスと同じ
- `resampling_factor:float` : マイクロメートルとピクセル単位系を変換するための値
- `z:int` : z座標指定用（デフォルト値は `None`）

---

`view_area(self, title:str = "", show_plot:bool = True, show_polygon:bool = True, image_keywords:List[str] = None)`
param:
- `self` : MerscopeDisplayAreaクラスのインスタンス
- `title:str`：描画画像のタイトルを指定する（デフォルトタイトルはなし）
- `show_plot:bool`：遺伝子発現点のプロットを描画するか否かを指示するフラグ（デフォルト値はTrue）
- `show_polygon:bool`：細胞涼気の多角形を〃
- `image_keywords:List[str]`：謎のリスト（デフォルト値は `None`）

-----------------------------------------------------------------------------

## viewer.py
事実上、merscopeとxeniumの共通処理関数をまとめた場所みたいになっている

### クラス：`Plot_Polygon_Viewwer`

`__init__(self, img_dir, fov, z, title, width, height, gene_list)`
- `self`：クラスのインスタンス
- `img_dir:str`：画像データのあるディレクトリのパス
- `fov:int`：対象視野の番号
- `z:int`：z座標
- `title:str`：タイトル
- `width:int, height:int`：画像の縦横幅（デフォルト値は2048, 2048）
- `gene_list:List[str]`：遺伝子リスト（デフォルト値は `[]`）

---

`load_image(self, png_path)`
param:
- `self` : クラスのインスタンス
- `png_path:str`：読み込む画像のパス

return:
`Image.open(png_path)`：読み込んだ後の画像オブジェクト

---

`draw_gene_plot(self, ax: plt.Axes, df: DataFrame, show_legend:bool, size: int = 2)`
遺伝子発現点のプロットを描画する関数
（内部ではseabornのscatterplotを呼ぶ）

param:
- `self`：クラスのインスタンス
- `ax: plt.Axes`：描画するaxオブジェクト
- `df: DataFrame`：発現データのデータフレーム
- `show_legend:bool`：凡例を描画するか否かのフラグ
- `size: int`：プロットサイズ（defaultでは2）

---

`draw_polygon(self, ax: plt.Axes, x_data: float, y_data: float)`
１点を描画する関数
（内部ではmatplotlibのplotを呼ぶ）

param:
- `self`：クラスのインスタンス
- `ax: plt.Axes`：描画するaxオブジェクト
- `x_data: float, y_data: float`：点を描画する座標

---

`draw_poly()`
draw_polygonと全く同じ関数

---

`draw_multipolygon(self, ax: plt.Axes, fov_cell_df: DataFrame)`
多角形を描画する関数

param:
- `self`：クラスのインスタンス
- `ax: plt.Axes`：描画するaxオブジェクト
- `cell_df: DataFrame`: 細胞領域の座標情報を持つデータフレーム

---

`show_whole_fov(self, ax: plt.Axes, fov_cell_df: DataFrame, selected_df: DataFrame, show_plot: bool = True, show_polygon: bool = True)` 
フラグに応じて必要な関数を呼び出す処理

param:
- `self`
- `fov_cell_df: DataFrame`：細胞領域の座標情報を持つデータフレーム
- `selected_df: DataFrame`：遺伝子発現点のデータフレーム
- `show_plot:bool`：遺伝子発現点のプロットを描画するか否かを指示するフラグ（デフォルト値はTrue）
- `show_polygon:bool`：細胞涼気の多角形を〃

---

`get_imagepath(self, input_dir: str, image_keywords)`
image_keywords内の文字列を含む画像パスのリストを返す

param:
- `self`
- `input_dir: str`
- `image_keywords`: 

return:
- `List[str]`

---

`display_area(self, fov_cell_df: DataFrame, selected_df: DataFrame ,show_plot: bool = True, show_polygon: bool = True, image_keywords:List[str]=[], is_save: bool = False, ax: Optional[plt.Axes] = None)`

param:
- `self`
- `fov_cell_df: DataFrame`：細胞領域の座標情報を持つデータフレーム
- `selected_df: DataFrame`：遺伝子発現点のデータフレーム
- `show_plot:bool`：遺伝子発現点のプロットを描画するか否かを指示するフラグ（デフォルト値はTrue）
- `show_polygon:bool`：細胞涼気の多角形を〃
- `image_keywords:List[str]`：謎のリスト（デフォルト値は `None`）
- `is_save: bool`：保存するか否かのフラグ、デフォルト値は`False`
- `ax: Optional[plt.Axes]`：axオブジェクト、デフォルト値は`None`

================================================================

# reader系統

## reader_abc.py

### 抽象クラス：`AbstractReadData`
`__init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int)`
param
- `self`：インスタンス
- `input_dir:str`：入力データのディレクトリパス
- `output_dir:str`：出力データのディレクトリパス
- `fov:int`：対象視野の番号
- `width:int, height:int`：小領域のサイズ

---

`get_path(self)`
サブセットのディレクトリとその中の画像ディレクトリへのパスを返す

param
- `self`：AbstractReadDataクラスのインスタンス

return
- `subset_dir:str`
- `img_dir:str`

---

`read_transcripts_data(self)`
遺伝子データの読み込みを行う関数

param
- `self`：AbstractReadDataクラスのインスタンス

return
- `gene_df:DataFrame`：遺伝子発現点のデータのデータフレーム
- `gene_name_list:List[str]`：遺伝子名のリスト

---

`read_cell_data(self)`
抽象メソッド

param
- `self`：AbstractReadDataクラスのインスタンス


## reader_merscope.py

### 子クラス：`ReadData(AbstractReadData)`

`__init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int, z:int = None)`
param:
- `self`から`height`：継承元と同じ
- `z:int`：z座標（デフォルト値は`None`）

---

`transfer_position(self, x, y, x_init, y_init)`
やっていることは分かるけど何に使うか分からない

param:
- `self`：MerscopeReadDataクラスのインスタンス
- `x, y, x_init, y_init`

return:
- `x - x_init, y - y_init`：実際の値を返す

---

`transform_multipolygon(self, multipolygon, x_init, y_init)`
グローバル系でのマルチポリゴンの頂点座標をローカルな座標系に変換

param
- `self`：MerscopeReadDataクラスのインスタンス
- `multipolygon`：マルチポリゴン（頂点座標）
- `x_init, y_init`：ローカルな座標に変換するための始点

return
- `multipolygon.apply(lambda geom: transform(lambda x, y, z=None: self.transfer_position(x, y, x_init, y_init), geom))`

---

`read_cell_data(self, mosaic_path:Optional[str]=None)`
細胞データの読み込み
※抽象メソッドの中身

param:
- `self`：MerscopeReadDataクラスのインスタンス
- `mosaic_path`：細胞領域のデータがあるファイルパス

return
- `selected_fov_cell_df: DataFrame`：細胞領域の頂点データののデータフレーム


## reader_xenium.py

### 子クラス：`XeniumReadData`

`__init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int,resampling_factor:float =0.2125, z:int = None)`
param:
- `self`から`height`：継承元と同じ
- `resampling_factor:float`：マイクロメートル単位系とピクセル単位系とを変換するための倍率
- `z:int`：z座標（デフォルト値は`None`）

---

`filter_and_create_polygons(self, filtered_df: pd.DataFrame)`

param:
- `self`:XeniumReadDataクラスのインスタンス
- `filtered_df: DataFrame`：指定したfov+座標変換済のデータ

return:
- `gpd.GeoDataFrame(polygons, geometry='Geometry_local', crs="EPSG:4326"):gpd.GeoDataFrame`：ポリゴンに変換した後のデータ（merscopeと同じ関数を適用するため）

---

`read_cell_data(self, parquet_path:Optionl[str]=None)`
細胞データの読み込み
※抽象メソッドの中身

param:
- `self`：MerscopeReadDataクラスのインスタンス
- `parquet_path`：細胞領域のデータがあるファイルパス

return
- `cell_df: DataFrame`：細胞領域の頂点データののデータフレーム

