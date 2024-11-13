# subset系のプログラム用設定書
## 要望
空間トランスクリプトームデータ全体を読み込み、小領域に分割する操作を行いたい

### 現状
膨大な空間トランスクリプトームデータ全体はそのまま扱うことが困難

### 理想
扱いやすいように小領域（subset）に分割したい
subsetごとにディレクトリを作成し、各データを格納する
空間トランスクリプトームデータの種類（技術）を問わず、以降の操作を汎用的に適用できるようにしたい

## 要求
- 空間トランスクリプトームデータ（MERSCOPE、Xenium）の読み込み
- 分割座標の決定
- 分割操作の実行（遺伝子発現データ、画像データ、その他）
- 分割座標の保存
- ピクセル単位系、ローカル座標系への変換

### input
must
- 入力データのあるディレクトリのパス
- 出力先ディレクトリのパス
- 測定方法（Xenium or MERSCOPE）が分かるような何か

option
- 分割領域のサイズ
- 分割する視野（特定の視野だけ抜き出す用）

### output
- 正常/異常終了用のメッセージ
- 下記ファイル構造への出力
指定した出力先ディレクトリ/<br>
| -- 分割座標を記録したファイル<br>
| -- subsetN/<br>
|     | -- 分割後のテーブルデータ（ローカルのピクセル座標系付き）<br>
|     | -- 分割後の画像データ<br>


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

# subset_creator_abc.py: 抽象クラス用
`class AbstractSubsetCreator(ABC)`
```{python}
def __init__(self, input_dir, output_dir, fov, width, height):
    self.input_dir = input_dir
    self.output_dir = output_dir
    self.fov = fov
    self.width = width
    self.height = height   
```
抽象メソッド一覧
- `make_directory()`
- `decide_crop_area()`
- `crop_table_data()`
- `crop_image()`
- `create_subset()`

# subset_creator_merscope.py : MERSCOPE用の子クラス


# subset_creator_xenium.py : Xenium用の子クラス

# subset_creator_common.py : 共通する関数をまとめたクラス
## `read_crop_info_json(json_path:str, fov:int)`
分割座標を読みだして返す関数

@params:
- `json_path:str`: 分割座標を記録したファイルのパス
- `fov:int`: 分割座標を知りたい視野
@return:
- `Tuple[int, int, int, int]`: 正常時は分割座標をピクセル単位系で返す [ x座標の視点、y座標の視点、x座標の終点、y座標の終点 ]
異常時は[ -1, -1, -1, -1 ] を返す

## `write_crop_info_json(json_path, data_to_write)`
分割座標を記録する関数

@params:
- `json_path:str`: 分割座標を記録するファイルのパス
- `data_to_write:dict`: ファイルに書き込みたい内容を辞書型で記述
@return: None

## `read_info_yaml(yaml_path:str)`
yamlファイルを読み出し、画像サイズを返す関数（Xenium用）

@params:
- `yaml_path:str`: yamlデータが保存されているパス
@return:
- `Tuple[int, int]`: 元データの画像のサイズ（width, height）

## `write_info_yaml(yaml_path:str, data_to_write)`
画像サイズをファイルに保存する関数

@params:
- `yaml_path:str`: yamlデータを保存するパス
- `data_to_write:dict`: ファイルに書き込みたい内容を辞書型で記述
@return: None

## `global_to_local(global_x:float, global_y:float, x_offset, y_offset)`
グローバル座標系→ローカル座標系に変換する関数（オフセットを取る）
※すべてピクセル単位系又はマイクロメートル単位系で統一すること

@params:
- `global_x:float, global_y:float`: グローバル座標系のxy座標の値
- `x_offset, y_offset`: オフセット（始点の値）
@returns:
- `Tuple[float, float]` : ローカル座標系に変換後の値(local_x, local_y)

## `microns_to_pixel(x_micron:float, y_micron:float, microns_per_pixel:float, bbox_x:float, bbox_y:float)` 
マイクロメートル単位系→ピクセル単位系に変換する関数

@params:
- `x_micron:float, y_micron:float`：マイクロメートル単位系のxyの値
- `microns_per_pixel:float`：倍率
- `bbox_x:float, bbox_y:float`：オフセット？
@return
- `Tuple[float, float]`：ピクセル単位系のxyの値

## `pixel_to_microns(x_pixel:float, y_pixel:float, microns_per_pixel:float, bbox_x:float, bbox_y:float)`
ピクセル単位系→マイクロメートル単位系に変換する関数

@params:
- `x_pixel:float, y_pixel:float`：ピクセル単位系のxyの値
- `microns_per_pixel:float`：倍率
- `bbox_x:float, bbox_y:float`：オフセット？
@return
- `Tuple[float, float]`：マイクロメートル単位系のxyの値



# subset.py : 雑多な処理のためのプログラムをまとめて置いておくところ
## `make_fov_info_json_xenium(input_json_path: str, output_json_path: str, resampling_factor: float)`
全てのfovに対して分割座標を一括計算し、書き込む関数

@params
- `input_json_path: str`：必要な情報が書かれたjsonファイルへのパス
- `output_json_path: str`：分割座標を保存するjsonパス
- `resampling_factor: float`：座標変換用の倍率
@return : none

## `image_size_xenium(tif_path:str, yaml_path:str)`
tifファイルを開き、画像サイズを取得し書き込む関数

@params
- `tif_path:str`：画像ファイルへのパス
- `yaml_path:str`：画像サイズを保存するyamlファイルへのパス
@return:
- `Tuple[int, int]`：画像サイズの幅と高さ（image_width, image_height）

## `decide_fov_hint(data_format:str, input_dir:str, output_dir:str, resampling_factor:Optional[float]=None)`
分割サイズを決定するためのヒントとして、画像サイズを表示する関数
"image size (width, height): {image_width}, {image_height}"

@params:
- `data_format:str`: データのフォーマットを指定（`xenium` or `merscope`）
- `input_dir:str`: 入力データのディレクトリ
- `output_dir:str`：xeniumの時、画像サイズを保存するためのファイルの置き場となるディレクトリへのパス
- `resampling_factor:Optional[float]=None`：xeniumの時、座標変換用の倍率
@return: none 


# subset_creator_merscope.py
## `error_handler()`
エラー処理用の関数

以下はクラス内`class MerscopeSubsetCreator(AbstractSubsetCreator)`
## `__init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int, image_keyword, z)`
- `image_keyword`: 切り抜く座標の文字列（空でもよい）
- `z`：z座標（使うかもしれないからおいてみた）

以下は抽象メソッドの具体的な実装
## `make_directory(self)`
サブセット用のディレクトリを作成しパスを返す関数

@params:
- `self`: MerscopeSubsetCreatorクラスのインスタンス
@return:
- `subset_dir`: 結果を保存するサブセット用ディレクトリのパス

## `decide_crop_area(self)`
1. 切り抜いた座標を保存するjsonパスを見に行く
2. 既に座標が保存されているならば切り抜く座標を返す
3. そうでないならば計算して保存する

@params:
- `self`: MerscopeSubsetCreatorクラスのインスタンス
@return:
- `Tuple[int, int, int, int]`: 切り抜く座標の始点と終点をピクセル単位系で返す（x_init, y_init, x_last, y_last）

## `crop_table_data(self, subset_dir:str, x_init:int, y_init:int, x_last:int, y_last:int)`
遺伝子発現データを切り抜いて保存する関数

@params:
- `self`: MerscopeSubsetCreatorクラスのインスタンス
- `subset_dir:str`: 保存先サブセットディレクトリのパス
- `x_init:int, y_init:int, x_last:int, y_last:int`: 切り抜く座標の始点と終点、ピクセル単位系
@return:none


以下はこのクラス独自の関数
## `read_manifest_json(self)`
manifest.jsonファイルを読み込み必要な情報を返す

@params
- `self`: MerscopeSubsetCreatorクラスのインスタンス
@return
- `Tuple[float, float, float, int, int]`: "microns_per_pixel, bbox_microns[0], bbox_microns[1], mosaic_width_pixels, mosaic_height_pixels"の値を返す


## `calculate_crop_coordinates(self, json_path:str, x_min_micron:float, y_min_micron:float)`
切り抜く始点をマイクロメートル単位系で与え、切り抜く座標を計算し保存する関数

@params:
- `self`: MerscopeSubsetCreatorクラスのインスタンス
- `json_path:str`: 切り抜く座標を保存するjsonパス
- `x_min_micron:float, y_min_micron:float`: マイクロメートル単位系での最小値（始点）
@return:
- `Tuple[int, int, int, int]`: 切り抜く座標の始点と終点をピクセル単位系で返す（x_init, y_init, x_last, y_last）

## `get_min_coordinates(self)`
発現データから指定fov内での最小のxyの値を抽出して返す

@params:
- `self`: MerscopeSubsetCreatorクラスのインスタンス
@return:
- `Tuple[float, float]`: xyの最小値（グローバル単位系かつマイクロメートル単位系：発現データそのまま）


## `get_imagepath(self, image_keywords)`
tifファイルのパスを取得する関数

@params:
- `self`: MerscopeSubsetCreatorクラスのインスタンス
- `image_keywords:Optional[none, str, List[str]]`: これに含まれる文字列を含むパスを選択する
@returns:
- `List[str]`: .tifファイルのパスのリスト



## subset_creator_xenium.py
### `error_handler()`
エラー処理用の関数

以下はクラス内`class XeniumSubsetCreator(AbstractSubsetCreator)`
## `__init__(self, input_dir:str, output_dir:str, fov:int, width:int, height:int, resampling_factor:float)`
- `resampling_factor:float`: 倍率


以下は抽象メソッドの具体的な実装
## `make_directory(self)`
サブセット用のディレクトリを作成しパスを返す関数

@params:
- `self`: XeniumSubsetCreatorクラスのインスタンス
@return:
- `subset_dir`: 結果を保存するサブセット用ディレクトリのパス

## `decide_crop_area(self)`
1. 切り抜いた座標を保存するjsonパスを見に行く
2. 既に座標が保存されているならば切り抜く座標を返す
3. そうでないならば計算して保存する

@params:
- `self`: XeniumSubsetCreatorクラスのインスタンス
@return:
- `Tuple[int, int, int, int]`: 切り抜く座標の始点と終点をピクセル単位系で返す（x_init, y_init, x_last, y_last）


## `crop_table_data(self, subset_dir:str, x_init:int, y_init:int, x_last:int, y_last:int)`
遺伝子発現データを切り抜いて保存する関数

@params:
- `self`: XeniumSubsetCreatorクラスのインスタンス
- `subset_dir:str`: 保存先サブセットディレクトリのパス
- `x_init:int, y_init:int, x_last:int, y_last:int`: 切り抜く座標の始点と終点、ピクセル単位系
@return:none

## `crop_image(self, subset_dir:str, x_init:int, y_init:int, x_last:int, y_last:int, is_show:bool=True)`
指定された座標で画像を切り抜き保存する関数
領域を提示することもある関数

@params:
- `self`: XeniumSubsetCreatorクラスのインスタンス
- `subset_dir:str`: 保存先サブセットディレクトリのパス
- `x_init:int, y_init:int, x_last:int, y_last:int`: 切り抜く座標の始点と終点、ピクセル単位系
- `is_show:bool=True`: 領域の提示が必要ならばTrue
@return: none


## `create_subset(self)`
最終的に全ての処理を実施する関数

@params:
- `self`: XeniumSubsetCreatorクラスのインスタンス
@return: none


以下はこのクラス独自の関数
## `calculate_crop_coordinates_Xenium(self, json_path:str, x_min_pixel:int, y_min_pixel:int, image_width:int, image_height:int)`
切り抜く座標を計算し保存する関数

@params:
- `self`: XeniumSubsetCreatorクラスのインスタンス
- `json_path:str`: 切り抜く座標を保存するファイルパス
- `x_min_pixel:int, y_min_pixel:int`：始点となるxyの最小値（ピクセル単位系）
- `image_width:int, image_height:int`：画像のサイズ
@return:
- `Tuple[int, int, int, int]`: 切り抜く座標の視点（x_init, y_init, x_last, y_last）






