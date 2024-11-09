import os
import numpy as np
import pandas as pd
from mip import *
import time
from typing import Optional

# Rから呼ぶ関数
def my_LP_test(sce:pd.DataFrame, spe:pd.DataFrame, output_path: Optional[str] = None):
    """
    @param
    - sce: シングルセル参照データ
    - spe: 空間トランスクリプトーム解析のデータ（共通する遺伝子のみであるということを想定）
    - output_path: 出力フォルダを作りたいところへのパス
    
    出力ファイルパス
    - output_path/LP_outputs -> pd.csv, md.csv, x.csv, calc.txt が保存される
    - output_path/LP_outputs/LPfiles -> .lp が保存される
    """
    tmp = list(range(len(sce.columns)))
    ref_list = [tmp for _ in range(len(spe.columns))]
    
    # sceで全ての細胞において発現していない遺伝子を削除
    # 警告を回避するためにインデックスを一致させたフィルタリング
    filter_condition = (sce.sum(axis=1) != 0)
    spe = spe.loc[filter_condition]  # 明示的にインデックスを一致させてフィルタリング
    sce = sce.loc[filter_condition]
    # spe = spe[sce.sum(axis=1) != 0]
    # sce = sce[sce.sum(axis=1) != 0]

    # 遺伝子が発現している細胞のインデックスリストを作成
    cell_indexlist = []
    for gene in sce.index:
        cell_indexlist += np.nonzero(np.array(np.isnan(sce.loc[gene,]) != True))

    # 出力ファイルパス
    if output_path is None:
        output_path = os.getcwd()
    
    output_dir = os.path.join(output_path, "LP_outputs")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    lp_path = os.path.join(output_dir, "LPfiles")
    if not os.path.exists(lp_path):
        os.makedirs(lp_path)

    # 線形計画法を実行
    pd_df, md_df, x_df, text = LP(sce, spe, cell_indexlist, ref_list, lp_path)

    # データを保存
    pd_df.to_csv(output_dir + "/pd.csv")
    md_df.to_csv(output_dir + "/md.csv")
    x_df.T.to_csv(output_dir + "/x.csv")
      
    f = open(output_dir+"/calc.txt", "w")
    f.write(text)
    f.close()

    return x_df.T, pd_df, md_df

# 最適化の目的関数を定義
def obj_def(s, spe, spot):
    s += "".join([" + pd,"+str(gene)+","+str(spot)+" + md,"+str(gene)+","+str(spot) for gene in range(np.shape(spe)[0])])
    return s

# 最適化の制約条件を定義
def const_def(s, spe, sce , cell_indexlist, spot, ref_list):
    for gene in range(np.shape(spe)[0]):
        tmp = " "
        tmp += "".join([" + "+str(sce.iloc[gene,cell]) + " x,"+str(cell)+","+str(spot) for cell in list(set(ref_list[spot]) & set(cell_indexlist[gene]))])
        s += "c,"+str(gene)+","+str(spot)+":"+ tmp + " + pd,"+str(gene)+","+str(spot)+" - md,"+str(gene)+","+str(spot)+" = "+ str(spe.iloc[gene,spot]) +"\n"
    return s

# 線形計画法を実行して結果を返す
def LP(sce, spe, cell_indexlist, ref_list, lp_path):
    pd_df = pd.DataFrame(index=spe.index, columns=spe.columns, dtype=float)
    md_df = pd.DataFrame(index=spe.index, columns=spe.columns, dtype=float)
    x_df = pd.DataFrame(index=sce.columns, columns=spe.columns, dtype=float)
    
    time_sta = time.time()
    for spot in range(np.shape(spe)[1]):
        lp_file = lp_path + spe.columns[spot] +".lp"

        s = "minimize\n"
        s = obj_def(s, spe, spot)
        s += "\n\nsubject to\n"
        with open(lp_file , mode='w') as f:
            f.write(s)
        s = ""

        s = const_def(s, spe, sce , cell_indexlist, spot, ref_list)
        s += "\nend"
        with open(lp_file, mode='a') as f:
            f.write(s)
        s = ""

        m = Model(solver_name=CBC)#GUROBI or CBC
        m.read(lp_file)
        m.write(lp_file)
        # ログの非表示化
        m.verbose = 0
        # 最適化
        m.optimize()

        for i in range(len(m.vars)):
            vars = m.vars[i].name.split(",")
            if vars[0] == "pd":
                pd_df.iloc[int(vars[1]),int(vars[2])] = float(m.vars[i].x)
            if vars[0] == "md":
                md_df.iloc[int(vars[1]),int(vars[2])] = float(m.vars[i].x)
            if vars[0] == "x":
                x_df.iloc[int(vars[1]),int(vars[2])] = float(m.vars[i].x)
        
    time_end = time.time()
    time_rec = time_end - time_sta
    text = "time: " + str(round(time_rec, 2)) + "\nspot, cell, gene: " + str(len(spe.columns)) + " " + str(len(sce.columns)) + " " + str(len(sce.index)) + "\ntime/spot: " + str(round(time_rec/(spot+1),3))
    #print("time: ",round(time_rec, 2))
    #print("spot, cell, gene: ",len(spe.columns),len(sce.columns),len(sce.index))
    #print("time/spot: ", round(time_rec/(spot+1),3))
    
    return pd_df, md_df, x_df, text
