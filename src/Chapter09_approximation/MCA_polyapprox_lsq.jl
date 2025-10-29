function MCA_polyapprox_lsq(x,y,m)
    V = x .^ (0:m-1)'   # ヴァンデルモンド行列の生成
    c = V \ y           # 最小二乗問題の求解
    
    return c
end