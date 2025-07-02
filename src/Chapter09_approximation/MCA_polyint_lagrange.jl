function MCA_polyint_lagrange(x,y)
    n = size(x,1)
    V = x .^ (0:n-1)'   # ヴァンデルモンド行列の生成
    c = V \ y           # 線形方程式の求解
    
    return c
end