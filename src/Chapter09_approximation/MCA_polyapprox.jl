function MCA_polyapprox(x,y,m)
    n = size(x,1)
    V = ones(n,m)
    for i = 2:m
        V[:,i] = x.^(i-1)
    end
    c = V \ y      # 最小二乗問題の求解
    return c
end