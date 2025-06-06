function MCA_polyint_lagrange(x,y)
    n = size(x,1)
    V = ones(n,n)
    for i = 2:n
        V[:,i] = x.^(i-1)
    end
    c = V \ y      # 線形方程式の求解
    return c
end