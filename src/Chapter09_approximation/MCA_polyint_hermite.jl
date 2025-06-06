function MCA_polyint_hermite(x,y,z)
    n = size(x,1)
    V1 = ones(n,2*n)
    V2 = zeros(n,2*n)
    for i = 2:2n
        V1[:,i] = x.^(i-1)
        V2[:,i] = (i-1) * x.^(i-2)
    end
    c = [V1; V2] \ [y; z]    # 線形方程式の求解
    return c
end