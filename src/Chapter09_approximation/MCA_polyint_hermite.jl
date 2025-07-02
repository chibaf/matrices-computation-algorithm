function MCA_polyint_hermite(x,y,z)
    n = size(x,1)
    V1 = x .^ (0:2*n-1)'
    V2 = [zeros(n) V1[:,1:2*n-1] .* (1:2*n-1)']
    c = [V1; V2] \ [y; z]    # 線形方程式の求解

    return c
end