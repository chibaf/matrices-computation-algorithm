using LinearAlgebra
function MCA_nlsolve_companion(p)
    n = size(p,1) - 1
    p = p ./ p[n+1]     # モニックであることを保証する
    C = [[zeros(1,n-1); I(n-1)] (-p[1:n])]
    F = eigen(C)        # LinearAlgebraパッケージの eigen を利用

    return F.values
end