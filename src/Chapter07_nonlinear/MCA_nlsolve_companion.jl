using LinearAlgebra
function MCA_nlsolve_companion(p)
    n = size(p,1)
    p = p ./ p[n]       # モニックであることを保証する
    C = [[zeros(1,n-2); I(n-2)] (-p[1:n-1])]
    F = eigen(C)        # LinearAlgebraパッケージの eigen を利用

    return F.values
end