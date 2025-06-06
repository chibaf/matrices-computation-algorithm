using LinearAlgebra
function MCA_linsolve_lu(A,b)
    n = size(A,1)
    x = zeros(n); y = zeros(n)
    F = lu(A)           # LU分解（LinearAlgebraパッケージを利用）
    for i = 1:n         # 前進代入による Ly=b の求解
        y[i] = b[F.p[i]] - F.L[i,1:i-1]' * y[1:i-1]
    end
    for i = n:-1:1      # 後退代入による Ux=y の求解
        x[i] = (y[i] - F.U[i,i+1:n]' * x[i+1:n]) / F.U[i,i]
    end

    return x
end