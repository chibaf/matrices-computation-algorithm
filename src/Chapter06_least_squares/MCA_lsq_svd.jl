using LinearAlgebra
function MCA_lsq_svd(A, b)
    F = svd(A)                         # LinearAlgebra の svd を利用
    r = sum(F.S / F.S[1] .> eps())     # (数値的な)ランクを計算
    x = F.V[:,1:r] * (F.S[1:r] .\ (F.U[:,1:r]' * b))

    return x
end