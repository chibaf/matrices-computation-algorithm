using LinearAlgebra
function MCA_qr_cholesky(A)
    F = cholesky(A' * A)
    Q = A / F.U

    return (Q=Q, R=F.U)
end