using LinearAlgebra
function MCA_qr_cholesky(A)
    F = cholesky(A' * A)
    R = F.U
    Q = A / R

    return (Q=Q, R=R)
end