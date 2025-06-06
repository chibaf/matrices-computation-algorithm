using LinearAlgebra
function MCA_lsq_qr(A, b)
    n = size(A,2)
    F = qr(A)     　# LinearAlgebra の qr を利用
    x = F.R \ (F.Q[:,1:n]' * b)

    return x
end