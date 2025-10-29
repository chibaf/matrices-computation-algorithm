using LinearAlgebra
function MCA_qr_cgs(A)
    m, n = size(A)
    Q = zeros(m,n)
    R = zeros(n,n)
    for j = 1:n
        aj = A[:,j]
        R[1:j-1,j] = Q[:,1:j-1]' * aj
        aj -= Q[:,1:j-1] * R[1:j-1,j]
        R[j,j] = norm(aj)
        Q[:,j] = aj / R[j,j]
    end

    return (Q=Q, R=R)
end