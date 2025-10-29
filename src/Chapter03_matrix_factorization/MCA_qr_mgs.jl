using LinearAlgebra
function MCA_qr_mgs(A)
    m, n = size(A)
    Q = zeros(m,n)
    R = zeros(n,n)
    for j = 1:n
        aj = A[:,j]
        for i = 1:j-1
            R[i,j] = Q[:,i]' * aj
            aj -= Q[:,i] * R[i,j]
        end
        R[j,j] = norm(aj)
        Q[:,j] = aj / R[j,j]
    end

    return (Q=Q, R=R)
end