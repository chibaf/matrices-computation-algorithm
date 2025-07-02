using LinearAlgebra
function MCA_sqrtm_schur(A)
    n = size(A,1)
    S, Q = schur(A)               # LinearAlgebraパッケージのシューア分解
    F = diagm(sqrt.(diag(S)))
    for j = 2:n
        for i = j-1:-1:1
            F[i,j] = S[i,j] - F[i,i+1:j-1]' * F[i+1:j-1,j]
            F[i,j] = F[i,j] / (F[i,i] + F[j,j])
        end
    end
    X = Q * F * Q'

    return X
end