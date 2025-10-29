using LinearAlgebra
function MCA_funcm_schur(A)
    n = size(A,1)
    S, Q = schur(A)               # LinearAlgebraパッケージのシューア分解
    F = diagm(MCA_setf.(diag(S)));# 別途関数 f を定義
    for j = 2:n
        for i = j-1:-1:1
            fs = F[i,i+1:j-1]' * S[i+1:j-1,j]
            sf = S[i,i+1:j-1]' * F[i+1:j-1,j]
            F[i,j] = S[i,j] * (F[i,i] - F[j,j]) + (fs - sf)
            F[i,j] = F[i,j] / (S[i,i] - S[j,j])
        end
    end
    X = Q * F * Q'

    return X
end