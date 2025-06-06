using LinearAlgebra
function MCA_cholesky(A)
    n = size(A,1)
    L = copy(LowerTriangular(float(A)))
    for j = 1:n
        L[j,j] = sqrt(L[j,j])
        L[j+1:n,j] = L[j+1:n,j] / L[j,j]
        for i = j+1:n
            L[i,j+1:i] -= L[i,j] * L[j+1:i,j]
        end
    end

    return L
end