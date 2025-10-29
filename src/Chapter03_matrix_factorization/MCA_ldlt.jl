using LinearAlgebra
function MCA_ldlt(A)
    n = size(A,1)
    L = copy(LowerTriangular(float(A)))
    D = zeros(n)
    for j = 1:n
        D[j] = L[j,j]
        L[j,j] = 1
        L[j+1:n,j] = L[j+1:n,j] / D[j]
        for i = j+1:n
            L[i,j+1:i] -= D[j] * L[i,j] * L[j+1:i,j]
        end
    end

    return (L=L, D=D)
end