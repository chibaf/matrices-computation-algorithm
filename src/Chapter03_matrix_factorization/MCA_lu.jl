function MCA_lu(A)
    n = size(A,1)
    L = zeros(n,n)
    U = copy(float(A))
    for j = 1:n-1
        L[j,j] = 1
        for i = j+1:n
            lij = U[i,j] / U[j,j]
            U[i,j] = 0
            U[i,j+1:n] -= lij * U[j,j+1:n]
            L[i,j] = lij
        end
    end
    L[n,n] = 1

    return (L=L, U=U)
end