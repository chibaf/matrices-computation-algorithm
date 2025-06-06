function MCA_lu_pivoting(A)
    n = size(A,1)
    L = zeros(n,n)
    U = copy(float(A))
    p = Vector(1:n)
    for j = 1:n-1
        L[j,j] = 1

        # （部分）ピボット選択
        pivot = argmax(abs.(U[j:n,j])) + (j-1)
        p[j], p[pivot] = p[pivot], p[j]
        U[j,j:n], U[pivot,j:n] = U[pivot,j:n], U[j,j:n]
        L[j,1:j-1], L[pivot,1:j-1] = L[pivot,1:j-1], L[j,1:j-1]
        
        for i = j+1:n
            lij = U[i,j] / U[j,j]
            U[i,j] = 0
            U[i,j+1:n] -= lij * U[j,j+1:n]
            L[i,j] = lij
        end
    end
    L[n,n] = 1

    return (L=L, U=U, p=p)
end