using LinearAlgebra
function MCA_svd_gkl(A, q; iter = 100, tol = 1e-10)
    m, n = size(A)
    P = zeros(m,iter)
    Q = zeros(n,iter)
    B = zeros(iter,iter)

    Q[:,1] = q / norm(q)
    p = A * Q[:,1]
    for k = 1:iter
        B[k,k] = norm(p)
        P[:,k] = p / B[k,k]

        U, S, V = svd(B[1:k,1:k])           # Bk の特異値分解
        
        q = A' * P[:,k] - B[k,k] * Q[:,k]
        B[k,k+1] = norm(q)
        Q[:,k+1] = q / B[k,k+1]
        
        norm_r = abs(B[k,k+1] * U[k,1])
        if norm_r < tol || k == iter
            global u = P[:,1:k] * U[:,1]
            global v = Q[:,1:k] * V[:,1]
            global sigma = S[1]
            break
        end
        
        p = A * Q[:,k+1] - B[k,k+1] * P[:,k]
    end
    
    return (sigma = sigma, u = u, v = v)
end