using LinearAlgebra
function MCA_expmv_arnoldi(A,b; iter = 20, tol = 1e-10)
    n = size(A,1)
    H = zeros(iter+1,iter)
    V = zeros(n,iter+1)
    e = zeros(iter); e[1] = norm(b)
    V[:,1] = b / e[1]
    for k = 1:iter
        w = A * V[:,k]                 # アーノルディ過程
        for i = 1:k                    # （修正グラム・シュミット）
            H[i,k] = V[:,i]' * w
            w -= H[i,k] * V[:,i]
        end
        H[k+1,k] = norm(w)
        V[:,k+1] = w / H[k+1,k]

        y = exp(H[1:k,1:k]) * e[1:k]   # H_k の行列関数ベクトル積
        norm_r = H[k+1,k] * abs(y[k])
        if norm_r / e[1] < tol
            global x = V[:,1:k] * y
            break
        end
    end
    
    return x
end