using LinearAlgebra
function MCA_eigen_arnoldi(A;
        m=20,v=randn(size(A,1)),iter=100,tol=1e-10)
    
    n = size(A,1)
    for k = 1:iter
        H = zeros(m+1,m)
        V = zeros(n,m+1)
        V[:,1] = v / norm(v)
        for j = 1:m
            w = A * V[:,j]                  # アーノルディ過程
            for i = 1:j                     # （修正グラム・シュミット）
                H[i,j] = V[:,i]' * w
                w -= H[i,j] * V[:,i]
            end
            
            d, Y = eigen(H[1:j,1:j])        # レイリー・リッツの技法
            id = argmax(abs.(d))            # 絶対値最大固有値を計算
            global lambda = d[id]
            global x = V[:,1:j] * Y[:,id]
            H[j+1,j] = norm(w)
            norm_r = abs(H[j+1,j] * Y[j,id])
            
            if norm_r < tol                 # 収束判定
                return (lambda = lambda, x = x)
            end
            
            V[:,j+1] = w / H[j+1,j]
        end
        v = x
    end

    return (lambda = lambda, x = x)
end