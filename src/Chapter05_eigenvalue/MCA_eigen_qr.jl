using LinearAlgebra
function MCA_eigen_qr(A; iter=20, tol=1e-15)  
    n = size(A,1)

    # ハウスホルダー法によるヘッセンベルク化
    W = zeros(n,n)
    H = copy(float(A))
    for j = 1:n-2
        W[j+1:n,j], r = MCA_householder_gen(H[j+1:n,j])
        H[j+1:n,j+1:n] = MCA_householder_op(W[j+1:n,j],H[j+1:n,j+1:n])
        H[j+1,j] = r; H[j+2:n,j] = zeros(n-j-1)
        H[:,j+1:n]= MCA_householder_op(W[j+1:n,j],H[:,j+1:n]')'
    end

    # ギブンス回転によるヘッセンベルク行列のQR法
    cc = zeros(n); ss = zeros(n)
    G = Matrix{Float64}(I,n,n)
    m = n
    for k = 1:iter
        s = H[m,m]                               # レイリー商シフト
        H = H - s * I(n)
        for i = 1:m-1
            cc[i], ss[i], r = MCA_givens_gen(H[i,i],H[i+1,i])
            H[i:i+1,i+1:n] = MCA_givens_op(cc[i],ss[i],H[i:i+1,i+1:n])
            H[i,i] = r; H[i+1,i] = 0
            G[i:i+1,:] = MCA_givens_op(cc[i],ss[i],G[i:i+1,:])
        end
        for i = 1:m-1
            H[1:i+1,i:i+1] = MCA_givens_op(cc[i],ss[i],H[1:i+1,i:i+1]')'
        end
        H = H + s * I(n)
        
        if abs(H[m,m-1] / H[m,m]) < tol          # 収束判定（個別）
            m = m - 1                            # デフレーション
            if m == 1                            # 収束判定（全体）
                break
            end
        end

    end

    # 固有ベクトル計算
    X = Matrix{Float64}(I,n,n)
    H = UpperTriangular(H);
    for i = 1:n
        X[1:i-1,i] = (H[i,i] * I(i-1) - H[1:i-1,1:i-1]) \ H[1:i-1,i]
    end
    X = G' * X
    for j = n-2:-1:1
        X[j+1:n,:] = MCA_householder_op(W[j+1:n,j],X[j+1:n,:])
    end
    X = X ./ sqrt.(sum(X.^2, dims=1))            # 固有ベクトルの正規化
    
    return (D = diag(H), X = X)
end