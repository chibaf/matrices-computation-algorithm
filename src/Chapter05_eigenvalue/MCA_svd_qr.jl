using LinearAlgebra
function MCA_svd_qr(A; iter=50, tol=1e-15)
    m, n = size(A)

    # ハウスホルダー法による上二重対角化
    P = zeros(m,n)
    Q = zeros(n,n)
    B = copy(float(A))
    for j = 1:n
        P[j:m,j] = MCA_householder_gen(B[j:m,j])
        B[j:m,j+1:n] = MCA_householder_op(P[j:m,j],B[j:m,j+1:n])
        B[j,j] = norm(B[j:m,j]); B[j+1:m,j] = zeros(m-j)

        if j <= n-2
            Q[j+1:n,j] = MCA_householder_gen(B[j,j+1:n])
            B[j+1:m,j+1:n] = MCA_householder_op(Q[j+1:n,j],B[j+1:m,j+1:n]')'
            B[j,j+1] = norm(B[j,j+1:n]); B[j,j+2:n] = zeros(n-j-1)
        end
    end
    B = B[1:n,:]

    # ギブンス回転による上二重対角QR法
    U = Matrix{Float64}(I,m,n)
    V = Matrix{Float64}(I,n,n)
    t = n
    for k = 1:iter
        # ウィルキンソンシフト
        a = B[t-1,t-1]; b = B[t-1,t]; d = B[t,t]
        e = eigen([a^2+b^2 b*d; b*d d^2])
        shift = minimum(e.values)
        
        c, s = MCA_givens_gen(B[1,1]^2-shift,B[1,1]*B[1,2])
        B[1:2,1:2] = MCA_givens_op(c,s,B[1:2,1:2]')'
        V[:,1:2] = MCA_givens_op(c,s,V[:,1:2]')'
        for i = 1:t-2
            c,s = MCA_givens_gen(B[i,i],B[i+1,i])
            B[i:i+1,i:i+2] = MCA_givens_op(c,s,B[i:i+1,i:i+2])
            U[:,i:i+1] = MCA_givens_op(c,s,U[:,i:i+1]')'
            c,s = MCA_givens_gen(B[i,i+1],B[i,i+2])
            B[i:i+2,i+1:i+2] = MCA_givens_op(c,s,B[i:i+2,i+1:i+2]')'
            V[:,i+1:i+2] = MCA_givens_op(c,s,V[:,i+1:i+2]')'
        end
        c,s = MCA_givens_gen(B[t-1,t-1],B[t,t-1])
        B[t-1:t,t-1:t] = MCA_givens_op(c,s,B[t-1:t,t-1:t])
        U[:,t-1:t] = MCA_givens_op(c,s,U[:,t-1:t]')'
        
        if abs(B[t-1,t] / B[t,t]) < tol          # 収束判定（個別）
            t = t - 1                            # デフレーション
            if t == 1                            # 収束判定（全体）
                break
            end
        end
        
    end
    S = diag(B)

    # 右特異ベクトル計算
    for j = n-2:-1:1
        V[j+1:n,:] = MCA_householder_op(Q[j+1:n,j],V[j+1:n,:])
    end

    # 左特異ベクトル計算
    for j = n:-1:1
        U[j:m,:] = MCA_householder_op(P[j:m,j],U[j:m,:])
    end
    
    return (U = U, S = S, V = V)
end