using LinearAlgebra
function MCA_qr_householder(A; Type="implicit")   
    m, n = size(A)
    W = zeros(m,n)
    R = copy(float(A))
    for j = 1:n
        W[j:m,j], r = MCA_householder_gen(R[j:m,j])
        R[j:m,j+1:n] = MCA_householder_op(W[j:m,j],R[j:m,j+1:n])
        R[j,j] = r; R[j+1:m,j] = zeros(m-j)
    end

    if Type == "thin"             # thin QR分解
        Q = Matrix{Float64}(I,m,n)
        for j = n:-1:1
            Q[j:m,j:n] = MCA_householder_op(W[j:m,j],Q[j:m,j:n])
        end
        R = R[1:n,:]
    elseif Type == "full"         # full QR分解
        Q = Matrix{Float64}(I,m,m)
        for j = n:-1:1
            Q[j:m,j:m] = MCA_householder_op(W[j:m,j],Q[j:m,j:m])
        end
    elseif Type == "implicit"     # Qを陽に構築しない
        Q = W
        R = R[1:n,:]
    end
    
    return (Q=Q, R=R)
end

function MCA_qr_householder_operate_Q(W,X,Type)
    m, n = size(W)
    Y = copy(float(X))
    if Type == "QX"      # Q*Xを計算
        for j = n:-1:1
            Y[j:m,:] = MCA_householder_op(W[j:m,j],Y[j:m,:])
        end
    elseif Type == "QTX" # Q^T*Xを計算
        for j = 1:n
            Y[j:m,:] = MCA_householder_op(W[j:m,j],Y[j:m,:])
        end
    end
    
    return Y
end