using LinearAlgebra
function MCA_qr_givens(A; Type="thin")
    m, n = size(A)
    C = zeros(m,n); S = zeros(m,n)
    R = copy(float(A))
    for j = 1:n
        for i = m-1:-1:j
            C[i,j], S[i,j], r = MCA_givens_gen(R[i,j],R[i+1,j])
            R[i:i+1,j+1:n] = MCA_givens_op(C[i,j],S[i,j],R[i:i+1,j+1:n])
            R[i,j] = r; R[i+1,j] = 0
        end
    end

    if Type == "thin"             # thin QR分解
        Q = Matrix{Float64}(I,m,n)
        for j = n:-1:1
            for i = j:m-1
                Q[i:i+1,j:n] = MCA_givens_op(C[i,j],-S[i,j],Q[i:i+1,j:n])
            end
        end
        R = R[1:n,:]
    elseif Type == "full"         # full QR分解
        Q = Matrix{Float64}(I,m,m)
        for j = n:-1:1
            for i = j:m-1
                Q[i:i+1,j:m] = MCA_givens_op(C[i,j],-S[i,j],Q[i:i+1,j:m])
            end
        end
    elseif Type == "implicit"     # Qを陽に構築しない
        Q = (C, S)
        R = R[1:n,:]
    end
    
    return (Q=Q, R=R)
end

function MCA_qr_givens_operate_Q(C,S,X,Type)
    m, n = size(C)
    Y = copy(float(X))
    if Type == "QX"      # Q*Xを計算
        for j = n:-1:1
            for i = j:m-1
                Y[i:i+1,:] = MCA_givens_op(C[i,j],-S[i,j],Y[i:i+1,:])
            end
        end
    elseif Type == "QTX" # Q^T*Xを計算
        for j = 1:n
            for i = m-1:-1:j
                Y[i:i+1,:] = MCA_givens_op(C[i,j],S[i,j],Y[i:i+1,:])
            end
        end
    end
    
    return Y
end