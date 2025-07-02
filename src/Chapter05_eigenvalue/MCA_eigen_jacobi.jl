using LinearAlgebra
function MCA_eigen_jacobi(A; iter=100,tol=1e-15)
    function MCA_inner_givens_jacobi(Aij)
        if abs(Aij[1,2]) < tol
            c = 1; s = 0 
            a = Aij[1,1]; b = Aij[1,2]; d = Aij[2,2]
        elseif Aij[1,1] == Aij[2,2]
            c = s = sqrt(2)/2
            a = Aij[1,1]-Aij[1,2]; b = 0; d = Aij[1,1]+Aij[1,2]
        else
            tau = (Aij[2,2] - Aij[1,1]) / (2*Aij[1,2])
            t = sign(tau) / (abs(tau) + sqrt(1 + tau^2))
            c = 1 / sqrt(1 + t^2); s = t*c
            a = Aij[1,1] * c^2 + Aij[2,2] * s^2 - 2 * Aij[1,2] * c * s
            b = 0
            d = Aij[1,1] * s^2 + Aij[2,2] * c^2 + 2 * Aij[1,2] * c * s
        end
        return c, s, [a b; b d]
    end
    
    n = size(A,1)
    D = copy(float(A))
    X = Matrix{Float64}(I,n,n)
    for k = 1:iter
        for j = 1:n
            for i = 1:j-1
                ij =[i,j]
                c, s, Bij = MCA_inner_givens_jacobi(D[ij,ij])
                D[ij,:] = MCA_givens_op(c,s,D[ij,:])
                X[ij,:] = MCA_givens_op(c,s,X[ij,:])
                D[:,ij] = D[ij,:]'
                D[ij,ij] = Bij
            end
        end
        
        dD = diag(D)
        if norm(D - diagm(dD)) / norm(dD) < tol       # 収束判定
            break
        end
    end
    X = X'; X = X ./ sqrt.(sum(X.^2, dims=1))         # 固有ベクトルの正規化
    
    return (D = diag(D), X = X)
end