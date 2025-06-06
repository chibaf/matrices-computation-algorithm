using LinearAlgebra
function MCA_eigen_shifted_inverse(A; 
        sigma=0,v=randn(size(A,1)),iter=100,tol=1e-10)
    
    F = lu(A - sigma * I(size(A,1)))   # (A-sigma I)のLU分解
    x = v / norm(v)
    for k = 1:iter
        v = F \ x                      # LU分解を再利用して方程式を求解
        mu = (x' * v)
        global lambda = 1 / mu + sigma
        norm_v = norm(v)
        r = (x + (sigma - lambda) * v) / norm_v
        x = v / norm_v
        
        if norm(r) < tol      # 収束判定
            break
        end
        
        
    end
    
    return (lambda = lambda, x = x)
end