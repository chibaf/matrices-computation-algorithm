using LinearAlgebra
function MCA_eigen_power(A; v=randn(size(A,1)),iter=100,tol=1e-10)
    x = v / norm(v)
    for k = 1:iter
        v = A * x
        global lambda = x' * v
        r = v - lambda * x
        
        if norm(r) < tol     　　　　　# 収束判定
            break
        end
        
        x = v / norm(v)
    end
    
    return (lambda = lambda, x = x)
end