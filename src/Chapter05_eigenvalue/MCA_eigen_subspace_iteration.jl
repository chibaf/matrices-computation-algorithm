using LinearAlgebra
function MCA_eigen_subspace_iteration(A; 
        m=2,V=randn(size(A,1),m),iter=100,tol=1e-10)

    F = qr(V); V = Array(F.Q)          # thin QR分解
    for i = 1:iter
        W = A * V
        global lambda, Y = eigen(V'*W) # レイリー・リッツの技法
        global X = V * Y
        R = W * Y - X * diagm(lambda)
        
        if norm(R) < tol               # 収束判定
            break
        end
        
        F = qr(W); V = Array(F.Q)      # thin QR分解
    end

    return (lambda = lambda, X = X)
end