using LinearAlgebra
function MCA_lsq_nrsor(A,b; x=zeros(size(A,2)),omega=1.0,iter=100,tol=1e-10)
    n = size(A,2)
    bnorm = norm(A'*b)
    Dinv = 1 ./ sum(A.^2, dims=1)
    r = b - A * x
    for k = 1:iter
        xm = copy(x);
        for i = 1:n
            d = omega * Dinv[i] * r' * A[:,i]
            x[i] = x[i] + d
            r = r - d * A[:,i]
        end

        if norm(x-xm) / bnorm < tol    # 収束判定
            break
        end
    end

    return x
end