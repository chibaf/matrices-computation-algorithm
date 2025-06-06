using LinearAlgebra
function MCA_lsq_cimminonr(A,b;
        x=zeros(size(A,2)),omega=1.0,iter=100,tol=1e-10)
    bnorm = norm(A'*b)
    Dinv = 1 ./ sum(A.^2, dims=1)
    r = b - A * x
    ATr = A' * r
    for k = 1:iter
        d = omega * diagm(Dinv[1,:]) * ATr
        x = x + d
        r = r - A * d
        ATr = A' * r;

        if norm(ATr) / bnorm < tol     # 収束判定
            break
        end
    end

    return x
end