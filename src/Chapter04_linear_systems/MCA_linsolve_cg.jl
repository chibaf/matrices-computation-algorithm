using LinearAlgebra
function MCA_linsolve_cg(A,b; x=zeros(size(b)),iter=100,tol=1e-10)
    bnorm = norm(b)
    resvec = zeros(iter+1)

    p = r = b - A*x
    rr = r' * r;
    resvec[1] = sqrt(rr) / bnorm
    for k = 1:iter
        ap = A * p
        pap = p' * ap
        alpha = rr / pap
        x = x + alpha * p
        r = r - alpha * ap
        rrm = rr
        rr = r' * r

        resvec[k+1] = sqrt(rr) / bnorm
        if resvec[k+1] < tol           # 収束判定
            resvec = resvec[1:k+1]
            break
        end

        beta = rr / rrm
        p = r + beta * p
    end

    return (x=x, resvec=resvec)
end