using LinearAlgebra
function MCA_lsq_cgls(A,b; x=zeros(size(A,2)),iter=100,tol=1e-10)
    bnorm = norm(A'*b)
    resvec = zeros(iter+1)

    r = b - A*x
    p = z = A'*r
    zz = z' * z
    resvec[1] = sqrt(zz) / bnorm
    for k = 1:iter
        ap = A * p
        apap = ap' * ap
        alpha = zz / apap
        x = x + alpha * p
        r = r - alpha * ap
        z = A'*r
        zzm = zz
        zz = z' * z

        resvec[k+1] = sqrt(zz) / bnorm
        if resvec[k+1] < tol      # 収束判定
            resvec=resvec[1:k+1]
            break
        end

        beta = zz / zzm
        p = z + beta * p
    end

    return (x=x, resvec=resvec)
end