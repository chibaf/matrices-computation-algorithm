using LinearAlgebra
function MCA_lsq_lsqr(A,b; x=zeros(size(A,2)),iter=100,tol=1e-10)
    bnorm = norm(A'*b)
    resvec = zeros(iter+1)

    r = b - A * x
    beta1 = norm(r)
    u = r / beta1
    for k = 1:iter
        
    end
    
    return (x=x, resvec=resvec)
end