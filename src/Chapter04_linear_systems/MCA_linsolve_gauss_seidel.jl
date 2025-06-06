using LinearAlgebra
function MCA_linsolve_gauss_seidel(A,b; x=zeros(size(b)),iter=100,tol=1e-10)
    n = size(A,1)
    bnorm = norm(b)
    for k = 1:iter
        xm = copy(x)
        for i = 1:n
            index = [1:i-1; i+1:n]
            x[i] = (b[i] - A[i,index]' * x[index]) / A[i,i]
        end
        
        if norm(x-xm) / bnorm < tol    # 収束判定
            break
        end
    end

    return x
end