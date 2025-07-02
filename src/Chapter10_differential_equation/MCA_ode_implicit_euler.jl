using LinearAlgebra
function MCA_ode_implicit_euler(A,f,x,h,tmax)
    n = size(x,1)
    kmax = Int(round(tmax/h))
    xvec = zeros(kmax,n);  tvec = zeros(kmax)
    
    F = lu(I(n) - h * A)     # 反復の外でLU分解を実施
    for k = 1:kmax
        t = k*h
        x = F \ (x + h*f(t)) # 線形方程式の求解（LU分解を繰り返し利用）
        xvec[k,:] = x
        tvec[k] = t
    end

    return xvec, tvec
end