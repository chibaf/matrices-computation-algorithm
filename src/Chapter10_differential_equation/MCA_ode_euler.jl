using LinearAlgebra
function MCA_ode_euler(A,f,x,h,tmax)
    n = size(x,1)
    kmax = Int(round(tmax/h))
    xvec = zeros(kmax,n); tvec = zeros(kmax)
    
    Ah = (I(n) + h * A)
    for k = 1:kmax
        t = k*h
        x = Ah*x + h*f(t-h)
        xvec[k,:] = x
        tvec[k] = t
    end

    return xvec, tvec
end