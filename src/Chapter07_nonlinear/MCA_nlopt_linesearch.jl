function MCA_nlopt_linesearch(x,d; c=[1e-4,0.1],iter=20)
    a = 0; b = Inf; alpha = 1
    fx = MCA_setf(x)
    gx = MCA_setg(x)
    gxd = gx'*d
    for i = 1:iter
        xp = x + alpha*d
        if MCA_setf(xp) > fx + alpha * c[1] * gxd
            b = alpha; alpha = (a+b)/2
        elseif MCA_setg(xp)'*d < c[2] * gxd
            a = alpha; alpha = (a+b)/2
            if b == Inf
                alpha = 2 * a
            end
        else
            break
        end
    end
    
    return alpha
end