function MCA_ode_expm(A,x0,t)
    n = size(A,1); m = size(t,1)
    xt = zeros(n,m)
    for i = 1:m
        xt[:,i] = exp(t[i]*A)*x0  # 行列指数関数の計算
    end

    return xt'
end