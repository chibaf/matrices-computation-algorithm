using LinearAlgebra
function MCA_nlopt_quasinewton(x; iter=50,tol=1e-8,c=[1e-4,0.1],liter=20)
    function MCA_inner_BFGS(Binv,s,y)
        ys = y' * s; By = Binv * y; Bys = By * s'; ss = s * s'
        Binv = Binv + ((1 + y'*By/ys)*ss - Bys - Bys') / ys
        return Binv
    end
    
    Binv = I(size(x,1))
    gx = MCA_setg(x)
    resvec = zeros(iter)
    for k = 1:iter
        resvec[k] = norm(gx)
        if resvec[k] < tol                  # 収束判定
            resvec = resvec[1:k]
            break
        end

        d = -Binv * gx
        # ウルフ条件のための二分法による直線探索
        alpha = MCA_nlopt_linesearch(x,d,c=c,iter=liter)
        s = alpha * d; x = x + s
        gxm = gx; gx = MCA_setg(x)          # 勾配の計算
        y = gx - gxm
        Binv = MCA_inner_BFGS(Binv,s,y)     # BFGS法によるB^-1の更新
    end

    return (x = x, resvec = resvec)
end

