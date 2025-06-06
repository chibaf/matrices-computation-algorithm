function MCA_nlopt_cg(x; iter=50,tol=1e-8,c=[1e-4,0.1],liter=20)
    gx = MCA_setg(x)
    p = -copy(gx)
    resvec = zeros(iter+1)
    resvec[1] = norm(gx)
    for k = 1:iter
        # ウルフ条件のための二分法による直線探索
        alpha = MCA_nlopt_linesearch(x,p,c=c,iter=liter)
        x = x + alpha * p

        resvec[k+1] = norm(gx)
        if resvec[k+1] < tol           # 収束判定
            resvec = resvec[1:k+1]
            break
        end
        
        gxm = gx
        gx = MCA_setg(x)               # 勾配の計算
        beta = (gx'*gx) / (gxm'*gxm)   # FRの公式に基づくbetaの設定
        p = -gx + beta * p
    end
    
    return (x = x, resvec = resvec)
end