using LinearAlgebra
function MCA_nlopt_newton(x; iter=50,tol=1e-8)
    resvec = zeros(iter)
    for k = 1:iter
        gx = MCA_setg(x)     # 勾配の計算（問題に応じて別途関数を定義）

        resvec[k] = norm(gx)
        if resvec[k] < tol   # 収束判定
            resvec = resvec[1:k]
            break
        end
        
        Hx = MCA_setH(x)     # ヘッセ行列の計算（問題に応じて別途関数を定義）
        x = x - Hx \ gx
    end

    return (x = x, resvec = resvec)
end