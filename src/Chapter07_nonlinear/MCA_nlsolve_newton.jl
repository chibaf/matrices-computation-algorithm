using LinearAlgebra
function MCA_nlsolve_newton(x; iter=30,tol=1e-8)
    resvec = zeros(iter)
    for k = 1:iter
        fx = MCA_setf(x)     # 関数値の計算（問題に応じて別途関数を定義）

        resvec[k] = norm(fx)
        if resvec[k] < tol   # 収束判定
            resvec = resvec[1:k]
            break
        end
        
        Jx = MCA_setJ(x)     # ヤコビ行列の計算（問題に応じて別途関数を定義）
        x = x - Jx \ fx
    end

    return (x = x, resvec = resvec)
end