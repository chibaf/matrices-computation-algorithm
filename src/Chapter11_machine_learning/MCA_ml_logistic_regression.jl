using LinearAlgebra
sigmoid(z) = 1 / (1 + exp(-z))
function MCA_setg(w)
    p = sigmoid.(MCA_X*w)
    return MCA_X' * (p-MCA_Y)
end
function MCA_setH(w)
    p = sigmoid.(MCA_X*w)
    return MCA_X' * diagm(p .* (1 .- p)) * MCA_X
end

function MCA_ml_logistic_regression(Xtrain,Ytrain,w; iter=50,tol=1e-8)
    global MCA_X = [ones(size(Xtrain,1)) Xtrain]
    global MCA_Y = Ytrain

    # ニュートン法により非線形最適化問題を求解
    w, resvec = MCA_nlopt_newton(w, iter=iter,tol=tol)

   return (w = w, resvec = resvec)
end