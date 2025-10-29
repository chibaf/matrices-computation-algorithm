using LinearAlgebra
function MCA_ml_ridge_regression(Xtrain,Ytrain,lambda)
    n, m = size(Xtrain)
    x = sum(Xtrain,dims=1); y = sum(Ytrain,dims=1)
    A = [n x; x' Xtrain'*Xtrain + lambda * I(m)]
    b = [y; Xtrain'*Ytrain]

    W = A \ b                     # 正規方程式の求解
    
    return W
end