using LinearAlgebra

sigmoid(z) = 1.0 / (1.0 + exp(-z))
function MCA_setf(w)
    p = sigmoid.(X*w);
    return -sum(Y .* log.(p) .+ (1 .- Y) .* log.(1 .- p))
end;
function MCA_setg(w)
    p = sigmoid.(X*w);
    return X' * (p-Y)
end;

function MCA_logistic_regression(Xtrain,Ytrain,w; iter=50,tol=1e-8,c=[1e-4,0.1],liter=20)
    global X = [ones(size(Xtrain,1)) Xtrain];
    global Y = Ytrain;
    w, res = MCA_nlopt_quasinewton(w; iter=iter,tol=tol,c=c,liter=liter);

   return (w = w, resvec = resvec)
end