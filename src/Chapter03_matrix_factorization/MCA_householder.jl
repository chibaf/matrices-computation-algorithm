using LinearAlgebra
function MCA_householder_gen(u)
    e = zeros(size(u)); e[1] = norm(u); w = u - e
    if norm(w) != 0
        w = w / norm(w)
    end
    
    return w, e[1]
end

function MCA_householder_op(w,x)
    return y = x - 2 * w * (w' * x)
end