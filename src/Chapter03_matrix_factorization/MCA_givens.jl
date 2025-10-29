using LinearAlgebra
function MCA_givens_gen(x1,x2)
    c = 1; s = 0; r = sqrt(x1^2 + x2^2)
    if r != 0
        c = x1 / r; s = -x2 / r
    end
    
    return c, s, r
end

function MCA_givens_op(c,s,x)
    return y = [c s; -s c]' * x
end 