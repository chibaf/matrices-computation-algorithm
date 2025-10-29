function MCA_bitstring(x; sign=1, exp=8, signif=23, bias=127)
    xbit = bitstring(x)
    
    k = sign + exp + 1
    x_signif = xbit[k:k+3]
    for i = 1:ceil(signif/4)-2
        k = k + 4
        x_signif = x_signif*" "*xbit[k:k+3]
    end
    x_signif = x_signif*" "*xbit[k+4:sign+exp+signif]    

    x_sign = (-1)^Int(xbit[1])
    x_exp = Int(parse(Int, xbit[2:exp+1], base=2)) - bias
    
    println(x_sign, ".", x_signif, " x 2^", x_exp)
end

function MCA_err(x,xast)
    println("絶対誤差 eA = ", abs(float(x) - xast))
    println("相対誤差 eR = ", abs(float(x) - xast) / abs(xast))
end