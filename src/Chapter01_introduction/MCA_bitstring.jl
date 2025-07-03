function MCA_bitstring(x; exp=8, sig=23, bias=127)
    xbit = bitstring(x)
    s = 1 + exp + 1
    xx = xbit[s:s+3]
    for i = 1:ceil(sig/4)-2
        s = s + 4
        xx = xx*" "*xbit[s:s+3]
    end
    xx = xx*" "*xbit[s+4:1+exp+sig]
    println((-1)^Int(xbit[1]), ".", xx, " x 2^", Int(parse(Int, xbit[2:exp+1], base=2)) - bias)
end