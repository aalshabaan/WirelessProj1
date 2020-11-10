function preamble=preamble_generate(length)

    preamble = zeros(length, 1);
    preamble(1:8) = ones(8,1);

    for i = 1: length-9
        k = preamble(8+i:-1:1+i);
        new = k(4) + k(5) + k(6) + k(8);
        new = mod(new, 2);
        preamble(i+9) = new;
    end    
end