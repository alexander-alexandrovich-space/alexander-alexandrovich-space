function bits = LoRa_Hamming_dec(coded, CR)

coded = coded(:).';

n = 4 + CR;
Nblk = length(coded)/n;

bits = [];

for k = 1:Nblk
    
    cw = coded((k-1)*n+1:k*n);
    
    d = cw(1:4);
    
    d1=d(1); d2=d(2); d3=d(3); d4=d(4);
    
    % пересчёт parity
    p1 = mod(d1 + d2 + d4,2);
    p2 = mod(d1 + d3 + d4,2);
    p3 = mod(d2 + d3 + d4,2);
    
    if CR >= 1
        s1 = mod(p1 + cw(5),2);
    else
        s1 = 0;
    end
    
    if CR >= 2
        s2 = mod(p2 + cw(6),2);
    else
        s2 = 0;
    end
    
    if CR >= 3
        s3 = mod(p3 + cw(7),2);
    else
       CR = 0;
    end
    
    syndrome = s1 + 2*s2 + 4*s3;
    
    % простая коррекция одной ошибки
    if syndrome ~= 0 && syndrome <= 4
        d(syndrome) = mod(d(syndrome)+1,2);
    end
    
    bits = [bits d];
    
end

end