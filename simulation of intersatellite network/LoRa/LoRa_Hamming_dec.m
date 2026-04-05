function bits = LoRa_Hamming_dec(coded, CR)
    coded = coded(:).';
    n = 4 + CR;
    Nblk = length(coded)/n;
    bits = zeros(1, Nblk * 4);
    
    for k = 1:Nblk
        cw = coded((k-1)*n+1:k*n);
        d = cw(1:4);
        d1=d(1); d2=d(2); d3=d(3); d4=d(4);
        
        % Пересчет parity
        p1 = mod(d1 + d2 + d4,2);
        p2 = mod(d1 + d3 + d4,2);
        p3 = mod(d2 + d3 + d4,2);
        
        s1 = 0; s2 = 0; s3 = 0;
        if CR >= 1, s1 = mod(p1 + cw(5),2); end
        if CR >= 2, s2 = mod(p2 + cw(6),2); end
        if CR >= 3, s3 = mod(p3 + cw(7),2); end
        
        syndrome = s1 + 2*s2 + 4*s3;
        
        % Коррекция ошибок (только если CR >= 3)
        if syndrome ~= 0
            switch syndrome
                case 3
                    d(1) = ~d(1);
                case 5
                    d(2) = ~d(2);
                case 6
                    d(3) = ~d(3);
                case 7
                    d(4) = ~d(4);
                % Синдромы 1, 2, 4 означают ошибку в самом parity-бите, данные целы
            end
        end
        
        bits((k-1)*4+1:k*4) = d;
    end
end