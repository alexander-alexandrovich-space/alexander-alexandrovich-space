function coded = LoRa_Hamming_enc(bits, CR)

% bits — вектор 0/1
% CR   — 1..4

bits = bits(:).';
Nb = length(bits);

% padding до кратности 4
pad = mod(4 - mod(Nb,4),4);
bits = [bits zeros(1,pad)];

Nblk = length(bits)/4;

coded = [];

for k = 1:Nblk
    
    d = bits((k-1)*4+1:k*4);
    
    d1=d(1); d2=d(2); d3=d(3); d4=d(4);
    
    % стандартные parity
    p1 = mod(d1 + d2 + d4,2);
    p2 = mod(d1 + d3 + d4,2);
    p3 = mod(d2 + d3 + d4,2);
    p4 = mod(d1 + d2 + d3 + d4,2); % используется только при CR=4
    
    switch CR
        
        case 1 % (5,4)
            cw = [d p1];
            
        case 2 % (6,4)
            cw = [d p1 p2];
            
        case 3 % (7,4)
            cw = [d p1 p2 p3];
            
        case 4 % (8,4)
            cw = [d p1 p2 p3 p4];
            
        otherwise
            error('CR must be 1..4');
    end
    
    coded = [coded cw];
    
end

end