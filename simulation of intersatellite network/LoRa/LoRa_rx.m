function [data_hat, BER, BLER] = LoRa_rx(rx, data_tx, cfg)

SF = cfg.SF;
BW = cfg.BW;

M  = 2^SF;
Ns = M;

Ts = M/BW;
Fs = 2*BW;

t = (0:Ns-1)/Fs;

base_chirp = exp(1j*2*pi*(BW/(2*Ts)*t.^2));
down = conj(base_chirp);

start_payload = cfg.Preamble*Ns + 1;

data_hat = zeros(1,cfg.Nsym);

bitErrors = 0;
totalBits = cfg.Nsym * SF;
blockErrors = 0;

for n = 1:cfg.Nsym
    
    idx = start_payload + (n-1)*Ns;
    rxsym = rx(idx:idx+Ns-1);
    
    dechirped = rxsym .* down;
    
    spectrum = fft(dechirped, M);
    
    [~, m_hat] = max(abs(spectrum));
    m_hat = m_hat - 1;
    
    data_hat(n) = m_hat;
    
    % BER
    err = bitxor(uint32(data_tx(n)), uint32(m_hat));
    bitErrors = bitErrors + sum(bitget(err,1:SF));
    
    if err ~= 0
        blockErrors = blockErrors + 1;
    end
    
end

BER  = bitErrors / totalBits;
BLER = blockErrors / cfg.Nsym;

end