function [data_hat, BER, BLER] = LoRa_rx(rx, data_tx, cfg)

SF = cfg.SF;
BW = cfg.BW;

M  = 2^SF;
Ns = M;
Fs = BW;

Ts = M/BW;

t = (0:Ns-1)/Fs;

base = exp(1j*2*pi*(BW/(2*Ts))*t.^2);
down = conj(base);

start = cfg.Preamble*Ns + 1;

data_hat = zeros(1,cfg.Nsym);

bitErrors = 0;
blockErrors = 0;
totalBits = cfg.Nsym*SF;

for n = 1:cfg.Nsym
    
    idx = start + (n-1)*Ns;
    r = rx(idx:idx+Ns-1);
    
    d = r .* down;
    
    S = fft(d,M);
    [~,k] = max(abs(S));
    
    m_hat = k-1;
    
    data_hat(n) = m_hat;
    
    err = bitxor(uint32(data_tx(n)),uint32(m_hat));
    bitErrors = bitErrors + sum(bitget(err,1:SF));
    
    if err~=0
        blockErrors = blockErrors+1;
    end
end

BER = bitErrors/totalBits;
BLER = blockErrors/cfg.Nsym;

end