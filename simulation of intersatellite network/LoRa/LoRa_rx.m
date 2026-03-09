function [ all_rx_bits, BER, BLER] = LoRa_rx(rx, data_tx, cfg)

SF = cfg.SF;
BW = cfg.BW;

M  = 2^SF;
Ns = M;

Ts = M/BW;
Fs = BW;

t = (0:Ns-1)/Fs;

base_chirp = exp(1j*2*pi*(-BW/2*t + BW/(2*Ts)*t.^2));
down = conj(base_chirp);

start_payload = cfg.Preamble*Ns + 1;

bitErrors = 0;
totalBits = cfg.Nsym * SF;
blockErrors = 0;

all_rx_bits = zeros(1, cfg.Nsym*SF);

for n = 1:cfg.Nsym
    
    idx = start_payload + (n-1)*Ns;
    rxsym = rx(idx:idx+Ns-1);
    
    dechirped = rxsym .* down;
    
    spectrum = fft(dechirped, M);
    
    [~, m_hat] = max(abs(spectrum));
    m_hat = m_hat - 1;
    
    bits = de2bi(m_hat,SF,'left-msb');

%    
    
end

BER  = bitErrors / totalBits;
BLER = blockErrors / cfg.Nsym;

end