clear;
clc;


BW = 125e3;
SF = 12;
rng(42);
M  = 2^SF;
Ns = M;

Ts = M/BW;
Fs = BW;

t = (0:Ns-1)/Fs;

base_chirp = exp(1j*2*pi*(-BW/2*t + BW/(2*Ts)*t.^2));
cfg.Nsym = 2;
cfg.Preamble = 8;
data_tx = randi([0 1],1,cfg.Nsym*SF);


data_sym = bi2de(reshape(data_tx,SF,[]).','left-msb');


tx = [];

for m = data_sym
symbol = circshift(base_chirp,m);
symbol = symbol(:).';
    tx = [tx symbol];
end

preamble = repmat(base_chirp,1,cfg.Preamble);
tx = [preamble tx];

rx=tx;

down = conj(base_chirp);
start_payload = cfg.Preamble*Ns + 1;
%z  start_payload = 1;

rx_bits = [];
data_hat = zeros(1,cfg.Nsym);
for n = 1:cfg.Nsym
    
    idx = start_payload + (n-1)*Ns;
    rxsym = rx(idx:idx+Ns-1);
    
    dechirped = rxsym .* down;

    spectrum = fft(dechirped, M);
    

    [~, m_hat] = max(abs(spectrum));
    m_hat = m_hat - 1;
    data_hat(n)=m_hat;
    bits = de2bi(m_hat,SF,'left-msb');
    bits = bits(:).';
   % all_rx_bits((n-1)*SF+1:n*SF) = bits;
  rx_bits = [rx_bits bits];
    
end

    huy = reshape(de2bi(data_hat,SF,'left-msb').',1,[]);


if isequal(rx_bits,data_tx)
    disp(100000)
end