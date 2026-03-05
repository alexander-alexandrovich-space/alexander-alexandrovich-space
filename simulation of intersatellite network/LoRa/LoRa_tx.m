function [tx, data] = LoRa_tx(cfg)

SF = cfg.SF;
BW = cfg.BW;

M  = 2^SF;
Ns = M;
Fs = BW;

Ts = M/BW;

t = (0:Ns-1)/Fs;

% ИДЕАЛЬНЫЙ upchirp
base = exp(1j*2*pi*(BW/(2*Ts))*t.^2);

data = randi([0 M-1],1,cfg.Nsym);

preamble = repmat(base,1,cfg.Preamble);
payload = [];

for m = data
    payload = [payload base .* exp(1j*2*pi*m/M*(0:Ns-1))];
end

tx = [preamble payload];

end