function [tx, data] = LoRa_tx(cfg)

SF = cfg.SF;
BW = cfg.BW;

M  = 2^SF;
Ns = M;              % КРИТИЧНО

Ts = M/BW;
Fs = BW;

t = (0:Ns-1)/Fs;

% Правильный upchirp (стандартная форма)
base_chirp = exp(1j*2*pi*(BW/(2*Ts)*t.^2));

data = randi([0 M-1],1,cfg.Nsym);

tx = [];

for m = data
    % Циклический сдвиг через умножение
    symbol = base_chirp .* exp(1j*2*pi*m/M*(0:Ns-1));
    tx = [tx symbol];
end

% Preamble
preamble = repmat(base_chirp,1,cfg.Preamble);
tx = [preamble tx];

end