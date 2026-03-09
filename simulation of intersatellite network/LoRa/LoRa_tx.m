function [tx, data_tx, data_sym] = LoRa_tx(cfg)

SF = cfg.SF;
BW = cfg.BW;

M  = 2^SF;
Ns = M;

Ts = M/BW;
Fs = BW;

t = (0:Ns-1)/Fs;

base_chirp = exp(1j*2*pi*(-BW/2*t + BW/(2*Ts)*t.^2));

data_tx = randi([0 1],1,cfg.Nsym*SF);

if cfg.FEC
    data_bits = LoRa_Hamming_enc(data_tx,cfg.CR);
else
    data_bits = data_tx;
end

data_sym = bi2de(reshape(data_bits,SF,[]).','left-msb');

tx = [];

for m = data_sym
    symbol = base_chirp .* exp(1j*2*pi*m/M*(0:Ns-1));
    symbol = symbol(:).';
    tx = [tx symbol];
end

preamble = repmat(base_chirp,1,cfg.Preamble);
tx = [preamble tx];

%%
if cfg.graph
    Fc = 868e6;
    t = (0:length(tx)-1)/Fs;
    rf = real(tx .* exp(1j*2*pi*Fc*t));
    
    figure;
    plot(t(1:1000), rf(1:1000));
    grid on;
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('LoRa RF Signal (Time Domain)');
    
    N = length(rf);
    f = (-N/2:N/2-1)*(Fs/N);
    
    RF_spectrum = fftshift(abs(fft(rf)));
    
    figure;
    plot(Fc + f, 20*log10(RF_spectrum));
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('LoRa Spectrum around Carrier');
    
    window = 1024;
    noverlap = 900;
    nfft = 2048;
    
    figure;
    spectrogram(tx, window, noverlap, nfft, cfg.BW, 'yaxis');
    colormap jet;
    colorbar;
    title('LoRa Baseband Chirp');
end

end