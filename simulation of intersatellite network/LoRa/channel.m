function rx = channel(tx, cfg, SNRdB)

Fs = cfg.BW;              % частота дискретизации

%% =========================
% 1) Допплеровский сдвиг
%% =========================
if isfield(cfg,'DopplerHz') && cfg.DopplerHz ~= 0
    t = (0:length(tx)-1)/Fs;
    tx = tx .* exp(1j*2*pi*cfg.DopplerHz*t);
end

%% =========================
% 2) Многолучевой канал
%% =========================
switch upper(cfg.ChannelType)

    case 'AWGN'
        sig_chan = tx;

    case 'EPA'
        % 3GPP LTE Extended Pedestrian A
        chan = comm.RayleighChannel( ...
            'SampleRate', Fs, ...
            'PathDelays', [0 30 70 90 110 190]*1e-9, ...
            'AveragePathGains', [0 -1 -2 -3 -8 -17], ...
            'MaximumDopplerShift', 0);

        sig_chan = chan(tx.').';
        
    case 'EVA'
        % Extended Vehicular A
        chan = comm.RayleighChannel( ...
            'SampleRate', Fs, ...
            'PathDelays', [0 30 150 310 370 710 1090 1730 2510]*1e-9, ...
            'AveragePathGains', [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9], ...
            'MaximumDopplerShift', 0);

        sig_chan = chan(tx.').';

    case 'ETU'
        % Extended Typical Urban
        chan = comm.RayleighChannel( ...
            'SampleRate', Fs, ...
            'PathDelays', [0 50 120 200 230 500 1600 2300 5000]*1e-9, ...
            'AveragePathGains', [-1 -1 -1 0 -3 -8 -17.2 -20.8 -23.9], ...
            'MaximumDopplerShift', 0);

        sig_chan = chan(tx.').';

    otherwise
        error('Unknown channel type');
end

%% =========================
% 3) Добавление AWGN
%% =========================
    Es = mean(abs(tx).^2);
    Eb = Es / cfg.SF;

    SNR_linear = 10^(SNRdB/10);
    N0 = Es / SNR_linear; % итерируюсь по SNR
    noise_var = N0 / 2;
    
    noise = sqrt(noise_var) * (randn(size(tx)) + 1j*randn(size(tx)));
    rx = sig_chan + noise;

end