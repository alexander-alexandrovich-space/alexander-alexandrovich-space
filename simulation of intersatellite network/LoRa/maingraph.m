clear; clc;

%% CONFIG
cfg.SF = 12;
cfg.BW = 125e3;
cfg.Fc = 868e6;      % безопасно для Fs=1e6
cfg.Nsym = 10;
cfg.Preamble = 8;
Niter = 10;          % количество итераций на каждую точку


cfg.DopplerHz = 23e3;


SNRdB = -40:2:10;    % диапазон SNR
cfg.ChannelType = 'AWGN';

BER_curve  = zeros(size(SNRdB));
BLER_curve = zeros(size(SNRdB));

%% LOOP OVER SNR

        
        % --- TX ---
        [tx, data_tx] = LoRa_tx(cfg);
        
        % --- CHANNEL ---
        rx = channel(tx,cfg,SNRdB(1));
        
        % --- RX ---
        [data_rx, BER, BLER] = LoRa_rx(rx, data_tx, cfg);
        
        
  

%% PLOT
