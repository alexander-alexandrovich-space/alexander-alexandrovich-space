clear; clc;close all;

%% CONFIG
cfg.SF = 12;
cfg.BW = 125e3;
cfg.Fc = 868e6;      % безопасно для Fs=1e6
cfg.Nsym = 10;
cfg.Preamble = 8;
cfg.graph = false;
cfg.FEC = false;
cfg.CR = 4;
Niter = 10;          % количество итераций на каждую точку


%cfg.DopplerHz = 5e3;


SNRdB = -40:2:10;    % диапазон SNR
cfg.ChannelType = 'AWGN';

BER_curve  = zeros(size(SNRdB));
BLER_curve = zeros(size(SNRdB));

% LOOP OVER SNR
for s = 1:length(SNRdB)
    
    BER_acc  = 0;
    BLER_acc = 0;
    
    for k = 1:Niter
        
        % --- TX ---
        [tx, data_tx, data_sym] = LoRa_tx(cfg);
        
        % --- CHANNEL ---
        %rx = channel(tx,cfg,SNRdB(s));
        rx = tx;
        % --- RX ---
        [data_decoded, BER, BLER] = LoRa_rx(rx, data_sym, cfg);
        
        BER_acc  = BER_acc  + BER;
        BLER_acc = BLER_acc + BLER;
        
    end
    
    BER_curve(s)  = BER_acc  / Niter;
    BLER_curve(s) = BLER_acc / Niter;
    
    fprintf('SNR = %d dB | BER = %.3e | BLER = %.3e\n', ...
            SNRdB(s), BER_curve(s), BLER_curve(s));
end

% PLOT
figure;
semilogy(SNRdB, BER_curve, '-o','LineWidth',2); hold on;
grid on;
xlabel('SNR (dB)');
ylabel('Error Rate');
legend('BER');
title(['LoRa BER curve | SF = ', num2str(cfg.SF)]);
xlim([-40,20]);
figure;
semilogy(SNRdB, BLER_curve, '-s','LineWidth',2);
grid on;
xlabel('SNR (dB)');
ylabel('Error Rate');
legend('BLER');
title(['LoRa BLER curve | SF = ', num2str(cfg.SF)]);
xlim([-40,20]);
%%
% close all;
% 
% cfg.SF = 12;
% cfg.BW = 125e3;
% cfg.Fc = 868e6;      % безопасно для Fs=1e6
% cfg.Nsym = 1000;
% cfg.Preamble = 8;
% cfg.graph = false;
% Niter = 10;          % количество итераций на каждую точку
% 
% cfg.graph = true;
% [tx, data_tx] = LoRa_tx(cfg);
