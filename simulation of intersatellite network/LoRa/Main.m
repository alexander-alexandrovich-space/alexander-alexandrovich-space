clear; clc;close all;

%% CONFIG
cfg.SF = 12;
cfg.BW = 125e3;
cfg.Fc = 868e6;      % безопасно для Fs=1e6
cfg.Nsym = 100;
cfg.Preamble = 8;
cfg.graph = false;
cfg.FEC = true;
Niter = 10;          % количество итераций на каждую точку

rng(42);

cfg.DopplerHz = 5;

CR_modes = [0 4];

legend_str = cell(1, length(CR_modes));
figure(1); hold on; grid on; title('LoRa BER curve'); xlabel('SNR (dB)'); ylabel('BER'); set(gca, 'YScale', 'log');
figure(2); hold on; grid on; title('LoRa BLER curve'); xlabel('SNR (dB)'); ylabel('BLER'); set(gca, 'YScale', 'log');

SNRdB = -20:1:15;    % диапазон SNR
cfg.ChannelType = 'EPA';





% LOOP OVER SNR

for CR_id = 1:length(CR_modes)
cfg.CR = CR_modes(CR_id);
    if cfg.CR == 0
        cfg.FEC = false;
        legend_str{CR_id} = 'No FEC';
    else
        cfg.FEC = true;
        legend_str{CR_id} = ['FEC, CR = ', num2str(cfg.CR)];
    end

    BER_curve  = zeros(size(SNRdB));
    BLER_curve = zeros(size(SNRdB));

    for s = 1:length(SNRdB)
        BER_acc = 0;
        BLER_acc = 0;
      
        for k = 1:Niter
            
            % --- TX ---
            [tx, data_tx, Nsym_tx] = LoRa_tx(cfg);
            
            % --- CHANNEL ---
            rx = channel(tx,cfg,SNRdB(s));
            %rx = tx;
            % --- RX ---
            [BER, BLER] = LoRa_rx(rx, data_tx, cfg, Nsym_tx);
            
            BER_acc=BER_acc+BER;
            BLER_acc=BLER_acc+BLER;
        
            
        end
        
       BER_curve(s)  = BER_acc  / Niter;
       BLER_curve(s) = BLER_acc / Niter;
    
        fprintf('SNR = %d dB | BER = %.3e | BLER = %.3e\n', ...
                SNRdB(s), BER_curve(s), BLER_curve(s));

    end
    figure(1); plot(SNRdB, BER_curve, 'Marker', 'o', 'LineWidth', 2);
    figure(2); plot(SNRdB, BLER_curve, 'Marker', 's', 'LineWidth', 2);
end

figure(1); legend(legend_str, 'Location', 'southwest'); ylim([1e-4, 1]);
figure(2); legend(legend_str, 'Location', 'southwest'); ylim([1e-4, 1]);

% PLOT
% figure;
% semilogy(SNRdB, BER_curve, '-o','LineWidth',2); hold on;
% grid on;
% xlabel('SNR (dB)');
% ylabel('Error Rate');
% legend('BER');
% title(['LoRa BER curve | SF = ', num2str(cfg.SF)]);
% xlim([-40,20]);
% figure;
% semilogy(SNRdB, BLER_curve, '-s','LineWidth',2);
% grid on;
% xlabel('SNR (dB)');
% ylabel('Error Rate');
% legend('BLER');
% title(['LoRa BLER curve | SF = ', num2str(cfg.SF)]);
% xlim([-40,20]);
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
