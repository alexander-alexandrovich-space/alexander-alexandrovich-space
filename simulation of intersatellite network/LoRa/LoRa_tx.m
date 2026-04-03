function [tx, data_tx, Nsym_tx] = LoRa_tx(cfg)
    SF = cfg.SF;
    BW = cfg.BW;
    M  = 2^SF;
    Ns = M;
    Ts = M/BW;
    Fs = BW;
    t = (0:Ns-1)/Fs;
    base_chirp = exp(1j*2*pi*(-BW/2*t + BW/(2*Ts)*t.^2));
    
    if cfg.FEC
        % Целевое количество данных (чтобы получить примерно Nsym символов)
        target_data_bits = floor(cfg.Nsym * 4 * SF / (4 + cfg.CR));
        data_tx = randi([0 1], 1, target_data_bits);
        
        % 1. ПАДДИНГ (ДО кодирования)
        block_data_size = 4 * SF;
        pad_len = mod(block_data_size - mod(length(data_tx), block_data_size), block_data_size);
        if pad_len == block_data_size, pad_len = 0; end
        
        data_tx_padded = [data_tx, zeros(1, pad_len)];
        
        % 2. КОДИРОВАНИЕ
        coded_bits = LoRa_Hamming_enc(data_tx_padded, cfg.CR);
        
        % 3. ИНТЕРЛИВЕР (Диагональное перемежение)
        N_blocks = length(coded_bits) / (SF * (4 + cfg.CR));
        interleaved_bits = zeros(size(coded_bits));
        bits_per_block = SF * (4 + cfg.CR);
        
        for b = 1:N_blocks
            idx = (b-1)*bits_per_block + 1 : b*bits_per_block;
            block = coded_bits(idx);
            
            cw_matrix = reshape(block, 4 + cfg.CR, SF).';
            
            for i = 1:SF
                cw_matrix(i, :) = circshift(cw_matrix(i, :), [0, -(i-1)]);
            end
            
            interleaved_bits(idx) = cw_matrix(:).';
        end
        
        symbols_bits = interleaved_bits;
    else
        % Без FEC просто делаем паддинг до кратности SF
        target_data_bits = cfg.Nsym * SF;
        data_tx = randi([0 1], 1, target_data_bits);
        
        pad_len = mod(SF - mod(length(data_tx), SF), SF);
        if pad_len == SF, pad_len = 0; end
        
        symbols_bits = [data_tx, zeros(1, pad_len)];
    end
    
    % Модуляция
    Nsym_tx = length(symbols_bits) / SF;
    data_sym = bi2de(reshape(symbols_bits, SF, []).', 'left-msb');
    
    % --- 4. КОДИРОВАНИЕ ГРЕЯ ---
    % Преобразуем бинарное значение символа в код Грея
    data_sym_int = uint32(data_sym);
    data_sym_gray = double(bitxor(data_sym_int, bitshift(data_sym_int, -1)));
    
    tx = [];
    for m = data_sym_gray.'
        symbol = circshift(base_chirp, -m);
        tx = [tx symbol(:).'];
    end
    
    preamble = repmat(base_chirp, 1, cfg.Preamble);
    tx = [preamble tx];
    
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