function [BER, BLER] = LoRa_rx(rx, data_tx, cfg, Nsym_tx)
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
    
    rx_bits = zeros(1, Nsym_tx * SF);
    
    for n = 1:Nsym_tx
        idx = start_payload + (n-1)*Ns;
        rxsym = rx(idx:idx+Ns-1); 
        
        dechirped = rxsym .* down;
        spectrum = fft(dechirped, M);
        
        [~, m_hat] = max(abs(spectrum));
        m_hat = m_hat - 1;
        
        % --- ДЕКОДИРОВАНИЕ ГРЕЯ ---
        % Восстанавливаем оригинальное бинарное значение символа
        m_bin = uint32(m_hat);
        mask = bitshift(m_bin, -1);
        while mask > 0
            m_bin = bitxor(m_bin, mask);
            mask = bitshift(mask, -1);
        end
        m_bin = double(m_bin);
        
        bits = de2bi(m_bin, SF, 'left-msb');
        rx_bits((n-1)*SF + 1 : n*SF) = bits;
    end
    
    if cfg.FEC
        % 1. ДЕИНТЕРЛИВЕР
        N_blocks = Nsym_tx / (4 + cfg.CR);
        bits_per_block = SF * (4 + cfg.CR);
        deinterleaved_bits = zeros(1, length(rx_bits));
        
        for b = 1:N_blocks
            idx = (b-1)*bits_per_block + 1 : b*bits_per_block;
            block = rx_bits(idx);
            
            cw_matrix_shifted = reshape(block, SF, 4 + cfg.CR);
            
            for i = 1:SF
                cw_matrix_shifted(i, :) = circshift(cw_matrix_shifted(i, :), [0, i-1]);
            end
            
            cw_matrix_orig = cw_matrix_shifted.';
            deinterleaved_bits(idx) = cw_matrix_orig(:).';
        end
        
        % 2. ДЕКОДИРОВАНИЕ
        rx_data_bits = LoRa_Hamming_dec(deinterleaved_bits, cfg.CR);
        
        % 3. ОТСЕЧЕНИЕ ПАДДИНГА
        rx_data_bits = rx_data_bits(1:length(data_tx));
    else
        rx_data_bits = rx_bits(1:length(data_tx));
    end
     
    % --- Вычисление ошибок ---
    bitErrors = sum(rx_data_bits ~= data_tx); 
    BER = bitErrors / length(data_tx);
    
    len = length(data_tx);
    len = len - mod(len, 4);
    rx_matrix = reshape(rx_data_bits(1:len), 4, []);
    tx_matrix = reshape(data_tx(1:len), 4, []);
    
    blockErrors = sum(any(rx_matrix ~= tx_matrix, 1));
    BLER = blockErrors / size(tx_matrix, 2);
end