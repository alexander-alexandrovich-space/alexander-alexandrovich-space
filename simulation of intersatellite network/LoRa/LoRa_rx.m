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
        
        bits = de2bi(m_hat, SF, 'left-msb');
        rx_bits((n-1)*SF + 1 : n*SF) = bits;
    end
    
    if cfg.FEC
        % Отбрасываем паддинг, чтобы длина была кратна блоку Хэмминга (4+CR)
        valid_coded_len = floor(length(rx_bits) / (4+cfg.CR)) * (4+cfg.CR);
        rx_coded_bits = rx_bits(1:valid_coded_len);
        
        rx_data_bits = LoRa_Hamming_dec(rx_coded_bits, cfg.CR);  
        % Отсекаем лишнее до оригинального размера data_tx
        rx_data_bits = rx_data_bits(1:length(data_tx));
    else
        rx_data_bits = rx_bits(1:length(data_tx));
    end
     
    % 1. Битовая ошибка (BER)
    bitErrors = sum(rx_data_bits ~= data_tx); 
    BER = bitErrors / length(data_tx);
    
    % 2. Блоковая ошибка (BLER) 
    % Группируем по 4 бита (информационное слово) для справедливого сравнения
    len = length(data_tx);
    len = len - mod(len, 4);
    rx_matrix = reshape(rx_data_bits(1:len), 4, []);
    tx_matrix = reshape(data_tx(1:len), 4, []);
    
    blockErrors = sum(any(rx_matrix ~= tx_matrix, 1));
    BLER = blockErrors / size(tx_matrix, 2);
end