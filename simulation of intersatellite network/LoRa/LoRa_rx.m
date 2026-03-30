function [BER, BLER] = LoRa_rx(rx, data_tx, cfg)
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
    
    rx_bits = [];
    data_hat = zeros(1,cfg.Nsym);
    
    for n = 1:cfg.Nsym
        
        idx = start_payload + (n-1)*Ns;
        rxsym = rx(idx:idx+Ns-1); 
        
        dechirped = rxsym .* down;
        spectrum = fft(dechirped, M);
        
        [~, m_hat] = max(abs(spectrum));
        m_hat = m_hat - 1;
        data_hat(n) = m_hat;
        bits = de2bi(m_hat,SF,'left-msb');
        bits = bits(:).';
        rx_bits = [rx_bits bits];
        
    end
    
    
 if cfg.FEC
      rx_bits = LoRa_Hamming_dec(rx_bits,cfg.CR);  
 end
     
    % --- Вычисление ошибок ---
    totalBits = cfg.Nsym * SF;
    
    % 1. Битовая ошибка (BER)
    % Считаем количество несовпадающих бит
    bitErrors = sum(rx_bits ~= data_tx); 
    BER = bitErrors / totalBits;
    
    % 2. Блоковая (символьная) ошибка (BLER)
    % Преобразуем массивы в матрицы, где каждый столбец = 1 символ (длиной SF)
    rx_matrix = reshape(rx_bits, SF, []);
    tx_matrix = reshape(data_tx, SF, []);

    % Если хотя бы один бит в столбце (символе) не совпадает (any(..., 1)),
    % то символ считается ошибочным. Суммируем ошибочные символы.
    blockErrors = sum(any(rx_matrix ~= tx_matrix, 1));
    BLER = blockErrors / cfg.Nsym;

 

end