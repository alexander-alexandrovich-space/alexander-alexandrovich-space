function [tx, data_tx, Nsym_tx] = LoRa_tx(cfg)
    SF = cfg.SF;
    BW = cfg.BW;
    M  = 2^SF;
    Ns = M;
    Ts = M/BW;
    Fs = BW;
    t = (0:Ns-1)/Fs;
    base_chirp = exp(1j*2*pi*(-BW/2*t + BW/(2*Ts)*t.^2));
    
    target_bits = cfg.Nsym * SF;
    
    if cfg.FEC
        % Рассчитываем количество бит так, чтобы на выходе получить примерно Nsym символов
        num_data_bits = floor(target_bits * 4 / (4 + cfg.CR));
        num_data_bits = num_data_bits - mod(num_data_bits, 4); % Кратность 4
        data_tx = randi([0 1], 1, num_data_bits);
        coded_bits = LoRa_Hamming_enc(data_tx, cfg.CR);
    else
        data_tx = randi([0 1], 1, target_bits);
        coded_bits = data_tx;
    end
    
    % Паддинг до кратности SF
    pad_len = mod(SF - mod(length(coded_bits), SF), SF);
    if pad_len == SF, pad_len = 0; end
    coded_bits_padded = [coded_bits zeros(1, pad_len)];
    
    Nsym_tx = length(coded_bits_padded) / SF; % Сообщаем приемнику, сколько символов ждать
    data_sym = bi2de(reshape(coded_bits_padded, SF, []).', 'left-msb');
    
    tx = [];
    for m = data_sym.'
        symbol = circshift(base_chirp, -m);
        tx = [tx symbol(:).'];
    end
    
    preamble = repmat(base_chirp, 1, cfg.Preamble);
    tx = [preamble tx];
end