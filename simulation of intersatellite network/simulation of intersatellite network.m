clc; clear; close all;

%% === Параметры сети LoRa ===
%numNodes = 100;               % Количество узлов в сети
txPower = 20;                % Передаточная мощность (dBm)
sensitivity = -137;          % Чувствительность приемника (dBm)
freq = 868e6;                % Частота передачи (868 MHz)

rng(42);
%% === Физические параметры и настройки орбитальной модели ===
mu = 3.986004418e14;    % Гравитационный параметр Земли, м^3/с^2
R_earth = 6400000;      % Радиус Земли, м

% Задаем группы орбит: выбираем несколько орбитальных высот
% (r - расстояние от центра Земли)
r_orbits = [7014000, 7628000];  % примерные орбитальные радиусы (м)

% Задаем время симуляции (секундах)
t = 1000;  % например, через 1000 секунд после начального момента

%% === Параметры размещения спутников ===
numPerOrbit = 20;  % число спутников на каждой орбите

% Для каждой группы орбит будем генерировать спутники для 3 типов:
% 1. Экваториальные: i = 0°, RAAN = 0.
% 2. Полярные: i = 90°, RAAN = 0.
% 3. Наклонные: i = 30° (можно добавить и i = 60°), RAAN = 0.
%
% Для каждой орбиты задаём начальный сдвиг (phaseShift) по фазе, чтобы обеспечить равномерное распределение.

% Зададим для каждой группы массив структур:
satellites = [];  % итоговый массив спутников
%--- Группа 1: Экваториальные орбиты (i = 0°) ---
inclination_eq = 0;
RAAN_eq = deg2rad(0*360);
phaseShift_eq = 0;
for idx = 1:length(r_orbits)
    r = r_orbits(idx);
    h = r - R_earth;
    T = 2*pi*sqrt(r^3/mu);
    omega = 2*pi/T;
    initPhases_deg = (0:(numPerOrbit-1))*(360/numPerOrbit) + phaseShift_eq;
    initPhases = deg2rad(initPhases_deg);
    for k = 1:numPerOrbit
        sat.id = (idx-1)*numPerOrbit + k;
        sat.r = r;
        sat.h = h;
        sat.inclination = inclination_eq;
        sat.RAAN = RAAN_eq;
        sat.initPhase = initPhases(k);
        sat.omega = omega;
        sat.trueAnomaly = mod(sat.initPhase + omega*t, 2*pi);
        sat.lat = 0;
        sat.lon = rad2deg(sat.trueAnomaly + RAAN_eq);
        sat.x = r * (cos(RAAN_eq)*cos(sat.trueAnomaly) - sin(RAAN_eq)*sin(sat.trueAnomaly));
        sat.y = r * (sin(RAAN_eq)*cos(sat.trueAnomaly) + cos(RAAN_eq)*sin(sat.trueAnomaly));
        sat.z = 0;
        sat.group = 'Equatorial';
        satellites = [satellites; sat];
    end
end

%--- Группа 2: Полярные орбиты (i = 90°) ---
inclination_pol = deg2rad(90);
RAAN_pol = deg2rad(0);
phaseShift_pol = 0;
for idx = 1:length(r_orbits)
    r = r_orbits(idx);
    h = r - R_earth;
    T = 2*pi*sqrt(r^3/mu);
    omega = 2*pi/T;
    initPhases_deg = (0:(numPerOrbit-1))*(360/numPerOrbit) + phaseShift_pol;
    initPhases = deg2rad(initPhases_deg);
    for k = 1:numPerOrbit
        sat.id = (length(satellites)+1);
        sat.r = r;
        sat.h = h;
        sat.inclination = inclination_pol;
        sat.RAAN = RAAN_pol;
        sat.initPhase = initPhases(k);
        sat.omega = omega;
        sat.trueAnomaly = mod(sat.initPhase + omega*t, 2*pi);
        sat.x = r * (cos(RAAN_pol)*cos(sat.trueAnomaly));
        sat.y = r * (sin(RAAN_pol)*cos(sat.trueAnomaly));
        sat.z = r * sin(sat.trueAnomaly);
        sat.lat = rad2deg(asin(sat.z/r));
        sat.lon = rad2deg(atan2(sat.y, sat.x));
        sat.group = 'Polar';
        satellites = [satellites; sat];
    end
end

% --- Группа 3: Наклонные орбиты (i = 30°) ---


for m=1:12%m=1:round(2*pi*6560000/1500000)
    inclination_inc = deg2rad(60);
    RAAN_inc = deg2rad(30*m); %RAAN_inc = deg2rad(rad2deg(acos((-(1500)^2+2*(6560)^2)/(2*(6560)^2)))*m);
    phaseShift_inc = 0;
    for idx = 1:length(r_orbits)
        r = r_orbits(idx);
        h = r - R_earth;
        T = 2*pi*sqrt(r^3/mu);
        omega = 2*pi/T;
        initPhases_deg = (0:(numPerOrbit-1))*(360/numPerOrbit) + phaseShift_inc;
        initPhases = deg2rad(initPhases_deg);
        for k = 1:numPerOrbit
            sat.id = (length(satellites)+1);
            sat.r = r;
            sat.h = h;
            sat.inclination = inclination_inc;
            sat.RAAN = RAAN_inc;
            sat.initPhase = initPhases(k);
            sat.omega = omega;
            sat.trueAnomaly = mod(sat.initPhase + omega*t, 2*pi);
            sat.x = r * (cos(RAAN_inc)*cos(sat.trueAnomaly) - sin(RAAN_inc)*sin(sat.trueAnomaly)*cos(inclination_inc));
            sat.y = r * (sin(RAAN_inc)*cos(sat.trueAnomaly) + cos(RAAN_inc)*sin(sat.trueAnomaly)*cos(inclination_inc));
            sat.z = r * sin(sat.trueAnomaly)*sin(inclination_inc);
            sat.lat = rad2deg(asin(sat.z/r));
            sat.lon = rad2deg(atan2(sat.y, sat.x));
            sat.group = 'Inclined';
            satellites = [satellites; sat];
        end
    end
end
%% ===============Задание количества узлов в сети=============
numNodes= length(satellites);

%% === Генерация расположения узлов (экваториальный пояс ±30°) (cтарая случайная)===
% r_min = 6560000;                  % Минимальный орбитальный радиус (м)
% r_max = 8400000;                  % Максимальный орбитальный радиус (м)
% 
% theta = rand(1, numNodes) * 2 * pi; % Равномерное распределение по долготе
% 
% % Концентрация узлов в экваториальном поясе (±30° широты)
% delta_deg = 30;                       % Ширина пояса в градусах
% delta_rad = deg2rad(delta_deg);        % Преобразование в радианы
% U_min = cos(pi/2 + delta_rad);        % cos(120°) = -0.5
% U_max = cos(pi/2 - delta_rad);        % cos(60°) = 0.5
% U = U_min + (U_max - U_min)*rand(1, numNodes); % Рандомные значения U
% phi = acos(U);                        % Угол phi для концентрации у экватора
% 
% 
% r_nodes = r_min + (r_max - r_min) * rand(1, numNodes);
% 
% % Преобразование в декартовы координаты
% x_nodes = r_nodes .* sin(phi) .* cos(theta);
% y_nodes = r_nodes .* sin(phi) .* sin(theta);
% z_nodes = r_nodes .* cos(phi);

%% ============= Присвоение позиций================

%nodes = [x_nodes' y_nodes' z_nodes']; %случайный способ

nodes = zeros(length(satellites), 3);
for i = 1:length(satellites)
    nodes(i, :) = [satellites(i).x, satellites(i).y, satellites(i).z];
end
% детерминированный способ

%% === Вычисление дальности связи и построение связей ===

PL_margin = 147.55;          % Константа FSPL
c=3e8;

meshRange = 10^((txPower - sensitivity - 20*log10(freq) + PL_margin)/20);
connections = zeros(numNodes);

realmeshranges=zeros(1,numNodes);

W = inf(numNodes);

noiseOffset=zeros(1,numNodes);

realrange=zeros(1,numNodes);

for i=1:numNodes
    for j = i+1:numNodes
        noiseOffset(i,j)=rand*5;
        noiseOffset(j,i)=noiseOffset(i,j);
    end
end

for i = 1:numNodes
    for j = i+1:numNodes
        dist = norm(nodes(i,:) - nodes(j,:));   
        realrange(i)=10^((txPower - sensitivity - 20*log10(freq) + PL_margin-noiseOffset(i,j))/20);
        if dist < realrange(i)
            connections(i,j)=1;
            connections(j,i)=1;
            totalDelay = dist / c;  % сек.
            W(i, j) = totalDelay;
            W(j, i) = totalDelay;
        end
    end
end

rr=mean(realrange);
       



fprintf('\nРасчетная дальность %.2f км:\n', meshRange/1000);
fprintf('\nРеальная дальность (+ шум) %.2f км:\n', rr/1000);

%% изолированные узлы и кластеры
isolatedNodes = [];
for i = 1:numNodes
    isIsolated = true;
    for j = 1:numNodes
        if i == j, continue; end % Пропускаем сравнение с собой
        if connections(i,j)==1
            isIsolated = false;
            break; % Прерываем цикл при первой найденной связи
        end
    end
    if isIsolated
        isolatedNodes(end+1) = i;
    end
end

fprintf('\nИзолированные узлы (расстояние > %.2f км):\n', rr/1000);
disp(isolatedNodes);




G = graph(connections);
components = conncomp(G);  % components(i) содержит номер компоненты для узла i
numComponents = max(components);


if numComponents == 1
    fprintf('Все спутники связаны в единую сеть.\n');
else
    fprintf('Обнаружено %d изолированных кластеров:\n', numComponents);
    for comp = 1:numComponents
        clusterNodes = find(components == comp);
        fprintf('Кластер %d: %s\n', comp, mat2str(clusterNodes'));
    end
end


%% === Поиск пары узлов с максимальным физическим расстоянием, для которых существует маршрут по сети ===
D_phys = pdist2(nodes, nodes);  % матрица физических расстояний
% Построим взвешенный граф для маршрутизации
% Создаем граф по матрице connections (только прямые связи)
G_conn = graph(connections);

maxD = -Inf;
best_i = 0;
best_j = 0;
for i = 1:numNodes
    for j = i+1:numNodes
        % Проверяем, что узлы i и j находятся в одной связной компоненте (существует маршрут)
        if ~isempty(shortestpath(G_conn, i, j))
            if D_phys(i,j) > maxD
                maxD = D_phys(i,j);
                best_i = i;
                best_j = j;
            end
        end
    end
end

%
best_i=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best_j=500;
%

if best_i == 0 || best_j == 0
    fprintf('Не найдено пары узлов с существующим маршрутом.\n');
else
    fprintf('\nПара узлов с максимальным физическим расстоянием и существующим маршрутом: %d и %d (%.2f км)\n', best_i, best_j, D_phys(best_i,best_j)/1000);
end





G_weight = graph(W);  % взвешенный граф, где ребра есть только для действующих соединений
[path, totalDelay] = shortestpath(G_weight, best_i, best_j);

if isempty(path)
    fprintf('Между узлами %d и %d нет связного пути.\n', best_i, best_j);
else
    fprintf('Наиболее выгодный маршрут между узлами %d и %d:\n', best_i, best_j);
    disp(path);
    fprintf('Общая накопленная задержка по лучшему маршруту: %.3f секунд\n', totalDelay);
end

%% === Поиск наихудшего маршрута (максимальная задержка) по задержке ===
% Построим матрицу кратчайших задержек между всеми парами узлов
allDelays = distances(G_weight);

% Извлекаем конечные значения задержек (исключая Inf и 0 для диагонали)
finiteDelays = allDelays(~isinf(allDelays) & allDelays > 0);

if isempty(finiteDelays)
    fprintf('Нет связных путей в сети.\n');
else
    worstDelay = max(finiteDelays);
    fprintf('Наихудшая (максимальная) накопленная задержка по сети: %.3f секунд\n', worstDelay);
    
    % Найдем первую пару узлов, для которых задержка равна worstDelay
    [worst_i, worst_j] = find(allDelays == worstDelay, 1);  
    % Получаем маршрут между ними с накоплением задержки
    worstPath = shortestpath(G_weight, worst_i, worst_j);
    disp('Наихудший маршрут (узлы):');
    disp(worstPath);
end

%% === Моделирование скорости передачи и оценка ошибок для наилучшего маршрута ===
% Предположим, что bestPath и totalDelay получены ранее из поиска маршрута по графу G_weight.
    % Параметры пакета и передачи
     % Параметры пакета
  % Симуляция передачи данных по лучшему маршруту с подсчетом ошибок для каждого узла

 bestPath=shortestpath(G_conn, best_i, best_j);
packetSize = 240; % размер пакета в битах
hopErrors = zeros(1, length(bestPath)-1);
hopBERs = zeros(1, length(bestPath)-1);
hopSNRs = zeros(1, length(bestPath)-1);

for k = 1:length(bestPath)-1
    i1 = bestPath(k);
    i2 = bestPath(k+1);
    dist = norm(nodes(i1,:) - nodes(i2,:));
    effectivePathLoss = 20*log10(dist) + 20*log10(freq) - PL_margin + noiseOffset(i1, i2);
    rxPower = txPower - effectivePathLoss;
    SNR_hop = rxPower - sensitivity;
    hopSNRs(k) = SNR_hop;
   % Расчет BER на основе SNR
    if SNR_hop < -20
        BER_hop = 1; % При SNR ниже -20 дБ считаем, что все биты ошибочны
    else
        % Примерная оценка BER на основе SNR
        BER_hop = 0.5 * erfc(sqrt(10^(SNR_hop / 10)));
    end

    hopBERs(k) = BER_hop;
    errors = binornd(packetSize, BER_hop);
    hopErrors(k) = errors;
end



totalErrors = sum(hopErrors);
fprintf('\nСимуляция передачи по лучшему маршруту:\n');
for k = 1:length(bestPath)-1
    fprintf('Хоп от узла %d к узлу %d: SNR = %.2f дБ, BER = %.4e, ошибок = %d из %d бит\n', ...
        bestPath(k), bestPath(k+1), hopSNRs(k), hopBERs(k), hopErrors(k), packetSize);
end
fprintf('\nОбщее число ошибок по маршруту: %d из %d бит (накоплено по хопам)\n', ...
    totalErrors, packetSize * (length(bestPath)-1));



% --------------------непосредственная модуляция через chirp----------------



% Параметры LoRa CSS
SF = 12;                % Spreading Factor (количество бит на символ)
M = 2^SF;              % Размер созвездия
BW = 125e3;            % Полоса частот (Гц)
symDuration = 2^SF / BW;  % Длительность символа (сек)
Fs = 125e1;              % Частота дискретизации (Гц)
t_sym = linspace(0, symDuration, round(symDuration*Fs));
% Генерация базового чирпа (up-chirp) от 0 до BW
baseChirp = chirp(t_sym, 0, symDuration, BW);

% Количество символов в пакете
packetSymbols = packetSize / SF; 

% Для каждого перехода по маршруту bestPath будем моделировать передачу пакета
Niter = 1;  % число итераций для усреднения
hopBER = zeros(1, length(bestPath)-1);

for k = 1:length(bestPath)-1
    i1 = bestPath(k);
    i2 = bestPath(k+1);

    % Используем ранее рассчитанный SNR для данного хопа
    % Если SNR < -20, ставим его равным -20 для модели (порог LoRa)
    SNR_hop = hopSNRs(k);
    % Перевод SNR (дБ) в линейное отношение
    snr_linear = 10^(SNR_hop/10);
    % Энергия базового чирпа
    Es = mean(abs(baseChirp).^2);
    % Для передачи символа мощность сигнала равна Es; дисперсия шума определяется как:
    noiseVar = Es / snr_linear;

    totalBitErrors = 0;
    totalBits = packetSize * Niter;
    for iter = 1:Niter
        % Генерация случайных битов
        data_bits = randi([0 1], packetSize, 1);
        % Группируем биты по SF и преобразуем в символы
        dataSymbols = bi2de(reshape(data_bits, SF, []).', 'left-msb');

        txSignal = [];
        for s = 1:length(dataSymbols)
            % Для LoRa CSS символ модулируется циклическим сдвигом базового чирпа
            shift = round(dataSymbols(s) * length(baseChirp) / M);
            modChirp = circshift(baseChirp, shift);
            txSignal = [txSignal, modChirp];  %#ok<AGROW>
        end
        txSignal = txSignal(:);

        % Генерация шума AWGN с дисперсией noiseVar
        noise = sqrt(noiseVar/2) * (randn(size(txSignal)) + 1j*randn(size(txSignal)));
        rxSignal = txSignal + noise;
        
        % Демодуляция: разбиваем сигнал на символы и для каждого выполняем корреляцию
        rxSymbols = zeros(length(dataSymbols), 1);
        L_chirp = length(baseChirp);
        
        % Преобразуем baseChirp в вектор-столбец для корректных операций
        baseChirp = baseChirp(:); 
        
        for s = 1:length(dataSymbols)
            segment = rxSignal((s-1)*L_chirp+1 : s*L_chirp);
            correlations = zeros(1, M);
            
            % Преобразуем segment в столбец
            segment = segment(:); 
            
            for m = 0:M-1
                % Генерируем опорный чирп как вектор-столбец
                refChirp = circshift(baseChirp, m);
                % Поэлементное умножение векторов-столбцов
                corrValue = abs(sum(segment .* conj(refChirp))); 
                correlations(m+1) = corrValue;
            end
            
            [~, detected] = max(correlations);
            rxSymbols(s) = detected - 1;
        end

        % Преобразуем полученные символы обратно в биты
        rxBits = de2bi(rxSymbols, SF, 'left-msb').';
        rxBits = rxBits(:);
        % Считаем ошибки
        totalBitErrors = totalBitErrors + sum(data_bits ~= rxBits(1:packetSize));
    end
    hopBER(k) = totalBitErrors / (packetSize * Niter);
    fprintf('\n Переход %d -> %d: SNR = %.2f dB, BER = %.4e\n', bestPath(k), bestPath(k+1), SNR_hop, hopBER(k));
end

overallBER = prod(1 - hopBER);
fprintf('\n Общая вероятность безошибочной передачи по маршруту: %.4e\n', overallBER);
fprintf('Общая вероятность ошибки (PER) по маршруту: %.4e\n', 1 - overallBER);


% --------------------калькулятор по среднему снр--------------------------

    dataRate=290;
   % Расчет времени передачи без FEC (только для пакета)
    txTime_noFEC = packetSize / dataRate;  % в секундах
    
    % Предположим, что для каждого ребра маршрута мы уже рассчитали SNR (в предыдущей части)
    % Здесь для простоты расчитаем усредненное значение SNR по маршруту.
 
   
    
    % Усреднённое SNR по маршруту
    avgSNR = mean( hopSNRs);
    % Без FEC: при среднем SNR, модель exp(-max(SNR,0)/10) дает:
    BER_noFEC = exp(-max(avgSNR,0)/10);  % max(-5.43, 0) = 0, поэтому BER_noFEC = exp(0)=1
    PER_noFEC = 1 - (1 - BER_noFEC)^packetSize;  % Практически PER = 1
    txTime_noFEC = packetSize / dataRate;  % Время передачи полезного пакета без FEC
    overallTime_noFEC = length(bestPath)*(txTime_noFEC+0.2) + totalDelay;  % Общая задержка маршрута + время передачи
    effectiveRate_noFEC = packetSize / overallTime_noFEC;  % Эффективная скорость передачи полезных бит
    
    % С применением FEC:
    % Предположим, что алгоритм FEC снижает эффективное BER до фиксированного уровня:
    effectiveBER_FEC = BER_noFEC/10000;  % Допустим, остаточный BER = 10^-5 с хорошо настроенным FEC
    PER_FEC = 1 - (1 - effectiveBER_FEC)^packetSize;  
    codeRate = 0.5;
    txTime_FEC = (packetSize / codeRate) / dataRate;  % Время передачи с учетом избыточности
    overallTime_FEC = length(bestPath)*(txTime_FEC+0.2) + totalDelay;  % Общая задержка маршрута + время передачи
    effectiveRate_FEC = packetSize / overallTime_FEC;  % Эффективная скорость полезной передачи
    
    fprintf('\n--- Без FEC ---\n');
    fprintf('Среднее SNR: %.2f дБ\n', avgSNR);
    fprintf('BER (без FEC): %.4e\n', BER_noFEC);
    fprintf('PER (без FEC): %.4f\n', PER_noFEC);
    fprintf('Время передачи пакета (без FEC): %.3f сек\n', overallTime_noFEC);
    fprintf('Эффективная скорость передачи (без FEC): %.2f бит/сек\n', effectiveRate_noFEC);
    
    fprintf('\n--- С FEC (кодовое отношение = %.2f, эффективное BER = %.4e) ---\n', codeRate, effectiveBER_FEC);
    fprintf('PER (с FEC): %.4f\n', PER_FEC);
    fprintf('Время передачи пакета (с FEC): %.3f сек\n', overallTime_FEC);
    fprintf('Эффективная скорость передачи (с FEC): %.2f бит/сек\n', effectiveRate_FEC);


    %% === Моделирование скорости передачи и оценка ошибок для наихудшего маршрута ===
% Предположим, что bestPath и totalDelay получены ранее из поиска маршрута по графу G_weight.
    % Параметры пакета и передачи
     % Параметры пакета
  % Симуляция передачи данных по лучшему маршруту с подсчетом ошибок для каждого узла

worstPath = shortestpath(G_weight, worst_i, worst_j);
packetSize = 240; % размер пакета в битах
hopErrors = zeros(1, length(worstPath)-1);
hopBERs = zeros(1, length(worstPath)-1);
hopSNRs = zeros(1, length(worstPath)-1);

for k = 1:length(worstPath)-1
    i1 = worstPath(k);
    i2 = worstPath(k+1);
    dist = norm(nodes(i1,:) - nodes(i2,:));
    effectivePathLoss = 20*log10(dist) + 20*log10(freq) - PL_margin + noiseOffset(i1, i2);
    rxPower = txPower - effectivePathLoss;
    SNR_hop = rxPower - sensitivity;
    hopSNRs(k) = SNR_hop;
   % Расчет BER на основе SNR
    if SNR_hop < -20
        BER_hop = 1; % При SNR ниже -20 дБ считаем, что все биты ошибочны
    else
        % Примерная оценка BER на основе SNR
        BER_hop = 0.5 * erfc(sqrt(10^(SNR_hop / 10)));
    end

    hopBERs(k) = BER_hop;
    errors = binornd(packetSize, BER_hop);
    hopErrors(k) = errors;
end



totalErrors = sum(hopErrors);
fprintf('\nСимуляция передачи по наихудшему маршруту:\n');
for k = 1:length(worstPath)-1
    fprintf('Хоп от узла %d к узлу %d: SNR = %.2f дБ, BER = %.4e, ошибок = %d из %d бит\n', ...
        worstPath(k), worstPath(k+1), hopSNRs(k), hopBERs(k), hopErrors(k), packetSize);
end
fprintf('\nОбщее число ошибок по маршруту: %d из %d бит (накоплено по хопам)\n', ...
    totalErrors, packetSize * (length(worstPath)-1));

%----------------- непосредственная модуляция через chirp-------------------


% Параметры LoRa CSS
SF = 12;                % Spreading Factor (количество бит на символ)
M = 2^SF;              % Размер созвездия
BW = 125e3;            % Полоса частот (Гц)
symDuration = 2^SF / BW;  % Длительность символа (сек)
Fs = 125e1;              % Частота дискретизации (Гц)
t_sym = linspace(0, symDuration, round(symDuration*Fs));
% Генерация базового чирпа (up-chirp) от 0 до BW
baseChirp = chirp(t_sym, 0, symDuration, BW);

% Количество символов в пакете
packetSymbols = packetSize / SF; 

% Для каждого перехода по маршруту bestPath будем моделировать передачу пакета
Niter = 1;  % число итераций для усреднения
hopBER = zeros(1, length(worstPath)-1);

for k = 1:length(worstPath)-1
    i1 = worstPath(k);
    i2 = worstPath(k+1);

    % Используем ранее рассчитанный SNR для данного хопа
    % Если SNR < -20, ставим его равным -20 для модели (порог LoRa)
    SNR_hop = hopSNRs(k);
    % Перевод SNR (дБ) в линейное отношение
    snr_linear = 10^(SNR_hop/10);
    % Энергия базового чирпа
    Es = mean(abs(baseChirp).^2);
    % Для передачи символа мощность сигнала равна Es; дисперсия шума определяется как:
    noiseVar = Es / snr_linear;

    totalBitErrors = 0;
    totalBits = packetSize * Niter;
    for iter = 1:Niter
        % Генерация случайных битов
        data_bits = randi([0 1], packetSize, 1);
        % Группируем биты по SF и преобразуем в символы
        dataSymbols = bi2de(reshape(data_bits, SF, []).', 'left-msb');

        txSignal = [];
        for s = 1:length(dataSymbols)
            % Для LoRa CSS символ модулируется циклическим сдвигом базового чирпа
            shift = round(dataSymbols(s) * length(baseChirp) / M);
            modChirp = circshift(baseChirp, shift);
            txSignal = [txSignal, modChirp];  %#ok<AGROW>
        end
        txSignal = txSignal(:);

        % Генерация шума AWGN с дисперсией noiseVar
        noise = sqrt(noiseVar/2) * (randn(size(txSignal)) + 1j*randn(size(txSignal)));
        rxSignal = txSignal + noise;

        % Демодуляция: разбиваем сигнал на символы и для каждого выполняем корреляцию
        rxSymbols = zeros(length(dataSymbols), 1);
        L_chirp = length(baseChirp);

        % Преобразуем baseChirp в вектор-столбец для корректных операций
        baseChirp = baseChirp(:); 

        for s = 1:length(dataSymbols)
            segment = rxSignal((s-1)*L_chirp+1 : s*L_chirp);
            correlations = zeros(1, M);

            % Преобразуем segment в столбец
            segment = segment(:); 

            for m = 0:M-1
                % Генерируем опорный чирп как вектор-столбец
                refChirp = circshift(baseChirp, m);
                % Поэлементное умножение векторов-столбцов
                corrValue = abs(sum(segment .* conj(refChirp))); 
                correlations(m+1) = corrValue;
            end

            [~, detected] = max(correlations);
            rxSymbols(s) = detected - 1;
        end

        % Преобразуем полученные символы обратно в биты
        rxBits = de2bi(rxSymbols, SF, 'left-msb').';
        rxBits = rxBits(:);
        % Считаем ошибки
        totalBitErrors = totalBitErrors + sum(data_bits ~= rxBits(1:packetSize));
    end
    hopBER(k) = totalBitErrors / (packetSize * Niter);
    fprintf('\n Переход %d -> %d: SNR = %.2f dB, BER = %.4e\n', worstPath(k), worstPath(k+1), SNR_hop, hopBER(k));
end

overallBER = prod(1 - hopBER);
fprintf('\n Общая вероятность безошибочной передачи по маршруту: %.4e\n', overallBER);
fprintf('Общая вероятность ошибки (PER) по маршруту: %.4e\n', 1 - overallBER);

% --------------------калькулятор по среднему снр--------------------------



    dataRate=290;
   % Расчет времени передачи без FEC (только для пакета)
    txTime_noFEC = packetSize / dataRate;  % в секундах

    % Предположим, что для каждого ребра маршрута мы уже рассчитали SNR (в предыдущей части)
    % Здесь для простоты расчитаем усредненное значение SNR по маршруту.



    % Усреднённое SNR по маршруту
    avgSNR = mean( hopSNRs);
    % Без FEC: при среднем SNR, модель exp(-max(SNR,0)/10) дает:
    BER_noFEC = exp(-max(avgSNR,0)/10);  % max(-5.43, 0) = 0, поэтому BER_noFEC = exp(0)=1
    PER_noFEC = 1 - (1 - BER_noFEC)^packetSize;  % Практически PER = 1
    txTime_noFEC = packetSize / dataRate;  % Время передачи полезного пакета без FEC
    overallTime_noFEC = length(worstPath)*(txTime_noFEC+0.2) + worstDelay;  % Общая задержка маршрута + время передачи
    effectiveRate_noFEC = packetSize / overallTime_noFEC;  % Эффективная скорость передачи полезных бит

    % С применением FEC:
    % Предположим, что алгоритм FEC снижает эффективное BER до фиксированного уровня:
    effectiveBER_FEC = BER_noFEC/10000;  % Допустим, остаточный BER = 10^-5 с хорошо настроенным FEC
    PER_FEC = 1 - (1 - effectiveBER_FEC)^packetSize;  
    codeRate = 0.5;
    txTime_FEC = (packetSize / codeRate) / dataRate;  % Время передачи с учетом избыточности
    overallTime_FEC = length(worstPath)*(txTime_FEC+0.2) + worstDelay;  % Общая задержка маршрута + время передачи
    effectiveRate_FEC = packetSize / overallTime_FEC;  % Эффективная скорость полезной передачи

    fprintf('\n--- Без FEC ---\n');
    fprintf('Среднее SNR: %.2f дБ\n', avgSNR);
    fprintf('BER (без FEC): %.4e\n', BER_noFEC);
    fprintf('PER (без FEC): %.4f\n', PER_noFEC);
    fprintf('Время передачи пакета (без FEC): %.3f сек\n', overallTime_noFEC);
    fprintf('Эффективная скорость передачи (без FEC): %.2f бит/сек\n', effectiveRate_noFEC);

    fprintf('\n--- С FEC (кодовое отношение = %.2f, эффективное BER = %.4e) ---\n', codeRate, effectiveBER_FEC);
    fprintf('PER (с FEC): %.4f\n', PER_FEC);
    fprintf('Время передачи пакета (с FEC): %.3f сек\n', overallTime_FEC);
    fprintf('Эффективная скорость передачи (с FEC): %.2f бит/сек\n', effectiveRate_FEC);

%% Визуализация

orbitParams = [];
for i = 1:length(satellites)
    if strcmp(satellites(i).group, 'Inclined')
        orbitParams = [orbitParams; satellites(i).r, satellites(i).inclination, satellites(i).RAAN]; %#ok<AGROW>
    end
end
[uniqueOrbits, ia, ~] = unique(orbitParams, 'rows');

% Определяем значения true anomaly от 0 до 2pi
theta_vals = linspace(0, 2*pi, 360);

% Определим разные цвета для разных групп
colors = struct('Equatorial', 'b', 'Polar', 'b', 'Inclined', 'b');

figure; hold on; grid on; axis equal;




% Отображение орбит (уникальных)
for i = 1:size(uniqueOrbits,1)
    r_val = uniqueOrbits(i,1);
    incl = uniqueOrbits(i,2);
    RAAN_val = uniqueOrbits(i,3);
    x_inc = r_val * (cos(RAAN_val)*cos(theta_vals) - sin(RAAN_val)*sin(theta_vals)*cos(incl));
    y_inc = r_val * (sin(RAAN_val)*cos(theta_vals) + cos(RAAN_val)*sin(theta_vals)*cos(incl));
    z_inc = r_val * sin(theta_vals)*sin(incl);
    plot3(x_inc, y_inc, z_inc, 'Color', colors.Inclined, 'LineWidth', 0.25);
end

% 
% 
% % Визуализируем орбиты для каждой высоты (для каждой группы отдельно)
for idx = 1:length(r_orbits)
    r = r_orbits(idx);
    theta_circle = linspace(0, 2*pi, 360);
    % Экваториальная орбита: z = 0
    plot3(r*cos(theta_circle), r*sin(theta_circle), zeros(size(theta_circle)), 'Color', colors.Equatorial, 'LineWidth', 0.25);
    % Полярная орбита: (при i=90, RAAN=0) x = r*cos(theta), z = r*sin(theta), y=0
    plot3(r*cos(theta_circle), zeros(size(theta_circle)), r*sin(theta_circle), 'Color', colors.Polar, 'LineWidth', 0.25);

end




% Отображение спутников
for i = 1:length(satellites)
    scatter3(satellites(i).x, satellites(i).y, satellites(i).z, 100, colors.(satellites(i).group), 'filled');
    % Подпись: номер спутника
    text(satellites(i).x, satellites(i).y, satellites(i).z, num2str(satellites(i).id),...
        'FontSize', 8, 'Color', 'w', 'HorizontalAlignment', 'center');
end


title(sprintf('Система спутников по орбитам (t = %d с)', t));
xlabel('X (м)'); ylabel('Y (м)'); zlabel('Z (м)');
view(3);
rotate3d on;


% Отображение связей
for i = 1:numNodes
    for j = i+1:numNodes
        if connections(i, j) == 1
            plot3([nodes(i,1), nodes(j,1)], [nodes(i,2), nodes(j,2)], [nodes(i,3), nodes(j,3)], 'g--');
        end
    end
end

% Выделяем лучший маршрут между наиболее удалёнными узлами (если найден)
if ~isempty(path)
    for k = 1:length(path)-1
        i1 = path(k);
        i2 = path(k+1);
        plot3([nodes(i1,1), nodes(i2,1)], [nodes(i1,2), nodes(i2,2)], [nodes(i1,3), nodes(i2,3)], 'm-', 'LineWidth', 2);
    end
end



% Выделяем маршрут с наихудшей задержкой, если он найден
if exist('worstPath','var') && ~isempty(worstPath)
    for k = 1:length(worstPath)-1
        i1 = worstPath(k);
        i2 = worstPath(k+1);
        % Выделяем маршрут более толстой красной линией
        plot3([nodes(i1,1), nodes(i2,1)], [nodes(i1,2), nodes(i2,2)], [nodes(i1,3), nodes(i2,3)], 'r-', 'LineWidth', 3);
    end
end

%% ЗЕМЛЯ
% Параметры Земли
radius = 6400000; % Радиус в метрах

% Создание сферы (30x30 граней для гладкости)
[X, Y, Z] = sphere(50);

% Масштабирование до заданного радиуса
X = X * radius;
Y = Y * radius;
Z = Z * radius;

% Загрузка текстуры Земли (требуется Image Processing Toolbox)
earthTexture = imread('earthmap.jpg');

surf(X, Y, Z, 'CData', flipud(earthTexture), 'FaceColor', 'texturemap',...
    'EdgeColor', 'none');
axis equal;

% Настройка освещения и материала
material shiny;
light('Style','infinite');
lighting phong;

% Настройка осей и вида
view(3);
rotate3d on;