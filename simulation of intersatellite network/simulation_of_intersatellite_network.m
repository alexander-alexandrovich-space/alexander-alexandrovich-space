clc; clear; close all;

%% === Параметры сети LoRa ===
txPower = 20;                
sensitivity = -137;          % Чувствительность приемника (dBm)
freq = 868e6;               
Gtx=2;
Grx=2;

rng(42); % сид
%% === Физические параметры и настройки орбитальной модели ===
mu = 3.986004418e14;    % Гравитационный параметр Земли, м^3/с^2
R_earth = 6400000;      

r_orbits = [(160000+2200000)/2+6400000];  %  высоты орбит в одной плоскости

t = 100000;  % время симуляции

%% === Параметры размещения спутников ===
numPerOrbit = 28;  % число спутников на каждой орбите

% Для каждой группы орбит будем генерировать спутники для 3 типов:
% 1. Экваториальные: i = 0, RAAN = 0.
% 2. Полярные: i = 90, RAAN = 0.
% 3. Наклонные: i = 60 
%

satellites = [];  
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
        sat.isIsolated=true;
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
        sat.isIsolated=true;
        satellites = [satellites; sat];
    end
end

% --- Группа 3: Наклонные орбиты  ---


for m=1:10%m=1:round(2*pi*6560000/1500000)
    inclination_inc = deg2rad(60);
    RAAN_inc = deg2rad(36*m); %RAAN_inc = deg2rad(rad2deg(acos((-(1500)^2+2*(6560)^2)/(2*(6560)^2)))*m);
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
            sat.isIsolated=true;
            satellites = [satellites; sat];
        end
    end
end
%% === Переход к узлам сети ===
numNodes= length(satellites);
nodes = zeros(length(satellites), 3);
for i = 1:length(satellites)
    nodes(i, :) = [satellites(i).x, satellites(i).y, satellites(i).z];
end
%% === Диаграмма направленности антенн и ослабление ===

min_gain_dB = -20;   
attenuation_dB = zeros(numNodes);


% Поворот вокруг радиус-вектора
beta_angle = deg2rad(90);

% Для каждого спутника вычисляем вектор ориентации (по направлению движения)
for i = 1:numNodes
    pos = [satellites(i).x; satellites(i).y; satellites(i).z];
    nodes(i,:) = pos';
    
    % Ось орбиты (нормаль к плоскости)
    incl = satellites(i).inclination;
    RAAN = satellites(i).RAAN;
    normal = [sin(RAAN)*sin(incl); -cos(RAAN)*sin(incl); cos(incl)];

    % Скорость как cross(omega, r)
    vel = cross(normal, pos);
    v_dir = vel / norm(vel);
    
    satellites(i).antennaOrient = v_dir;
    
     %  Находим вектор в плоскости орбиты, перпендикулярный радиальному,
     %  поворачиваем базис
        perp_vec = cross(normal, v_dir);
        perp_vec = perp_vec / norm(perp_vec);
        
    
        R = [cos(beta_angle) -sin(beta_angle) 0;
             sin(beta_angle)  cos(beta_angle) 0;
             0               0              1];
        

        M = [ v_dir perp_vec normal];
        
     
        rotated_orient_local = R * [0; 0; 1];  

        % Запоминаем ориентацию антенны
        satellites(i).antennaOrient = M * rotated_orient_local;
end


% Расчет ослабления для каждого соединения (дипольная модель)
for i = 1:numNodes
    ori = satellites(i).antennaOrient/norm(satellites(i).antennaOrient);
    for j = 1:numNodes
        if i ~= j
            vec_ij = nodes(j,:)' - nodes(i,:)';
            d_dir = vec_ij / norm(vec_ij);
            
            % Косинус угла между направлением антенны и вектором к цели
            cos_angle = dot(ori, d_dir);
            angle = acos(cos_angle);
            
        
            gain_linear = sin(angle);
            
          
            angle_loss_dB = 20*log10(max(10^(min_gain_dB/20), gain_linear));
            attenuation_dB(i,j) = angle_loss_dB;
        else
            attenuation_dB(i,j) = 0; 
        end
    end
end

%% === Вычисление дальности связи и сбор канальных метрик ===
c=3e8;
SF =12;
BW = 125000;

CodeRateFec = 4/8;
CodeRateNoFec = 1;

CODERATE=CodeRateFec;

CASE=1; % 1- эмпирическая инженерная модель, 2 - q-функция; A*exp(margin*B)
% margin - на сколько дб выше порогового с-ш, А - пороговая BERмощностит апмпапаот, B: BER=A*e^-margin*B

connections = zeros(numNodes);
W = inf(numNodes);
realrange=zeros(numNodes);
rxPower=zeros(numNodes);
SNR = zeros(numNodes);
BER = zeros(numNodes);
RATE = zeros(numNodes);
LoRalimit = zeros(numNodes);

lambda = c/freq;
FSPLc = 20*log10(4*pi/lambda);

noiseOffset=5;

for i = 1:numNodes
    for j = i+1:numNodes
        
        dist = norm(nodes(i,:) - nodes(j,:));
        
        rxPower(i,j) = txPower + Gtx + Grx + attenuation_dB(i,j) + attenuation_dB(j,i) - 20*log10(dist) - FSPLc;
        rxPower(j,i) = txPower + Gtx + Grx + attenuation_dB(i,j) + attenuation_dB(j,i) - 20*log10(dist) - FSPLc;

        % рассчет метрик
        [SNR(i,j), BER(i,j), RATE(i,j), LoRalimit(i,j)] = linkMetrics_LoRa_BER(rxPower(i,j),BW,SF, CASE, CODERATE);
        [SNR(j,i), BER(j,i), RATE(j,i), LoRalimit(j,i)] = linkMetrics_LoRa_BER(rxPower(j,i),BW,SF, CASE, CODERATE);

        realrange(i,j)=10^((txPower - sensitivity + Gtx + Grx + attenuation_dB(i,j) + attenuation_dB(j,i)  - FSPLc - noiseOffset)/20);              
        realrange(j,i)=10^((txPower - sensitivity + Gtx + Grx + attenuation_dB(i,j) + attenuation_dB(j,i)  - FSPLc - noiseOffset)/20);
        if SNR(i,j) >= min_gain_dB && SNR(j,i) >= min_gain_dB && dist<realrange(i,j) && dist<realrange(j,i)
            connections(i,j)=1;
            connections(j,i)=1;
            totalDelay = dist / c; 
            W(i, j) = totalDelay;
            W(j, i) = totalDelay;
        end
    end
end

rr=10^((txPower - sensitivity + Gtx + Grx - FSPLc - noiseOffset)/20);
       

fprintf('\nРеальная дальность (дополнительный шум 5дБ без диаграмм направленности) %.2f км:\n', rr/1000);


%% === Изолированные узлы и кластеры ===
isolatedNodes = [];
for i = 1:numNodes
    for j = 1:numNodes
        if i == j, continue; end 
        if connections(i,j)==1
            satellites(i).isIsolated = false;
            break; 
        end
    end
    if satellites(i).isIsolated == true
        satellites(i).group = 'Isolated';
        isolatedNodes(end+1) = i;
    end
end

fprintf('\nИзолированные узлы (расстояние > %.2f км):\n', rr/1000);
disp(isolatedNodes);


G = graph(connections);
components = conncomp(G); 
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
G_conn = graph(connections);

maxD = -Inf;
best_i = 0;
best_j = 0;
for i = 1:numNodes
    for j = i+1:numNodes
        if ~isempty(shortestpath(G_conn, i, j))
            if D_phys(i,j) > maxD
                maxD = D_phys(i,j);
                best_i = i;
                best_j = j;
            end
        end
    end
end

% Произвольные узлы
%best_i=1;
%best_j=500;

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

%% === Поиск наихудшего маршрута (максимальная задержка) ===
allDelays = distances(G_weight);

finiteDelays = allDelays(~isinf(allDelays) & allDelays > 0);

if isempty(finiteDelays)
    fprintf('Нет связных путей в сети.\n');
else
    worstDelay = max(finiteDelays);
    fprintf('Наихудшая (максимальная) накопленная задержка по сети: %.3f секунд\n', worstDelay);
    
    % Найдем первую пару узлов, для которых задержка равна worstDelay
    [worst_i, worst_j] = find(allDelays == worstDelay, 1);  
    worstPath = shortestpath(G_weight, worst_i, worst_j);
    fprintf('\nПара узлов с максимальной задержкой: %d и %d (%.2f км)\n', worst_i, worst_j, D_phys(worst_i,worst_j)/1000);
    disp('Наихудший маршрут (узлы):');
    disp(worstPath);

end

%% === Моделирование скорости передачи и оценка ошибок для наилучшего маршрута ===

bestpath = shortestpath(G_weight, best_i, best_j);

packetSize = 240; % размер пакета в битах

hopBERs = zeros(1, length(bestpath)-1);
hopSNRs = zeros(1, length(bestpath)-1);
hopErrors = zeros(1, length(bestpath)-1);
dataRate = zeros(1, length(bestpath)-1);
times = zeros(1, length(bestpath)-1);

for k = 1:length(bestpath)-1
    i1 = bestpath(k);
    i2 = bestpath(k+1);
    hopSNRs(k) = SNR(i1,i2);
    hopBERs(k) = BER(i1,i2);
    dataRate(k) = LoRalimit(i1,i2);
    if dataRate(k) > 293
        dataRate(k)=293;
    end
    times(k)=packetSize/dataRate(k);
    errors = binornd(packetSize, BER(i1,i2));
    hopErrors(k) = errors;
end

totalErrors = sum(hopErrors);
overallTime = sum(times) + 0.2 * length(bestpath) + totalDelay; 
relativeRate = packetSize/overallTime;

fprintf('\nСимуляция передачи по "наилучшему" маршруту:\n');
for k = 1:length(bestpath)-1
    fprintf('Хоп от узла %d к узлу %d: SNR = %.2f дБ, BER = %.4e%%, ошибок = %d из %d бит\n', ...
        bestpath(k), bestpath(k+1), hopSNRs(k), hopBERs(k)*100, hopErrors(k), packetSize);
end
fprintf('\nОбщее число ошибок по маршруту: %d из %d бит (накоплено по хопам)\n', ...
    totalErrors, packetSize * (length(bestpath)-1));

avgSNR1 = mean (hopSNRs);

 fprintf('Время передачи пакета: %.3f сек\n', overallTime);
 fprintf('Эффективная скорость передачи: %.2f бит/сек\n', relativeRate);
    %% === Моделирование скорости передачи и оценка ошибок для наихудшего маршрута ===

worstpath = shortestpath(G_weight, worst_i, worst_j);

packetSize = 240; % размер пакета в битах

hopBERs = zeros(1, length(worstpath)-1);
hopSNRs = zeros(1, length(worstpath)-1);
hopErrors = zeros(1, length(worstpath)-1);
dataRate = zeros(1, length(worstpath)-1);
times = zeros(1, length(worstpath)-1);
worstlen = zeros(1, length(worstPath)-1);

for k = 1:length(worstpath)-1
    i1 = worstpath(k);
    i2 = worstpath(k+1);
    worstlen(k)=D_phys(i1,i2);
    hopSNRs(k) = SNR(i1,i2);
    hopBERs(k) = BER(i1,i2);
    dataRate(k) = LoRalimit(i1,i2);
    if dataRate(k) > 293
        dataRate(k)=293;
    end
    times(k)=packetSize/dataRate(k);
    errors = binornd(packetSize, BER(i1,i2));
    hopErrors(k) = errors;
end

totalErrors = sum(hopErrors);
overallTime = sum(times) + 0.2 * length(worstpath) + worstDelay; 
relativeRate = packetSize/overallTime;

fprintf('\nСимуляция передачи по "наихудшему" маршруту:\n');
for k = 1:length(worstpath)-1
    fprintf('Хоп от узла %d к узлу %d: SNR = %.2f дБ, BER = %.4e%%, ошибок = %d из %d бит\n', ...
        worstpath(k), worstpath(k+1), hopSNRs(k), hopBERs(k)*100, hopErrors(k), packetSize);
end
fprintf('\nОбщее число ошибок по маршруту: %d из %d бит (накоплено по хопам)\n', ...
    totalErrors, packetSize * (length(worstpath)-1));

 fprintf('Время передачи пакета: %.3f сек\n', overallTime);
 fprintf('Эффективная скорость передачи: %.2f бит/сек\n', relativeRate);
avgSNR2 = mean (hopSNRs);

%% === Визуализация ===


figure; hold on; grid on; axis equal;

 scale = 1000000; % Масштаб для видимости стрелки антенны
 
orbitParams = [];
for i = 1:length(satellites)
    quiver3(satellites(i).x, satellites(i).y, satellites(i).z,scale*satellites(i).antennaOrient(1),scale*satellites(i).antennaOrient(2),scale*satellites(i).antennaOrient(3),'cyan', 'LineWidth',1);
    if strcmp(satellites(i).group, 'Inclined')
        orbitParams = [orbitParams; satellites(i).r, satellites(i).inclination, satellites(i).RAAN]; %#ok<AGROW>
    end
end
[uniqueOrbits, ia, ~] = unique(orbitParams, 'rows');



% Определяем значения true anomaly от 0 до 2pi
theta_vals = linspace(0, 2*pi, 360);

% Определим разные цвета для разных групп
colors = struct('Equatorial', 'b', 'Polar', 'b', 'Inclined', 'b', 'Isolated', 'y');

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



% Визуализация диаграммы направленности (тороидальная форма)
sphere_res = 10;
for i = 1:numNodes
    pos = nodes(i,:);
    v_dir = satellites(i).antennaOrient;
    
    % Создаем тороидальную поверхность
    [X, Y, Z] = create_toroidal_pattern(pos, v_dir, 700000, sphere_res);
    [x1, y1, z1] = create_toroidal_pattern(pos, -v_dir, 700000, sphere_res);

    % Отрисовка полупрозрачной поверхности
    surf(X, Y, Z, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'b');
    surf(x1, y1, z1, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'b');

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



%% === Земля ===
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



axis off                   % отключить оси
set(gca, 'Color', 'none') % сделать фон области графика прозрачным
set(gcf, 'Color', 'none') % сделать фон всей фигуры прозрачным
box off                   % отключить рамку вокруг графика

%% === Вспомогательная функция для создания тороидальной диаграммы ===
function [X, Y, Z] = create_toroidal_pattern(center, orientation, radius, res)
    [theta, phi] = meshgrid(linspace(0, 2*pi, res), linspace(0, pi, res));
    
    % Параметры тора для диаграммы направленности
    R_torus = radius * 0.7;  
    r_torus = radius * 0.3;  
    
    % Координаты тора в локальной системе (ось Z - ориентация антенны)
    x_loc = (R_torus + r_torus*cos(phi)) .* cos(theta);
    y_loc = (R_torus + r_torus*cos(phi)) .* sin(theta);
    z_loc = r_torus * sin(phi);
    
    % Поворачиваем локальную систему в соответствии с ориентацией антенны
    z_axis = orientation(:);
    z_axis = z_axis / norm(z_axis);
    
    % Создаем ортогональную систему координат
    if abs(z_axis(3)) > 0.8
        x_axis = [1; 0; 0];
    else
        x_axis = [0; 0; 1];
    end
    y_axis = cross(z_axis, x_axis);
    y_axis = y_axis / norm(y_axis);
    x_axis = cross(y_axis, z_axis);
    
    X = center(1) + x_axis(1)*x_loc + y_axis(1)*y_loc + z_axis(1)*z_loc;
    Y = center(2) + x_axis(2)*x_loc + y_axis(2)*y_loc + z_axis(2)*z_loc;
    Z = center(3) + x_axis(3)*x_loc + y_axis(3)*y_loc + z_axis(3)*z_loc;
end
%% === Функция для шума и BER ===
function [SNR_dB, BER, C_shannon, Rb] = linkMetrics_LoRa_BER(P_rx_dBm, BW_Hz, SF, CASE, R)
     NF_dB = 6;
switch CASE
    case 1
        k = 1.38064852e-23;
        T = 290;
        F = 10^(NF_dB/10);
    
        % Шумовая мощность
        Pn = k * T * BW_Hz * F;
        Pn_dBm = 10 * log10(Pn * 1e3);  
    
        % SNR
        SNR_dB = P_rx_dBm(:) - Pn_dBm;
        SNR_lin = 10.^(SNR_dB/10);
        N = numel(SNR_dB);
        BER = zeros(N,1);
    
        % Табличные минимальные SNR для успешной демодуляции (Semtech AN1200.22)
        SNR_thresh = [-7.5, -10, -12.5, -15, -17.5, -20];  % для SF7…SF12
        sf_vals = 7:12;
        
        % Защита от выхода за пределы таблицы
        idx = find(sf_vals == SF, 1);
        if isempty(idx)
            error('Unsupported SF value');
        end
        SNR_min = SNR_thresh(idx);
    
        % Модель BER: экспоненциальное снижение при превышении порога
        for i = 1:N
            margin = SNR_dB(i) - SNR_min;
            if margin < 0
                BER(i) = 0.5;  % связь невозможна
            else
                BER(i) = 0.01 * exp(-0.6 * margin / R);  % эмпирическая зависимость
            end
        end
    case 2
       
        % P_rx_dBm — вектор мощностей на входе, дБм (Nx1)
        N = numel(P_rx_dBm);
        % 1) средняя шумовая мощность
        k = 1.38064852e-23;    % Дж/К
        T = 290;               % К
        F = 10^(NF_dB/10);     % линеаризованный NF
        Pn_mean = k * T * BW_Hz * F;  % средняя в ваттах
    
        % 2) генерируем AWGN (комплексный) для каждого узла
        sigma = sqrt(Pn_mean/2);
        noise_I = sigma * randn(N,1);
        noise_Q = sigma * randn(N,1);
        noise_W = noise_I.^2 + noise_Q.^2;
    
        % 3) переводим в дБм и считаем SNR
        Pn_dBm = 10*log10(noise_W*1e3);
        SNR_dB = P_rx_dBm(:) - Pn_dBm;
        SNR_lin = 10.^(SNR_dB/10);
    
        % 4) SER и BER для CSS (прибл.)
        M = 2^SF;
        SER = qfunc( sqrt(2 * SNR_lin * SF / M) );
        BER = SER ./ (SF * R);
end

    C_shannon  = BW_Hz * log2(1 + SNR_lin);  % [бит/с]
    Rb = (BW_Hz/2^SF) * SF * R;  % [бит/с]


end