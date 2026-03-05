clear; close all; clc;

%% =========================
% 1. ПАРАМЕТРЫ
%% =========================

fc  = 868e6;       % несущая
BW  = 125e3;       % полоса
SF  = 10;          
v   = 8000;        
c   = 3e8;

Fs  = 10*BW;       % частота дискретизации (увеличена для точности)
beta = v/c;

fprintf('beta = %.3e\n', beta);

%% =========================
% 2. ВРЕМЯ И ЧИРП
%% =========================

Ts = (2^SF)/BW;
t  = 0:1/Fs:Ts-1/Fs;

k = BW/Ts;

% baseband LoRa up-chirp
s = exp(1j*pi*k*t.^2);

%% =========================
% 3. ТОЧНАЯ МОДЕЛЬ ДОППЛЕРА
%% =========================

% масштабированное время
t_scaled = t*(1 - beta);

s_scaled = interp1(t, s, t_scaled, 'spline', 0);

% RF допплер после смешения вниз остаётся CFO = beta*fc
fD = beta*fc;

r_exact = s_scaled .* exp(1j*2*pi*fD*t);

%% =========================
% 4. ПРОСТАЯ CFO МОДЕЛЬ (для сравнения)
%% =========================

r_cfo = s .* exp(1j*2*pi*fD*t);

%% =========================
% 5. СПЕКТРЫ
%% =========================

Nfft = 2^16;

S        = fftshift(fft(s, Nfft));
R_exact  = fftshift(fft(r_exact, Nfft));
R_cfo    = fftshift(fft(r_cfo, Nfft));

f_axis = linspace(-Fs/2, Fs/2, Nfft)/1e3;

%% =========================
% 6. ГРАФИКИ
%% =========================

figure;
plot(f_axis,20*log10(abs(S)/max(abs(S))),'LineWidth',1.2);
grid on;
title('Исходный LoRa спектр');
xlabel('Частота (кГц)');
ylabel('Амплитуда (дБ)');

figure;
plot(f_axis,20*log10(abs(R_cfo)/max(abs(R_cfo))),'LineWidth',1.2); hold on;
plot(f_axis,20*log10(abs(R_exact)/max(abs(R_exact))),'--','LineWidth',1.2);
grid on;
legend('Простая CFO модель','Точная модель');
title('Сравнение моделей допплера');
xlabel('Частота (кГц)');
ylabel('Амплитуда (дБ)');