import matlab.engine
eng = matlab.engine.start_matlab()

%clear; close all; clc;
%Matlab excercise 1 : Hamed Talebian - 150360360

%3.1: Signal generation-time domain
Fs=16000; % sampling frequency  
Ts=1/Fs; % sampling interval 
t=0:Ts:0.1; % sampling time instants [s] 
x=sin(2*pi*800*t); % signal vector 
figure % Opens a new figure 
plot(t,x) % plot in time domain;
title('Time domain Plot of x(t)') 
xlabel('t [s]')  
ylabel('Amplitude') 
axis([0 0.02 -1.1 1.1]);

% Discrete Fourier Transform (DFT)
F_x=fft(x);                         % DFT of X, saved to Fx  
N=length(x);  
Fo=1/(Ts*N);                       % frequency resolution. 
freq1=0:Fo:(N-1)*Fo;               % One-sided frequency Axis  
figure 
plot(freq1,abs(F_x)/N)             % One-sided amplitude Spectrum
title('One side spectrum (0-F)') 
xlabel('f')  
ylabel('X [f]')
freq2=-N/2*Fo:Fo:(N/2-1)*Fo;      % Two-sided frequency Axis  
figure 
plot(freq2,fftshift(abs(F_x)/N))   % Two-sided amplitude Spectrum
title('Two sideed spectrum (-Fs/2 - Fs/2)') 
xlabel('f')  
ylabel('X [f] ')

%% 3.2 new signal and multiplication 750MHz sinusoidal 
m=sin(2*pi*750*t); %making sinusodial signal m
s=m.*x; %multiplication
figure 
plot(t,s)
title('Time domain Plot of s(t)')
xlabel('t [s]')  
ylabel('S(t) : Amplitude')
axis([0 0.01 -1.1 1.1]);
%DFT : lenght s = lenght x----> just plot the F(s)
F_s=fft(s);                         % DFT of X, saved to Fx  
figure 
plot(freq2,fftshift(abs(F_s)/N))   % Two-sided amplitude Spectrum
title('Two sideed spectrum (-Fs/2 - Fs/2)') 
xlabel('f')  
ylabel('F [s] ')

%% 3.3 Adding a Noise signal 
n = randn(size(x));
y = 10*s + n;
figure
plot(t,y)
axis([0 0.02 -15 15]);
title('Time domain Plot of y(t)')
xlabel('t [y]')  
ylabel('y(t) : Amplitude')
axis([0 0.01 -1.1 1.1]);
%DFT 
F_y=fft(y);                         % DFT of Y, saved to Fy  
plot(freq2/1e6,fftshift(abs(F_y)/N))   % Two-sided spectrum in 1 MHz range 
title('Two sideed spectrum (-Fy/2 - Fy/2)') 
xlabel('f in MHZ')  
ylabel('F(y)');

%% 4. Linear filtering 
% 4.1: Low pass butter worth filter
f_cut = 200; 
order = 10;
fr = f_cut/(Fs/2); % Cut-off frequency normalized to 1. 
[b,a]= butter(order,fr); % Coefficients of the filter
freqz(b,a,N,Fs)  
title('Frequency response of the Butterworth filter');

%4.1.1. filtering y(t) with butter filter 
y_filtered_butter = filter(b,a,y);

%plot time domain
figure 
plot(t,y)
hold on
plot(t,y_filtered_butter)%time-domain filtered signal y_filtered_butter(t)
legend('Input Data y(t)','Filtered Data')
title('Time domain Plot of y(t) and Butter filtered signal')
xlabel('t [s]')  
ylabel('S(t) : Amplitude')
hold off;


% DFT + plot two side amplitude spectrum + F(y)
F_y_Filtered_butter =fft(y_filtered_butter);       % DFT of Filtered Y
figure
plot(freq2/1e6,fftshift(abs(F_y)/N))          % two-sided amplitude Spectrum
hold on
plot(freq2/1e6,fftshift(abs(F_y_Filtered_butter)/N)) %original f(y) signal before filterning
title('Two sided spectrum and Butter filtered of F(y)')
xlabel('f in MHZ') 
ylabel('F(y)');
legend('Input Data F(y)','Filtered Data')
hold off;


%% 4.2. BAND PASS FIR (FINITE IMPULSE RESPONSE) FILTER
order = 60; % filter order 
f_filter = [0 0.8e9 1.3e9 1.8e9 2.3e9 1.6e10/2]/(1.6e10/2); 
a_filter = [0 0 1 1 0 0];
b = firpm(order,f_filter,a_filter);
figure
stem(-order/2:order/2,b)
title('Impulse response')
ylabel('Desired amplitude')
xlabel('filter order');
F_b = fft(b,N); % We use same length with the frequency vector for the fft
y_filtered_FIR = filter(b,1,y);
F_y_filtered_FIR = fft(y_filtered_FIR);
%plot
figure 
plot(t,y)
hold on
plot(t,y_filtered_FIR)%time-domain filtered signal y_filtered_butter(t)
legend('Input Data y(t)','Filtered Data')
title('Time domain Plot of y(t) and Butter FIR signal')
xlabel('t [s]')  
ylabel('S(t) : Amplitude')
hold off
%axis([0 0.01 -7 7])
% plot two side amplitude spectrum + F(y)
figure
plot(freq2/1e6,fftshift(abs(F_y)/N))          % two-sided amplitude Spectrum
xlabel('f in MHZ')  
hold on 
plot(freq2/1e6,fftshift(abs(F_y_filtered_FIR)/N)); %original f(y) signal before filterning
title('Two sided spectrum and FIR of F(y)')
xlabel('f in MHZ') 
ylabel('F(y)');
legend('Input Data F(y)','Filtered Data')
hold off

