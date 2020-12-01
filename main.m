clear
close all

%% VAR

N = 5000;%pow2(10);
p=[0.95*exp(i*pi/5),0.9*exp(i*3*pi/4),0.95*exp(i*pi/2),0.95*exp(-i*pi/5),0.9*exp(-i*3*pi/4),0.95*exp(-i*pi/2)];
A = poly(p);
B = [1];
sigma2=10;
bb=sqrt(sigma2)*randn(1,N); %Bruit pour generer les bails
x_k = filter(B,A,bb); %signal filtré
f=-1/2:1/N:1/2-1/N;
powerspectre=fftshift((abs(fft(x_k)))).^2/N;
% p_welch = pwelch(x_k,256);

Y = freqz(B,A,2*pi*f); %RI du filtre
DSP = Y.*conj(Y);

%signal powerspectre dsp
subplot 311
plot(x_k)
title('x')
subplot 312
plot(f,powerspectre)
title('Spectre de puissance')
subplot 313
plot(f,DSP)
title('DSP')

figure
plot(f,abs(Y))

% figure
% plot(p_welch)