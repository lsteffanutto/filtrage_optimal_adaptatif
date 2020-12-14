clear
close all

%% A RENDRE LE 22 JANVIER %%

%% VAR

disp1=1;

N = 5000;%pow2(10);
f=-1/2:1/N:1/2-1/N;
sigma2=10;
ordre_p=6;

%% PART1: Processus Auto-Regressif
p=[0.95*exp(i*pi/5),0.9*exp(i*3*pi/4),0.95*exp(i*pi/2),0.95*exp(-i*pi/5),0.9*exp(-i*3*pi/4),0.95*exp(-i*pi/2)];

A = poly(p);
B = [1];

bb=sqrt(sigma2)*randn(1,N); %Bruit pour generer les bails
x_k = filter(B,A,bb); %signal filtré
powerspectre=fftshift((abs(fft(x_k)))).^2/N;
% p_welch = pwelch(x_k,256);

Y = freqz(B,A,2*pi*f); %RI du filtre
DSP = Y.*conj(Y);

%% PART2: Yule-Walker et Filtrage Adaptatif

%%% 1.YULE-WALKER (voir cahier) %%%

%Estimation des paramètres AR

%vecteur de corrélation
R_xx_tot=xcorr(x_k);
% figure,plot(R_xx_tot)

%vecteur de corrélation d'ordre p=6 => [Rxx(0),...,Rxx(p)]
r_xx_p=R_xx_tot(N+1:N+ordre_p)';

%matrice R_X de covariance (Matrice Toeplitz)
c = R_xx_tot(N:N+ordre_p-1);    % 1ère colonne
r = R_xx_tot(N:-1:N-ordre_p+1); % 1ère ligne
R_X=toeplitz(c,r);              % Matrice finale

AR_estime=(-(R_X^-1)*r_xx_p)'  % Paramètres AR ESTIMES

AR_reel=A(2:6)                % Paramètres AR REELS

%DSP avec paramètres AR estimés
Y_reel = freqz(B,[1 AR_estime],2*pi*f); %RI du filtre
DSP_reelle = Y_reel.*conj(Y_reel);

figure,plot(f,DSP,'LineWidth',2)
hold on; plot(f,DSP_reelle,'LineWidth',1);
title("Distance spectrale")
legend("Paramètres AR réels","Paramètres AR estimés")

%Filtrage Adaptatif (regarder les pdf de Grivel)










%% PLOTS
if disp1==1
    figure,
%signal powerspectre dsp
subplot 311
plot(x_k)
title('x')
subplot 312
% plot(f,20*log10(powerspectre/pi))
plot(f,powerspectre/N*100)

title('Spectre de puissance')
subplot 313
plot(f,DSP)
title('DSP')

figure
plot(f,abs(Y))
end
