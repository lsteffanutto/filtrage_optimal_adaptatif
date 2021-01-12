clear all
close all
clc;
beep off;
%% A RENDRE LE 22 JANVIER %%

%% VAR

disp_part_1=1;
disp_part_estimation_parametre_AR=0;
disp_part_LMS=0;
disp_part_Kalman=0;

% nombre de parametre AR a ajouter pour une surestimation
surestimation=0;

N = 5000;%pow2(10);
f=-1/2:1/N:1/2-1/N;
sigma2=10;
ordre_p=6;
p1=0.95*exp(i*pi/5);
p2=0.9*exp(i*3*pi/4);
p3=0.95*exp(i*pi/2);

%% PART1: Processus Auto-Regressif (faire fonction)

% Pôles du filtre
% p=[0.95*exp(i*pi/5),0.9*exp(i*3*pi/4),0.95*exp(i*pi/2),0.95*exp(-i*pi/5),0.9*exp(-i*3*pi/4),0.95*exp(-i*pi/2)];
p=[p1,p2,p3,conj(p1),conj(p2),conj(p3)];

% retourne les coeff AR à partir des pôles
A = poly(p)';   

% retourne les pôles à partir des coeff AR
%poles=roots(A);

B = 1;

bb=sqrt(sigma2)*randn(1,N); %Bruit pour generer le processus AR
x_k = filter(B,A,bb);       %signal filtré => Processus AR d'ordre 6

% Spectre de Puissance
powerspectre=fftshift((abs(fft(x_k)))).^2/N;
% p_welch = pwelch(x_k,256);

%RI et DSP du filtre
Y = freqz(B,A,2*pi*f); 
DSP = Y.*conj(Y);

%% PART 2
%% 1. Estimation des paramètres AR - YULE-WALKER (voir cahier) %%%

%% faire une fonction pour ça

% Vecteur de corrélation
R_xx_tot=xcorr(x_k);
% figure,plot(R_xx_tot)

% Vecteur rx de corrélation d'ordre p=6 => [Rxx(0),...,Rxx(p)]
r_xx_p=R_xx_tot(N+1:N+ordre_p)';
% figure,plot(r_xx_p)

% Matrice R_X de covariance     (Matrice Toeplitz (1.85) p.16 cours PDF1)
c = R_xx_tot(N:N+ordre_p-1);    % 1ère colonne
r = R_xx_tot(N:-1:N-ordre_p+1); % 1ère ligne
R_X=toeplitz(c,r);              % Matrice finale

% Paramètres AR ESTIMES (Forme matricielle (1.87) p.17 PDF1)
AR_estime=(-(R_X^-1)*r_xx_p)  

% Paramètres AR REELS
AR_reel=A(2:7)                 

% DSP avec paramètres AR estimés
Y_reel = freqz(B,[1 AR_estime'],2*pi*f); %RI du filtre avec paramètres estimés
DSP_reelle = Y_reel.*conj(Y_reel);       %DSP à partir RI du filtre

% Comparer selon le nombre d'echantillon N (tend vers inf , apporche mieux les paramètres)
% N=500; N=5000 ; N=50000

%% 2.Filtrage Adaptatif (p.5 cours PDF2)
if disp_part_LMS==1

% Si surestimation de l'ordre
if surestimation~=0
    surestimation_LMS=surestimation
    ordre_p=ordre_p+surestimation_LMS;
end

% On initialise au hasard paramètres AR
H_coeff=zeros(ordre_p,1);

% nombre aleatoire entre -0.5 et 0.5 pour init les coeff 
% for i = 1:ordre_p
%     H_coeff(i)=rand-0.5;
% end

% Alpha optimal avec critère vp_max
% [alpha] = alpha_optimal_LMS(x_k);

% Vaut mieux prendre un alpha + petit pour que ça converge bien
alpha=0.0002; 

% Graphique rapport
% 1. NMC=1; alpha=0.0002
% 2. NMC=10; alpha=0.00015; N=5000

% Monte Carlo
NMC=1;

% Iterations LMS
boucle_filtre=N;

%stocker les valeurs des paramètres AR et moyenner les réalisations
H_coeff_tot=zeros(ordre_p,boucle_filtre); 

% Filtre LMS + MC (pour moyenner réalisations)
for i=1:NMC
    
% iterations LMS: X_n=input=bb, d(k)=x_k=signal désiré (on prend les 6 derniers échantillons de chaque)
for k = ordre_p+1:boucle_filtre
    
    [H_coeff] = f_LMS(H_coeff,alpha,-x_k(k-1:-1:k-ordre_p),x_k(k)); %MAJ coeff filtre

%     H_coeff_real(:,k)=H_coeff; % juste une réalisation pour voir
    
    H_coeff_tot(:,k)= H_coeff_tot(:,k)+H_coeff; % Stocker toutes les réalisations
    
end
    
end

% On moyenne les réalisations
for i =1:N
    H_coeff_tot(:,i)=H_coeff_tot(:,i)/NMC; 
end
% AR estimés LMS
AR_LMS=H_coeff 

end % FIN LMS



%% 3.Filtrage de Kalman (p.19 cours PDF3)
if disp_part_Kalman==1

% Si surestimation de l'ordre
if surestimation~=0
    surestimation_AR=surestimation
    ordre_p=ordre_p+surestimation_AR;
end

% Init paramètres de Kalman

% Au début beaucoup d'erreur car paramètres init au hasard
alpha = 1e4; 

% Matrice covariance de l'erreur
P = alpha*eye(ordre_p);       

% Phi
Phi=eye(ordre_p);

% H contient les [x(k)...x(k-6)] échantillons
H = [1 zeros(1,ordre_p-1)];
% Q=eye(ordre_p);

% Bruit de modèle = 0 car parametre a estimer ne bouge pas et reste
% constant
Q=0;

% Bruit de mesure
sigma2 = 10;

% Iterations Kalman
boucle_filtre=N/5;

% On init le vecteur d'état
x=zeros(ordre_p,1);

%stocker les valeurs des paramètres AR à chaque itération
H_coeff_tot_Kalman=zeros(ordre_p,boucle_filtre);

for k = ordre_p+1:boucle_filtre
    [x,P] = Kalman_processus_AR(x,P,x_k(k),Q,sigma2,Phi,-x_k(k-1:-1:k-ordre_p));
    H_coeff_tot_Kalman(:,k)= x;
end

AR_Kalman=x
    
end


%% PLOTS
if disp_part_1==1           %% PART 1 %%
figure,
%signal powerspectre dsp
subplot 211
plot(x_k)
title("Processus AR d'ordre 6 généré")
xlabel('Nombre échantillons');
ylabel('Amplitude');

subplot 212
% plot(f,20*log10(powerspectre/pi))     %Truc de Wedji pour plot en log
plot(f,powerspectre/N*100)
xlabel('Nombre échantillons');
ylabel('Amplitude');

title('Spectre de puissance et DSP');
hold on; plot(f,DSP)
legend('Spectre de puissance','DSP');

figure,plot(f,abs(Y))
xlabel('Fréquence normalisée (Hz)');
ylabel('Amplitude');
title('RII du filtre');
end

if disp_part_estimation_parametre_AR==1
figure,

%DSP avec paramètres AR estimés vs Reelle
plot(f,DSP,'LineWidth',2)
hold on; plot(f,DSP_reelle,'LineWidth',1);
title("Distance spectrale")
xlabel('Fréquence normalisée (Hz)');
ylabel('Amplitude');
legend("RII du filtre estimé","RII du filtre réel")
    
end

if disp_part_LMS==1         %% LMS %%
x=1:boucle_filtre;    
figure,
for i =1:ordre_p-surestimation
    hold on;plot([0 boucle_filtre],[AR_reel(i) AR_reel(i)],'LineWidth',2);
end
xlabel('itération');
ylabel('évolution estimation');
% ylim([-0.6 1]);
title("estimation des coeff AR par filtre LMS à chaque itération, alpha="+alpha);
hold on;scatter(x,H_coeff_tot(1,:),'b.');
hold on;scatter(x,H_coeff_tot(2,:),'r.');
hold on;scatter(x,H_coeff_tot(3,:),'y.');
hold on;scatter(x,H_coeff_tot(4,:),'m.');
hold on;scatter(x,H_coeff_tot(5,:),'g.');
hold on;scatter(x,H_coeff_tot(6,:),'c.');
hold on;scatter(x,H_coeff_tot(7,:),'k+');
legend('AR 1 réel','AR 2 réel','AR 3 réel','AR 4 réel','AR 5 réel','AR 6 réel','AR 1 estimé','AR 2 estimé','AR 3 estimé','AR 4 estimé','AR 5 estimé','AR 6 estimé','AR 7 sur-estimé');

end

if disp_part_Kalman==1         %% Kalman %%
x=1:boucle_filtre;    
figure,
for i =1:ordre_p-surestimation
    hold on;plot([0 boucle_filtre],[AR_reel(i) AR_reel(i)],'LineWidth',2);
end
xlabel('itération');
ylabel('évolution estimation');
ylim([-2 2]);
title("estimation des coeff AR par filtre Kalman à chaque itération");
hold on;scatter(x,H_coeff_tot_Kalman(1,:),'b.');
hold on;scatter(x,H_coeff_tot_Kalman(2,:),'r.');
hold on;scatter(x,H_coeff_tot_Kalman(3,:),'y.');
hold on;scatter(x,H_coeff_tot_Kalman(4,:),'m.');
hold on;scatter(x,H_coeff_tot_Kalman(5,:),'g.');
hold on;scatter(x,H_coeff_tot_Kalman(6,:),'c.');
hold on;scatter(x,H_coeff_tot_Kalman(7,:),'r*');
hold on;scatter(x,H_coeff_tot_Kalman(8,:),'k+');
legend('AR 1 réel','AR 2 réel','AR 3 réel','AR 4 réel','AR 5 réel','AR 6 réel','AR 1 estimé','AR 2 estimé','AR 3 estimé','AR 4 estimé','AR 5 estimé','AR 6 estimé','AR 7 sur-estimé','AR 8 sur-estimé');

end