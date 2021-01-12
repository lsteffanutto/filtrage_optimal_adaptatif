clear all
close all
clc;
beep off;
%% A RENDRE LE 22 JANVIER %%

%% VAR

disp=1;

N = 10000;%pow2(10);
f=-1/2:1/N:1/2-1/N;
sigma2=10;
ordre_p=4;
p1=0.95*exp(i*pi/5);
p2=0.9*exp(i*3*pi/4);
p3=0.95*exp(i*pi/2);

%% ALGO

%pôles obtenus
poles=[p1,p2,conj(p1),conj(p2)];

%nombre de fois ou on va changer l'argument des pôles
nb_chgmt=20;

%Création du processus AR variant dans le temps
[x_k,AR_tot]=processus_AR_var(poles,sigma2,N,nb_chgmt);

% ! Pourquoi un parametre AR est constant au cours du temps meme en faisant varier l'argument ! ?

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

% Bruit de modèle != 0 car parametre a estimer bouge
% constant
Q=eye(ordre_p)*0.000005;

% Bruit de mesure
sigma2 = 10;

% Iterations Kalman
boucle_filtre=N;

% On init le vecteur d'état
x=zeros(ordre_p,1);

%stocker les valeurs des paramètres AR à chaque itération
H_coeff_tot_Kalman=zeros(ordre_p,boucle_filtre);

for k = ordre_p+1:boucle_filtre
    [x,P] = Kalman_processus_AR(x,P,x_k(k),Q,sigma2,Phi,-x_k(k-1:-1:k-ordre_p));
    H_coeff_tot_Kalman(:,k)= x;
end

AR_Kalman=x




%% Figures
if disp==1      

% signal AR évoluant dans le temps
figure,
subplot 211,plot(x_k);
title('processus AR evoluant au cours du temps');
xlabel('échantillon');
ylabel('amplitude');
%spectro
fech=8000;
win  = 8;
spectro_noverlap = 0.5*win; 
subplot 212,spectrogram(x_k,win,spectro_noverlap,[],fech,'yaxis')
title('spectrogramme du processus AR généré');  

% paramètres AR évoluant au cours du temps et leurs estimation
x=1:N;    
figure,
for i=1:ordre_p
    hold on;plot(x,AR_tot(i,:),'LineWidth',2);
end

xlabel('itération');
ylabel('évolution estimation');
title('parametre estimé');

hold on;scatter(x,H_coeff_tot_Kalman(1,:),'b.');
hold on;scatter(x,H_coeff_tot_Kalman(2,:),'r.');
hold on;scatter(x,H_coeff_tot_Kalman(3,:),'y.');
hold on;scatter(x,H_coeff_tot_Kalman(4,:),'m.');
legend('AR 1 réel','AR 2 réel','AR 3 réel','AR 4 réel','AR 1 estimé','AR 2 estimé','AR 3 estimé','AR 4 estimé'); 
 
    
% x=1:boucle_filtre;    
% figure,
% for i =1:ordre_p-surestimation
%     hold on;plot([0 boucle_filtre],[AR_reel(i) AR_reel(i)],'LineWidth',2);
% end
% xlabel('itération');
% ylabel('évolution estimation');
% ylim([-2 2]);
% title("estimation des coeff AR par filtre Kalman à chaque itération");
% hold on;scatter(x,H_coeff_tot_Kalman(1,:),'b.');
% hold on;scatter(x,H_coeff_tot_Kalman(2,:),'r.');
% hold on;scatter(x,H_coeff_tot_Kalman(3,:),'y.');
% hold on;scatter(x,H_coeff_tot_Kalman(4,:),'m.');
% hold on;scatter(x,H_coeff_tot_Kalman(5,:),'g.');
% hold on;scatter(x,H_coeff_tot_Kalman(6,:),'c.');
% hold on;scatter(x,H_coeff_tot_Kalman(7,:),'r*');
% hold on;scatter(x,H_coeff_tot_Kalman(8,:),'k+');
% legend('AR 1 réel','AR 2 réel','AR 3 réel','AR 4 réel','AR 5 réel','AR 6 réel','AR 1 estimé','AR 2 estimé','AR 3 estimé','AR 4 estimé','AR 5 estimé','AR 6 estimé','AR 7 sur-estimé','AR 8 sur-estimé');

end