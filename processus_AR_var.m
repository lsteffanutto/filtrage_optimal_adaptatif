function [x_k,AR_tot] = processus_AR_var(poles,sigma2,N,nb_chgmt)
% INPUT=pôles, ecart type, nb_ech, nb_changmt (10000 => 1000)
% OUPUT= processus AR variant dans le temps

%init le signal à fabriquer
x_k=zeros(1,N);

%nb paramètres AR
nb_pole=length(poles)

%taille d'une trame
taille_trames=N/nb_chgmt

%stock les paramètres AR qui varient au cours du temps
AR_tot=zeros(nb_pole,N);

% on va changer 20 fois par marche aleatoire les parametres AR
for i = 1:nb_chgmt 
    
    % Signal généré avec avec les pôles actuels
    
    % retourne les coeff AR à partir des pôles
    A = poly(poles)'   
    
    B = 1;
    
    %processus générateur
    bb=sqrt(sigma2)*randn(1,taille_trames); 
    
    %processus AR à partir des poles
    x_k(((i-1)*taille_trames)+1:i*taille_trames) = filter(B,A,bb); 
    
    %stockage des paramètres AR qui ne varient pas sur la longueur de la trame
    AR_tot(:,((i-1)*taille_trames)+1:i*taille_trames)=A(2:5).*ones(nb_pole,taille_trames);
    
    %% Changement arguments des pôles
    
    % créer la marche aléatoire à appliquer sur les arguments des pôles
    a = -0.05;
    b = 0.05;
    delta_arg = (b-a).*rand(1,nb_pole/2) + a;
    delta_arg=delta_arg
    

    % récupère les arguments actuels des pôles et ajout d'un incrément
    arg1=angle(poles(1))+delta_arg(1);
    arg2=angle(poles(2))+delta_arg(2);
    
    % mets à jour les pôles
    p1=abs(poles(1)).*exp(1i*arg1);
    p2=abs(poles(2)).*exp(1i*arg2);
    
    poles=[p1,p2,conj(p1),conj(p2)]'
    
    
end


end

