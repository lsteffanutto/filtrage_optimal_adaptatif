function [H] = f_LMS(H,alpha,x,d) 
%% INPUT => Notations cours p.1 PDF2
% H = coefficient AR du filtre
% alpha = paramètre contrôlant la convergence de l'algo ( CV si 0<alpha<2/tr(Rxx))
% x = vecteur avec les p last échantillons du processus AR en input (~ x(i-ordre_p+1:i) )
% d = sortie désirée à l'instant k = le vrai processus AR à l'instant k

%% OUTPUT
% H_maj = Filtre adaptatif mis à jour à l'instant k en fonction de l'erreur calculée
% => même nom de variable pour pouvoir boucler

%% ALGO


% On calcul la sortie réelle
y=H'*x';

% Erreur entre sortie réelle et désirée
e = d - y;

% Mise à jour des coefficients du filtre
H = H + alpha*x'*e;

%% Explications

% chaque instant on prend les 6 derniers échantillons du signal d'entrée
% i.e le bb, est x_k(i) qui est le signal désiré à cette instant
end

