function [alpha] = alpha_optimal_LMS(signal)
%% calcul les vp d'un signal, la vp_max pour trouver le pas optimal du LMS
%% ALGO

% Rxx du signal
z=autocorr(signal);
auto_corr_matrix=toeplitz(z);       

% SVD et vp max de Rxx
[~,lambda] = eig(auto_corr_matrix); 
lambda_max=max(max(lambda));       

% 0<alpha<2/lambda_max
borne_sup_alpha = 2/lambda_max;

% On retourne soit le alpha max soit un plus petit
alpha=borne_sup_alpha/2;

end

