function [Phi] = compagnon_matrix(AR_reel)
%% Retourne la matrice compagnon avec les paramètres AR souhaité pour
% obtenir une matrice de transition Phi souhaité PDF3 p.17
ordre_p=size(AR_reel,1);
u = [1 AR_reel'];
Phi = compan(u)';
Phi(:,1)=zeros(ordre_p,1);
Phi(end,:)= -(AR_reel(end:-1:1)');

end

