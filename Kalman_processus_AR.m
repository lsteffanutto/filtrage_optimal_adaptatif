function [x,P] = Kalman_processus_AR(x,P,y,Q,R,Phi,H)
        
%ajouter G comme dans pdf ? maisbalek
%     [Phi] = compagnon_matrix(Phi)
    
    %% Etape de prediction
    x = Phi*x;          % pas d'étape de prédiction que de crorrection on lui dit juste c les valeurs d'avant
    P = Phi*P*Phi' + Q;
%     Mat_P(:,1:2) = P;
    
    %% Etape de Correction
    K = (H*P*H' + R)\(P*H'); % gain     
    x = x + K*(y-H*x);
    P = P - K*H*P; 
%     Mat_P(:,3:4) = P;
    

       
end