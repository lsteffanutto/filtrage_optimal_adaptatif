function [x,P,Mat_X,Mat_P] = Kalman_processus_AR(x,P,y,Q,R,Phi,H)
        
%     Mat_X = zeros(ordre,2);
      Mat_X=1;
      Mat_P=1;
%     Mat_P = zeros(2,4);
    
%ajouter G comme dans pdf ? maisbalek

    %% Etape de prediction
    x = Phi*x; 
    Mat_X(:,1) = x;
    P = Phi*P*Phi' + Q;
%     Mat_P(:,1:2) = P;
    
    %% Etape de Correction
    K = (H*P*H' + R)\(P*H'); % gain     
    x = x + K*(y-H*x);
    Mat_X(:,2) = x;
    P = P - K*H*P; 
%     Mat_P(:,3:4) = P;
    

       
end