clc;
clear;
clear global;
global Q R;
% load('HPC_Sensorimotor_1e6_1.mat');
load('HPC_DoubleInt_1e4_1.mat');
% set the weighting matrices
Q = eye(n);
R = eye(m);
% Qp = [2000 -4000;-400 1000];
% Qv = [20 -1;-1 20];
% Qa = [0.01 0;0 0.01];
% Q = [Qp zeros(2) zeros(2);
%     zeros(2) Qv zeros(2);
%     zeros(2) zeros(2) Qa;];
% R = 0.01*eye(m);

% set parameters for learning
I = 10;
Ipe = 10;
P_tr = zeros(n,n,I);
theta_tr = zeros(n2,n2,I);
P_hat = zeros(n,n,I);
theta = zeros(n2,n2,I);
K_tr = zeros(m,n,I);
K_hat = zeros(m,n,I);
K_tr(:,:,1) = K;
K_hat(:,:,1) = K;
P_tilde(:,:,I) = zeros(n);
for i = 1:I-1
    P_tr(:,:,i) = lyap((A-B*K_tr(:,:,i))', Q + K_tr(:,:,i)'*R*K_tr(:,:,i));
    K_tr(:,:,i+1) = R\B'*P_tr(:,:,i);
    [tP, P_tmp] = ode45(@(t,y) PE_dyn(t,y,coef_est,K_hat(:,:,i))...
        ,[0,Ipe],sm2vec(zeros(n)));
    P_hat(:,:,i) = vec2sm(P_tmp(end,:),n);
    theta_tmp = vec2sm(coef_est*P_tmp(end,:)',n2);
    K_hat(:,:,i+1) = R\theta_tmp(n+1:end,1:n);
    P_tilde(:,:,i) = lyap((A-B*K_hat(:,:,i))', Q...
        + K_hat(:,:,i)'*R*K_hat(:,:,i));
end
i = i+1;
P_tr(:,:,i) = lyap((A-B*K_tr(:,:,i))', Q + K_tr(:,:,i)'*R*K_tr(:,:,i));
P_tilde(:,:,i) = lyap((A-B*K_hat(:,:,i))', Q...
    + K_hat(:,:,i)'*R*K_hat(:,:,i));

% optimal solution
[Pstar,~,Kstar] = care(A,B,Q,R);

save('result_HPC_DoubleInt_1e4_1.mat');