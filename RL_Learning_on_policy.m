clc;
clear;
clear global;
global A B G Q R n2 n3 n4 m K p q n ww mag;
% mass = 1;
% b = 10;
% tau = 0.05;
% c1 = 0.075;
% c2 = 0.025;
% 
% ff = 0.7*[13 -18;18 13];
% A = [zeros(2) eye(2) zeros(2);
%     zeros(2) ff/mass-b/mass*eye(2) eye(2)/mass;
%     zeros(2) zeros(2) -eye(2)/tau;];
% B = [zeros(2);zeros(2);eye(2)/tau];
% G{1} = [c1 0;
%     c2 0;];
% G{2} = [0 c2;
%     0 c1];
A = [0 1;-1 2];
B = [0;1];
G{1} = [0.1];
% rng('shuffle');
rng(1);

q = length(G);
[n, m] = size(B);
n2 = m+n;            % dimension of z
n3 = n2*(n2+1)/2;      % dimension of z_tilt
n4 = n3*(n3+1)/2;      % dimension of z_tilt*z_tilt

% exploration noise
ww = -250 + 500*randn(m,100);
mag = 1;

% initial control gain
% K = place(A,B,[-1, -2, -3, -4, -5,-6]);
K = [11,9];
% set the weighting matrices
Q = eye(n);
R = eye(m);

% set parameters for learning
I = 5;
Ipe = 100;
P_tr = zeros(n,n,I);
theta_tr = zeros(n2,n2,I);
P_hat = zeros(n,n,I);
theta = zeros(n2,n2,I);
K_tr = zeros(m,n,I);
K_hat = zeros(m,n,I);
K_tr(:,:,1) = K;
K_hat(:,:,1) = K;
P_tilde(:,:,I) = zeros(n);
% parameters for policy evaluation
M = 7;
tstep = 0.05;
x0 = ones(n,1);
N = 1e4;
dt = 0.001;
for i = 1:I-1
    P_tr(:,:,i) = lyap((A-B*K_tr(:,:,i))', Q + K_tr(:,:,i)'*R*K_tr(:,:,i));
    K_tr(:,:,i+1) = R\B'*P_tr(:,:,i);
    [coef_est, err_coef] = collect_data(K_hat(:,:,i),x0,tstep,dt,N,M);
    [tP, P_tmp] = ode45(@(t,y) PE_dyn(t,y,coef_est,K_hat(:,:,i))...
        ,[0,Ipe],sm2vec(zeros(n)));
    P_hat(:,:,i) = vec2sm(P_tmp(end,:),n);
    theta_tmp = vec2sm(coef_est*P_tmp(end,:)',n2);
    K_hat(:,:,i+1) = R\theta_tmp(n+1:end,1:n);
    P_tilde(:,:,i) = lyap((A-B*K_hat(:,:,i))', Q...
        + K_hat(:,:,i)'*R*K_hat(:,:,i));
    disp(['Iteration ',num2str(i)]);
end
i = i+1;
P_tr(:,:,i) = lyap((A-B*K_tr(:,:,i))', Q + K_tr(:,:,i)'*R*K_tr(:,:,i));
P_tilde(:,:,i) = lyap((A-B*K_hat(:,:,i))', Q...
    + K_hat(:,:,i)'*R*K_hat(:,:,i));

% optimal solution
[Pstar,~,Kstar] = care(A,B,Q,R);