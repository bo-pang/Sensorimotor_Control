clc;
clear;
clear global;
global A B G Q R n2 n3 n4 m K p q n ww mag;
% mass = 1.3;
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
% A = [0 1;-1 2];
% B = [0;1];
% G{1} = [0.1];
rng('shuffle');
% rng(1);
load('Model_DF');
% load('Model_VF');
% load('Model_NF');

q = length(G);
[n, m] = size(B);
n2 = m+n;            % dimension of z
n3 = n2*(n2+1)/2;      % dimension of z_tilt
n4 = n3*(n3+1)/2;      % dimension of z_tilt*z_tilt

% exploration noise
ww = -100 + 200*randn(m,2000);
mag = 2000;
% mag = 10;

% initial control gain
% K = place(A,B,[-6, -9, -20, -21, -22,-23]);
% K = [11,9];
% set the weighting matrices
% Q = eye(n);
% R = eye(m);

% Weight matrices inspired by Tao, 2016, TAC
% Qp = [200 -40;-40 1000];
% Qv = [20 -1;-1 20];
% Qa = [0.01 0;0 0.01];
% Q = [Qp zeros(2) zeros(2);
%     zeros(2) Qv zeros(2);
%     zeros(2) zeros(2) Qa;];
% R = 0.01*eye(m);

% Weight matrices in by Tao, 2016, TAC
Qp = [2000 -40;-40 1000];
Qv = [20 -1;-1 20];
Qa = [0.01 0;0 0.01];
Q = [Qp zeros(2) zeros(2);
    zeros(2) Qv zeros(2);
    zeros(2) zeros(2) Qa;];
R = 0.01*eye(m);

% Weight matrices in Yu
% dm = 0.15;
% a0 = 5e4;
% a1 = 1e5;
% a2 = 1e6;
% % a2 = 0;
% angle = [cos(deg2rad(144));sin(deg2rad(144))];
% Q0 = [a0 0 ; 0 a1;] + a2*dm*(angle*angle');
% Q = blkdiag(Q0,1e-2*Q0,5e-5*Q0);
% R = eye(2);


% optimal solution
[Pstar,~,Kstar] = care(A,B,Q,R);
% K = Kstar;

% set parameters for learning
I = 30;
Ipe = 200;
P_tr = zeros(n,n,I);
theta_tr = zeros(n2,n2,I);
P_hat = zeros(n,n,I);
theta = zeros(n2,n2,I);
K_tr = zeros(m,n,I);
K_hat = zeros(m,n,I);
P_tilde = zeros(n,n,I);
DeltaG = zeros(I,1);
G_tilde = zeros(m+n,m+n,I);
err_G = zeros(I,1);
% K = [100 0 10 0 10 0;
%     0 100 0 10 0 10;];
K = load('Kstar_approx_NF');
% K = 3*K.Kstar_approx;       % initial admissible gain for VF
K = K.Kstar_approx;
K(1,1) = K(1,1) + 600;
K_tr(:,:,1) = K;
K_hat(:,:,1) = K;
% parameters for policy evaluation
M = 50;
tstep = 0.05;
tspan = [0:M]*tstep;
x0 = ones(n,1);
xx0 = x0*x0';
z0 = [x0;xx0(:);zeros(n3,1)];
dt = 0.001;
mag_perturb_1 = 1;
mag_perturb_2 = 1;
% mag_perturb = 0.00001;
for i = 1:I-1
    P_tr(:,:,i) = lyap((A-B*K_tr(:,:,i))', Q + K_tr(:,:,i)'*R*K_tr(:,:,i));
    K_tr(:,:,i+1) = R\B'*P_tr(:,:,i);
    % Check stability
    if ~all(eig(A-B*K_hat(:,:,i))<0)
       disp('Control gain is not stabilizing!'); 
    end
    % Check the rank condition
    [tz,z] = ode45(@(t,z) exp_dyn(t,z,K_hat(:,:,i)),tspan,z0);
    Xtr = zeros(n*(n+1)/2,M+1);
    Ztr = zeros(n3,M+1);
    for j=1:M+1
        Xtr(:,j) = sm2vec(reshape(z(j,n+1:n^2+n),[n,n]));
        Ztr(:,j) = z(j,n+n^2+1:end)';
    end
    Phi = Ztr(:,2:end)'-Ztr(:,1:end-1)';
%     delta_Phi = reshape(normrnd(0,mag_perturb,[n3*M,1]),[M,n3]);
%     Phi_est = Phi + delta_Phi;
    if rank(Phi) ~= n3
        disp('Rank defficient');
    end
    Psi = Xtr(:,2:end)'-Xtr(:,1:end-1)';
%     delta_Psi = reshape(normrnd(0,mag_perturb,[n*(n+1)/2*M,1]),[M,n*(n+1)/2]);
%     Psi_est = Psi + delta_Psi;
    P = randn([n,n]);
    P = (P+P')/2;
    theta_P_v = Phi\Psi*sm2vec(P);
    Pi = zeros(m);
    for j=1:q
        Pi = Pi + G{j}'*B'*P*B*G{j};
    end
    theta_P_v_tr = sm2vec([A'*P + P*A, P*B;
        B'*P, Pi;]);
%     coef = Phi\Psi;
%     delta_coef = reshape(normrnd(0,mag_perturb,[n3*n*(n+1)/2,1])...
%         ,[n3,n*(n+1)/2]);
%     delta_coef = -mag_perturb + 2*mag_perturb*rand([n3,n*(n+1)/2]);
%     coef_est = coef + delta_coef;
    Phi_perturb = -mag_perturb_1 + 2*mag_perturb_1*rand([n3,M]);
    Psi_perturb = -mag_perturb_2 + 2*mag_perturb_2*rand([n*(n+1)/2,M]);
    Phi_hat = Phi + Phi_perturb';
    Psi_hat = Psi + Psi_perturb';
    coef = Phi\Psi;
    coef_est = Phi_hat\Psi_hat;
%     coef_est = coef;
    [tP, P_tmp] = ode45(@(t,y) PE_dyn(t,y,coef_est,K_hat(:,:,i))...
        ,[0,Ipe],sm2vec(zeros(n)));
    P_hat(:,:,i) = vec2sm(P_tmp(end,:),n);
    theta_tmp = vec2sm(coef_est*P_tmp(end,:)',n2);
    K_hat(:,:,i+1) = R\theta_tmp(n+1:end,1:n);
    P_tilde(:,:,i) = lyap((A-B*K_hat(:,:,i))', Q...
        + K_hat(:,:,i)'*R*K_hat(:,:,i));
    DeltaG(i) = norm(theta_tmp(1:n,1:n) - A'*P_tilde(:,:,i)...
        - P_tilde(:,:,i)*A,'fro') + 2*norm(theta_tmp(n+1:end,1:n)...
        -B'*P_tilde(:,:,i),'fro');
    G_tilde(:,:,i) = [ A'*P_tilde(:,:,i)...
        - P_tilde(:,:,i)*A + Q, (B'*P_tilde(:,:,i))';
        B'*P_tilde(:,:,i), R;];
    err_G(i) = DeltaG(i)/norm(G_tilde(:,:,i));
    disp(['Iteration ',num2str(i)]);
end
i = i+1;
P_tr(:,:,i) = lyap((A-B*K_tr(:,:,i))', Q + K_tr(:,:,i)'*R*K_tr(:,:,i));
P_tilde(:,:,i) = lyap((A-B*K_hat(:,:,i))', Q...
    + K_hat(:,:,i)'*R*K_hat(:,:,i));
DeltaG(i) = norm(theta_tmp(1:n,1:n) - A'*P_tilde(:,:,i)...
    - P_tilde(:,:,i)*A,'fro') + 2*norm(theta_tmp(n+1:end,1:n)...
    -B'*P_tilde(:,:,i),'fro');
G_tilde(:,:,i) = [ A'*P_tilde(:,:,i)...
    - P_tilde(:,:,i)*A + Q, (B'*P_tilde(:,:,i))';
    B'*P_tilde(:,:,i), R;];
err_G(i) = DeltaG(i)/norm(G_tilde(:,:,i));
save('Sensorimotor_mimic_DF.mat');