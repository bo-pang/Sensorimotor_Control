% Try small gain
clear;
clc;
global A B G C Q R n2 n3 n4 n m K p q c3;
mass = 1.3;
b = 10;
tau = 0.05;
c1 = 0.00075;
c2 = 0.00075;
c3 = 0.01;

A = [zeros(2) eye(2) zeros(2);
    zeros(2) -b/mass*eye(2) eye(2)/mass;
    zeros(2) zeros(2) -eye(2)/tau;];
B = [zeros(2);zeros(2);eye(2)/tau];
G{1} = [c1 0;
    c2 0;];
G{2} = [0 c2;
    0 c1];

q = length(G);
[n, m] = size(B);

C = c3*eye(n);
[~,p] = size(C);
Q = 1.1*eye(n);
R = 1.1*eye(m);

K = place(A,B,[-1, -2, -3, -4, -5,-6]);

n2 = m+n+1;            % dimension of z
n3 = n2*(n2+1)/2;      % dimension of z_tilt
n4 = n3*(n3+1)/2;      % dimension of z_tilt*z_tilt

Kest = {K};
Gest = {};
kmax = 1e4;
k = 1;
t = 10;
dt = 0.001;
tstep = 0.05;
x_his = [];

% Collect data
while k<kmax
    if k==1
        xaug0 = [0.1;0.2;0.3;0.4;0.5;0.6;0;0;];
        y0 = zeros(n3+n4,1);
        X0 = [xaug0;y0];
        [xaug, I_Psi_xx, Psi_Psi, Psi_r] = EM(@mydrift,@mydiffusion,X0,dt,0,t);
    else
        t = t + tstep;
        [xaug, d_I_Psi_xx, d_Psi_Psi, d_Psi_r] = EM(@mydrift,@mydiffusion,X,dt,t-tstep,t);
        I_Psi_xx = I_Psi_xx + d_I_Psi_xx;
        Psi_Psi = Psi_Psi + d_Psi_Psi;
        Psi_r = Psi_r + d_Psi_r;
    end

    x_his = [xaug,x_his];
    X = [xaug;y0];
    rk = rank(Psi_Psi/t);
%     if rk ~= n3
%         disp('Rank defficient');
%         break;
%     end
    k = k + 1;
end
 save('/scratch/bp1471/mytest1/Sensorimotor_hpc_data_7.mat','A','B','C','G',...
        'Q','R','kmax','t','tstep','K','dt','I_Psi_xx','Psi_Psi','Psi_r','x_his');

