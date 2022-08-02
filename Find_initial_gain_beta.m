clear;
clc;

nominal = load('Model_DF.mat');

m = ureal('m',2,'PlusMinus',[-1.5,3]);
b = ureal('b',15,'PlusMinus',[-13,5]);
tau = ureal('tau',0.3,'PlusMinus',[-0.28,0.2]);
beta = ureal('beta',400,'PlusMinus',[-200,100]);
% m = 1.3;
% b = 10;
% tau = 0.05;
% beta = 300;

A = [0 0 1 0 0 0;
    0 0 0 1 0 0;
    beta/m 0 -b/m 0 1/m 0;
    0 0 0 -b/m 0 1/m;
    0 0 0 0 -1/tau 0;
    0 0 0 0 0 -1/tau;];

B = [0 0;
    0 0;
    0 0;
    0 0;
    1/tau 0;
    0 1/tau;];

C = eye(6);

D = zeros(6);

k11 = realp('k11',0);
k12 = realp('k12',0);
k13 = realp('k13',0);
k14 = realp('k14',0);
k15 = realp('k15',0);
k16 = realp('k16',0);

k21 = realp('k21',0);
k22 = realp('k22',0);
k23 = realp('k23',0);
k24 = realp('k24',0);
k25 = realp('k25',0);
k26 = realp('k26',0);

K = [k11 k12 k13 k14 k15 k16;
    k21 k22 k23 k24 k25 k26];
    
aug_C = [C;-K];
aug_D = [D;zeros(2,6)];

closed_sys = ss(A-B*K,C,aug_C,aug_D);
closed_sys.name = 'hand';

closed_sys.InputName = {'w1','w2','w3','w4','w5','w6'};

closed_sys.OutputName = {'x1','x2','x3','x4','x5','x6','u1','u2'};

Qz = blkdiag(eye(6),1e2*eye(2));
R = TuningGoal.LQG(closed_sys.InputName,{'x1','x2','x3','x4','x5','x6','u1','u2'},1,Qz);
[rs_sys,~] = systune(closed_sys,R);

k11_val = getBlockValue(rs_sys,'k11');
k12_val = getBlockValue(rs_sys,'k12');
k13_val = getBlockValue(rs_sys,'k13');
k14_val = getBlockValue(rs_sys,'k14');
k15_val = getBlockValue(rs_sys,'k15');
k16_val = getBlockValue(rs_sys,'k16');

k21_val = getBlockValue(rs_sys,'k21');
k22_val = getBlockValue(rs_sys,'k22');
k23_val = getBlockValue(rs_sys,'k23');
k24_val = getBlockValue(rs_sys,'k24');
k25_val = getBlockValue(rs_sys,'k25');
k26_val = getBlockValue(rs_sys,'k26');

Kval = [k11_val k12_val k13_val k14_val k15_val k16_val;
    k21_val k22_val k23_val k24_val k25_val k26_val;];

EIGEN = eig(nominal.A - nominal.B*Kval);
EIGEN_open = eig(nominal.A);
if all(EIGEN_open<0)
    disp('System stable before search!');
else
    disp('System unstable before search!');
end
if all(EIGEN<0)
    disp('Stabilizing gain found!');
else
    disp('Failed to find stabilizing gain!');
end