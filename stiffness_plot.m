clc;
clear;
load('Kstar_approx_NF');
Kopt_NF = Kstar_approx(1:2,1:2);
load('Sensorimotor_mimic_VF');
Kopt_VF = K_hat(1:2,1:2,end);
load('Sensorimotor_mimic_DF');
Kopt_DF = K_hat(1:2,1:2,end);

figure;
stiffness(Kopt_NF,'g');
stiffness(Kopt_VF,'r');
xlim([-800,800]);
ylim([-800,800]);
xlabel('$x$ component of stiffness $(N\cdot m^{-1})$','Interpreter','latex');
ylabel('$y$ component of stiffness $(N\cdot m^{-1})$','Interpreter','latex');

figure;
stiffness(Kopt_NF,'g');
stiffness(Kopt_DF,'r');
xlim([-2000,2000]);
ylim([-800,800]);
xlabel('$x$ component of stiffness $(N\cdot m^{-1})$','Interpreter','latex');
ylabel('$y$ component of stiffness $(N\cdot m^{-1})$','Interpreter','latex');