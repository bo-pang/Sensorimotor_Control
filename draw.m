clc;
clear;

load('result_HPC_DoubleInt_1e4_1.mat');
rel_err_K = zeros(I,1);
rel_err_P = zeros(I,1);
for i = 1:I
    rel_err_K(i) = norm(K_hat(:,:,i)-Kstar,'fro')/norm(Kstar,'fro');
    rel_err_P(i) = norm(P_tilde(:,:,i)-Pstar,'fro')/norm(Pstar,'fro');
end
figure(1);
plot(1:I,rel_err_K,'--x');hold on;
plot(1:I,rel_err_P,'--o');hold on;
xlabel('Iteration Index');
ylabel('Relative Error');
legend({'$\Vert \hat{K}_i-K^*\Vert_F/\Vert K^*\Vert_F$',...
    '$\Vert \hat{P}_i-P^*\Vert_F/\Vert P^*\Vert_F$'},'Interpreter','latex');