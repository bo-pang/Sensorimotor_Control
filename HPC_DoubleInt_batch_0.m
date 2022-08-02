clear;
clear global;
clc;

NN = 1;

% tstart = tic;

for ii=1:NN
    run('HPC_DoubleInt_0.m');
%     save(['HPC_DoubleInt_1e',num2str(log10(N)),'_',num2str(ii),'.mat'],'A','B','G',...
%         'N','dt','tstep','M','tspan','Psi','Phi','Psi_est','Phi_est',...
%         'rel_err','theta_P_v_tr','theta_P_v','Xdis','X','ww','Xtr');
    save(['HPC_DoubleInt_1e',num2str(log10(N)),'_',num2str(ii),'.mat']);
end