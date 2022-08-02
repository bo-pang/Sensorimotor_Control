clc;
clear;
clear global;
global A B G C n m;

load("Sensorimotor_mimic_VF.mat");
% G{1} = G{1}/tau;
% G{2} = G{2}/tau;
% K0 = [289.93 79.02 95.11 -17.78 2.56 -0.55;
% -243.53 376.30 -6.49 135.07 -0.55 2.86];
% K0 = [293.88 79.91 95.92 -17.92 2.59 -0.56;
%     -248.62 381.36 -7.17 136.41 -0.56 2.89];
% K0 = Kstar;
% K_yu = [289.93 79.02 95.11 -17.78 2.56 -0.55;
%     -243.53 376.3 -6.49 135.07 -0.55 2.86];
load('Kstar_approx_NF');

[n, m] = size(B);

q = length(G);
dt = 0.0001;
N = 18000;
t = dt*[0:N-1];
color = ['mcbgk'];

fig1 = figure;
fig2 = figure;

% 5 trials in NF
load('Model_NF');
for j=1:5
    
    % check stability
    Mat_stab = stochastic_sys_mat(Kstar_approx);
    if ~all(eig(Mat_stab)<0)
        disp("Instability");
    end
    
    % Generate data
    rv1 = randn(q,N);
    x0 = [0 -0.25 0 0 0 0]';
    x = x0;
    for i=1:N-1
        u = -Kstar_approx*x(:,end);
        dx = (A*x(:,end)+B*u)*dt +...
            B*(G{1}*u*rv1(1,i)+G{2}*u*rv1(2,i))*sqrt(dt);
        x = [x, x(:,end)+dx];
    end
    
    % plot 2-D trajectory
    figure(fig1);
    subplot(1,4,1);
    plot(0, 0, '.r', 'MarkerSize',50); hold on;
    axis([-0.1 0.1 -0.3 0.1]);
    lg(j) = plot(x(1,:),x(2,:),'b','LineWidth',1);hold on;
    title('A');
    xlabel('X-position (m)');
    ylabel('Y-position (m)');
    
    % plot trajectoris of components
    figure(fig2);
    
    subplot(4,4,1);
    plot(t,x(3,:),'b');hold on;
    xlabel('time (s)');
    ylabel('X-velocity (m/s)');
    title('A');
    xlim([0,0.6]);
    ylim([-0.5,0.5]);
    
    subplot(4,4,5);
    plot(t,x(4,:),'b');hold on;
    xlabel('time (s)');
    ylabel('Y-velocity (m/s)');
    xlim([0,0.6]);
    ylim([-0.5,1.5]);
    
    subplot(4,4,9);
    plot(t,x(5,:),'b');hold on;
    xlabel('time (s)');
    ylabel('X-endpoint force (N)');
    xlim([0,0.6]);
    
    subplot(4,4,13);
    plot(t,x(6,:),'b');hold on;
    xlabel('time (s)');
    ylabel('Y-endpoint force (N)');
    xlim([0,0.6]);
end


% learning in VF
load('Model_VF');

% check stability
Mat_stab = stochastic_sys_mat(Kstar_approx);
if ~all(eig(Mat_stab)<0)
    disp("Instability");
end

% Generate data
rv1 = randn(q,N);
x0 = [0 -0.25 0 0 0 0]';
x = x0;
for i=1:N-1
    u = -Kstar_approx*x(:,end);
    dx = (A*x(:,end)+B*u)*dt +...
        B*(G{1}*u*rv1(1,i)+G{2}*u*rv1(2,i))*sqrt(dt);
    x = [x, x(:,end)+dx];
end

% plot 2-D trajectory
figure(fig1);
subplot(1,4,2);
plot(0, 0, '.r', 'MarkerSize',50); hold on;
axis([-0.12 0.07 -0.3 0.1]);
lg(j) = plot(x(1,:),x(2,:),color(1),'LineWidth',1);hold on;
title('B');
xlabel('X-position (m)');

% plot trajectoris of components
figure(fig2);

subplot(4,4,2);
plot(t,x(3,:),'b');hold on;
xlabel('time (s)');
title('B');
xlim([0,0.6]);
ylim([-1,0.5]);

subplot(4,4,6);
plot(t,x(4,:),'b');hold on;
xlabel('time (s)');
xlim([0,0.6]);
ylim([-1,1.5]);

subplot(4,4,10);
plot(t,x(5,:),'b');hold on;
xlabel('time (s)');
xlim([0,0.6]);

subplot(4,4,14);
plot(t,x(6,:),'b');hold on;
xlabel('time (s)');
xlim([0,0.6]);
    
for j=1:4
    
    % check stability
    Mat_stab = stochastic_sys_mat(K_hat(:,:,j));
    if ~all(eig(Mat_stab)<0)
        disp("Instability");
    end

    % Generate data
    rv1 = randn(q,N);
    x0 = [0 -0.25 0 0 0 0]';
    x = x0;
    for i=1:N-1
        u = -K_hat(:,:,j)*x(:,end);
        dx = (A*x(:,end)+B*u)*dt +...
            B*(G{1}*u*rv1(1,i)+G{2}*u*rv1(2,i))*sqrt(dt);
        x = [x, x(:,end)+dx];
    end

    % plot 2-D trajectory
    figure(fig1);
    subplot(1,4,2);
    plot(0, 0, '.r', 'MarkerSize',50); hold on;
    axis([-0.13 0.07 -0.3 0.1]);
    lg(j) = plot(x(1,:),x(2,:),color(j+1),'LineWidth',1);hold on;

    % plot trajectoris of components
    figure(fig2);
    
    subplot(4,4,2);
    plot(t,x(3,:),'b');hold on;
    xlim([0,0.6]);
    ylim([-1,0.5]);
    
    subplot(4,4,6);
    plot(t,x(4,:),'b');hold on;
    xlim([0,0.6]);
    ylim([-0.5,1.5]);
    
    subplot(4,4,10);
    plot(t,x(5,:),'b');hold on;
    xlim([0,0.6]);
    
    subplot(4,4,14);
    plot(t,x(6,:),'b');hold on;
    xlim([0,0.6]);

end


% Trials after learning in VF
for j=1:5
    
    % check stability
    Mat_stab = stochastic_sys_mat(K_hat(:,:,end));
    if ~all(eig(Mat_stab)<0)
        disp("Instability");
    end

    % Generate data
    rv1 = randn(q,N);
    x0 = [0 -0.25 0 0 0 0]';
    x = x0;
    for i=1:N-1
        u = -K_hat(:,:,end)*x(:,end);
        dx = (A*x(:,end)+B*u)*dt +...
            B*(G{1}*u*rv1(1,i)+G{2}*u*rv1(2,i))*sqrt(dt);
        x = [x, x(:,end)+dx];
    end

    % plot 2-D trajectory
    figure(fig1);
    subplot(1,4,3);
    plot(0, 0, '.r', 'MarkerSize',50); hold on;
    axis([-0.1 0.1 -0.3 0.1]);
    lg(j) = plot(x(1,:),x(2,:),'b','LineWidth',1);hold on;
    title('C');
    xlabel('X-position (m)');

    % plot trajectoris of components
    figure(fig2);
    
    subplot(4,4,3);
    plot(t,x(3,:),'b');hold on;
    xlabel('time (s)');
    title('C');
    xlim([0,0.6]);
    ylim([-0.5,0.5]);
    
    subplot(4,4,7);
    plot(t,x(4,:),'b');hold on;
    xlabel('time (s)');
    xlim([0,0.6]);
    ylim([-0.5,1.5]);
    
    subplot(4,4,11);
    plot(t,x(5,:),'b');hold on;
    xlabel('time (s)');
    xlim([0,0.6]);
    
    subplot(4,4,15);
    plot(t,x(6,:),'b');hold on;
    xlabel('time (s)');
    xlim([0,0.6]);

end

% After effects
load('Model_NF');
for j=1:6
    
    % check stability
    Mat_stab = stochastic_sys_mat(K_hat(:,:,end));
    if ~all(eig(Mat_stab)<0)
        disp("Instability");
    end

    % Generate data
    rv1 = randn(q,N);
    x0 = [0 -0.25 0 0 0 0]';
    x = x0;
    for i=1:N-1
        u = -K_hat(:,:,end)*x(:,end);
        dx = (A*x(:,end)+B*u)*dt +...
            B*(G{1}*u*rv1(1,i)+G{2}*u*rv1(2,i))*sqrt(dt);
        x = [x, x(:,end)+dx];
    end

    % plot 2-D trajectory
    figure(fig1);
    subplot(1,4,4);
    plot(0, 0, '.r', 'MarkerSize',50); hold on;
    axis([-0.1 0.1 -0.3 0.1]);
    lg(j) = plot(x(1,:),x(2,:),'b','LineWidth',1);hold on;
    title('D');
    xlabel('X-position (m)');

    % plot trajectoris of components
    figure(fig2);
    
    subplot(4,4,4);
    plot(t,x(3,:),'b');hold on;
    xlabel('time (s)');
    title('D');
    xlim([0,0.6]);
    ylim([-0.5,0.5]);
    
    subplot(4,4,8);
    plot(t,x(4,:),'b');hold on;
    xlabel('time (s)');
    xlim([0,0.6]);
    ylim([-0.5,1.5]);
    
    subplot(4,4,12);
    plot(t,x(5,:),'b');hold on;
    xlabel('time (s)');
    xlim([0,0.6]);
    
    subplot(4,4,16);
    plot(t,x(6,:),'b');hold on;
    xlabel('time (s)');
    xlim([0,0.6]);
    
end

figure;
plot(1:I,err_G,'b--o');
xlabel('Iteration');
ylabel('$\Vert \Delta G_i\Vert_F/\Vert G(\hat{P}_i)\Vert_F$',...
    'Interpreter','latex');
xlim([1,30]);
ylim([0,0.2]);
