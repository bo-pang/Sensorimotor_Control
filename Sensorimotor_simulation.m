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
K_yu = [289.93 79.02 95.11 -17.78 2.56 -0.55;
    -243.53 376.3 -6.49 135.07 -0.55 2.86];

[n, m] = size(B);

q = length(G);
dt = 0.0001;
N = 12000;
t = dt*[0:N-1];
color = ['ymcbgk'];

fig1 = figure;
plot(0, 0, '.r', 'MarkerSize',50); hold on;
axis([-0.1 0.1 -0.3 0.1]);
fig2 = figure;
for j=1:6
    
    % check stability
    Mat_stab = stochastic_sys_mat(K_hat(:,:,j));
    if ~all(eig(Mat_stab)<0)
        disp("Instability");
    end
    
    rv1 = randn(q,N);
    x0 = [0 -0.25 0 0 0 0]';
    x = x0;
    for i=1:N-1
        u = -K_hat(:,:,j)*x(:,end);
        dx = (A*x(:,end)+B*u)*dt +...
            B*(G{1}*u*rv1(1,i)+G{2}*u*rv1(2,i))*sqrt(dt);
        x = [x, x(:,end)+dx];
    end
    figure(fig1);
    lg(j) = plot(x(1,:),x(2,:),color(j),'LineWidth',1);hold on;
    figure(fig2);
    subplot(2,2,1);
    plot(t,x(3,:),'b');hold on;
    title('X-velocity');
    xlim([0,0.6]);
    subplot(2,2,2);
    plot(t,x(4,:),'b');hold on;
    title('Y-velocity');
    xlim([0,0.6]);
    subplot(2,2,3);
    plot(t,x(5,:),'b');hold on;
    title('X-accelaration');
    xlim([0,0.6]);
    subplot(2,2,4);
    plot(t,x(6,:),'b');hold on;
    title('Y-accelaration');
    xlim([0,0.6]);
end
legend(lg,'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5','Trial 6');

% figure(4);
% for j=1:6
%     subplot(2,2,1);
%     plot(t,x(3,:));
% end

figure;
plot(0, 0, '.r', 'MarkerSize',50); hold on;
axis([-0.1 0.1 -0.3 0.1]);
dt = 0.0001;
for j=1:6
    
    rv1 = randn(q,N);
    x0 = [0 -0.25 0 0 0 0]';
    x = x0;
    for i=1:N-1
        u = -K_hat(:,:,end-j+1)*x(:,end);
%         u = -K_yu*x(:,end);
%         u = -Kstar*x(:,end);
        dx = (A*x(:,end)+B*u)*dt +...
            B*(G{1}*u*rv1(1,i)+G{2}*u*rv1(2,i))*sqrt(dt);
        x = [x, x(:,end)+dx];
    end
    lg(j) = plot(x(1,:),x(2,:),'b','LineWidth',1);hold on;
end
% legend(lg);

% After effects
figure;
% load('Model_NF');
K_NF = load('Kstar_approx_NF.mat');
K_NF = K_NF.Kstar_approx;
% Kstar = place(A,B,[1,-1,-2,-3,-4,-5]);
% Kstar = 100*Kstar;
plot(0, 0, '.r', 'MarkerSize',50); hold on;
% axis([-0.1 0.1 -0.3 0.1]);
axis([-0.05 0.05 -0.3 0.1]);
line([-0.025,-0.025],[-0.2,0.01],'color','k','LineWidth',1);
line([0.025,0.025],[-0.2,0.01],'color','k','LineWidth',1);
for j=1:6
    rv1 = randn(q,N);
    x0 = [0 -0.25 0 0 0 0]';
    x = x0;
    for i=1:N-1
%         u = -K_hat(:,:,end)*x(:,end);
        u = -K_NF*x(:,end);
        dx = (A*x(:,end)+B*u)*dt +...
            B*(G{1}*u*rv1(1,i)+G{2}*u*rv1(2,i))*sqrt(dt);
        xnext = x(:,end)+dx;
        if abs(xnext(1))>0.025 || xnext(2)>0.01
            break;
        else
            x = [x, xnext];
        end
    end
    lg(j) = plot(x(1,:),x(2,:),'b','LineWidth',1);hold on;
end
% legend(lg);
