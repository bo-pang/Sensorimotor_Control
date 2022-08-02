function [coef_est, err_coef] = PE_func(K,x0,tstep,dt,N,M)
global A B G Q R n2 n3 n4 m p q n ww mag;
L1 = stochastic_sys_mat(K);
if ~all(eig(L1)<0)
   disp('Initial control gain is not stabilizing!'); 
end

% Check the rank condition
tspan = [0:M]*tstep;
xx0 = x0*x0';
z0 = [x0;xx0(:);zeros(n3,1)];
[tz,z] = ode45(@(t,z) exp_dyn(t,z,K),tspan,z0);
Xtr = zeros(n*(n+1)/2,M+1);
Ztr = zeros(n3,M+1);
for i=1:M+1
    Xtr(:,i) = sm2vec(reshape(z(i,n+1:n^2+n),[n,n]));
    Ztr(:,i) = z(i,n+n^2+1:end)';
end
Phi = Ztr(:,2:end)'-Ztr(:,1:end-1)';
if rank(Phi) ~= n3
    disp('Rank defficient');
end
Psi = Xtr(:,2:end)'-Xtr(:,1:end-1)';
P = randn([n,n]);
P = (P+P')/2;
theta_P_v = Phi\Psi*sm2vec(P);
Pi = zeros(m);
for i=1:q
    Pi = Pi + G{i}'*B'*P*B*G{i};
end
theta_P_v_tr = sm2vec([A'*P + P*A, P*B;
    B'*P, Pi;]);

% Qp = [2000 -4000;-400 1000];
% Qv = [20 -1;-1 20];
% Qa = [0.01 0;0 0.01];
% Q = [Qp zeros(2) zeros(2);
%     zeros(2) Qv zeros(2);
%     zeros(2) zeros(2) Qa;];
% R = 0.01*eye(m);


p = 1;
t = 0;
X = zeros(n*(n+1)/2,M+1);
Z = zeros(n3,M+1);
Nt = floor((tspan(end)-t)/dt);
sNt = floor(tstep/dt);

% Collect data
for p=1:N
    x = x0;
    intZ = zeros(n3,1);
    j = 1;
    t = 0;
    for k=1:Nt
        rv = normrnd(0,1,[q,1]);
        y = exp_noise(t);
        u = -K*x + y;
        g = [];
        for i=1:q
            g = [g, B*G{i}*u];
        end
        if mod(k-1,sNt)==0
            X(:,j) = X(:,j)*((p-1)/p) + sm2vec(x*x')/p;
            Z(:,j) = Z(:,j)*((p-1)/p) + intZ/p;
            j = j + 1;
            intZ = zeros(n3,1);
        end
        xnext = x + (A*x +B*u)*dt + g*rv*dt^0.5;
        intZ = intZ + sm2vec([x;u]*[x',u'])*dt;
        x = xnext;
        t = t + dt;
    end
    X(:,j) = X(:,j)*((p-1)/p) + sm2vec(x*x')/p;
    Z(:,j) = Z(:,j)*((p-1)/p) + intZ/p;
end
Phi_est = Z(:,2:end)';
Psi_est = X(:,2:end)'-X(:,1:end-1)';

% compute true statistics of the discretized process
Xdis = zeros(n*(n+1)/2,M+1);
xm = x0;
xxm = xx0(:);
Adis = eye(n) + dt*A;
Bdis = dt*B;
coef_xxm = kron((Adis-Bdis*K),(Adis-Bdis*K));
for i=1:q
    Gdis{i} = G{i}*dt^0.5;
    coef_xxm = coef_xxm + kron(B*Gdis{i}*K,B*Gdis{i}*K);
end
j = 1;
intZm = zeros(n3,1);
Zdis = zeros(n3,M+1);
for k=1:Nt
    y = exp_noise((k-1)*dt);
    xmnext = (Adis-Bdis*K)*xm + Bdis*y;
    coef_xm = kron(Adis-Bdis*K,Bdis*y) + kron(Bdis*y,Adis-Bdis*K);
    coef = Bdis*(y*y')*Bdis';
    for i=1:q
        coef_xm = coef_xm - kron(B*Gdis{i}*K,B*Gdis{i}*y)...
            - kron(B*Gdis{i}*y,B*Gdis{i}*K);
        coef = coef + B*Gdis{i}*(y*y')*Gdis{i}'*B';
    end
    xxmnext = coef_xxm*xxm + coef_xm*xm + coef(:);
    XXM = reshape(xxm,[n,n]);
    if mod(k-1,sNt)==0
        Xdis(:,j) = sm2vec(XXM);
        Zdis(:,j) = intZm;
        j = j + 1;
        intZm = zeros(n3,1);
    end
    xum = -XXM*K' + xm*y';
    uum = K*XXM*K' - y*xm'*K' - K*xm*y' + y*y';
    intZm = intZm + sm2vec([XXM, xum; xum', uum;])*dt;
    xxm = xxmnext;
    xm = xmnext;
end
Xdis(:,j) = sm2vec(reshape(xxm,[n,n]));
Zdis(:,j) = intZm;
Phi_dis = Zdis(:,2:end)';
Psi_dis = Xdis(:,2:end)'-Xdis(:,1:end-1)';

err_coef = zeros(3,1);
coef_tr = Phi\Psi;
coef_est = Phi_est\Psi_est;
coef_dis = Phi_dis\Psi_dis;
err_coef(1) = norm(coef_est - coef_tr)/norm(coef_tr);
err_coef(2) = norm(coef_est - coef_dis)/norm(coef_dis);
err_coef(3) = norm(coef_dis - coef_tr)/norm(coef_tr);
end
