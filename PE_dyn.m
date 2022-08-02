function dx = PE_dyn(t,x,model,K)
global Q R;
[m,n] = size(K);
theta_x = vec2sm(model*x,m+n);
dx = sm2vec(theta_x(1:n,1:n) - K'*theta_x(n+1:end,1:n)...
    - theta_x(n+1:end,1:n)'*K + Q + K'*R*K);
end