function dz = exp_dyn(t,z,K)
global A B G n n3 q;
y = exp_noise(t);
dz = zeros(n+n^2+n3,1);
dz(1:n) = (A-B*K)*z(1:n) + B*y;
coef_avg = kron(eye(n),B*y) + kron(B*y,eye(n));
coef = zeros(n);
for i=1:q
    coef_avg = coef_avg - kron(B*G{i}*K,B*G{i}*y)...
        - kron(B*G{i}*y,B*G{i}*K);
    coef = coef + B*G{i}*(y*y')*G{i}'*B';
end
dz(n+1:n^2+n) = stochastic_sys_mat(K)*z(n+1:n^2+n) + coef_avg*z(1:n)...
    + coef(:);
xx = reshape(z(n+1:n^2+n),[n,n]);
xu = -xx*K' + z(1:n)*y';
uu = K*xx*K' - K*z(1:n)*y' - y*z(1:n)'*K' +y*y';
dz(n+n^2+1:end) = sm2vec([xx, xu; xu', uu]);
end