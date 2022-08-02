function X = stochastic_sys_mat(K)
global A B G n q;
X = kron((A-B*K),eye(n)) + kron(eye(n),(A-B*K));
for i=1:q
    X = X + kron(B*G{i}*K,B*G{i}*K);
end