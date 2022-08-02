function P = PI(K)
global Q R n;
    p = stochastic_sys_mat(K)'...
        \vec(-Q-K'*R*K);
    Pr = reshape(p,[n,n]);
    P = (Pr'+Pr)/2;
end