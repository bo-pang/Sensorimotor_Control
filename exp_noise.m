function y = exp_noise(t)
global ww mag;
y = mag*sum(sin(ww*t),2);
end