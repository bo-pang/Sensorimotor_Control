function [x,y] = stiffness(K,color)
angle = [0:100]/100*2*pi;
x0 = sin(angle);
y0 = cos(angle);
ellipse = K*[x0;y0];
plot(ellipse(1,:),ellipse(2,:),color);hold on;
% axis([-1100,1100,-1000,1000])
end