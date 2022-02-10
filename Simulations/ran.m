n=1000
radius=1e-6;
r = (radius-1.0e-9) * (((rand(1,n))').^3);
% r = rand(n,1).^(1/3);
costheta = 2*rand(n,1)-1;
sintheta = sqrt(1-costheta.^2);
phi=rand(n,1)*(2*pi);
x=r.*sintheta.*cos(phi);
y=r.*sintheta.*sin(phi);
z=r.*costheta;
figure;
plot3(x,y,z,'.')