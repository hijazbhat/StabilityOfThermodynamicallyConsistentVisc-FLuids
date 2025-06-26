N = 200;
[D_up,y_up] = cheb(N);
B = zeros(1, N+1);
B(1,1) = 1;

y=y_up(2:end);
R=1000;
Wi=2;
y0=y;
y=2*R*Wi*y;

Z = null(B);
A = (108+8*y.^6+12*(12*y.^6+81).^(1/2)).^(1/3);
integ = -y.*(4*y.^4+2*A.*y.^2+A.^2)./(6*A);

C_proj=Z'*D_up*Z;

u=linsolve(C_proj,integ);
u = [0;u];
u = u/(4*R*Wi^2);
uprime = D_up*u;
udoubleprime = D_up*uprime;
residual = D_up*(uprime./(1+(Wi*uprime).^2).^(1/3))+2*Re;
% residual = udoubleprime./(1+(Wi*uprime).^2).^1/3 - 2*Wi^2*uprime.^2.*udoubleprime./(1+(Wi*uprime).^2).^4/3;
plot(y_up, residual);
hold on;
%plot(u/max(u),y_up)

