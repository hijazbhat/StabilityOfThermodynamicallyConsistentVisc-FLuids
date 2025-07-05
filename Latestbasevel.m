Wi_values = 12:1:14;
figure()
hold on
Re_values = 1000:500:3500;
for i=1:length(Wi_values)
    Wi = Wi_values(i)
    N = 100;
    [D_up,y_up] = cheb(N);
    B = zeros(1, N+1);
    B(1,1) = 1;
    Wi =Wi_values(i);
    y=y_up(2:end);
    Re=3000;
    
    y0=y;
    y=2*Re*Wi*y;
    
    Z = null(B);
    P = (12.*(y.^6)+81).^(1/2);
    % A = 108 + 8*(y.^6)+12.*P;
    % integ = -(y.*(4*(y.^4)+2*(y.^2).*(A.^(1/3))+A.^(2/3)))./(6*(A.^(1/3)));
    D_proj = Z'*D_up*Z;
    D_proj = D_proj*2*Wi*Re;
    A = (108+8.*y.^6+12.*(12.*y.^6+81).^(1/2)).^(1/3);
    integ =-2*Re*y.*(4.*y.^4+2.*A.*y.^2+A.^2)./(6*A);
    u = linsolve(D_proj, integ);
    u = [0;u];
    %
    % 
    %u =u/(4*Wi^2*Re);
    %u = u/(max(abs(u)));
    uprime = D_up*u;
    udoubleprime = D_up*uprime;
    residue = D_up*(uprime./(1+(Wi*uprime).^2).^(1/3))+2*Re;
    %residue = (udoubleprime./(1+Wi^2*uprime.^2).^(1/3) -(2*Wi^2*(uprime.^2).*udoubleprime)./(3*(1+Wi^2*uprime.^2).^(4/3)));% +2*Re;
    plot(y_up,residue/(max(u)));
    
 
    %plot(u,y_up);
    disp(norm(residue/max(u)))
   
    % figure(2)
    % plot(u,y_up);
    % figure(3)
    % plot(uprime, y_up)
    % figure(4)
    % plot(udoubleprime, y_up)
end
