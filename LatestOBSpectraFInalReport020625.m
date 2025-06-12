clear;clc;
N= 500;
Re =800;
E=1e-10;
Wi = E*Re;
Wi = Wi;
beta =0.8;

alpha =1.5;
[D, y] = cheb(N);
% D = D(2:N, 2:N);
% y = y(2:N);
u = 1 - y.^2;
uprime = -2 .* y;
udoubleprime = -2;
Txx = 8 * ((1 - beta) / Re) * Wi * y.^2;
Txxprime = 16* ((1 - beta) / Re) * Wi .* y;
Txy = 2 * ((beta - 1) / Re) .* y;
Txyprime = 2 * ((beta - 1) / Re);

% Higher-order differentiation matrices
D1= D;
D2 = D1^2;
D3 = D1^3;
D4 = D1^4;

I = eye(N+1);

%LHS Matrix(System Matrix)
firstRowmatrix = [1i*alpha*diag(u)*(D2-alpha^2*I)-1i*alpha*diag(udoubleprime)*I-(beta/Re)*(D4-2*alpha^2*D2 + alpha^4*I), -1i*(alpha)*D1, -(alpha^2*I + D2),1i*(alpha)*D1];
secondRowmatrix = [Wi*(-1i*alpha*diag(Txxprime)-2*diag(Txy)*D2 -2*1i*alpha*diag(Txx)*D1)-2*((1-beta)/Re)*1i*alpha*D1, (I+(Wi*alpha*1i*diag(u))) , -2*Wi*diag(uprime), 0*I];
thirdRowmatrix = [Wi*(-1i*alpha*diag(Txyprime)*I-alpha^2*diag(Txx))-((1-beta)/Re)*(D2+alpha^2*I), 0*I, I+(1i*alpha*Wi*diag(u)), -Wi*diag(uprime)];
fourthRowmatrix = [-2*(alpha^2)*Wi*diag(Txy)+2*((1-beta)/Re)*1i*alpha*D1, 0*I, 0*I,I+1i*Wi*alpha*diag(u)];

% RHS Matrix(Mass Mattrix)
RHSfirstrow = [1i*(D2-alpha^2*I), 0*I, 0*I, 0*I];
RHSsecondrow = [0*I, 1i*Wi*I, 0*I, 0*I];
RHSthirdrow = [ 0*I, 0*I, 1i*Wi*I, 0*I];
RHSfourthrow = [ 0*I, 0*I, 0*I, 1i*Wi*I];
% thirdRowmatrix(1,1:N-1) = D1(1, :);
% thirdRowmatrix(N-1,1:N-1) = D1(N-1, :);
% RHSthirdrow([1,N-1], :) = 0;

%Boundary Conditions
d1 = [D1(1,:) , zeros(1,3*(N+1))];          %psi'(-1)=0
d2 = [D1(end,:) , zeros(1,3*(N+1))];        %psi'(1)=0
d3 = [1, zeros(1,N) , zeros(1,3*(N+1))];    %psi(-1)=0
d4 = [zeros(1,N) , 1, zeros(1,3*(N+1))];    %psi(1)=0

d5 = [zeros(1,(N+1)), zeros(1,N) , 1, zeros(1,2*(N+1))];    %tau(-1)=0
d6 = [zeros(1,2*(N+1)), zeros(1,N) , 1, zeros(1,(N+1))];
d7 = [zeros(1,3*(N+1)), zeros(1,N) , 1];

d8 = [zeros(1,(N+1)), 1, zeros(1,N)  , zeros(1,2*(N+1))];   %tau(1)=0
d9 = [zeros(1,2*(N+1)),  1,zeros(1,N) , zeros(1,(N+1))];
d10 =[zeros(1,3*(N+1)), 1, zeros(1,N) ];

d11 = [zeros(1,(N+1)), D1(1,:), zeros(1,2*(N+1))];    %tau'(-1)=0
d12 = [zeros(1,2*(N+1)), D1(1,:), zeros(1,(N+1))];
d13 = [zeros(1,3*(N+1)), D1(1,:)];

d14 = [zeros(1,(N+1)), D1(end,:), zeros(1,2*(N+1))];    %tau'(-1)=0
d15 = [zeros(1,2*(N+1)), D1(end,:), zeros(1,(N+1))];
d16 = [zeros(1,3*(N+1)), D1(end,:)];

NB = null([d1 ; d2; d3; d4; d5; d6; d7; d8; d9; d10]); %; d11; d12; d13; d14; d15; d16]);
%Block Matrix Assembly
A = [firstRowmatrix;
    secondRowmatrix;
    thirdRowmatrix;
    fourthRowmatrix];
B = [RHSfirstrow;
    RHSsecondrow;
    RHSthirdrow;
    RHSfourthrow];


AN=NB'*A*NB; %Null Space Projection
BN=NB'*B*NB; %Null Space Projection

[T1, T2] = balance2(AN,BN);
Abalanced = T1*AN*T2;
Bbalanced = T1*BN*T2;
[EV,evs]= eig(Abalanced, Bbalanced);
eeOB = diag(evs);
%ix = real(eeOB)>=-0.5 & real(eeOB)<=2.5;
evals = eeOB;
[~, idx] = max(imag(evals));

ee = evals(idx);
figure(1)
plot(evals,'.', MarkerSize=45);
axis([0 1.5 -0.8 0.2]);

hold on;
plot(ee, '*r', MarkerSize=45);
max_imag = imag(ee);
msg = sprintf('Max Im(\\omega): %.4f', max_imag);
title(['Eigenvalue Spectrum for Re=', num2str(Re), ' ,Wi =',num2str(Wi), ' and beta=',num2str(beta)],FontSize=22, FontWeight='bold');
xlabel('Real(\omega)',FontSize=22,FontWeight='bold');
ylabel('Imag(\omega)',FontSize=22,FontWeight='bold');
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
xLimits = xlim;
yLimits = ylim;
text(xLimits(2) - 0.2*range(xLimits),yLimits(1) + 0.05*range(yLimits),msg,'FontSize', 24, 'FontWeight', 'bold', 'BackgroundColor', 'w');
hold off;
ev_unstable_bal = EV(:, idx);
ev_unstable = NB * T2 * ev_unstable_bal;
psi = ev_unstable(1:N+1);
%psi = psi / max(abs(psi));
psi =imag(psi);
psi = psi/max(psi);
figure(2)
plot(y,psi, 'LineWidth', 1.5);
axis([0 1 -1 1])
xlabel('Imag(\psi)');
ylabel('y');
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.XAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 12;
title(['Imaginary Part of Streamfunction, for \omega_{max} = ', num2str(ee)]);
