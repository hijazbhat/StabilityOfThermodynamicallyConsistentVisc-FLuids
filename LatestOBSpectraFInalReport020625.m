clear;clc;
N= 500;
Re =800;
%E=0.05;
%Wi = E*Re;
Wi = 0.8;
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
ix = real(eeOB)>=-2 & real(eeOB)<=2;
evals = eeOB(ix);
%evals = eeOB;
[~, idx] = max(imag(evals));

ee = evals(idx);
figure(1)
set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 600]);
plot(evals, '.', 'MarkerSize', 25, 'HandleVisibility','off');
hold on;
plot(ee, '*r', 'MarkerSize', 25,  'HandleVisibility','off');
hold on;
yline(0,'k', LineWidth=1.5,HandleVisibility = 'off');
xlim([0.5 1.5]);
ylim([-1.2 0.2]);
xlabel('$\mathbf{Re(\omega)}$', 'Interpreter', 'latex', 'FontSize', 24, FontWeight='bold');
ylabel('$\mathbf{Im(\omega)}$', 'Interpreter', 'latex', 'FontSize', 24, FontWeight='bold');

grid off;
box on;
ax = gca;
ax.FontSize = 22;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.2;
max_imag = imag(ee);
msg = sprintf('$\\omega_{cr} = %.4f$', max_imag);

h_legend = plot(nan, nan, '*r', MarkerSize=25);  % red asterisk as placeholder
legend(h_legend, {sprintf('$Im(\\omega_{cr}) = %.4f$', imag(ee))}, ...
       'Interpreter', 'latex', 'FontSize', 26, 'Location', 'southeast');
% ax_inset = axes('Position', [0.20, 0.15, 0.15, 0.15]);
% distances = abs(evals - ee);
% [~, sorted_indices] = sort(distances);
% nearest_indices = sorted_indices(1:100);  % 30 closest eigenvalues
% x_min = real(ee)-0.7;
% x_max = real(ee)+0.10;
% y_min = imag(ee)-0.15;
% y_max =imag(ee)+0.15;
% 
% real_part = real(evals);
% imag_part = imag(evals);
% 
% in_region = (real_part >= x_min) & (real_part <= x_max) & (imag_part >= y_min) & (imag_part <= y_max);
% 
% evals_region = evals(in_region);
% %evals_zoom = evals(evals_region);
% 
% plot(evals_region, '.', 'MarkerSize', 15, 'HandleVisibility','off');
% hold on;
% plot(ee, '*r', 'MarkerSize', 10, 'HandleVisibility','off');
% yline(0, 'k', 'LineWidth', 1.2);
% x_center = real(ee);
% y_center = imag(ee);
% xlim([x_center - 0.05, x_center + 0.05]);
% ylim([y_center - 0.05, y_center + 0.05]);
% set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', 'LineWidth', 1.0);

ev_unstable_bal = EV(:, idx);
ev_unstable = NB * T2 * ev_unstable_bal;
psi = ev_unstable(1:N+1);
%psi = psi / max(abs(psi));
psi =real(1i*alpha*psi);
psi = psi/max(abs(psi));
figure(2)
plot(y,psi, 'LineWidth', 1.5);
axis([-1 1 -1 1])
ylabel('Im($\psi$)','Interpreter','latex', 'FontSize',24,'FontWeight','bold');
xlabel('y', 'FontSize',24,'FontWeight','bold');
ax = gca;
% ax.XAxis.FontWeight = 'bold';
ax.XAxis.FontSize = 22;
%ax.YAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 22;
%title(['Imaginary Part of Streamfunction, for $\omega_{max}$ = ', num2str(ee)],'Interpreter','latex', 'FontSize', 22);
