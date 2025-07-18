clear; clc;
N = 500;
Re = 800;
Wi = 80;
beta = 0.8;
alpha = 1.5;

[D, y] = cheb(N);
u = 1 - y.^2;
uprime = -2 .* y;
udoubleprime = -2;
Txx = 8 * ((1 - beta) / Re) * Wi * y.^2;
Txxprime = 16 * ((1 - beta) / Re) * Wi .* y;
Txy = 2 * ((beta - 1) / Re) .* y;
Txyprime = 2 * ((beta - 1) / Re);

D1 = D;
D2 = D1^2;
D3 = D1^3;
D4 = D1^4;
I = eye(N + 1);

firstRowmatrix = [1i*alpha*diag(u)*(D2 - alpha^2 * I) - 1i*alpha*diag(udoubleprime) * I - (beta/Re)*(D4 - 2*alpha^2*D2 + alpha^4*I), -1i*alpha*D1, -(alpha^2*I + D2), 1i*alpha*D1];
secondRowmatrix = [Wi*(-1i*alpha*diag(Txxprime) - 2*diag(Txy)*D2 - 2i*alpha*diag(Txx)*D1) - 2*((1 - beta)/Re)*1i*alpha*D1, (I + Wi * alpha * 1i * diag(u)), -2*Wi*diag(uprime), 0*I];
thirdRowmatrix = [Wi*(-1i*alpha*diag(Txyprime) - alpha^2*diag(Txx)) - ((1 - beta)/Re)*(D2 + alpha^2 * I), 0*I, I + 1i*alpha*Wi*diag(u), -Wi*diag(uprime)];
fourthRowmatrix = [-2*alpha^2*Wi*diag(Txy) + 2*((1 - beta)/Re)*1i*alpha*D1, 0*I, 0*I, I + 1i*Wi*alpha*diag(u)];

RHSfirstrow = [1i*(D2 - alpha^2 * I), 0*I, 0*I, 0*I];
RHSsecondrow = [0*I, 1i*Wi*I, 0*I, 0*I];
RHSthirdrow = [0*I, 0*I, 1i*Wi*I, 0*I];
RHSfourthrow = [0*I, 0*I, 0*I, 1i*Wi*I];

d1 = [D1(1,:), zeros(1,3*(N+1))];
d2 = [D1(end,:), zeros(1,3*(N+1))];
d3 = [1, zeros(1,N), zeros(1,3*(N+1))];
d4 = [zeros(1,N), 1, zeros(1,3*(N+1))];
d5 = [zeros(1,(N+1)), zeros(1,N), 1, zeros(1,2*(N+1))];
d6 = [zeros(1,2*(N+1)), zeros(1,N), 1, zeros(1,(N+1))];
d7 = [zeros(1,3*(N+1)), zeros(1,N), 1];
d8 = [zeros(1,(N+1)), 1, zeros(1,N), zeros(1,2*(N+1))];
d9 = [zeros(1,2*(N+1)), 1, zeros(1,N), zeros(1,(N+1))];
d10 = [zeros(1,3*(N+1)), 1, zeros(1,N)];

NB = null([d1 ; d2; d3; d4; d5; d6; d7; d8; d9; d10]);
A = [firstRowmatrix; secondRowmatrix; thirdRowmatrix; fourthRowmatrix];
B = [RHSfirstrow; RHSsecondrow; RHSthirdrow; RHSfourthrow];

AN = NB' * A * NB;
BN = NB' * B * NB;

[T1, T2] = balance2(AN, BN);
Abalanced = T1 * AN * T2;
Bbalanced = T1 * BN * T2;
[EV, evs] = eig(Abalanced, Bbalanced);
eeOB = diag(evs);
evals = eeOB;
[~, idx] = max(imag(evals));
ee = evals(idx);

figure(1)
set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 600]);
plot(evals, '.', 'MarkerSize', 40, 'HandleVisibility','off');
hold on;
plot(ee, '*r', 'MarkerSize', 40,  'HandleVisibility','off');
hold on;
yline(0,'k', LineWidth=1.5,HandleVisibility = 'off');
xlim([0.5 1.5]);
ylim([-1.2 0.2]);
xlabel('$\mathbf{Re(\omega)}$', 'Interpreter', 'latex', 'FontSize', 30, FontWeight='bold');
ylabel('$\mathbf{Im(\omega)}$', 'Interpreter', 'latex', 'FontSize', 30, FontWeight='bold');

grid off;
box on;
ax = gca;
ax.FontSize = 26;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.2;
max_imag = imag(ee);
msg = sprintf('$\\omega_{cr} = %.4f$', max_imag);

h_legend = plot(nan, nan, '*r', MarkerSize=40);  % red asterisk as placeholder
legend(h_legend, {sprintf('$Im(\\omega_{cr}) = %.4f$', imag(ee))}, ...
       'Interpreter', 'latex', 'FontSize', 30, 'Location', 'northeast');
ax_inset = axes('Position', [0.6, 0.25, 0.25, 0.25]);
distances = abs(evals - ee);
[~, sorted_indices] = sort(distances);
nearest_indices = sorted_indices(1:30);  % 30 closest eigenvalues
evals_zoom = evals(nearest_indices);

plot(evals_zoom, '.', 'MarkerSize', 25, 'HandleVisibility','off');
hold on;
plot(ee, '*r', 'MarkerSize', 20, 'HandleVisibility','off');
x_center = real(ee);
y_center = imag(ee);
xlim([x_center - 0.05, x_center + 0.05]);
ylim([y_center - 0.05, y_center + 0.05]);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', 'LineWidth', 1.0);


ev_unstable_bal = EV(:, idx);
ev_unstable = NB * T2 * ev_unstable_bal;
psi = ev_unstable(1:N+1);
psi = imag(psi);
psi = psi / max(psi);

figure(2)
set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 600]);
plot(y, psi, 'LineWidth', 2);
axis([0 1 -1 1]);
xlabel('$y$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\mathrm{Imag}(\psi)$', 'Interpreter', 'latex', 'FontSize', 18);
title(sprintf('Imaginary Part of Streamfunction\n$\\omega_{\\max} = %.5f$', ee), 'Interpreter', 'latex', 'FontSize', 20);
grid on;
box on;
ax = gca;
ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.2;
