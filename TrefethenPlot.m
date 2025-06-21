R = 5900; clf, [ay,ax] = meshgrid([.56 .04],[.1 .5]);
for N = 300
[D,x] = cheb(N); D2 = D^2; D2 = D2(2:N,2:N);
S = diag([0; 1 ./(1-x(2:N).^2); 0]);
D4 = (diag(1-x.^2)*D^4 - 8*diag(x)*D^3 - 12*D^2)*S;
D4 = D4(2:N,2:N);
I = eye(N-1);
A = (D4-2*D2+I)/R - 2i*I - 1i*diag(1-x(2:N).^2)*(D2-I);
B = -1i*(D2-I);
[EV, evals] = eig(A,B);
evals = diag(evals);
[~, idx] = max(imag(evals));
ee = evals(idx);
figure(1)
set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 600]);
%plot(evals, '.', 'MarkerSize', 25, 'HandleVisibility','off');
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
h_legend = plot(nan, nan, '*r', MarkerSize=25);
legend(h_legend, {sprintf('$Im(\\omega_{cr}) = %.4f$', imag(ee))}, ...
       'Interpreter', 'latex', 'FontSize', 26, 'Location', 'southeast');

Re_vals = real(evals);
Im_vals = imag(evals);
A_idx = find(Re_vals < 0.66);
P_idx = find(Re_vals > 0.7);
S_idx = find(Re_vals >= 0.66 & Re_vals <= 0.7);
ix = P_idx(5);
plot(evals(A_idx), 'ob', 'MarkerSize', 10, 'DisplayName', 'A-branch');
plot(evals(S_idx), 'dg', 'MarkerSize', 10, 'DisplayName', 'S-branch');
plot(evals(P_idx), 'sm', 'MarkerSize', 10, 'DisplayName', 'P-branch');
legend('show', 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeast');
ev_unstable = EV(:, idx);
%ev_unstable = NB * T2 * ev_unstable_bal;
psi = 1i*alpha*ev_unstable(1:N-1);
psi = psi / max(abs(psi));
psi =real(psi);
%psi = psi/max(psi);
figure(2)
plot(x(2:N),psi, 'LineWidth', 1.5);
xlabel('y');
ylabel('Imag(\psi)');
end
