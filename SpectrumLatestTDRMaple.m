clc; clear; close all;
n = 1;
iter = 1:1:n;
N =500;
eigenvalues_all = {};
figure(1)
%mp.Digits(30)
tic
for i=1:numel(iter)
    N = N;
    [D, y] = cheb(N);
    %y = mp(y);
    %D=mp(D);
    D2 = D^2;
    D3 = D^3;
    D4 = D^4;
    Re = 1400;
    Wi=12;
    I = eye(N+1);
    alpha = 1;
    B = zeros(1, N+1);
    B(1,1) = 1;

    y=y(2:end);

    y0=y;
    y=2*Re*Wi*y;

    Z = null(B);
    A = (108+8*y.^6+12*(12*y.^6+81).^(1/2)).^(1/3);
    integ = -Re*y.*(4*y.^4+2*A.*y.^2+A.^2)./(6*A);

    C_proj=Z'*D*Z;
    %C_proj = 2*Wi*Re*C_proj;

    u=linsolve(C_proj,integ);
    u_n = [0;u];
    %u_n = u/(Wi^2*Re);
    u_n = u_n/max(abs(u_n));
    %u_n = 1-y.^2;
    uprime = D*u_n;
    % y_plot = [1;y/(Wi*Re)];
    udoubleprime = D2*u_n;
    lambda = (1./(1+(Wi*uprime).^2)).^(1/3); %base state lagrange mult

    %plot(u_n, [-1;y]);

    %residual = uprime./(1+(Wi*uprime).^2).^(1/3) + Re*y_plot;
    %residual = udoubleprime./(1+(Wi*uprime).^2).^(1/3) - 2*Wi^2*uprime.^2.*udoubleprime./(1+(Wi*uprime).^2).^(4/3);
    %plot(y_plot,residual);
    %lambda = ones(1,3);
    %lambda = diag(lambdavec); %lambda matrix
    Bxx = lambda.*(1 + 2*(uprime*Wi).^2);
    Byy = lambda;
    Bxy = lambda.*((Wi*uprime));
    Bzz = lambda;

    Bxxprime = D*diag(Bxx);
    Byyprime = D*diag(Byy);
    Bxyprime = D*diag(Bxy);
    Bzzprime = D*diag(Bzz);
    q = diag(lambda.^2/(3*Wi));

    sfLHS =[diag(u_n)*1i*alpha*(D2-I*alpha^2)-1i*alpha*I*diag(udoubleprime), -(1/(Wi*Re))*1i*alpha*D, -(1/(Wi*Re))*(D2+alpha^2*I), (1/(Wi*Re))*1i*alpha*D,0*I];
    xxstLHS= [-2*diag(Bxy)*D2-1i*alpha*Bxxprime-2i*alpha*diag(Bxx)*D,1i*alpha*diag(u_n)+(1/Wi)*I+q*(diag(Bzz)+diag(Byy)),-2*diag(uprime)-2*q*diag(Bxy),q*diag(Bzz+Bxx), q*diag(Byy+Bxx)];
    xystLHS = [-1i*alpha*Bxyprime-alpha^2*diag(Bxx)-diag(Byy)*D2, 0*I,1i*alpha*diag(u_n)+(1/Wi)*I,-diag(uprime),0*I];
    yystLHS = [-1i*alpha*Byyprime + 2i*alpha*diag(Byy)*D-2*alpha^2*diag(Bxy), q*diag(Byy+Bzz), -2*q*diag(Bxy), 1i*alpha*diag(u_n) + (1/Wi)*I + q*diag(Bzz+Bxx), q*diag(Bxx+Byy)];
    zzstLHS = [-1i*alpha*Bzzprime, q*diag(Byy+Bzz), -2*q*diag(Bxy), q*diag(Bxx+Bzz),1i*alpha*diag(u_n)+(1/Wi)*I+q*diag(Bxx+Byy)];


    sfRHS = [1i*(D2-I*alpha^2), 0*I, 0*I, 0*I, 0*I];
    xxStRHS = [0*I, 1i*I, 0*I, 0*I, 0*I];
    xyStRHS = [0*I, 0*I, 1i*I, 0*I, 0*I];
    yyStRHS = [0*I, 0*I, 0*I, 1i*I, 0*I];
    zzStRHS = [0*I, 0*I, 0*I, 0*I, 1i*I];

    d1 = [D(1,:) , zeros(1,4*(N+1))];
    d2 = [D(end,:) , zeros(1,4*(N+1))];
    d3 = [1, zeros(1,N) , zeros(1,4*(N+1))];
    d4 = [zeros(1,N) , 1, zeros(1,4*(N+1))];
    d5 = [zeros(1,(N+1)), zeros(1,N) , 1, zeros(1,3*(N+1))];
    d6 = [zeros(1,2*(N+1)), zeros(1,N) , 1, zeros(1,2*(N+1))];
    d7 = [zeros(1,3*(N+1)), zeros(1,N) , 1, zeros(1, N+1)];
    d8 = [zeros(1,4*(N+1)), zeros(1,N) , 1];
    d9 = [zeros(1,(N+1)), 1, zeros(1,N)  , zeros(1,3*(N+1))];
    d10 = [zeros(1,2*(N+1)),  1,zeros(1,N) , zeros(1,2*(N+1))];
    d11 =[zeros(1,3*(N+1)), 1, zeros(1,N), zeros(1, N+1)];
    d12 =[zeros(1,4*(N+1)),1, zeros(1,N)];

    A = [sfLHS;
        xxstLHS;
        xystLHS;
        yystLHS;
        zzstLHS];

    B = [sfRHS;
        xxStRHS;
        xyStRHS;
        yyStRHS;
        zzStRHS];

    BCM = null([d1; d2; d3; d4; d5; d6; d7; d8; d9; d10; d11; d12]);
    AN = BCM'*A*BCM;
    BN = BCM'*B*BCM;
    [EVS, ~] = eig(AN,BN);

    [T1, T2] = balance2(AN, BN);
    Abalanced = T1*AN*T2;
    Bbalanced = T1*BN*T2;


    [EV, eigss] = eig(Abalanced, Bbalanced);
    eigss = diag(eigss);
    ix = real(eigss)>=-0.1 & real(eigss)<=1.5;
    %evals = eigss(ix);
    EV = EV(:,ix);
    EVS = EVS(:, ix);
    evals = eigss;
    [~, idx] = max(imag(evals));
    ee = evals(idx)
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

    ev_unstable_bal = EVS(:, idx);
    ev_unstable = BCM*ev_unstable_bal;
    psi = ev_unstable(1:N+1);
    %psi = psi / max(abs(psi));
    psi_imag =imag(psi);
    psi_imag = psi_imag/max(abs(psi_imag));
    psi_real =real(psi);
    psi_real = psi_real/max(abs(psi_real));
    psi_mag =abs(psi);
    psi_mag = psi_mag/max(abs(psi_mag));
    y_tot = y/(2*Wi*Re);
    y_tot = [1;y_tot];
    figure(2)
    plot(y_tot,psi_imag,'r', 'LineWidth', 1.5);
    hold on;
    plot(y_tot,psi_real,'c', 'LineWidth', 1.5);
    hold on;
    plot(y_tot,psi_mag,'b', 'LineWidth', 1.5);
    axis([-1 1 -1 1])
    legend({'Imag($\psi$)', 'Real($\psi$)', 'abs($\psi$)'}, ...
       'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
    ylabel('$\psi$','Interpreter','latex', 'FontSize',24,'FontWeight','bold');
    xlabel('y', 'FontSize',24,'FontWeight','bold');
    ax = gca;
    % ax.XAxis.FontWeight = 'bold';
    ax.XAxis.FontSize = 22;
    %ax.YAxis.FontWeight = 'bold';
    ax.YAxis.FontSize = 22;
end
toc
eigvec = EV(:, idx);
eigvec =T2*eigvec;
residue_vec = AN*eigvec - ee*BN*eigvec;
%residue_vec = BCM*residue_vec;
%plot(y_tot(2:end-1),abs((residue_vec(2:N))));
plot(y_tot(2:end-1), residue_vec(1:N-1));
norm_sf = norm(residue_vec)
% norm_res = norm(abs(residue_vec))
normAv = norm(AN*eigvec)
rel_norm = norm_sf/normAv;
disp(rel_norm)