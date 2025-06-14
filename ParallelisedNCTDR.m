clear; clc; clf;
N = 500;
alpha = 1;
Re = 5000;
Wi_values = 2:2:8;  
Re_critical = zeros(size(Wi_values));
sensits = zeros(size(Wi_values));


tic
parpool(4);  % parallel pool with 4 workers


parfor ii = 1:numel(Wi_values)
    eemax = 1i;
    Wi = Wi_values(ii);
    alpha_local = alpha; 
    Re_local = Re; 

    while abs(imag(eemax)) > 1e-9 && Re_local<50000
        fun = @(alpha) max_imag(Re_local, Wi, alpha, N);
        options = optimoptions('fminunc', 'StepTolerance', 1e-2, 'OptimalityTolerance', 1e-2, ...
                               'SpecifyObjectiveGradient', true, 'Display', 'off');
        alpha_local = fminunc(fun, alpha_local, options);
        [eemax, dwda, dwdRe] = OB(Re_local, Wi,  alpha_local, N);
        dRe = imag(eemax) / imag(dwdRe);
        relax = 1;
        Re_local = min(Re_local + relax * dRe, 50000);

    end
    
    Re_critical(ii) = Re_local;  
end

toc

figure;
plot(Re_critical,Wi_values ,'d-','LineWidth', 2);
xlabel('Critical Reynolds Number (Re_{cr})');
ylabel('Weissenberg Number (Wi)');
title('Neutral Stability Curve (TDR Model)');
grid on;
output_data = [Wi_values(:), Re_critical(:)];
writematrix(output_data, 'neutral_curve.csv');  % CSV export
save('neutral_curve.mat', 'Wi_values', 'Re_critical');  % Optional MAT file

function [max_im, dwdalpha_imag, dwdRe] = max_imag(Re, Wi,alpha, N)
    [eigg, dwdalphacalc, dwdRecalc] = OB(Re, Wi, alpha, N);
    max_im = -imag(eigg);
    dwdRe = imag(dwdRecalc); 
    dwdalpha_imag = imag(dwdalphacalc);
end

function [ee,dwdalpha, dwdRe] = OB(Re, Wi, alpha, N)    
    [D, y] = cheb(N);
    D2 = D^2;
    D3 = D^3;
    D4 = D^4;
    B = zeros(1, N+1);
    B(1,1) = 1;
    
    y=y(2:end);
    
    y0=y;
    y=Re*Wi*y;
    
    Z = null(B);
    A = (108+8*y.^6+12*(12*y.^6+81).^(1/2)).^(1/3);
    integ = -y.*(4*y.^4+2*A.*y.^2+A.^2)./(6*A);
    
    C_proj=Z'*D*Z;
    
    u=linsolve(C_proj,integ);
    u = [0;u];
    u_n = u/(Wi^2*Re);
    u_n = u_n/max(abs(u_n));
    
    %u_n = 1-y.^2;
    uprime = D*u_n;
    udoubleprime = D2*u_n;
    lambda = (1./(1+(Wi*uprime).^2)).^(1/3); %base state lagrange mult
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
    
    % Higher-order differentiation matrices
    D1= D;
    D2 = D1^2;
    D3 = D1^3;
    D4 = D1^4;
    
    I = eye(N+1);
    
    %LHS Matrix(System Matrix)
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
    
    [T1, T2] = balance2(AN, BN);
    Abalanced = T1*AN*T2;
    Bbalanced = T1*BN*T2;
    [EV,evs]= eig(Abalanced, Bbalanced);
    eeOB = diag(evs);
    ix = real(eeOB)>=-1 & real(eeOB)<=1;
    evals = eeOB(ix);
    EV_filtered = EV(:, ix);
    [~, w] = clencurt(N);

    [~, idx] = max(imag(evals));
    ee = evals(idx);
    AA = EV_filtered(:, idx);
    eeT = conj(ee);
    OB = Abalanced - ee*Bbalanced;
    %Projection of weights onto the null space
    W1 = BCM'*kron(eye(5),diag(1./w))*BCM;
    W2 =BCM'*kron(eye(5),diag(w))*BCM;

    %Adjoint of the OrrSommerfeld system for Sensitivity analysis 
    OBAdj =W1*ctranspose(OB)*W2;


    %Adjoint Eigenvector corresponding to max eigenvalue
    OBT = null(OBAdj);
    [~, cols] = size(OBT);
    
    for ii= 1:cols
        dA1dRe = [0*I,(1/(Wi*Re^2))*1i*alpha*D, (1/(Wi*Re^2))*(D2+alpha^2*I), -(1/(Wi*Re))*1i*alpha*D, 0*I]; %d/dRe of momentum eqn
        dA2dRe = [0*I, 0*I,0*I, 0*I, 0*I]; 
        dA3dRe = [0*I, 0*I,0*I, 0*I, 0*I]; 
        dA4dRe = [0*I, 0*I,0*I, 0*I, 0*I];
        dA5dRe = [0*I, 0*I,0*I, 0*I, 0*I];
        dAdRe= [dA1dRe;dA2dRe;dA3dRe;dA4dRe;dA5dRe];        
        dOBdRe = dAdRe;
        dOBdRe = BCM'*dOBdRe*BCM;
        dOBdRe = T1*dOBdRe*T2;
        dwdRe = -dot(W2*OBT(:, ii), dOBdRe*EV_filtered(:, idx))/dot(W2*OBT(:, ii), Bbalanced*EV_filtered(:,idx)); %sensitivity calculation dwdRe
        dA1dalpha = [diag(u_n)*1i*(D2-I*alpha^2)-2*I*alpha^2*1i*diag(u_n)-1i*diag(udoubleprime), -(1/(Wi*Re))*1i*D, -(1/(Wi*Re))*(2*alpha*I), (1/(Wi*Re))*1i*D,0*I];
        dA2dalpha = [-1i*Bxxprime-2i*diag(Bxx)*D,1i*diag(u_n),0*I,0*I, 0*I];
        dA3dalpha = [-1i*Bxyprime-2*alpha*diag(Bxx), 0*I,1i*diag(u_n),0*I,0*I];
        dA4dalpha = [-1i*Byyprime + 2i*diag(Byy)*D-4*alpha*diag(Bxy), 0*I, 0*I, 1i*diag(u_n),0*I];
        dA5dalpha = [-1i*Bzzprime, 0*I, 0*I, 0*I,1i*diag(u_n)];
        dB1dalpha  = [-2i*alpha*I,0*I 0*I, 0*I, 0*I];
        dB2dalpha  = [ 0*I,0*I, 0*I, 0*I, 0*I];
        dB3dalpha  = [ 0*I,0*I, 0*I, 0*I, 0*I];
        dB4dalpha  = [ 0*I,0*I, 0*I, 0*I, 0*I];
        dB5dalpha  = [ 0*I,0*I, 0*I, 0*I, 0*I];
        dAdalpha = [dA1dalpha;dA2dalpha;dA3dalpha;dA4dalpha;dA5dalpha];
        dBdalpha = [dB1dalpha;dB2dalpha;dB3dalpha;dB4dalpha;dB5dalpha];
        dOBdalpha = dAdalpha - ee*dBdalpha;
        dOBdalphaBC = BCM'*dOBdalpha*BCM;
        dOBdalphafinal = T1*dOBdalphaBC*T2;
        dwdalpha = -dot(W2*OBT(:, ii), dOBdalphafinal*EV_filtered(:, idx))/dot(W2*OBT(:, ii), Bbalanced*EV_filtered(:,idx));
        end
end