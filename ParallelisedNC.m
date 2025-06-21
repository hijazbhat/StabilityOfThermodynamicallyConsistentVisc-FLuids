clear; clc; clf;
N = 100;
alpha = 1;
Re = 5000;
beta = 0.9;
Wi_values = 1:1:15;  
Re_critical = zeros(size(Wi_values));
sensits = zeros(size(Wi_values));


tic
parpool(4);  % parallel pool with 4 workers


parfor ii = 1:numel(Wi_values)
    eemax = 1i;
    Wi = Wi_values(ii);
    alpha_local = alpha; 
    Re_local = Re; 

    while abs(imag(eemax)) > 1e-9
        figure(1), hold on;
        fun = @(alpha) max_imag(Re_local, Wi, beta, alpha, N);
        options = optimoptions('fminunc', 'StepTolerance', 1e-2, 'OptimalityTolerance', 1e-2, ...
                               'SpecifyObjectiveGradient', true, 'Display', 'off');
        alpha_local = fminunc(fun, alpha_local, options);
        [eemax, dwda, dwdRe] = OB(Re_local, Wi, beta, alpha_local, N);
        dRe = imag(eemax) / imag(dwdRe);
        relax = 1;
        Re_local = Re_local + relax * dRe;

        quiver(Re_local, imag(eemax), 0, 0.001 * sign(imag(dwdRe)), 'MaxHeadSize', 1);
    end
    
    Re_critical(ii) = Re_local;  
end

toc

figure;
plot(Re_critical,Wi_values ,'d-', 'MarkerSize',14, 'LineWidth', 4);
xlabel('Critical Reynolds Number ($Re_{cr}$)','FontSize',22,'Interpreter', 'latex');
ylabel('Weissenberg Number (Wi)','FontSize',22);
grid off;

function [max_im, dwdalpha_imag, dwdRe] = max_imag(Re, Wi, beta, alpha, N)
    [eigg, dwdalphacalc, dwdRecalc] = OB(Re, Wi, beta, alpha, N);
    max_im = -imag(eigg);
    dwdRe = imag(dwdRecalc); 
    dwdalpha_imag = imag(dwdalphacalc);
end

function [ee,dwdalpha, dwdRe] = OB(Re, Wi, beta, alpha, N)    
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
    evals = eeOB;
    [~, w] = clencurt(N);

    [~, idx] = max(imag(eeOB));
    ee = eeOB(idx);
    AA = EV(:, idx);
    eeT = conj(ee);
    OB = Abalanced - ee*Bbalanced;
    %Projection of weights onto the null space
    W1 = NB'*kron(eye(4),diag(1./w))*NB;
    W2 =NB'*kron(eye(4),diag(w))*NB;

    %Adjoint of the OrrSommerfeld system for Sensitivity analysis 
    OBAdj =W1*ctranspose(OB)*W2;


    %Adjoint Eigenvector corresponding to max eigenvalue
    OBT = null(OBAdj);
    [~, cols] = size(OBT);
    
    for ii= 1:cols
        dfirstrowdRe = [(beta/Re^2)*(D4-2*alpha^2*D2 + alpha^4*I), 0*I, 0*I, 0*I]; %d/dRe of momentum eqn
        dsecondrowdRe = [2*1i*((1-beta)/Re^2)*alpha*D1,0*I, 0*I, 0*I]; %d/dRe for Txx
        dthirdrowdRe = [((1-beta)/Re^2)*(D2+alpha^2*I), 0*I, 0*I, 0*I]; %d/dRe for Txy
        dfourthrowdRe = [-2*((1-beta)/Re^2)*1i*alpha*D1, 0*I, 0*I, 0*I];%d/dRe for Tyy
        dAdRe= [dfirstrowdRe;dsecondrowdRe;dthirdrowdRe;dfourthrowdRe];        
        dOBdRe = dAdRe;
        dOBdRe = NB'*dOBdRe*NB;
        dOBdRe = T1*dOBdRe*T2;
        dwdRe = -dot(W2*OBT(:, ii), dOBdRe*EV(:, idx))/dot(W2*OBT(:, ii), Bbalanced*EV(:,idx)); %sensitivity calculation dwdRe
        dA1dalpha = [-1i*diag(u)*D2-3i*alpha^2*diag(u)-1i*diag(udoubleprime)*I-(beta/Re)*(-4*alpha*D2+4*alpha^3*I),-1i*D1, -2*alpha*I, 1i*D1];
        dA2dalpha = [Wi*(-1i*diag(Txxprime) -2*1i*diag(Txx)*D1)-2*((1-beta)/Re)*1i*D1, (Wi*1i*diag(u)) , 0*I, 0*I];
        dA3dalpha = [Wi*(-1i*diag(Txyprime)*I-2*alpha*diag(Txx))-((1-beta)/Re)*(2*alpha*I), 0*I, (1i*Wi*diag(u)), 0*I];
        dA4dalpha = [-4*(alpha)*Wi*diag(Txy)+2*((1-beta)/Re)*D1, 0*I, 0*I,1i*Wi*diag(u)];
        dB1dalpha = [1i*(-2*alpha*I), 0*I, 0*I, 0*I];
        dB2dalpha  = [0*I, 0*I, 0*I, 0*I];
        dB3dalpha  = [ 0*I, 0*I, 0*I, 0*I];
        dB4dalpha  = [ 0*I, 0*I, 0*I, 0*I];
        dAdalpha = [dA1dalpha;dA2dalpha;dA3dalpha;dA4dalpha];
        dBdalpha = [dB1dalpha;dB2dalpha;dB3dalpha;dB4dalpha];
        dOBdalpha = dAdalpha - ee*dBdalpha;
        dOBdalphaBC = NB'*dOBdalpha*NB;
        dOBdalphafinal = T1*dOBdalphaBC*T2;
        dwdalpha = -dot(W2*OBT(:, ii), dOBdalphafinal*EV(:, idx))/dot(W2*OBT(:, ii), Bbalanced*EV(:,idx));
        end

end