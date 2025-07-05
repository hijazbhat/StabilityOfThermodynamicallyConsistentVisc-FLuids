
function [T1 T2,ilo,ihi] = balance2(A,B,rho)
mp.digits(30)
if (nargin < 2)
    error(message('Gcg:NotEnoughInputs'));
end
if ~isequal(size(A),size(B)),
    error('Balance2:nonsymmetric',...
        'Matrices A and B must both be square and of equal size');
end
% Default values
if nargin < 3 || isempty(rho),
    rho = mp(2);
end
% -------------------------------------------------------------------------
% Step 1: Reduce the order
n = size(A,1);                          % matrix size
% Incidence and permutaiton matrices
W = (A ~= 0 | B ~= 0);                  % incidence matrix, see [War1981, p.144]
p1 = 1:n;                               % incidence vector for row permuation matrix
p2 = 1:n;                               % incidence vector for column permuation matrix
% Row scanning, see [War1981, p.144]
ihi = n;
for i = n:-1:1,
    r = find(sum(W(1:i,1:i),2) == 1,1); % find row with one nonzero element
    if isempty(r),
        break;
    else
        c = find(W(r,1:i));
        p1([r i]) = p1([i r]);          % swap rows
        p2([c i]) = p2([i c]);          % swap columns
        W([r i],:) = W([i r],:);        % swap rows
        W(:,[c i]) = W(:,[i c]);        % swap columns
        ihi = i - 1;                    % update high index
    end
end
% Column scanning, see [War1981, p.144]
ilo = 1;
for i = 1:ihi,
    c = find(sum(W(i:ihi,i:ihi),1) == 1,1) + i - 1; % find column with one nonzero element
    if isempty(c),
        break;
    else
        r = find(W(i:ihi,c)) + i - 1;
        p1([r i]) = p1([i r]);          % swap rows
        p2([c i]) = p2([i c]);          % swap columns
        W([r i],:) = W([i r],:);        % swap rows
        W(:,[c i]) = W(:,[i c]);        % swap columns
        ilo = i + 1;                    % update low index
    end
end
P1 = sparse(1:n,p1,mp(ones(1,n)),n,n);      % row permutation matrix
P2 = sparse(p2,1:n,mp(ones(1,n)),n,n);      % column permuation matrix
% Return permutation matrices if no scaling is necessary
if ihi == ilo,
    T1 = P1; T2 = P2;
    ihi = 1; ilo = n;
    return;
end
% -------------------------------------------------------------------------
% Step 2: Scale A and B
% Non-triangularizable submatrices
m = ihi-ilo+1;                          % submatrix size
A22 = P1(ilo:ihi,:)*A*P2(:,ilo:ihi);    % submatrix of A
B22 = P1(ilo:ihi,:)*B*P2(:,ilo:ihi);    % submatrix of B
% Incidence matrices of nonzero elements
IA = (A22 ~= 0);                        % incidence matrix for A
IB = (B22 ~= 0);                        % incidence matrix for B
% Diagonal matrices whoes elements are the number of nonzeros, see
% [War1981, p.144]
F1 = diag(sum(IA+IB,2));                % number of nonzeros in the rows
F2 = diag(sum(IA+IB,1));                % number of nonzeros in the columns
% Sum of the incidence matrices, see [War1981, p.144]
E = IA + IB;
% Vectors of row and column sums of logarithms of nonzero elements, see
% [War1981, p.144]
logA = log(abs(A22))/log(rho); logA(~IA) = mp(0);
logB = log(abs(B22))/log(rho); logB(~IB) = mp(0);
g = sum(logA+logB,2);                   % row sums
h = sum(logA+logB,1)';                  % column sums
% Coefficient matrix to solve for powers of rho for scaling factors, see
% [War1981, p.144]
MN = [ F1  mp(E)
       mp(E')  F2 ];
% Right hand side vector, see [War1981, p.144]
b = [ -g
      -h ];
% Inverse of auxiliary system [War1981, p.145]
I = mp(eye(m));
e = mp(ones(m,1));
Minv = mp(1)/mp(2)*[mp(1/m)*I - mp(3)/(mp(4)*m^2)*(e*e')   mp(1)/(mp(4)*m^2)*(e*e')
            mp(1)/(mp(4)*m^2)*(e*e')           mp(1/m)*I - 3/(4*m^2)*(e*e')];
% Find powers of rho by generalized conjugate gradient method using Minv
% as auxiliary matrix
rc = round(gcg(MN,b,.5,m,Minv));
% Powers of rho for row and columns scaling
r = rc(1:m);
c = rc(m+1:end);
% Global scaling matrix
D1 = mp(speye(n)); D1(ilo:ihi,ilo:ihi) = diag(rho.^r);
D2 = mp(speye(n)); D2(ilo:ihi,ilo:ihi) = diag(rho.^c);
% -------------------------------------------------------------------------
% Step 3: Grade A/B
A22 = diag(rho.^r)*A22*diag(rho.^c);    % scaled submatrix of A
B22 = diag(rho.^r)*B22*diag(rho.^c);    % scaled submatrix of B
% ratio of row and columns norms, see [War1981, p.146]
rnorm = sum(abs(A22),2)./sum(abs(B22),2);   % ratio of row norms
cnorm = sum(abs(A22),1)./sum(abs(B22),1);   % ratio of column norms
[rnorm ir] = sort(rnorm,'descend');
[cnorm ic] = sort(cnorm,'descend');
% Permutation matrices
G1 = mp(speye(n)); G1(ilo:ihi,:) = G1(ilo+ir-1,:);
G2 = mp(speye(n)); G2(:,ilo:ihi) = G2(:,ilo+ic-1);
% -------------------------------------------------------------------------
% Output
T1 = G1*D1*P1;                          % left multiplication matrix
T2 = P2*D2*G2;                          % right multiplication matrix