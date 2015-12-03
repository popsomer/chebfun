function [A, B, U, V, E] = computeCoeffs(scheme, L, LR, dt)
%COMPUTECOEFFS   Compute coefficients for KROGSTAD.

% Author: Hadrien Montanelli.

% Set-up:
s = scheme.internalStages;
q = scheme.steps; 
A = cell(s);
B = cell(s, 1);
C = zeros(s, 1);
U = cell(s, q-1);
V = cell(q-1, 1);
phit = cell(s);

% Compute C:
C(1) = 0;
C(2) = 1/2;
C(3) = 1/2;
C(4) = 1;

% Compute the phi functions:
phi1 = mean(feval(computePhi(1), LR), 2);
phi2 = mean(feval(computePhi(2), LR), 2);
phi3 = mean(feval(computePhi(3), LR), 2);
phit{1,2} = mean(C(2)^1*feval(computePhi(1), C(2)*LR), 2);
phit{1,3} = mean(C(3)^1*feval(computePhi(1), C(3)*LR), 2);
phit{1,4} = mean(C(4)^1*feval(computePhi(1), C(4)*LR), 2);
phit{2,2} = mean(C(2)^2*feval(computePhi(2), C(2)*LR), 2);

% Take real part fo diffusive problems (real eigenvalues): 
if ( isreal(L) == 1 )
    phi1 = real(phi1);
    phi2 = real(phi2);
    phi3 = real(phi3);
    phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
end

% Compute A:
A{3,2} = 4*phit{2,2};
A{4,3} = 2*phi2;

% Compute B:
B{2} = 2*phi2 - 4*phi3;
B{3} = 2*phi2 - 4*phi3;
B{4} = -phi2 + 4*phi3;

% Compute the missing oefficients using the summation properties of the 
% coefficients:
[A, B, U, V, E] = spinscheme.computeMissingCoeffs(A, B, C, U, V, L, phi1, ...
    phit, dt);

end