function [A, B, U, V, E] = computeCoeffs(scheme, L, LR, dt)
%COMPUTECOEFFS   Compute coefficients for EGLM433.

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
C(3) = 1;

% Compute the phi functions:
phi1 = mean(feval(computePhi(1), LR), 2);
phi2 = mean(feval(computePhi(2), LR), 2);
phi3 = mean(feval(computePhi(3), LR), 2);
phi4 = mean(feval(computePhi(4), LR), 2);
phi5 = mean(feval(computePhi(5), LR), 2);
phit{1,2} = mean(C(2)^1*feval(computePhi(1), C(2)*LR), 2);
phit{1,3} = mean(C(3)^1*feval(computePhi(1), C(3)*LR), 2);
phit{2,2} = mean(C(2)^2*feval(computePhi(2), C(2)*LR), 2);
phit{3,2} = mean(C(2)^3*feval(computePhi(3), C(2)*LR), 2);

% Take real part fo diffusive problems (real eigenvalues): 
if ( isreal(L) == 1 )
    phi1 = real(phi1);
    phi2 = real(phi2);
    phi3 = real(phi3);
    phi4 = real(phi4);
    phi5 = real(phi5);
    phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
end

% Compute A:
A{3,2} = 16/15*phi2 + 16/5*phi3 + 16/5*phi4;

% Compute B:
B{2} = 32/15*(phi2 + phi3) - 64/5*phi4 - 128/5*phi5; 
B{3} = -1/3*phi2 + 1/3*phi3 + 5*phi4 + 8*phi5;

% Compute U:
U{2,1} = -2*(phit{2,2} + phit{3,2});
U{3,1} = -2/3*phi2 + 2*phi3 + 4*phi4;
U{2,2} = 1/2*phit{2,2} + phit{3,2};
U{3,2} = 1/10*phi2 - 1/5*phi3 - 6/5*phi4;

% Compute V:
V{1} = -1/3*phi2 + 5/3*phi3 - phi4 - 8*phi5;
V{2} = 1/30*phi2 - 2/15*phi3 - 1/5*phi4 + 8/5*phi5;

% Compute the missing oefficients using the summation properties of the 
% coefficients:
[A, B, U, V, E] = spinscheme.computeMissingCoeffs(A, B, C, U, V, L, phi1, ...
    phit, dt);

end