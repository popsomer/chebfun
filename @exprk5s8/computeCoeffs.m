function [A, B, U, V, E] = computeCoeffs(scheme, L, LR, dt)
%COMPUTECOEFFS   Compute coefficients for EXPRK5S8.

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
C(4) = 1/4;
C(5) = 1/2;
C(6) = 1/5;
C(7) = 2/3;
C(8) = 1;

% Compute the phi functions:
phi1 = mean(feval(computePhi(1), LR), 2);
phi2 = mean(feval(computePhi(2), LR), 2);
phi3 = mean(feval(computePhi(3), LR), 2);
phi4 = mean(feval(computePhi(4), LR), 2);
phit{1,2} = mean(C(2)^1*feval(computePhi(1), C(2)*LR), 2);
phit{1,3} = phit{1,2};
phit{1,4} = mean(C(4)^1*feval(computePhi(1), C(4)*LR), 2);
phit{1,5} = phit{1,2};
phit{1,6} = mean(C(6)^1*feval(computePhi(1), C(6)*LR), 2);
phit{1,7} = mean(C(7)^1*feval(computePhi(1), C(7)*LR), 2);
phit{1,8} = phi1;
phit{2,2} = mean(C(2)^2*feval(computePhi(2), C(2)*LR), 2);
phit{2,4} = mean(C(4)^2*feval(computePhi(2), C(4)*LR), 2);
phit{2,6} = mean(C(6)^2*feval(computePhi(2), C(6)*LR), 2);
phit{2,7} = mean(C(7)^2*feval(computePhi(2), C(7)*LR), 2);
phit{3,2} = mean(C(2)^3*feval(computePhi(3), C(2)*LR), 2);
phit{3,6} = mean(C(6)^3*feval(computePhi(3), C(6)*LR), 2);
phit{3,7} = mean(C(7)^3*feval(computePhi(3), C(7)*LR), 2);
phit{4,6} = mean(C(6)^4*feval(computePhi(4), C(6)*LR), 2);
phit{4,7} = mean(C(7)^4*feval(computePhi(4), C(7)*LR), 2);

% Take real part fo diffusive problems (real eigenvalues): 
if ( isreal(L) == 1 )
    phi1 = real(phi1);
    phi2 = real(phi2);
    phi3 = real(phi3);
    phi4 = real(phi4);
    phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
end

% Compute A:
A{3,2} = 2*phit{2,2};
A{4,3} = 2*phit{2,4};
A{5,3} = -2*phit{2,2} + 16*phit{3,2};
A{5,4} = 8*phit{2,2} - 32*phit{3,2};
A{6,4} = 8*phit{2,6} - 32*phit{3,6};
A{7,4} = -(125/162)*A{6,4};
A{6,5} = -2*phit{2,6} + 16*phit{3,6};
A{7,5} = (125/1944)*A{6,4} - (4/3)*phit{2,7} + (40/3)*phit{3,7};
Phi = (5/32)*A{6,4} - (25/28)*phit{2,6} + (81/175)*phit{2,7} - (162/25)*phit{3,7} + ...
    (150/7)*phit{4,6} + (972/35)*phit{4,7} + 6*phi4;
A{8,5} = -(16/3)*phi2 + (208/3)*phi3 - 40*Phi;
A{7,6} = (3125/3888)*A{6,4} + (25/3)*phit{2,7} - (100/3)*phit{3,7};
A{8,6} = (250/21)*phi2 - (250/3)*phi3 + (250/7)*Phi;
A{8,7} = (27/14)*phi2 - 27*phi3 + (135/7)*Phi;

% Compute B:
B{6} = (125/14)*phi2 - (625/14)*phi3 + (1125/14)*phi4;
B{7} = -(27/14)*phi2 + (162/7)*phi3 - (405/7)*phi4;
B{8} = (1/2)*phi2 - (13/2)*phi3 + (45/2)*phi4;

% Compute the missing oefficients using the summation properties of the 
% coefficients:
[A, B, U, V, E] = spinscheme.computeMissingCoeffs(A, B, C, U, V, L, phi1, ...
    phit, dt);

end