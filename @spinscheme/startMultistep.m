function c = startMultistep(c, Nc, Nv, L, LR, dt, q)
%STARTMULTISTEP  Get enough initial data when using a multistep scheme.
%    C = STARTMULTISTEP(C, NC, NV, L, LR, DT, Q) does Q-1 steps of ETDRK4 to get
%    enough initial data to start a multistep scheme of order Q, using the
%    nonlinear part of the operator in coefficient space NC and in value space 
%    NV, the linear part L, the linear part for complex means LR, and the 
%    timestep DT.

% Author: Hadrien Montanelli.

% Number of points to discretize the interval:
N = length(c);

% Create a NxQ matrix to store the coefficients at the Q steps:
coeffs = zeros(N, q);

% Store the initial conidition in the last column:
coeffs(:,q) = c;

% Do Q-1 steps of ETRDK4:
scheme = etdrk4();
[A, B, U, V, E] = computeCoeffs(scheme, L, LR, dt);
for i = 1:q-1
    c = spinscheme.oneStep(c, Nc, Nv, A, B, U, V, E);
    coeffs(:,q-i) = c;
end
c = coeffs;

end