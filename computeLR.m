function LR = computeLR(L, N, M, dt)
%COMPUTELR  Create a contour around each e-value of dt*L.
%   LR = COMPUTELR(L, N, M, dt) outputs a NxM matrix to be used for the complex
%   means. L is the operator (diagonal matrix stored as a vector), N is the 
%   number points to discretize the interval, M is the number of points to 
%   discretize the contour, and dt is the timestep.

% Author: Hadrien Montanelli.

% Roots of unity:
if ( isreal(L) == 1 )
    % Use only upper-half circle when eigenvalues are real:
    r = exp(1i*pi*((1:M) - .5)/M);
else
    r = exp(2i*pi*((1:M) - .5)/M);
end

% Move each root of unity around each entry of dt*L:
LR = dt*repmat(L, 1, M) + repmat(r, N, 1);

end