function [t, r, r_p, r_m] = trigpade(f, m, n, varargin)
%TRIGPADE   Trigonometric Pade approximation.
%   [P, Q, R_HANDLE] = CHEBPADE(F, M, N) computes polynomials P and Q of degree
%   M and N, respectively, such that the rational function P/Q is the type (M,
%   N) Chebyshev-Pade approximation of type Clenshaw-Lord to the CHEBFUN F. That
%   is, the Chebyshev series of P/Q coincides with that for the CHEBFUN F up to
%   the maximum possible order for the polynomial degrees permitted. R_HANDLE is
%   a function handle for evaluating the rational function.
%
%   [P, Q, R_HANDLE] = CHEBPADE(F, M, N, TYPE) allows one to additionally
%   specify the type of Chebyshev-Pade approximation sought. If TYPE is set to
%   'clenshawlord', the Clenshaw-Lord approximation as described above is used.
%   Alternatively, setting TYPE to 'maehly' computes a Maehly-type
%   approximation, which satisfies a linearized version of the Chebyshev-Pade
%   conditions.
%
%   [P, Q, R_HANDLE] = CHEBPADE(F, M, N, TYPE, K) uses only the K-th partial sum
%   in the Chebyshev expansion of F when computing the approximation. CHEPADE(F,
%   M, N, K) is shorthand for CHEBPADE(F, M, N, 'clenshawlord', K).
%
% See also PADEAPPROX.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Form the outputs.

% F is assumed on the default trig interval [-1, 1]:

% Extract the coefficients
coeffs  = f.coeffs;
c = (length(coeffs)-1)/2;

% Separate the +ve and -ve coefficients
coeffs_p = coeffs(c+1:end);
coeffs_m = coeffs(c+1:-1:1);

% Half the zeroth order coefficient
coeffs_p(1) = 1/2*coeffs_p(1);
coeffs_m(1) = 1/2*coeffs_m(1);

% Solve the two Pade approximation problems
[r_p, a_p, b_p] = padeapprox(coeffs_p, m, n, 1e-12);
[r_m, a_m, b_m] = padeapprox(coeffs_m, m, n, 1e-12);

% Construct the four trigonometric polynomials:
aa_p = [zeros(length(a_p)-1, 1); a_p];
bb_p = [zeros(length(b_p)-1, 1); b_p];
tdp = chebfun(aa_p, 'coeffs', 'trig' );
tnp = chebfun(bb_p, 'coeffs', 'trig' );

aa_m = [zeros(length(a_m)-1, 1); a_m];
aa_m = flipud(aa_m);
bb_m = [zeros(length(b_m)-1, 1); b_m];
bb_m = flipud(bb_m);

tdm = chebfun(aa_m, 'coeffs', 'trig' );
tnm = chebfun(bb_m, 'coeffs', 'trig' );


% Construct the full approximation:
t = tdp./tnp + tdm./tnm;

% Discard the imaginary rounding errors:
if ( norm(imag(t)) > norm(f) * 100 * eps )
	warning('Chebfun:trigpade:imag', 'imaginary part not negligible.');
else
	t = real(t);	
end

r = chebfun(@(x) t(x))
end
