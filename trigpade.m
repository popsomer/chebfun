function [t, r_p, r_m] = trigpade(f, m, n, varargin)
%TRIGPADE   Trigonometric Pade approximation.
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
[r_p, a_p, b_p] = padeapprox(coeffs_p, m, n);
[r_m, a_m, b_m] = padeapprox(coeffs_m, m, n);

%% Construct the four trigonometric polynomials:

% padd coefficients with zeros:
aa_p = [zeros(length(a_p)-1, 1); a_p];
bb_p = [zeros(length(b_p)-1, 1); b_p];

% denonminator and numerator for the +ve part
tdp = chebfun(aa_p, 'coeffs', 'trig' );
tnp = chebfun(bb_p, 'coeffs', 'trig' );

% padd coefficients with zeros
aa_m = [zeros(length(a_m)-1, 1); a_m];
aa_m = flipud(aa_m);
bb_m = [zeros(length(b_m)-1, 1); b_m];
bb_m = flipud(bb_m);

% denonminator and numerator for the -ve part
tdm = chebfun(aa_m, 'coeffs', 'trig' );
tnm = chebfun(bb_m, 'coeffs', 'trig' );


% Construct the full approximation:
t = tdp./tnp + conj(tdp./tnp); %tdm./tnm;

% Discard the imaginary rounding errors:
if ( norm(imag(t)) > norm(f) * 100 * eps )
	warning('Chebfun:trigpade:imag', 'imaginary part not negligible.');
else
	t = real(t);	
end

end
