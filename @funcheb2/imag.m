function f = imag(f)
%IMAG	Imaginary part of a FUNCHEB2.
%   IMAG(F) is the imaginary part of F.
%
%   See also REAL, ISREAL, CONJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the imaginary part of the values:
f.values = imag(f.values);

if ( ~any(f.values(:)) )
    
    % Input was real, so output a zero FUNCHEB2:
    f = funcheb2(zeros(1, size(f.values, 2)), f.vscale, f.epslevel);
    
else
    
    % Compute imaginary part of the coefficients:
    f.coeffs = imag(f.coeffs);
    
end
