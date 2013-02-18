function f = rdivide(f, c, pref)
%./   Right array divide for a FUNCHEB2
%   F ./ C Divides the FUNCHEB2 F by a an array C. If F is a vectorised FUNCHEB2
%   with M columns, then C must be either a scalar or a 1xM array. 

%   Alternatively if C is a FUNCHEB2 with the same number of columns as F, or if
%   F is a scalar, the resulting division is returned if C is found to have no
%   roots in [-1,1]. The division is performed columnwise.
%
% See also MRDIVIDE, TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(c, 'double') )
    % Dividing by a constant is easy.
    
    % This can never work (as size(f, 1) == inf):
    if ( size(c, 1) > 1 || ...
           ( numel(c) > 1 && size(f.values, 2) ~= size(c, 2) ) )
        error('CHEBUFN:FUNCHEB2:rdivide:size', ...
            'Matrix dimensions must agree.');
    end
    
    if ( ~any(c(:)) )  
        f = funcheb2(NaN(1, size(f, 2)));
        
    elseif ( numel(c) == 1 )
        % Scalar
        f.values = f.values/c;      % Divide values
        f.coeffs = f.coeffs/c;      % Divide coeffs
        f.vscale = f.vscale/abs(c); % Divide vscale
        
    else
        % Vectorised 
        n = size(f.values, 1);   
        f.values = f.values./repmat(c, n, 1);   % Divide values
        f.coeffs = f.coeffs./repmat(c, n, 1);   % Divide coeffs
        f.vscale = f.vscale/norm(c, inf);       % Divide vscale
        
        f.values(:, c == 0) = NaN;
        f.coeffs(:, c == 0) = NaN;
        
    end
    
else
    % Dividing by another FUNCHEB2 is harder. Call COMPOSE.
    
    % Obtain preferencess:
    if ( nargin < 3 )
        pref = funcheb2.pref;
    end
    
    % Check if c (the denominator) has any roots.
    if ( ~isempty(roots(c)) )
        error('CHEBFUN:FUNCHEB2:rdivide:DivideByZeros', ...
        'Cannot divide by a FUNCHEB2 with roots in [-1, 1].');
    else
%         warning('CHEBFUN:FUNCHEB2:RDIVIDE:DivideByFun2', ...
%         'Division by a FUNCHEB2');
    end
    
    % Call COMPOSE.
    if ( isa(f, 'funcheb2') )  % funcheb2 / funcheb2
        f = compose(f, @rdivide, c, pref);
    else                       % double / funcheb2
        op = @(x) f./x;
        f = compose(c, op, pref);
    end
end

end
