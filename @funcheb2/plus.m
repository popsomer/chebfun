function f = plus(f, g)
%+	Addition of two FUNCHEB2 objects.
%   F + G adds F and G, where F and G may be FUNCHEB2 objects or scalars.
%
% See also MINUS, UPLUS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % FUNCHEB2 + [] = []
    
    f = [];

elseif ( isa(g,'double') ) % FUNCHEB2 + double
    
    % Update values:
    f.values = f.values + g; 
    % Update coeffs:
    f.coeffs(end,:) = f.coeffs(end,:) + g;
    % Update scale:
    f.vscale = max(f.vscale, norm(f.values(:), inf)); 
    
elseif ( isa(f,'double') ) % double + FUNCHEB2
    
    % Switch argument order and call FUNCHEB2/plus again:
    f = plus(g, f);
    
else % FUNCHEB2 + FUNCHEB2
    
    % Make both FUNCHEB2 objects have the same length:
    nf1 = size(f.values, 1);
    nf2 = size(g.values, 1);
    if ( nf1 > nf2 )
        % Increase the length of f2 (via prolong):
        g = prolong(g, nf1);
    elseif ( nf1 < nf2 )
        % Increase the length of f1 (via prolong):
        f = prolong(f, nf2);
    end
    
    % Update values and coefficients:
    f.values = f.values + g.values;
    f.coeffs = f.coeffs + g.coeffs;
    
    % Update epslevel:
    f.epslevel = max(f.epslevel, g.epslevel);
    
    % Update scales:
    f.vscale = max([f.vscale, g.vscale, norm(f.values(:), inf)]);
    
    % Look for a zero output:
    if ( ~any(f.values(:)) || ~any(f.coeffs(:)) )
        % Creates a zero funcheb2:
        f = funcheb2(zeros(1, size(f.values, 2)), f.vscale, f.epslevel);
    end
  
end

end
