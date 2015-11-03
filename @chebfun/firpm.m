function [H, h, err, res] = firpm(order, freqs, f, varargin)
%FIRPM   FIR minimax filter design for real valued chebfuns.
%   H = FIRPM(M, FREQS, F) is a periodic chebfun of lenght M+1 representing
%   an order M filter.
%   ORDER corresponds to a filter length of ORDER+1.
%   FREQS is a vector with numbers in [0, 1]


% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Check arguments:
if ( order <= 0 || order ~= round(order) )
    error( 'CHEBFUN:CHEBFUN:firpm', ...
            'order must be a positive integer.');
else
    n = order;
end

if ( isa(f, 'chebfun') )
    dom = domain(f);
    if ( norm([dom(1), dom(end)] - [0, 1], inf) > 100*eps )
        error( 'CHEBFUN:CHEBFUN:firpm', ...
            'F must live on [0, 1]');
    end
    % Map the problem from the circle onto [-1, 1]:
    g = chebfun(@(x) feval(f, 1/pi*acos(x)), 'splitting', 'on');
end

if ( isa(f, 'function_handle') )    
    % Map the problem from the circle onto [-1, 1]:
    g = chebfun(@(x) feval(f, 1/pi*acos(x)), 'splitting', 'on');
end

if ( isvector(f) )
    f = double(f);
    f = chebfun.interp1(freqs, f, 'linear', [0, 1]);
    g = chebfun(@(x) feval(f, 1/pi*acos(x)), 'splitting', 'on');
end

% Remove the first and the last point from freqs which are
% just the end-points 0 and 1:
freqs = freqs(2:end-1);
freqs = freqs(:);
intervals = sort(cos(pi*freqs));
[p, err, status] = remez(g, n/2, 'ignore', intervals, varargin{:});
status.xk = sort(1/pi*acos(status.xk));
a = chebcoeffs(p);
h = [1/2*a(end:-1:2); a(1); 1/2*a(2:end);];
H = chebfun(h, 'trig', 'coeffs');
res = status;
end
