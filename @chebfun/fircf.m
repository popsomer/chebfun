function [H, h, err, res] = fircf(order, freqs, f, varargin)
%FIRCF   FIR minimax filter design using CF for real valued chebfuns.
%   H = FIRCF(M, FREQS, F) is a periodic chebfun of length M+1 representing
%   an order M filter.
%   ORDER corresponds to a filter length of ORDER+1.
%   FREQS is a vector with numbers in [0, 1]


% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Check arguments:
if ( order <= 0 || order ~= round(order) )
    error( 'CHEBFUN:CHEBFUN:fircf', ...
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

%%
freqs = freqs(:);
intervals = sort(cos(pi*freqs));
nPieces = length(intervals)-1;
nDer = 2;
funPieces = cell(1, nPieces);
for i = 1:nPieces
    if ( rem(i, 2) == 1 )
        funPieces{i} = chebfun(g.funs{i});
    else
        a = intervals(i);
        b = intervals(i+1);
        s = sigmoid(nDer, a, b);      
        gl = feval(g, intervals(i));        
        gr = feval(g, intervals(i+1));
        funPieces{i} = gl + ( gr - gl ) * s;         
    end
end
H = chebfun(funPieces, intervals.' );
H = chebfun(H);
%%
[pcf, err] = cf(H, n/2);
%[p, err, status] = remez(g, n/2, 'ignore', intervals);
%status.xk = sort(1/pi*acos(status.xk));
a = chebcoeffs(pcf);
h = [1/2*a(end:-1:2); a(1); 1/2*a(2:end);];
H = chebfun(h, 'trig', 'coeffs');
res = [];
end


function f = sigmoid(k, a, b)
% Creates a C^k sigmoidal function on [a, b] and returns it
% as a chebufun.
% INPUT:
%    k: A non-negative integer or inf
% OUTPUT:
%    f is the chebfun sigmoidal function on [a, b] and the function is C^k

if ( isinf(k) )
    % Rescale and remap on the bump on the domain: 
    fh = @(t) CinfbumpHandle(t, a, b);    
    f = chebfun(@(t) fh(t), [a, b]);
else
    fh = @(x) (x-a).^k.*(b-x).^k;
    % Construct the degree 2*k polynomial on [-a,a]:
    f = chebfun(@(t) fh(t), [a, b], 2*k+1);
end
% Integrate to create the sigmoidal:
f = cumsum(f)/sum(f);
end


% Function to create the standard C-inf bump on [-1, 1];
function y = CinfbumpHandle(t, a, b)
y = exp(-1./((b-t).*(t-a)));
y(isnan(y)) = 0;
end
