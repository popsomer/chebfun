function normA = norm(A, n)
%NORM   Norm of a CHEBMATRIX object.
%   NORM(A) computes the Frobenius norm of the CHEBMATRIX object A, defined as
%   the sum of the squares of the 2-norms of each of the blocks.
%
%   NORM(A, 2) or NORM(A, 'fro') is the same as above.
%
%   NORM(A, INF) computes the infinity norm of the CHEBMATRIX A, defined as the
%   maximum infinity norm of each of the blocks.
%
%   See also CHEBMATRIX, CHEBFUN/NORM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Add support for norms of operators (inf x inf blocks).

% Empty CHEBMATRIX has norm 0.
if ( isempty(A) )
    normA = 0;
    return
end

if ( nargin == 1 )
    n = 'fro'; 	% Frobenius norm is the default.
end

% The norm of a chebmatrix with inf x inf block(s) is not supported.
s = cellfun(@(b) min(size(b)), A.blocks);
if ( ~all(isfinite(s(:))) )
    error('CHEBFUN:chebmatrix:norm', ...
    'Norm of a chebmatrix with inf x inf block(s) is not supported.')
end

% Initialise.
normA = 0;

% Deal with different cases.
switch n
    
    case {'fro', 2}
        for k = 1:numel(A.blocks)
            normA = normA + norm(A.blocks{k}, 2).^2;
        end
        normA = sqrt(normA);
        
    case {inf, 'inf'}
        for k = 1:numel(A.blocks)
            normA = max(normA, norm(A.blocks{k}, inf));
        end
        
end

end