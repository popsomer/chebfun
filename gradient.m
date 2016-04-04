function varargout = gradient(G,N,f,prefs)
%GRADIENT  Compute functional gradient using the adjiont method.
%   DG = GRADIENT(G,N,F) returns the derivative with respect to F of the
%   following constrained minimization problem:
%     min_(u,f) sum(G(u,f)) subject to N(u) = f,
%   where G and N are CHEBOPS.
%
%   [G,DG] = GRADIENT(G,N,F) also returns the value of SUM(G(U,F)).
%
%   [G,DG,LAMBDA] = GRADIENT(G,N,F) also returns the adjoint variable
%   that was computed using the adjiont method.
%
%   [...] = GRADIENT(G,N,F,PREFS) takes an optional CHEBOPPREF object
%   that can be used during the internal system solves.
%
% Example:
%   G = chebop(@(x,u,f) u.^2 + f.^2, [0 1]);  
%   N = chebop(@(u) diff(u, 2) + sin(u), [0 1]);
%   N.lbc = @(u) [ u-1; diff(u) ];
%   f = chebfun('x', [0,1]);
%   dg = gradient(G,N,f);
%
% See also CHEBOP/ADJOINT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% check for prefs argument
if nargin < 4
    prefs = cheboppref();
end

% solve for u
u = mldivide(N,f,prefs);

% compute cost 
g = sum(G(u,f));

% compute Frechet derivatives
Guf = linearize(G,[u;f]);
Gu = linop(Guf(1,1));
Gf = linop(Guf(1,2));
Nu = linearize(N,u);

% compute adjoints of Frechet derivatives of G
bcG = getBCType(G);
Gu_star = adjoint(Gu,bcG);
Gf_star = adjoint(Gf,bcG);

% compute adjoint of Frechet derivative of N
bcN = getBCType(N);
[~,op,bcOpL,bcOpR,bcOpM] = adjoint(Nu,bcN);
Nu_star = chebop(N.domain);
Nu_star.op = op;
Nu_star.lbc = bcOpL;
Nu_star.rbc = bcOpR;
Nu_star.bc = bcOpM;

% compute lambda
one = chebfun('1',G.domain);
lambda = mldivide(Nu_star,Gu_star*one,prefs);

% compute dg
dg = chebfun(lambda + Gf_star*one);

% outputs
if nargout <= 1
    % return just the derivative
    varargout = {dg};
elseif nargout == 2
    % return functional value and derivative
    varargout = {g,dg};
    elseif nargout == 3
    %  return functional value, derivative and adjiont variable
    varargout = {g,dg,lambda};
else
    % give error if more than 3 outputs
    error('CHEBFUN:gradient','Maximum of 3 outputs allowed.');
end

end
