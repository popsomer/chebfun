function varargout = gradient(G,N,f,prefs)





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

% compute adjoints
bcG = getBCType(G);
Gu_star = adjoint(Gu,bcG);
Gf_star = adjoint(Gf,bcG);
bcN = getBCType(N);
Nu_star = adjoint(Nu,bcN);

% compute lambda
one = chebfun('1',G.domain);
prefs.discretization = @chebcolloc2;
lambda = linsolve(Nu_star,Gu_star*one,prefs);
lambda = chebfun(lambda);

% compute dg
dg = lambda + Gf_star*one;
dg = chebfun(dg);

% outputs
if nargout <= 1
    varargout = {dg};
elseif nargout == 2
    varargout = {g,dg};
    elseif nargout == 3
    varargout = {g,dg,lambda};
else
% error
end

end
