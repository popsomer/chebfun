function varargout = remez(varargin)
%REMEZ   Best polynomial or rational approximation for real valued chebfuns.
%   P = REMEZ(F, M) computes the best polynomial approximation of degree M to
%   the real CHEBFUN F in the infinity norm using the Remez algorithm.
%
%   [P, Q] = REMEZ(F, M, N) computes the best rational approximation P/Q of type
%   (M, N) to the real CHEBFUN F using the Remez algorithm.
%
%   [P, Q, R_HANDLE] = REMEZ(F, M, N) does the same but additionally returns a
%   function handle R_HANDLE for evaluating the rational function P/Q.
%
%   [...] = REMEZ(..., 'tol', TOL) uses the value TOL as the termination
%   tolerance on the increase of the levelled error.
%
%   [...] = REMEZ(..., 'display', 'iter') displays output at each iteration.
%
%   [...] = REMEZ(..., 'maxiter', MAXITER) sets the maximum number of allowable
%   iterations to MAXITER.
%
%   [...] = REMEZ(..., 'plotfcns', 'error') plots the error after each iteration
%   while the algorithm executes.
%
%   [P, ERR] = REMEZ(...) and [P, Q, R_HANDLE, ERR] = REMEZ(...) also returns
%   the maximum error ERR.
%
%   [P, ERR, STATUS] = REMEZ(...) and [P, Q, R_HANDLE, ERR, STATUS] = REMEZ(...)
%   also return a structure array STATUS with the following fields:
%      STATUS.DELTA  - Obtained tolerance.
%      STATUS.ITER   - Number of iterations performed.
%      STATUS.DIFFX  - Maximum correction in last trial reference.
%      STATUS.XK     - Last trial reference on which the error equioscillates.
%
%   This code is quite reliable for polynomial approximations but rather
%   fragile for rational approximations.  Better results can often be obtained
%   with CF(), especially if f is smooth.
%
% References:
%
%   [1] Pachon, R. and Trefethen, L. N.  "Barycentric-Remez algorithms for best
%   polynomial approximation in the chebfun system", BIT Numerical Mathematics,
%   49:721-742, 2009.
%
%   [2] Pachon, R.  "Algorithms for Polynomial and Rational Approximation".
%   D. Phil. Thesis, University of Oxford, 2010 (Chapter 6).
%
% See also CF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 ) 
    error('CHEBFUN:CHEBFUN:remez:nargin', ...
        'Not enough input arguments.');
else
    [f, m, n, opts] = parseInput(varargin{:});
end

if ( isempty(f) ) 
    varargout = {f};
    return;
end

if ( opts.rationalMode )
    [m, n] = adjustDegreesForSymmetries(f, m, n);
end

dom = f.domain([1, end]);
normf = norm(f);

% With zero denominator degree, the denominator polynomial is trivial.
if ( n == 0 ) 
    q = chebfun(1, dom);
    qmin = q;
end

% Initial values for some parameters.
iter = 0;                 % Iteration count.
deltaLevelError = max(normf, eps);  % Value for stopping criterion.
deltamin = inf;           % Minimum error encountered.
diffx = 1;                % Maximum correction to trial reference.

N = m + n;
% Compute an initial reference set to start the algorithm.
[xk, pmin, qmin] = getInitialReference(f, m, n, N);
if ( isempty(qmin) )
    qmin = chebfun(1, f.domain([1, end]));
end
if ( n > 0 )
%    h = norm(feval(f, xk) - feval(pmin,xk)./feval(qmin,xk));
end
xo = xk;

% Print header for text output display if requested.
if ( opts.displayIter )
    disp('It.     Max(|Error|)       |ErrorRef|      Delta ErrorRef      Delta Ref')
end



% Run the main algorithm.
while ( (deltaLevelError/normf > opts.tol) && (iter < opts.maxIter) && (diffx > 0) )
    fk = feval(f, xk);     % Evaluate on the exchange set.
    w = baryWeights(xk);   % Barycentric weights for exchange set.    
        
    % Compute trial function and levelled reference error.
    if ( n == 0 )
        [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, dom);
        opts.rh = [];
    else
        [p, q, h, rh] = computeTrialFunctionRational(fk, xk, w, m, n, N, dom);
        opts.rh = rh;
        opts.f = f;
    end
    
    % Perturb exactly-zero values of the levelled error.
    if ( h == 0 )
        h = 1e-19;
    end

    % Update the exchange set using the Remez algorithm with full exchange.   
    [xk, err, err_handle] = exchange(xk, h, 2, f, p, q, N + 2, opts);

    % If overshoot, recompute with one-point exchange.
    if ( err/normf > 1e5 )
        [xk, err, err_handle] = exchange(xo, h, 1, f, p, q, N + 2, opts);
    end

    % Update max. correction to trial reference and stopping criterion.
    diffx = max(abs(xo - xk));
    deltaLevelError = err - abs(h);

    % Store approximation with minimum norm.
    if ( deltaLevelError < deltamin )
        pmin = p;
        if ( n > 0 )
            qmin = q;
        end

        errmin = err;
        xkmin = xk;
        deltamin = deltaLevelError;
    end

    % Display diagnostic information as requested.
    if ( opts.plotIter )
        doPlotIter(xo, xk, err_handle, dom);
        pause();
    end

    if ( opts.displayIter )
        doDisplayIter(iter, err, h, deltaLevelError, normf, diffx);
    end
    
    if ( opts.demoMode )
        pause
    end

    xo = xk;
    iter = iter + 1;
end

% Take best results of all the iterations we ran.
p = pmin;
err = errmin;
xk = xkmin;
deltaLevelError = deltamin;

% Warn the user if we failed to converge.
if ( deltaLevelError/normf > opts.tol )
    warning('CHEBFUN:CHEBFUN:remez:convergence', ...
        ['Remez algorithm did not converge after ', num2str(iter), ...
         ' iterations to the tolerance ', num2str(opts.tol), '.']);
end

% Form the outputs.
status.delta = deltaLevelError/normf;
status.iter = iter;
status.diffx = diffx;
status.xk = xk;

p = simplify(p);
if ( opts.rationalMode )
    q = simplify(qmin);
    varargout = {p, q, @(x) feval(p, x)./feval(q, x), err, status};
else
    varargout = {p, err, status};
end

end

function [m, n] = adjustDegreesForSymmetries(f, m, n)
%ADJUSTDEGREESFORSYMMETRIES   Adjust rational approximation degrees to account
%   for function symmetries.
%
%   [M, N] = ADJUSTDEGREESFORSYMMETRIES(F, M, N) returns new degrees M and N to
%   correct the defect of the rational approximation if the target function is
%   even or odd.  In either case, the Walsh table is covered with blocks of
%   size 2x2, e.g.  for even function the best rational approximant is the same
%   for types [m/n], [m+1/n], [m/n+1] and [m+1/n+1], with m and n even. This
%   strategy is similar to the one proposed by van Deun and Trefethen for CF
%   approximation in Chebfun (see @chebfun/cf.m).

% Sample piecewise-smooth CHEBFUNs.
if ( (numel(f.funs) > 1) || (length(f) > 128) )
  f = chebfun(f, f.domain([1, end]), 128);
end

% Compute the Chebyshev coefficients.
c = chebcoeffs(f, length(f));
c(1) = 2*c(1);

% Check for symmetries and reduce degrees accordingly.
if ( max(abs(c(2:2:end)))/vscale(f) < eps )   % f is even.
    if ( mod(m, 2) == 1 )
        m = max(m - 1, 0);
    end
    if ( mod(n, 2) == 1 )
        n = max(n - 1, 0);
    end
elseif ( max(abs(c(1:2:end)))/vscale(f) < eps ) % f is odd.
    if ( mod(m, 2) == 0 )
        m = max(m - 1, 0);
    end
    if ( mod(n, 2) == 1 )
        n = max(n - 1, 0);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implementing the core part of the algorithm.

function [xk, p, q] = getInitialReference(f, m, n, N)

% If doing rational Remez, get initial reference from trial function generated
% by CF or Chebyshev-Pade.
flag = 0;
a = f.domain(1);
b = f.domain(end);

if ( n > 0 )
    if ( numel(f.funs) == 1 )
        %[p, q] = chebpade(f, m, n);
        [p, q] = cf(f, m, n);
    else
        %[p, q] = chebpade(f, m, n, 5*N);
        [p, q] = cf(f, m, n, 100*N);
    end
    opts = struct;
    opts.rh = [];
    [xk, err, e, flag] = exchange([], 0, 2, f, p, q, N + 2, opts);
end
% REMOVE THIS:
%flag = 0;

% In the polynomial case or if the above procedure failed to produce a reference
% with enough equioscillation points, just use the Chebyshev points.
if ( flag == 0 )
    xk = chebpts(N + 2, f.domain([1, end]), 1);
    %mid = round(length(xk)/2);
    %xk = [xk(1:mid)+1; xk(mid+1:end)-1];
    %xk = sort(xk);
    p = [];
    q = [];        
end


end

function [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, dom)

% Vector of alternating signs.
sigma = ones(N + 2, 1);
sigma(2:2:end) = -1;

h = (w'*fk) / (w'*sigma);                          % Levelled reference error.
pk = (fk - h*sigma);                               % Vals. of r*q in reference.

% Trial polynomial.
p = chebfun(@(x) bary(x, pk, xk, w), dom, m + 1);

end

function [p, q, h, rh] = computeTrialFunctionRational(fk, xk, w, m, n, N, dom)

% Vector of alternating signs.
sigma = ones(N + 2, 1);
sigma(2:2:end) = -1;

% Orthogonal matrix with respect to <,>_{xk}.
[C, ignored] = qr(fliplr(vander(xk)));

% Left and right rational interpolation matrices.
ZL = C(:,m+2:N+2).'*diag(fk)*C(:,1:n+1);
ZR = C(:,m+2:N+2).'*diag(sigma)*C(:,1:n+1);

% Solve generalized eigenvalue problem.
[v, d] = eig(ZL, ZR);

% Compute all possible qk and and look for ones with unchanged sign.
qk_all = C(:,1:n+1)*v;
pos =  find(abs(sum(sign(qk_all))) == N + 2);  % Sign changes of each qk.

if ( isempty(pos) || (length(pos) > 1) )
    error('CHEBFUN:CHEBFUN:remez:badGuess', ...
        'Trial interpolant too far from optimal');
end

qk = qk_all(:,pos);       % Keep qk with unchanged sign.
h = d(pos, pos);          % Levelled reference error.
pk = (fk - h*sigma).*qk;  % Vals. of r*q in reference.

% Trial numerator and denominator.
[xk_leja, idx] = leja(xk, 1, m+1);
pk_leja = pk(idx);
w_leja = baryWeights(xk_leja);
p = chebfun(@(x) bary(x, pk_leja, xk_leja, w_leja), dom, m + 1);

[xk_leja, idx] = leja(xk, 1, n+1);
qk_leja = qk(idx);
w_leja = baryWeights(xk_leja);
q = chebfun(@(x) bary(x, qk_leja, xk_leja, w_leja), dom, n + 1);

nn = round(length(xk)/2);
fvals = fk - h*sigma;
xx = xk; xx(nn) = [];
fx = fvals; fx(nn) = [];
A = berrut(xx, fx, m, n);
v = null(A);
rh = @(t) bary(t, fx, xx, v);
%r = chebfun(fh, dom, 'splitting', 'on');

end

function [xk, norme, err_handle, flag] = exchange(xk, h, method, f, p, q, Npts, opts)
%EXCHANGE   Modify an equioscillation reference using the Remez algorithm.
%   EXCHANGE(XK, H, METHOD, F, P, Q, W) performs one step of the Remez algorithm
%   for the best rational approximation of the CHEBFUN F of the target function
%   according to the first method (METHOD = 1), i.e. exchanges only one point,
%   or the second method (METHOD = 2), i.e. exchanges all the reference points.
%   XK is a column vector with the reference, H is the levelled error, P is the
%   numerator, and Q is the denominator of the trial
%   rational function P/Q and W is the weight function.
%
%   [XK, NORME, E_HANDLE, FLAG] = EXCHANGE(...) returns the modified reference
%   XK, the supremum norm of the error NORME (included as an output argument,
%   since it is readily computed in EXCHANGE and is used later in REMEZ), a
%   function handle E_HANDLE for the error, and a FLAG indicating whether there
%   were at least N+2 alternating extrema of the error to form the next
%   reference (FLAG = 1) or not (FLAG = 0).
%
%   [XK, ...] = EXCHANGE([], 0, METHOD, F, P, Q, N + 2) returns a grid of N + 2
%   points XK where the error F - P/Q alternates in sign (but not necessarily
%   equioscillates). This feature of EXCHANGE is useful to start REMEZ from an
%   initial trial function rather than an initial trial reference.

% Compute extrema of the error.
if ( ~isempty(opts.rh) )
    % Rational case:
    r = opts.rh;    
    err_handle = @(x) feval(f, x) - r(x);    
    rts = [];
    %doms = unique(sort([f.domain(:); xk]));
    doms = unique([f.domain(1); xk; f.domain(end)]);
    %doms = sort([doms; 0]);
    for k = 1:length(doms)-1
        ek = chebfun(@(x) err_handle(x), [doms(k), doms(k+1)], 'splitting', 'on' );
        %plot(ek)
        rts = [rts; roots(diff(ek), 'nobreaks')];  %#ok<AGROW>
    end    
else
    e_num = (q.^2).*diff(f) - q.*diff(p) + p.*diff(q);
    rts = roots(e_num, 'nobreaks');
    % Function handle output for evaluating the error.
    err_handle = @(x) feval(f, x) - feval(p, x)./feval(q, x);
end

rr = [f.domain(1) ; rts; f.domain(end)];

% Select exchange method.
if ( method == 1 )                           % One-point exchange.
    [ignored, pos] = max(abs(feval(err_handle, rr)));
    pos = pos(1);
else                                           % Full exchange.
    pos = find(abs(err_handle(rr)) >= abs(h)); % Values above levelled error
end

% Add extrema nearest to those which are candidates for exchange to the
% existing exchange set.
[r, m] = sort([rr(pos) ; xk]);
v = ones(Npts, 1);
v(2:2:end) = -1;
er = [feval(err_handle, rr(pos)) ; v*h];
er = er(m);

% Delete repeated points.
repeated = diff(r) == 0;
r(repeated) = [];
er(repeated) = [];



% Determine points and values to be kept for the reference set.
s = r(1);    % Points to be kept.
es = er(1);  % Values to be kept.
for i = 2:length(r)
    if ( (sign(er(i)) == sign(es(end))) && (abs(er(i)) > abs(es(end))) )
        % Given adjacent points with the same sign, keep one with largest value.
        s(end) = r(i);
        es(end) = er(i);
    elseif ( sign(er(i)) ~= sign(es(end)) )
        % Keep points which alternate in sign.
        s = [s ; r(i)];    %#ok<AGROW>
        es = [es ; er(i)]; %#ok<AGROW>
    end
end

% Of the points we kept, choose n + 2 consecutive ones 
% that include the maximum of the error.
[norme, index] = max(abs(es));
d = max(index - Npts + 1, 1);
if ( Npts <= length(s) )
    xk = s(d:d+Npts-1);
    flag = 1;

else
    xk = s;
    flag = 0;
end

% if (extraPts > 0)
    
% else
%     xk = s;
%     flag = 0;
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for displaying diagnostic information.

% Function called when opts.plotIter is set.
function doPlotIter(xo, xk, err_handle, dom)

xxk = linspace(dom(1), dom(end), max(3000, 10*length(xk)));
plot(xo, 0*xo, '.r', 'MarkerSize', 6)   % Old reference.
holdState = ishold;
hold on
plot(xk, 0*xk, 'ok', 'MarkerSize', 3)   % New reference.

plot(xxk, err_handle(xxk))               % Error function.
if ( ~holdState )                        % Return to previous hold state.
    hold off
end
xlim(dom)
legend('Current Ref.', 'Next Ref.', 'Error')
drawnow

end

% Function called when opts.displayIter is set.
function doDisplayIter(iter, err, h, delta, normf, diffx)

disp([num2str(iter), '        ', num2str(err, '%5.4e'), '        ', ...
    num2str(abs(h), '%5.4e'), '        ', ...
    num2str(delta/normf, '%5.4e'), '        ', num2str(diffx, '%5.4e')])

end

function status = remezParseFunction(f)
% Parse a Chebfun to see if Remez
% can be applied to it

% Add conditions needed on f:
if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN:remez:real', ...
        'REMEZ only supports real valued functions.');    
end

if ( numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:remez:quasi', ...
        'REMEZ does not currently support quasimatrices.');    
end

if ( issing(f) )
    error('CHEBFUN:CHEBFUN:remez:singularFunction', ...
        'REMEZ does not currently support functions with singularities.');
end

% If all are satisifed, we can go ahead:
status = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [f, m, n, opts] = parseInput(varargin)

%% Handle the first two arguments:
f = varargin{1};
remezParseFunction(f);
varargin(1) = [];

m = varargin{1};
varargin(1) = [];

n = 0;
opts.rationalMode = false;
if ( ~isempty(varargin) > 0 )
    if ( isa(varargin{1}, 'double') )
        n = varargin{1};
        varargin(1) = [];
        opts.rationalMode = true;
    end
end

if ( isa(m, 'vector') || isa(n, 'vector') || ...
        m < 0 || n < 0 || m ~= round(m) || n ~= round(n) )
   error('CHEBFUN:CHEBFUN:remez:degree', ...
        'Approximant degree must be a non-negative integer.');
end


% Parse name-value option pairs.
N = m + n;
opts.tol = 1e-16*(N^2 + 10); % Relative tolerance for deciding convergence.
opts.maxIter = 40;           % Maximum number of allowable iterations.
opts.displayIter = false;    % Print output after each iteration.
opts.plotIter = false;       % Plot approximation at each iteration.
opts.demoMode = false;

while ( ~isempty(varargin) )
    if ( strcmpi('tol', varargin{1} ) )
        opts.tol = varargin{2};
        varargin(1) =[];
        varargin(1) =[];        
    elseif ( strcmpi('maxiter', varargin{1}) )
        opts.maxIter = varargin{2};        
        varargin(1) =[];
        varargin(1) =[];        
    elseif ( strcmpi('display', varargin{1}) )
        varargin(1) = [];
        if ( strcmpi('iter', varargin{1}) )
            varargin(1) = [];
        end
        opts.displayIter = true;        
    elseif ( strcmpi('demo', varargin{1}) )        
        varargin(1) = [];
        opts.demoMode = true;
        opts.displayIter = true;
        opts.plotIter = true;        
    elseif ( strcmpi('plotfcns', varargin{1}) )
        varargin(1) = [];
        if ( strcmpi('error', varargin{1}) )
            varargin(1) = [];
        end            
        opts.plotIter = true;        
    else
        error('CHEBFUN:CHEBFUN:remez:badInput', ...
            'Unrecognized sequence of input parameters.')
    end
        
end

end


function [xx, pos] = leja(x, startIndex, nPts) 
% put NPTS from X in a Leja sequence
% starting from x(startIndex)
n = length(x);
p = zeros(n,1);
pos = zeros(nPts, 1);
xx = zeros(nPts, 1);
xx(1) = x(startIndex); 
pos(1) = startIndex;

for j = 2:nPts
    % we want to pick the jth point now:
    for i = 1:n
        %p(i) = prod(abs(x(i) - xx(1:j-1)));
        p(i) = sum(log(abs(x(i) - xx(1:j-1)))); % no overflow
    end  
    [val,pos(j)] = max(p);
    xx(j) = x(pos(j));
end

end

function r = mergePoints(rLeja, rOther, erOther, Npts)
rLeja   = rLeja(:);
rOther  = rOther(:);
erOther = erOther(:);

idx = rOther < rLeja(1);
rTemp = rOther(idx);
erTemp = erOther(idx);
[~, pos] = max(abs(erTemp));
r = rTemp(pos);
i = 1;
while i < length(rLeja)
    r = [r; rLeja(i)]; 
    k = i+1;
    while ( ~any((rOther > rLeja(i)) & (rOther < rLeja(k))) )
        k = k + 1;
    end
    idx = (rOther > rLeja(i)) & (rOther < rLeja(k));    
    rTemp = rOther(idx);
    erTemp = erOther(idx);
    [~, pos] = max(abs(erTemp));
    r = [r; rTemp(pos)];               
    i = k;
end
r = [r; rLeja(end)];
idx = rOther > rLeja(end);
rTemp = rOther(idx);
erTemp = erOther(idx);
[~, pos] = max(abs(erTemp));
r = [r; rTemp(pos)];
if ( length(r) ~= Npts )
    warning('You are likely to fail my friend.')
end

end

function A = berrut(x, f, m, n)
x = x(:); x = x.';
f = f(:); f = f.';
A = zeros(m+n, m+n+1);
for i = 1:m
    A(i, :) = x.^(i-1);
end

for i = 1:n
    A(m+i,:) = f.*x.^(i-1);
end
end

function L = lowner(y, x, r, N)

L = zeros(r, N-r);

for i = 1:r
    for j = r+1:N
        L(i,j-r) = (y(i) - y(j))/(x(i)-x(j));
    end
end

end