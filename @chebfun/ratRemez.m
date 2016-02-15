function varargout = ratRemez(varargin)
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
if ( n <= 0 ) 
    error('CHEBFUN:CHEBFUN:ratRemez:ratOnly', ...
        'Use remez for polynomial case' );
end

% Initial values for some parameters.
iter = 0;                 % Iteration count.
deltaLevelError = max(normf, eps);  % Value for stopping criterion.
deltamin = inf;           % Minimum error encountered.
deltaReference = 1;                % Maximum correction to trial reference.

N = m + n;
% Compute an initial reference set to start the algorithm.
[xk, pmin, qmin] = getInitialReference(f, m, n, N);

if ( isempty(qmin) )
    qmin = chebfun(1, f.domain([1, end]));
end



% Print header for text output display if requested.
if ( opts.displayIter )
    disp('It.     Max(|Error|)       |ErrorRef|      Delta ErrorRef      Delta Ref')
end


% Old reference and levelled error
x0 = xk;
h0 = 0;
% Run the main algorithm.
while ( (deltaLevelError/normf > opts.tol) && (iter < opts.maxIter) && (deltaReference > 0) )

    [p, q, rh, pqh, h, interpSuccess] = computeTrialFunctionRational(f, xk, m, n);      
    
    if ( abs(h) <= abs(h0) )
        % The levelled error has not increased
        disp('level error decreased' )
        xk = makeNewReference(xkPrev, xk);
    end
    % Perturb exactly-zero values of the levelled error.
    if ( h == 0 )
        h = 1e-19;
    end
    xkPrev = xk;
    % Update the exchange set using the Remez algorithm with full exchange.   
    [xk, err, err_handle] = exchange(xk, h, 2, f, p, q, rh, N + 2, opts);

    % Update max. correction to trial reference and stopping criterion.
    deltaReference = max(abs(x0 - xk));
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
        doPlotIter(x0, xk, err_handle, dom);
        pause();
    end

    if ( opts.displayIter )
        doDisplayIter(iter, err, h, deltaLevelError, normf, deltaReference);
    end
    
    if ( opts.demoMode )
        pause
    end

    % Save the old reference and
    % the old levelled error.
    x0 = xk;
    h0 = h;
    iter = iter + 1
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
status.deltaReference = deltaReference;
status.xk = xk;

p = simplify(p);
if ( opts.rationalMode )
    q = simplify(qmin);
    varargout = {p, q, @(x) feval(p, x)./feval(q, x), err, status};
else
    varargout = {p, err, status};
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

if ( numel(f.funs) == 1 )
    %[p, q] = chebpade(f, m, n);
    [p, q] = cf(f, m, n);
else
    %[p, q] = chebpade(f, m, n, 5*N);
    [p, q] = cf(f, m, n, 100*N);
end
pqh = @(x) feval(p, x)./feval(q, x);
[xk, err, e, flag] = exchange([], 0, 2, f, p, q, pqh, N + 2);

% If the above procedure failed to produce a reference
% with enough equioscillation points, just use the Chebyshev points.
if ( flag == 0 )
    xk = chebpts(N + 2, f.domain([1, end]), 1);
    xk = [xk(round(length(xk)/2)+1:length(xk)) - 1; xk(1:round(length(xk)/2))+1];
    xk(round(length(xk)/2))=-1;
    xk = sort(xk);    
end

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
function doDisplayIter(iter, err, h, delta, normf, deltaReference)

disp([num2str(iter), '        ', num2str(err, '%5.4e'), '        ', ...
    num2str(abs(h), '%5.4e'), '        ', ...
    num2str(delta/normf, '%5.4e'), '        ', num2str(deltaReference, '%5.4e')])
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


function [xk, norme, err_handle, flag] = exchange(xk, h, method, f, p, q, rh, Npts, opts)
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
% Rational case:

rr = findExtrema(f, p, q, rh, h, xk);
err_handle = @(x) feval(f, x) - rh(x);

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

end

function xk = makeNewReference(x0, x1)

xk = x0;
if ( length(x0) ~= length(x1) )
    error('makeNewReference:not equal points in both references')    
end

for i = 1:length(x0)
    % Find the index of the nearest point:
    [~, idx] = min(abs(x1-x0(i)));
    xk(i) = x0(i) + (x1(idx) - x0(i))/2;
end
end



