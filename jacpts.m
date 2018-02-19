%n = 161; alpha = -0.7; beta = 1/sqrt(2); [xStd, wStd, vStd] = jacpts(n,alpha,beta);
%figure;for T =1:2:7,[xRH,wRH,vRH] =jacpts(n,alpha,beta,[-1,1], 'exp', [T,0,n+1]); semilogy(abs(xStd-xRH));hold on; end; legend(num2str((1:2:7)'));


function [x, w, v] = jacpts(n, a, b, int, meth, param)
% function [x, w, v] = jacpts(n, a, b, int, meth)
%JACPTS  Gauss-Jacobi quadrature nodes and weights.
%   X = JACPTS(N, ALPHA, BETA) returns the N roots of the degree N Jacobi
%   polynomial with parameters ALPHA and BETA (which must both be > -1)
%   where the Jacobi weight function is w(x) = (1-x)^ALPHA*(1+x)^BETA.
%
%   [X, W] = JACPTS(N, ALPHA, BETA) returns also a row vector W of weights.
%
%   [X, W, V] = JACPTS(N, ALPHA, BETA) returns additionally a column vector V of
%   weights in the barycentric formula corresponding to the points X.
%
%   JACPTS(N, ALPHA, BETA, INTERVAL, METHOD) or JACPTS(N, ALPHA, BETA, METHOD)
%   allows the user to select which method to use.
%    METHOD = 'REC' uses the recurrence relation for the Jacobi polynomials
%     and their derivatives to perform Newton iteration on the WKB approximation
%     to the roots. Default for N < 100.
%    METHOD = 'ASY' uses the Hale-Townsend fast algorithm based upon asymptotic
%     formulae, which is fast and accurate. Default for N >= 100.
%    METHOD = 'GW' uses the traditional Golub-Welsch eigenvalue method,
%     which is maintained mostly for historical reasons.
%
%   [X, W, V] = JACPTS(N, ALPHA, BETA, [A, B]) scales the nodes and weights for
%       the finite interval [A,B].
%
%   The cases ALPHA = BETA = -.5 and ALPHA = BETA = .5 correspond to
%   Gauss-Chebyshev nodes and quadrature, and are treated specially (as a closed
%   form expression for the nodes and weights is available). ALPHA = BETA = 0
%   calls LEGPTS, which is more efficient. The other cases with ALPHA = BETA call
%   ULTRAPTS, which is also faster.
% 
%   When ALPHA ~= BETA and MAX(ALPHA, BETA) > 5 the results may not be accurate. 
%
% See also CHEBPTS, LEGPTS, LOBPTS, RADAUPTS, HERMPTS, LAGPTS, TRIGPTS, and
% ULTRAPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'REC' by Nick Hale, July 2011
% 'ASY' by Nick Hale & Alex Townsend, May 2012 - see [2].
% 'RH' by Peter Opsomer, February 2015 - see [3].
% 'EXP' by Peter Opsomer, February 2018 - see [4].
%
% NOTE: The subroutines DO NOT SCALE the weights (with the exception of GW).
% This is done in the main code to avoid duplication.
%
% References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature
%       rules", Math. Comp. 23:221-230, 1969.
%   [2] N. Hale and A. Townsend, "Fast computation of Gauss-Jacobi 
%       quadrature nodes and weights", SISC, 2012.
%   [3] A. Deano, D. Huybrechs and P. Opsomer, Construction and implementation of asymptotic
%       expansions for Jacobi-type orthogonal polynomials", Adv. Comput. Math. 42 (4), 2016
%   [4] D. Huybrechs and P. Opsomer, "Arbitrary-order asymptotic expansions of
%       generalized Gaussian quadrature rules", (in preparation).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
interval = [-1, 1];
method = 'default';
method_set = 0;

if ( a <= -1 || b <= -1 )
    error('CHEBFUN:jacpts:sizeAB', 'Alpha and beta must be greater than -1.')
elseif (a~=b && max(a, b) > 5 )
    warning('CHEBFUN:jacpts:largeAB',...
        'ALPHA~=BETA and MAX(ALPHA, BETA) > 5. Results may not be accurate.')
end


% Check inputs:
if ( nargin > 3 )
%     if ( nargin == 5 )
    if ( nargin >= 5 ) % TODO remove temp extra args param etc
        % Calling sequence = JACPTS(N, INTERVAL, METHOD)
        interval = int;
        method = meth;
        method_set = 1;
    elseif ( nargin == 4 )
        if ( ischar(int) )
            % Calling sequence = JACPTS(N, METHOD)
            method = int;
            method_set = true;
        else
            % Calling sequence = JACPTS(N, INTERVAL)
            interval = int;
        end
    end
%     validStrings = {'default', 'GW', 'ASY', 'REC'};
    validStrings = {'default', 'GW', 'ASY', 'REC', 'RH', 'EXP'};
    if ( ~any(strcmpi(method, validStrings)) )
        if ( strcmpi(method, 'GLR') )
            error('CHEBFUN:jacpts:glr', ...
                'The GLR algorithm is no longer supported.');
        end
        error('CHEBFUN:jacpts:inputs', ['Unrecognised input string: ', method]);
    end
end

if ( any(isinf(interval)) )  % Inf intervals not yet supported.
    % TODO: How do we scale the weights?
    error('CHEBFUN:jacpts:infinterval', ... 
          'JACPTS() does not yet support infinite intervals');
elseif ( numel(interval) > 2 )
    warning('CHEBFUN:legpts:domain',...
        'Piecewise intervals are not supported and will be ignored.');
    interval = interval([1, end]);
end

% Deal with trivial cases:
if ( n < 0 )
    error('CHEBFUN:jacpts:n', 'First input should be a positive number.');
elseif ( n == 0 )   % Return empty vectors if n == 0:
    x = []; 
    w = []; 
    v = []; 
    return
elseif ( n == 1 )
    x0 = (b-a)/(a+b+2);
    x = diff(interval)/2 * (x0+1) + interval(1); % map from [-1,1] to interval. 
    w = 2^(a+b+1)*beta(a+1, b+1) * diff(interval)/2;
    v = 1;
    return
end

% Special cases:
% if ( a == b )
if ( a == b ) && ~strcmpi(method, 'rh') && ~strcmpi(method, 'exp')
    if ( a == 0 )  % Gauss-Legendre: alpha = beta = 0
        [x, w, v] = legpts(n, method);
        [x, w] = rescale(x, w, interval, a, b);
        return
    elseif ( a == -.5 )  % Gauss-Chebyshev: alpha = beta = -.5
        [x, ignored, v] = chebpts(n, interval, 1);
        w = repmat(pi/n,1,n);
        [ignored, w] = rescale(x, w, interval, a, b);
        return
    elseif ( a == .5 )   % Gauss-Chebyshev2: alpha = beta = .5
        x = chebpts(n+2, 2);
        x = x(2:n+1);
        w = pi/(n+1)*(1-x.^2).';
        v = (1-x.^2);
        v(2:2:end) = -v(2:2:end);
        [x, w] = rescale(x,w,interval,a,b);
        return
    else % Gauss-Gegenbauer: alpha = beta
        lambda = a + .5;
        [x, w, v] = ultrapts(n, lambda, interval);
        return
    end
end

% Choose an algorithm:
if strcmpi(method, 'rh')
    [x, w, v] = newton(n, a, b, 1); % RH (Newton on Riemann-Hilbert asy exp)
elseif strcmpi(method, 'exp')
    [x, w] = jacobiExp(n, a, b, param); % EXP (Explicit asy exp)
%     [x, w] = jacobiExp(n, a, b); % EXP (Explicit asy exp)
    v = sqrt(w'.*x);  % TODO FIX % Barycentric weights
%     v = (-1).^(0:n-1)'.*sqrt(w'.*x);  % Barycentric weights
%     v = v./max(abs(v));

elseif ( n < 20 || (n < 100 && ~method_set) || strcmpi(method, 'rec') )
%     [x, w, v] = rec(n, a, b); % REC (Recurrence relation)
    [x, w, v] = newton(n, a, b, 0); % REC (Recurrence relation)
    
elseif ( strcmpi(method, 'GW') )
    [x, w, v] = gw(n, a, b);  % GW  see [1]
    
else
    [x, w, v] = asy(n, a, b); % HT  see [2]
    
end

% Compute the constant for the weights:
if ( ~strcmpi(method,'GW') && ~strcmpi(method,'exp') )
%     if ( ~strcmpi(method,'GW') )
    if ( n >= 100 )
        cte1 = ratioOfGammaFunctions(n+a+1, b);
        cte2 = ratioOfGammaFunctions(n+1, b);
        C = 2^(a+b+1)*(cte2/cte1);
    else
        C = 2^(a+b+1) * exp( gammaln(n+a+1) - gammaln(n+a+b+1) + ...
            gammaln(n+b+1) - gammaln(n+1) );   
        % An alternative approach could be used: 
        %   C = 2^(a+b+1)*beta(a+1, b+1)/sum(w), 
        % but we prefer compute the weights independently.
    end
    w = C*w;
end

% Scale the nodes and quadrature weights:
[x, w] = rescale(x, w, interval, a, b);

% Scale the barycentric weights:
v = abs(v); 
v(2:2:end) = -v(2:2:end);
v = v./max(abs(v)); 

end

function [x, w] = rescale(x, w, interval, a, b)
%RESCALE   Rescale nodes and weights to an arbitrary finite interval.
    if ( ~all(interval == [-1, 1]) )
        c1 = .5*sum(interval); 
        c2 = .5*diff(interval);
        w = c2^(a+b+1)*w;
        x = c1 + c2*x;    
    end
end

function cte = ratioOfGammaFunctions(m,delta)
%RATIOGAMMA Compute the ratio gamma(m+delta)/gamma(m). See [2].
    % cte = gamma(m+delta)/gamma(m)
    ds = .5*delta^2/(m-1);
    s = ds;
    j = 1;
    while ( abs(ds/s) > eps/100 ) % Taylor series in expansion 
        j = j+1;
        ds = -delta*(j-1)/(j+1)/(m-1)*ds;
        s = s + ds;
    end
    p2 = exp(s)*sqrt(1+delta/(m-1))*(m-1)^(delta);
    % Stirling's series:
    g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
        5246819/75246796800, -534703531/902961561600, ...
        -4483131259/86684309913600, 432261921612371/514904800886784000];
    f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);
    cte = p2*(f(m+delta-1)/f(m-1));
end


%% ------------------------- Routines for GW ----------------------------
    
function [x, w, v] = gw(n, a, b)
    ab = a + b;
    ii = (2:n-1)';
    abi = 2*ii + ab;
    aa = [(b - a)/(2 + ab)
          (b^2 - a^2)./((abi - 2).*abi)
          (b^2 - a^2)./((2*n - 2+ab).*(2*n+ab))];
    bb = [2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ; 
          2*sqrt(ii.*(ii + a).*(ii + b).*(ii + ab)./(abi.^2 - 1))./abi];
    TT = diag(bb,1) + diag(aa) + diag(bb,-1); % Jacobi matrix.
    [V, x] = eig( TT );                       % Eigenvalue decomposition.
    x = diag(x);                              % Jacobi points.
    w = V(1,:).^2*2^(ab+1)*beta(a+1, b+1);    % Quadrature weights.
    v = sqrt(1-x.^2).*abs(V(1,:))';           % Barycentric weights.
end

%% ------------------------- Routines for REC ---------------------------


% function [x, w, v] = rec(n, a, b)
%    [x1, ders1] = rec_main(n, a, b, 1); % Nodes and P_n'(x)
%    [x2, ders2] = rec_main(n, b, a, 0); % Nodes and P_n'(x)
function [x, w, v] = newton(n, a, b, rh)
%REC   Compute nodes and weights using recurrrence relation.

   [x1, ders1] = rec_main(n, a, b, 1, rh); % Nodes and P_n'(x)
   [x2, ders2] = rec_main(n, b, a, 0, rh); % Nodes and P_n'(x)
   x = [-x2(end:-1:1) ; x1];
   ders = [ders2(end:-1:1) ; ders1];
   w = 1./((1-x.^2).*ders.^2)';        % Quadrature weights
   v = 1./ders;                        % Barycentric weights
   
end

function [x, PP] = rec_main(n, a, b, flag, rh)
%REC_MAIN   Jacobi polynomial recurrence relation.
if rh
    T = ceil(50/log(n) );
    Dinf = 2^(-a/2-b/2);
    Uright = zeros(2,2,T-1,ceil((T-1)/2));
    Uleft = zeros(2,2,T-1,ceil((T-1)/2));
    % About half of this tensor will not be used
    Uright(:,:,1,1) = (4*a^2-1)/16*[Dinf,0;0,Dinf^(-1)]*[-1,1i; 1i,1]*[Dinf^(-1),0;0,Dinf]; %A1
    Uleft(:,:,1,1) = (4*b^2-1)/16*[Dinf,0;0,Dinf^(-1)]*[1,1i; 1i,-1]*[Dinf^(-1),0;0,Dinf]; %B1
% %     if T <= 2,        break;    end
    Uright(:,:,2,1) = (4*a^2-1)/256*[Dinf,0;0,Dinf^(-1)]*[8*a+8*b-4*b^2+1, ...
        1i*(-8*a-8*b+4*a^2+4*b^2-10) ; 1i*(-8*a-8*b-4*a^2 ...
        -4*b^2+10), -8*a-8*b-4*b^2+1]*[Dinf^(-1),0;0,Dinf]; %A2
    Uleft(:,:,2,1) = (4*b^2-1)/256*[Dinf,0;0,Dinf^(-1)]*[-(8*b+8*a-4*a^2+1), ...
        1i*(-8*b-8*a+4*b^2+4*a^2-10) ; 1i*(-8*b-8*a-4*b^2 ...
        -4*a^2+10), -(-8*b-8*a-4*a^2+1)]*[Dinf^(-1),0;0,Dinf]; %B2
    
    a = a+1; b = b+1;
    Dinf = 2^(-a/2-b/2);
    Ursh = zeros(2,2,T-1,ceil((T-1)/2));
    Ulsh = zeros(2,2,T-1,ceil((T-1)/2));
    Ursh(:,:,1,1) = (4*a^2-1)/16*[Dinf,0;0,Dinf^(-1)]*[-1,1i; 1i,1]*[Dinf^(-1),0;0,Dinf]; %A1
    Ulsh(:,:,1,1) = (4*b^2-1)/16*[Dinf,0;0,Dinf^(-1)]*[1,1i; 1i,-1]*[Dinf^(-1),0;0,Dinf]; %B1
    Ursh(:,:,2,1) = (4*a^2-1)/256*[Dinf,0;0,Dinf^(-1)]*[8*a+8*b-4*b^2+1, ...
        1i*(-8*a-8*b+4*a^2+4*b^2-10) ; 1i*(-8*a-8*b-4*a^2 ...
        -4*b^2+10), -8*a-8*b-4*b^2+1]*[Dinf^(-1),0;0,Dinf]; %A2
    Ulsh(:,:,2,1) = (4*b^2-1)/256*[Dinf,0;0,Dinf^(-1)]*[-(8*b+8*a-4*a^2+1), ...
        1i*(-8*b-8*a+4*b^2+4*a^2-10) ; 1i*(-8*b-8*a-4*b^2 ...
        -4*a^2+10), -(-8*b-8*a-4*a^2+1)]*[Dinf^(-1),0;0,Dinf]; %B2
    a = a-1; b = b-1;
end
% Asymptotic formula (WKB) - only positive x.
if ( flag )
    r = ceil(n/2):-1:1;
else
    r = floor(n/2):-1:1;  
end
C = (2*r+a-.5)*pi/(2*n+a+b+1);
T = C + 1/(2*n+a+b+1)^2 * ((.25-a^2)*cot(.5*C) - (.25-b^2)*tan(.5*C));
x = cos(T).';

% Initialise:
dx = inf; 
l = 0;
% Loop until convergence:
while ( (norm(dx,inf) > sqrt(eps)/1000) && (l < 10) )
    l = l + 1;
    if rh
        P = polyAsyRH(n, x, a, b, Uright, Uleft);
%         error('der not impl');
%         warning('Check weights and possibly multiply with factor ratioOfGammaFunctions(m,delta)?');
        PP = n*polyAsyRH(n-1, x, a+1, b+1, Ursh, Ulsh); %Monic OP
%         PP = n*polyAsyRH(n-1, x+1, a+1, b+1, Ursh, Ulsh); %Wrong Monic OP but converged because stopped at first iter with der too high
        debug = 1;
    else
        [P, PP] = eval_Jac(x, n, a, b);
    end
    dx = -P./PP; 
    x = x + dx;
end
% Once more for derivatives:
[ignored, PP] = eval_Jac(x, n, a, b);

end

function [P, Pp] = eval_Jac(x, n, a, b)
%EVALJAC   Evaluate Jacobi polynomial and derivative via recurrence relation.

% Initialise:
ab = a + b;
P = .5*(a-b+(ab+2)*x);  
Pm1 = 1; 
Pp = .5*(ab+2);         
Ppm1 = 0; 

% n = 0 case:
if ( n == 0 )
    P = Pm1; 
    Pp = Ppm1; 
end

for k = 1:n-1
    % Useful values:
    A = 2*(k + 1)*(k + ab + 1)*(2*k + ab);
    B = (2*k + ab + 1)*(a^2 - b^2);
    C = prod(2*k + ab + (0:2)');
    D = 2*(k + a)*(k + b)*(2*k + ab + 2);

    % Recurrence:
    Pa1 = ( (B+C*x).*P - D*Pm1 ) / A;
    Ppa1 = ( (B+C*x).*Pp + C*P - D*Ppm1 ) / A;

    % Update:
    Pm1 = P; 
    P = Pa1;  
    Ppm1 =  Pp; 
    Pp = Ppa1;
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------- Routines for ASY algorithm ------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = asy(n, a, b)
%ASY   Compute nodes and weights using asymptotic formulae.

    if ( n <= 20 ) % Use only boundary formula:
        [xbdy, wbdy, vbdy] = asy2(n, a, b, ceil(n/2));  
        [xbdy2, wbdy2, vbdy2] = asy2(n, b, a, floor(n/2));  
        x = [-xbdy2(end:-1:1) ; xbdy];
        w = [wbdy2(end:-1:1), wbdy];
        v = [vbdy2(end:-1:1) ; vbdy];
        return
    end

    % Determine switch between interior and boundary regions:
    nbdy = 10;
    bdyidx1 = n-(nbdy-1):n; 
    bdyidx2 = nbdy:-1:1;

    % Interior formula:
    [x, w, v] = asy1(n,a,b,nbdy);   

    % Boundary formula (right):
    [xbdy, wbdy, vbdy] = asy2(n, a, b, nbdy);  
    x(bdyidx1) = xbdy;  
    w(bdyidx1) = wbdy; 
    v(bdyidx1) = vbdy;
    
    % Boundary formula (left):
    if ( a ~= b )
        [xbdy, wbdy, vbdy] = asy2(n, b, a, nbdy);  
    end
    x(bdyidx2) = -xbdy; 
    w(bdyidx2) = wbdy; 
    v(bdyidx2) = vbdy;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ASY1 (Interior)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = asy1(n, a, b, nbdy)
% Algorithm for computing nodes and weights in the interior.

    % Approximate roots via asymptotic formula: (Gatteschi and Pittaluga, 1985)
    K = (2*(n:-1:1)+a-.5)*pi/(2*n+a+b+1);
    tt = K + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*K)-(.25-b^2)*tan(.5*K));

    % First half (x > 0):
    t = tt(tt <= pi/2);
    mint = t(end-nbdy+1);
    idx = 1:max(find(t < mint,1)-1, 1);

    dt = inf; j = 0;
    % Newton iteration
    while ( norm(dt,inf) > sqrt(eps)/100 && j < 10 )
        [vals, ders] = feval_asy1(n, a, b, t, idx, 0);  % Evaluate
        dt = vals./ders;                                % Newton update
        t = t + dt;                                     % Next iterate
        j = j + 1;
        dt = dt(idx);
    end
    [vals, ders] = feval_asy1(n, a, b, t, idx, 1);      % Once more for luck
    t = t + vals./ders;                                 % Newton update.

    % Store:
    x = cos(t);
    w = 1./ders.^2;
    v = (sin(t)./ders);

    % Second half (x < 0):
    tmp = a; 
    a = b; 
    b = tmp;
    t = pi - tt(1:(n-length(x)));
    mint = t(nbdy);
    idx = max(find(t > mint, 1), 1):numel(t);

    dt = inf; j = 0;
    % Newton iteration
    while ( norm(dt,inf) > sqrt(eps)/100 && j < 10 )
        [vals, ders] = feval_asy1(n, a, b, t, idx, 0);  % Evaluate.
        dt = vals./ders;                                % Newton update.
        t = t + dt;                                     % Next iterate.
        j = j + 1;
        dt = dt(idx);
    end
    [vals, ders] = feval_asy1(n, a, b, t, idx, 1);      % Once more for luck.
    t = t + vals./ders;                                 % Newton update.

    % Store:
    x = [-cos(t) x].';
    w = [1./ders.^2 w];
    v = [sin(t)./ders v].';

end

% -------------------------------------------------------------------------

function [vals, ders] = feval_asy1(n, a, b, t, idx, flag)
% Evaluate the interior asymptotic formula at x = cos(t).
    
    % Number of terms in the expansion:
    M = 20;

    % Some often used vectors/matrices:
    onesT = ones(1,length(t));
    onesM = ones(M,1);
    MM = transpose(0:M-1);

    % The sine and cosine terms:
    alpha = (.5*(2*n+a+b+1+MM))*onesT .* (onesM*t) - .5*(a+.5)*pi;
    cosA = cos(alpha);
    sinA = sin(alpha);

    if ( flag ) % Evaluate cos(alpha) using Taylor series.
        k = 1:numel(t);
        if ( idx(1) == 1 )
            k = fliplr(k);
        end
        % Hi-lo computation to squeeze an extra digit in the computation.
        ta = double(single(t));    tb = t - ta;
        hi = n*ta;                 lo = n*tb + (a+b+1)*.5*t; 
        pia = 3.1415927410125732;  pib = -8.742278000372485e-08; % pib = pi-pia;
        % Doing this means that pi - pia - pib = 10^-24
        dh = ( hi - (k-.25)*pia ) + lo - .5*a*pia - (k-.25+.5*a)*pib;
        tmp = 0; sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;       % Initialise.
        for j = 0:20
            dc = sgn*DH/fact;
            tmp = tmp + dc;
            sgn = -sgn;
            fact = fact*(2*j+3)*(2*j+2);
            DH = DH.*dh2;
            if ( norm(dc, inf) < eps/2000 )
                break
            end
        end
        tmp(2:2:end) = -tmp(2:2:end);          % }
        [~, loc] = max(abs(tmp));              %  } Fix up the sign.
        tmp = sign(cosA(1,loc)*tmp(loc))*tmp;  % }
        cosA(1,:) = tmp;
    end

    sinT = onesM*sin(t);
    cosT = onesM*cos(t);
    cosA2 = cosA.*cosT + sinA.*sinT;
    sinA2 = sinA.*cosT - cosA.*sinT;

    one = ones(1,length(t));
    sinT = [one ; cumprod(onesM(2:end)*(.5*csc(.5*t)))];
    cosT = .5*sec(.5*t);

    j = 0:M-2;
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*n+a+b+j+2);
    P1 = [1  cumprod(vec)];
    P1(3:4:end) = -P1(3:4:end);
    P1(4:4:end) = -P1(4:4:end);
    P2 = eye(M);
    for l = 0:M-1
        j = 0:(M-l-2);
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*n+a+b+j+l+2);
        P2(l+1+(1:length(j)),l+1) = cumprod(vec);
    end
    PHI = repmat(P1,M,1).*P2;

    j = 0:M-2;
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*(n-1)+a+b+j+2);
    P1 = [1  cumprod(vec)];
    P1(3:4:end) = -P1(3:4:end);
    P1(4:4:end) = -P1(4:4:end);
    P2 = eye(M);
    for l = 0:M-1
        j = 0:(M-l-2);
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*(n-1)+a+b+j+l+2);
        P2(l+1+(1:length(j)),l+1) = cumprod(vec);
    end
    PHI2 = repmat(P1,M,1).*P2;

    S = 0; S2 = 0;
    SC = sinT;
    for m = 0:M-1

        l = 0:2:m;
        phi = PHI(m+1,l+1);
        dS1 = phi*SC(l+1,:).*cosA(m+1,:);

        phi2 = PHI2(m+1,l+1);
        dS12 = phi2*SC(l+1,:).*cosA2(m+1,:);

        l = 1:2:m;
        phi = PHI(m+1,l+1);
        dS2 = phi*SC(l+1,:).*sinA(m+1,:);

        phi2 = PHI2(m+1,l+1);
        dS22 = phi2*SC(l+1,:).*sinA2(m+1,:);

        if m > 10 && norm(dS1(idx) + dS2(idx),inf) < eps/100, break, end

        S = S + dS1 + dS2;
        S2 = S2 + dS12 + dS22;

        SC(1:m+1,:) = bsxfun(@times,SC(1:m+1,:),cosT);
    end

    % Constant out the front:
    dsa = .5*(a^2)/n; dsb = .5*(b^2)/n; dsab = .25*(a+b)^2/n;
    ds = dsa + dsb - dsab; s = ds; j = 1; 
    dsold = ds; % to fix a = -b bug.
    while ( (abs(ds/s) + dsold) > eps/10 )
        dsold = abs(ds/s);
        j = j+1;
        tmp = -(j-1)/(j+1)/n;
        dsa = tmp*dsa*a;
        dsb = tmp*dsb*b;
        dsab = .5*tmp*dsab*(a+b);
        ds = dsa + dsb - dsab;
        s = s + ds;
    end
    p2 = exp(s)*sqrt(2*pi)*sqrt((n+a)*(n+b)/(2*n+a+b))/(2*n+a+b+1);
    g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
         5246819/75246796800 -534703531/902961561600 ...
         -4483131259/86684309913600 432261921612371/514904800886784000];
    f = @(z) sum(g.*[1 cumprod(ones(1,9)./z)]);
    C = p2*(f(n+a)*f(n+b)/f(2*n+a+b))*2/pi;
    C2 = C*(a+b+2*n).*(a+b+1+2*n)./(4*(a+n).*(b+n));

    vals = C*S;
    S2 = C2*S2;

    % Use relation for derivative:
    ders = (n*(a-b-(2*n+a+b)*cos(t)).*vals + 2*(n+a)*(n+b)*S2)/(2*n+a+b)./sin(t);
    denom = 1./real(sin(t/2).^(a+.5).*cos(t/2).^(b+.5));
    vals = vals.*denom;
    ders = ders.*denom;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ASY2 (Boundary)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = asy2(n, a, b, npts)
% Algorithm for computing nodes and weights near the boundary.

% Use Newton iterations to find the first few Bessel roots:
smallK = min(30, npts);
jk = besselroots(a, min(npts, smallK));
% Use asy formula for larger ones (See NIST 10.21.19, Olver 1974 p247)
if ( npts > smallK )
    mu = 4*a^2;
    a8 = 8*((length(jk)+1:npts).'+.5*a-.25)*pi;
    jk2 = .125*a8-(mu-1)./a8 - 4*(mu-1)*(7*mu-31)/3./a8.^3 - ...
          32*(mu-1)*(83*mu.^2-982*mu+3779)/15./a8.^5 - ...
          64*(mu-1)*(6949*mu^3-153855*mu^2+1585743*mu-6277237)/105./a8.^7;
    jk = [jk ; jk2];
end
jk = real(jk(1:npts));

% Approximate roots via asymptotic formula: (see Olver 1974, NIST, 18.16.8)
rho = n + .5*(a + b + 1); 
phik = jk/rho;
t = phik + ((a^2-.25)*(1-phik.*cot(phik))./(2*phik) - ...
    .25*(a^2-b^2)*tan(.5*phik))/rho^2;

% Only first half (x > 0):
if ( any(t > 1.1*pi/2) )
    warning('CHEBFUN:jacpts:asy2:theta', ...
        'Theta > pi/2. Result may be inaccurate.'); 
end

% Compute higher terms:
[tB1, A2, tB2, A3] = asy2_higherterms(a, b, t, n);

dt = inf; j = 0;
% Newton iteration:
while ( norm(dt,inf) > sqrt(eps)/200 && j < 10)
    [vals, ders] = feval_asy2(n, t, 1); % Evaluate via asymptotic formula.
    dt = vals./ders;                    % Newton update.
    t = t + dt;                         % Next iterate.
    j = j + 1;
end
[vals, ders] = feval_asy2(n, t, 1);     % Evaluate via asymptotic formula.
dt = vals./ders;                        % Newton update
t = t + dt;    
    
% flip:
t = t(npts:-1:1); 
ders = ders(npts:-1:1);
% vals = vals(npts:-1:1);

% Revert to x-space:
x = cos(t);      
w = (1./ders.^2).';   
v = sin(t)./ders;

    function [vals, ders] = feval_asy2(n, t, flag)
    % Evaluate the boundary asymptotic formula at x = cos(t).
    
        % Useful constants:
        rho2 = n + .5*(a + b - 1);
        A = (.25 - a^2);       
        B = (.25 - b^2);
        
        % Evaluate the Bessel functions:
        Ja = besselj(a, rho*t, 0);
        Jb = besselj(a + 1, rho*t, 0);
        Jbb = besselj(a + 1, rho2*t, 0);
        if ( ~flag )
            Jab = besselj(a, rho2*t, 0);
        else
            % In the final step, perform accurate evaluation
            Jab = besselTaylor(-t, rho*t, a);
        end

        % Evaluate functions for recurrsive definition of coefficients:
        gt = A*(cot(t/2)-(2./t)) - B*tan(t/2);
        gtdx = A*(2./t.^2-.5*csc(t/2).^2) - .5*B*sec(t/2).^2;
        tB0 = .25*gt;
        A10 = a*(A+3*B)/24;
        A1 = gtdx/8 - (1+2*a)/8*gt./t - gt.^2/32 - A10;
        % Higher terms:
        tB1t = tB1(t); 
        A2t = A2(t); 

        % VALS:
        vals = Ja + Jb.*tB0/rho + Ja.*A1/rho^2 + Jb.*tB1t/rho^3 + Ja.*A2t/rho^4;
        % DERS:
        vals2 = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + Jbb.*tB1t/rho2^3 + Jab.*A2t/rho2^4;
        
        % Higher terms (not needed for n > 1000).
        tB2t = tB2(t); A3t = A3(t);
        vals = vals + Jb.*tB2t/rho^5 + Ja.*A3t/rho^6;
        vals2 = vals2 + Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6;
        
        % Constant out the front (Computed accurately!)
        ds = .5*(a^2)/n;
        s = ds; jj = 1;
        while abs(ds/s) > eps/10
            jj = jj+1;
            ds = -(jj-1)/(jj+1)/n*(ds*a);
            s = s + ds;
        end
        p2 = exp(s)*sqrt((n+a)/n)*(n/rho)^a;
        g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
             5246819/75246796800 -534703531/902961561600 ...
             -4483131259/86684309913600 432261921612371/514904800886784000];
        f = @(z) sum(g.*[1 cumprod(ones(1,9)./z)]);
        C = p2*(f(n+a)/f(n))/sqrt(2);

        % Scaling:
        valstmp = C*vals;
        denom = sin(t/2).^(a+.5).*cos(t/2).^(b+.5);
        vals = sqrt(t).*valstmp./denom;

        % Relation for derivative:
        C2 = C*n/(n+a)*(rho/rho2)^a;
        ders = (n*(a-b-(2*n+a+b)*cos(t)).*valstmp + 2*(n+a)*(n+b)*C2*vals2)/(2*n+a+b);
        ders = ders.*(sqrt(t)./(denom.*sin(t)));
        
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TODO]: The codes below are only here temporarily.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ja = besselTaylor(t, z, a)
%BESSELTAYLOR    Accurate evaluation of Bessel function J_A for asy2. (See [2].)
% BESSELTAYLOR(T, Z, A) evaluates J_A(Z+T) by a Taylor series expansion about Z. 

npts = numel(t);
kmax = min(ceil(abs(log(eps)/log(norm(t, inf)))), 30);
H = bsxfun(@power, t, 0:kmax).';
% Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
[nu, JK] = meshgrid(-kmax:kmax, z);
Bjk = besselj(a + nu, JK, 0);
nck = abs(pascal(floor(1.25*kmax), 1)); nck(1,:) = []; % nchoosek
AA = [Bjk(:,kmax+1), zeros(npts, kmax)];
fact = 1;
for k = 1:kmax
    sgn = 1;
    for l = 0:k
        AA(:,k+1) = AA(:,k+1) + sgn*nck(k,l+1)*Bjk(:,kmax+2*l-k+1);
        sgn = -sgn;
    end
    fact = k*fact;
    AA(:,k+1) = AA(:,k+1)/2^k/fact;
end
% Evaluate Taylor series:
Ja = zeros(npts, 1);
for k = 1:npts
    Ja(k,1) = AA(k,:)*H(:,k);
end
end

function [tB1, A2, tB2, A3, tB3, A4] = asy2_higherterms(alph, bet, theta, n)
%ASY2_HIGHERTERMS   Higher-order terms for boundary asymptotic series.
% Compute the higher order terms in asy2 boundary formula. See [2]. 

% These constants are more useful than alph and bet:
A = (0.25 - alph^2);
B = (0.25 - bet^2);

% For now, just work on half of the domain:
c = max(max(theta), .5);
if ( n < 30 )
    N = ceil(40 - n);
elseif ( n >= 30 && c > pi/2-.5)
    N = 15;
else
    N = 10;
end
Nm1 = N - 1;

% Scaled 2nd-kind Chebyshev points and barycentric weights:
t = .5*c*( sin(pi*(-Nm1:2:Nm1)/(2*Nm1)).' + 1 );
v = [.5 ; ones(Nm1,1)];
v(2:2:end) = -1;
v(end) = .5*v(end);

% The g's:
g = A*(cot(t/2)  -2./t) - B*tan(t/2);
gp = A*(2./t.^2 - .5*csc(t/2).^2) - .5*(.25-bet^2)*sec(t/2).^2;
gpp = A*(-4./t.^3 + .25*sin(t).*csc(t/2).^4) - 4*B*sin(t/2).^4.*csc(t).^3;
g(1) = 0; gp(1) = -A/6-.5*B; gpp(1) = 0;

% B0:
B0 = .25*g./t;
B0p = .25*(gp./t-g./t.^2);
B0(1) = .25*(-A/6-.5*B);
B0p(1) = 0;

% A1:
A10 = alph*(A+3*B)/24;
A1 = .125*gp - (1+2*alph)/2*B0 - g.^2/32 - A10;
A1p = .125*gpp - (1+2*alph)/2*B0p - gp.*g/16;
A1p_t = A1p./t;
A1p_t(1) = -A/720 - A^2/576 - A*B/96 - B^2/64 - B/48 + alph*(A/720 + B/48);

% Make f accurately: (Taylor series approx for small t)
fcos = B./(2*cos(t/2)).^2;
f = -A*(1/12 + t.^2/240+t.^4/6048 + t.^6/172800 + t.^8/5322240 + ...
    691*t.^10/118879488000 + t.^12/5748019200 + ...
    3617*t.^14/711374856192000 + 43867*t.^16/300534953951232000);
idx = t > .5;
ti = t(idx);
f(idx) = A.*(1./ti.^2 - 1./(2*sin(ti/2)).^2);
f = f - fcos;

% Integrals for B1: (Note that N isn't large, so we don't need to be fancy).
C = chebcolloc2.cumsummat(N)*(.5*c);
D = chebcolloc2.diffmat(N)*(2/c);
I = (C*A1p_t);
J = (C*(f.*A1));

% B1:
tB1 = -.5*A1p - (.5+alph)*I + .5*J;
tB1(1) = 0;
B1 = tB1./t;
B1(1) = A/720 + A^2/576 + A*B/96 + B^2/64 + B/48 + ...
    alph*(A^2/576 + B^2/64 + A*B/96) - alph^2*(A/720 + B/48);

% A2:
K = C*(f.*tB1);
A2 = .5*(D*tB1) - (.5+alph)*B1 - .5*K;
A2 = A2 - A2(1);

if ( nargout < 3 )
    % Make function for output
    tB1 = @(theta) bary(theta,tB1,t,v);
    A2 = @(theta) bary(theta,A2,t,v);
    return
end

% A2p:
A2p = D*A2;
A2p = A2p - A2p(1);
A2p_t = A2p./t;
% Extrapolate point at t = 0:
w = pi/2-t(2:end);
w(2:2:end) = -w(2:2:end);
w(end) = .5*w(end);
A2p_t(1) = sum(w.*A2p_t(2:end))/sum(w);

% B2:
tB2 = -.5*A2p - (.5+alph)*(C*A2p_t) + .5*C*(f.*A2);
B2 = tB2./t;
% Extrapolate point at t = 0:
B2(1) = sum(w.*B2(2:end))/sum(w);

% A3:
K = C*(f.*tB2);
A3 = .5*(D*tB2) - (.5+alph)*B2 - .5*K;
A3 = A3 - A3(1);

tB1 = @(theta) bary(theta, tB1, t, v);
A2 = @(theta) bary(theta, A2, t, v);
tB2 = @(theta) bary(theta, tB2, t, v);
A3 = @(theta) bary(theta, A3, t, v);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% Routines for RH algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the expansion of the orthonormal polynomial based on some heuristics
% np = n or n-1, ...
function p = polyAsyRH(np, x, alpha, beta, Uright, Uleft) %, Dinf, psi)
% We could avoid these tests by splitting the loop k=1:mn into three parts
% with heuristics for the bounding indices.
%T = ceil(50/log(np) );
if ( (alpha^2+ beta^2)/np > 1 )
    warning('CHEBFUN:jacpts:inputs',['A large alpha or beta may lead to inaccurate ' ...
        'results because the weight is low and R(z) is not close to identity.']);
end
p = nan*x;
T = size(Uright,3);
Dinf = 2^(-alpha/2-beta/2);
for k = 1:length(x)
    psiz = (alpha+beta)/2*acos(x(k)) -alpha*pi/2; %psi(x(k));
    R = zeros(2,2, T-1);
    for kt = 1:T-1
        for m = 1:round(kt/2)
            R(:,:,kt) = R(:,:,kt) + Uright(:,:,kt,m)./(x(k)-1).^m + Uleft(:,:,kt,m)./(x(k)+1).^m;
        end
    end
    if x(k) < -0.8 % Does not occur due to splitting
        % The fixed delta in the RHP would mean k ~ sqrt(n) in x_k = -1+jbk^2/n^2
        p(k) = asyLeft(np, x(k), beta, T, psiz, Dinf, R);
    elseif x(k) > 0.8
        % The fixed delta in the RHP would mean k ~ sqrt(n) in x_k = 1 -jak^2/n^2
        p(k) = asyRight(np, x(k), alpha, T, psiz, Dinf, R);
    else
        RI = eye(2);
        for kt = 1:T-1
            RI = RI + R(:,:,kt)/np^kt;
        end
%         p(k) = [2.^(1/2-n)./(sqrt(w(z) ).*(1+z).^(1/4).*(1-z).^(1/4) ), 0]*RI*...
        p(k) = real([sqrt(2), 0]*RI*[Dinf*cos((np +1/2)*acos(x(k)) +psiz -pi/4); ...
            -1i/Dinf*cos((np -1/2)*acos(x(k)) +psiz -pi/4)]);
%         p(k) = asyBulk(np, x(k), alpha, beta, T);
    end
%     p(k) = p(k)*(1 - x(k))^(1/4 + alpha/2)*(1 + x(k))^(1/4 + beta/2)/2^np;
    p(k) = p(k)*(1 - x(k))^(-1/4 -alpha/2)*(1 + x(k))^(-1/4 -beta/2)/2^np;
end

end

% NOT NEEDED DUE TO SPLITTING: in newton(...)
% % Compute the expansion of the monic polynomial in the left disk without the common factor before
% function p = asyLeft(n, z, beta, T, psiz, Dinf, R)
% 
% if abs(z+1) < eps^(1/3)
%     warning('z is close to -1 so use Q-s and series expansions in "method"')
% end
% brac = @(k) (k == 0) + (k~= 0)*prod(4*beta^2 -(2*(1:k)-1).^2)/(2^(2*k)*factorial(k) ); % = (beta,k)
% % k = reshape(1:T-1, [1,1,T-1]);
% % w = nan;
% % s = arrayfun(@(y) brac(y-1), k)./2^k./(log(-z -sqrt(z-1).*sqrt(z+1)).^k)./(2*sqrt(z+1)*sqrt(z-1))*[Dinf 0; 0 Dinf^(-1)]*...
% DeltaL = @(k) brac(k-1)./2^k./(log(-z -sqrt(z-1).*sqrt(z+1)).^k)./(2*sqrt(z+1)*sqrt(z-1))*[Dinf 0; 0 Dinf^(-1)]*...
%     [sqrt(z + sqrt(z-1).*sqrt(z+1) ), 1i/sqrt(z + sqrt(z-1).*sqrt(z+1) ); ...
%     -1i/sqrt(z + sqrt(z-1).*sqrt(z+1) ), sqrt(z + sqrt(z-1).*sqrt(z+1) ) ]*...
%     [exp( 1i*(psiz -beta*pi/2) ) 0; 0 exp(-1i*(psiz -beta*pi/2))]*...
%     [ ((-1).^k)./k.*(beta^2+k/2-1/4), 1i*(k-1/2) ; (-1).^(k+1).*1i.*(k-1/2), (beta^2+k/2-1/4)./k]*...
%     [exp(-1i*(psiz -beta*pi/2)) 0; 0 exp(1i*(psiz -beta*pi/2) )]*...
%     [sqrt(z + sqrt(z-1).*sqrt(z+1) ), -1i/sqrt(z + sqrt(z-1).*sqrt(z+1) ); ...
%     1i/sqrt(z + sqrt(z-1).*sqrt(z+1) ), sqrt(z + sqrt(z-1).*sqrt(z+1) ) ]*[Dinf^(-1) 0; 0 Dinf];
% % Actually [exp( 1i*(-1)^(angle(z-1) <= 0)*(psiz -beta*pi/2) ) 0; 0 exp(-1i*(-1)^(angle(z-1) <= 0)*(psiz -beta*pi/2))]*...
% RL = eye(2);
% s = zeros(2,2,T-1);
% for k = 1:T-1
%     s(:,:,k) = DeltaL(k);
%     if mod(k,2)==0 % k is even: add extra matrix
%         s(:,:,k) = s(:,:,k) -brac(k-1)*(4*beta^2+2*k-1)/2^(k+1)/k...
%             /(log(-z-sqrt(z-1).*sqrt(z+1) ) )^k*eye(2);
%     end
%     RL = RL + (R(:,:,k)-s(:,:,k) )/n^k;
%     for j = 1:(k-1)
%         RL = RL - R(:,:,k-j)*s(:,:,j)/n^k;
%     end
% end
% % p = real( [sqrt(pi*n*acos(-z))./((-2)^n.*sqrt(w(z)).*(1+z).^(1/4).*(1-z).^(1/4) ) ,  0]*RL*...
% p = real( [sqrt(pi*n*acos(-z)).*(-1)^n ,  0]*RL*...
%         [Dinf*sin(psiz -beta*pi/2 + acos(z)/2)*besselj(beta,n*acos(-z) ) + cos(psiz -beta*pi/2 + acos(z)/2)*...
% 		(besselj(beta-1,n*acos(-z) ) - besselj(beta+1,n*acos(-z) ) )/2; ...
%         -1i/Dinf*sin(psiz -beta*pi/2 -acos(z)/2)*besselj(beta,n*acos(-z) ) + cos(psiz -beta*pi/2 -acos(z)/2)*...
%         (besselj(beta-1,n*acos(-z) ) - besselj(beta+1,n*acos(-z) ) )/2] );
% end

% Compute the expansion of the monic polynomial in the right disk without common
% TODO: merge with asyLeft??
function p = asyRight(n, z, alpha, T, psiz, Dinf, R)
if abs(z-1) < eps^(1/3)
    warning('z is close to -1 so use Q-s and series expansions in "method"')
end
brac = @(k) (k == 0) + (k~= 0)*prod(4*alpha^2 -(2*(1:k)-1).^2)/(2^(2*k)*factorial(k) ); % = (alpha,k)
% k = reshape(1:T-1, [1,1,T-1]);
% psi = nan;
% s = arrayfun(@(y) brac(y-1), k)./2^k./(log(z + sqrt(z-1).*sqrt(z+1)).^k)./(2*sqrt(z+1)*sqrt(z-1))*[Dinf 0; 0 Dinf^(-1)]*...
DeltaR = @(k) brac(k-1)./2^k./(log(z + sqrt(z-1).*sqrt(z+1)).^k)./(2*sqrt(z+1)*sqrt(z-1))*[Dinf 0; 0 Dinf^(-1)]*...
    [sqrt(z + sqrt(z-1).*sqrt(z+1) ), 1i/sqrt(z + sqrt(z-1).*sqrt(z+1) ); ...
    -1i/sqrt(z + sqrt(z-1).*sqrt(z+1) ), sqrt(z + sqrt(z-1).*sqrt(z+1) ) ]*...
    [exp( 1i*(psiz +alpha*pi/2) ) 0; 0 exp(-1i*(psiz +alpha*pi/2))]*...
    [ ((-1).^k)./k.*(alpha^2+k/2-1/4), -1i*(k-1/2) ; (-1).^k.*1i.*(k-1/2), (alpha^2+k/2-1/4)./k]*...
    [exp(-1i*(psiz +alpha*pi/2)) 0; 0 exp(1i*(psiz +alpha*pi/2) )]*...
    [sqrt(z + sqrt(z-1).*sqrt(z+1) ), -1i/sqrt(z + sqrt(z-1).*sqrt(z+1) ); ...
    1i/sqrt(z + sqrt(z-1).*sqrt(z+1) ), sqrt(z + sqrt(z-1).*sqrt(z+1) ) ]*[Dinf^(-1) 0; 0 Dinf];
% Actually [exp( 1i*(-1)^(angle(z-1) <= 0)*(psiz ...
RR = eye(2);
s = zeros(2,2,T-1);
for k = 1:T-1
    s(:,:,k) = DeltaR(k);
    if mod(k,2)==0 % k is even: add extra matrix
        s(:,:,k) = s(:,:,k) -brac(k-1)*(4*alpha^2+2*k-1)/2^(k+1)/k...
            /(log(z+sqrt(z-1).*sqrt(z+1) ) )^k*eye(2);
    end
    RR = RR + (R(:,:,k)-s(:,:,k) )/n^k; % Avoids defining R_0^{right} = I
    for j = 1:(k-1)
        RR = RR - R(:,:,k-j)*s(:,:,j)/n^k;
    end
end
% p = real([sqrt(pi*n*acos(z))./(2^n.*sqrt(w(z) ).*(1+z).^(1/4).*(1-z).^(1/4) ),  0]*RR*...
p = real([sqrt(pi*n*acos(z)),  0]*RR*...
    [Dinf*cos(psiz +alpha*pi/2 + acos(z)/2)*besselj(alpha,n*acos(z) ) + sin(psiz +alpha*pi/2 + acos(z)/2)*...
    (besselj(alpha-1,n*acos(z) ) - besselj(alpha+1,n*acos(z) ) )/2;  ...
    -1i/Dinf*cos(psiz +alpha*pi/2 -acos(z)/2)*besselj(alpha,n*acos(z) ) + sin(psiz +alpha*pi/2 -acos(z)/2)*...
    (besselj(alpha-1,n*acos(z) ) - besselj(alpha+1,n*acos(z) ) )/2] );
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%% Routine for explicit expansion %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = jacobiExp(n, alpha, beta, param) %T, il, ir)

% w = zeros(1, n);

if 1
    il = param(2); ir = param(3);
% %     il = 0; ir = n+1;
else
% %     il = max(round(sqrt(n)), 7);
    % Heuristics to switch between Bessel regions and bulk
% %     ir = n - max(round(sqrt(n)), 7);
end

% This is a heuristic for the number of terms in the expansions that follow.
% % T = ceil(50/log(n) );
T = param(1);
if ( alpha^2/n > 1 )
    warning('CHEBFUN:lagpts:inputs',['A large alpha may lead to inaccurate ' ...
        'results because the weight is low and R(z) is not close to identity.']);
end
d = 1/(2*n+alpha+beta+1);

% t = cos((4*(il+1:ir-1) +2*alpha +3)/(4*n +2*alpha +2*beta +2) )';
t = cos((4*n -4*(il+1:ir-1) +2*alpha +3)/(4*n +2*alpha +2*beta +2)*pi )';

% jak = besselroots(alpha, il); % To separate fct to allow using symmetry
% jbk = besselroots(beta, n-ir+1);

bulk = t*0;
wbulk = t*0;

if T >= 7
    bulk = bulk+ 1/240*(576*alpha^6 - 576*beta^6 + (96*alpha^6 + 96*beta^6 + 80*(8*alpha^2 - 3)*beta^4 ...
        - 240*alpha^4 + 2*(320*alpha^4 - 440*alpha^2 + 101)*beta^2 + 202*alpha^2 - 39)*t.^5 ...
        -320*(alpha^2 - 6)*beta^4 + 240*(alpha^2 - beta^2)*t.^4 - 1920*alpha^4 ...
        -10*(32*(5*alpha^2 + 3)*beta^4 + 96*alpha^4 + 2*(80*alpha^4 - 152*alpha^2 - 97)*beta^2 ...
        - 194*alpha^2 + 99)*t.^3 + 16*(20*alpha^4 - 127)*beta^2 + 160*(6*alpha^6 - 6*beta^6 +...
        2*(alpha^2 + 15)*beta^4 - 30*alpha^4 - (2*alpha^4 + 41)*beta^2 + 41*alpha^2)*t.^2 + ...
        2032*alpha^2 + 15*(96*alpha^6 + 96*beta^6 + 16*(4*alpha^2 - 23)*beta^4 - 368*alpha^4 ...
        + 2*(32*alpha^4 - 72*alpha^2 + 223)*beta^2 + 446*alpha^2 - 173)*t)*d^6./(t.^4 - 2*t.^2 +1);
end
if T >= 5
    bulk = bulk - 1/24*(32*alpha^4 - 32*beta^4 - (16*alpha^4 + 16*beta^4 + 4*(12*alpha^2 -5)*beta^2 - 20*alpha^2 + 5)*t.^3 ...
        - 24*(alpha^2 - beta^2)*t.^2 - 40*alpha^2 + 40*beta^2+ 3*(16*alpha^4 + 16*beta^4 + 4*(4*alpha^2 - 7)*beta^2 ...
        - 28*alpha^2 + 11)*t)*d^4./(t.^2- 1) ;
end
if T >= 3
    bulk = bulk + 1/2*(2*alpha^2 - 2*beta^2 + (2*alpha^2 + 2*beta^2 - 1)*t)*d^2;
end

[xL, wL] = jacobiDisk(n,alpha,beta,T,il);
[xR, wR] = jacobiDisk(n,beta,alpha,T,n-ir+1);

x = [xL; bulk+t; flipud(-xR)];
w = transpose((1-x).^alpha.*(1+x).^beta.*[wL; 4*pi*d^2*sqrt(1 -t.^2).*(1+wbulk); flipud(wR)]);

% serWeiLensGen = 2*(pi + 4*pi*alpha^2 + 4*pi*beta^2 - (pi + 4*pi*alpha^2 + 4*pi*beta^2 + 4*pi*alpha +4*(pi + pi*alpha)*beta)*t^2 +...
%     4*pi*alpha + 4*(pi + pi*alpha)*beta)*d^4/(sqrt(t +1)*sqrt(-t + 1)) + 4*(pi + pi*alpha + pi*beta ...
%     - (pi + pi*alpha +pi*beta)*t)*d^3*sqrt(t + 1)/sqrt(-t + 1) + 4*pi*d^2*sqrt(t + 1)*sqrt(-t + 1);

% serNodLensGen = 1/240*(576*alpha^6 - 576*beta^6 + (96*alpha^6 + 96*beta^6 + 80*(8*alpha^2 - 3)*beta^4
% - 240*alpha^4 + 2*(320*alpha^4 - 440*alpha^2 + 101)*beta^2 + 202*alpha^2 - 39)*t^5 -
% 320*(alpha^2 - 6)*beta^4 + 240*(alpha^2 - beta^2)*t^4 - 1920*alpha^4 -
% 10*(32*(5*alpha^2 + 3)*beta^4 + 96*alpha^4 + 2*(80*alpha^4 - 152*alpha^2 - 97)*beta^2
% - 194*alpha^2 + 99)*t^3 + 16*(20*alpha^4 - 127)*beta^2 + 160*(6*alpha^6 - 6*beta^6 +
% 2*(alpha^2 + 15)*beta^4 - 30*alpha^4 - (2*alpha^4 + 41)*beta^2 + 41*alpha^2)*t^2 +
% 2032*alpha^2 + 15*(96*alpha^6 + 96*beta^6 + 16*(4*alpha^2 - 23)*beta^4 - 368*alpha^4
% + 2*(32*alpha^4 - 72*alpha^2 + 223)*beta^2 + 446*alpha^2 - 173)*t)*d^6/(t^4 - 2*t^2 +1) 
% - 1/24*(32*alpha^4 - 32*beta^4 - (16*alpha^4 + 16*beta^4 + 4*(12*alpha^2 -5)*beta^2 - 20*alpha^2 + 5)*t^3 - 24*(alpha^2 - beta^2)*t^2 ...
%     - 40*alpha^2 + 40*beta^2+ 3*(16*alpha^4 + 16*beta^4 + 4*(4*alpha^2 - 7)*beta^2 - 28*alpha^2 + 11)*t)*d^4/(t^2- 1) 
% + 1/2*(2*alpha^2 - 2*beta^2 + (2*alpha^2 + 2*beta^2 - 1)*t)*d^2 + t

end

% Evaluate the asymptotic expansions of ix nodes and weights in one of the disks: default left disk 
function [x,w] = jacobiDisk(n, alpha, beta, T, ix)
jbk = besselroots(beta, ix);
x = 0*jbk;
w = 0*jbk;

d = 1/(2*n+alpha+beta+1);
% z = nan
if T >= 7
    x = x -1/2835*(9*jbk.^6 - 18*(7*alpha^2 + 5*beta^2 - 3)*jbk.^4 + (328*beta^4 + (1512*alpha^2- 575)*beta^2 + 567*alpha^2 - 113)*jbk.^2 ...
        - (2835*alpha^6 + 247*beta^6 +1407*(3*alpha^2 - 1)*beta^4 - 8505*alpha^4 + 21*(405*alpha^4 - 600*alpha^2 +...
        133)*beta^2 + 8379*alpha^2 - 1633))*d^6;
    w = w + 1/2835*(2835*alpha^6 + 247*beta^6 - 36*jbk.^6 + 1407*(3*alpha^2 - 1)*beta^4 + ...
        54*(7*alpha^2 + 5*beta^2 - 3)*jbk.^4 - 8505*alpha^4 + 21*(405*alpha^4 - 600*alpha^2 +133)*beta^2 ...
        - 2*(328*beta^4 + (1512*alpha^2 - 575)*beta^2 + 567*alpha^2 - 113)*jbk.^2+ 8379*alpha^2 - 1633)*d^6;
end
if T >= 5
    x = x + 1/45*(2*jbk.^4 - 3*(5*alpha^2 + 3*beta^2 - 2)*jbk.^2 + 45*alpha^4 + 7*beta^4 + 20*(3*alpha^2 - 1)*beta^2 - 60*alpha^2 + 13)*d^4;
    w = w + 1/45*(45*alpha^4 + 7*beta^4 + 6*jbk.^4 + 20*(3*alpha^2 - 1)*beta^2 - 6*(5*alpha^2 + 3*beta^2 - 2)*jbk.^2 - 60*alpha^2 + 13)*d^4;
end
if T >= 3
    x = x - 1/3*(jbk.^2 - (3*alpha^2 + beta^2 - 1))*d^2;
    w = w + 1/3*(3*alpha^2 + beta^2 - 2*jbk.^2 - 1)*d^2;
end

x = -1 + 2*d^2*jbk.^2.*(1 +x);
w = 8*d^2*besselj(beta-1,jbk).^(-2).*(1 +w);

% serNodBounStd =  -2/2835*(9*jbk^8 - 18*(7*alpha^2 + 5*beta^2 - 3)*jbk^6 + (328*beta^4 + (1512*alpha^2 ...
% - 575)*beta^2 + 567*alpha^2 - 113)*jbk^4 - (2835*alpha^6 + 247*beta^6 + ...
% 1407*(3*alpha^2 - 1)*beta^4 - 8505*alpha^4 + 21*(405*alpha^4 - 600*alpha^2 + ...
% 133)*beta^2 + 8379*alpha^2 - 1633)*jbk^2)*z^8 + 2/45*(2*jbk^6 - 3*(5*alpha^2 + ...
% 3*beta^2 - 2)*jbk^4 + (45*alpha^4 + 7*beta^4 + 20*(3*alpha^2 - 1)*beta^2 - 60*alpha^2 ...
% + 13)*jbk^2)*z^6 - 2/3*(jbk^4 - (3*alpha^2 + beta^2 - 1)*jbk^2)*z^4 + 2*jbk^2*z^2 - 1;
% 
% serWeiBounStd = 8/2835*(2835*alpha^6 + 247*beta^6 - 36*jbk^6 + 1407*(3*alpha^2 - 1)*beta^4 + ...
% 54*(7*alpha^2 + 5*beta^2 - 3)*jbk^4 - 8505*alpha^4 + 21*(405*alpha^4 - 600*alpha^2 + ...
% 133)*beta^2 - 2*(328*beta^4 + (1512*alpha^2 - 575)*beta^2 + 567*alpha^2 - 113)*jbk^2 ...
% + 8379*alpha^2 - 1633)*z^8/jbm^2 + 8/45*(45*alpha^4 + 7*beta^4 + 6*jbk^4 + ...
% 20*(3*alpha^2 - 1)*beta^2 - 6*(5*alpha^2 + 3*beta^2 - 2)*jbk^2 - 60*alpha^2 + ...
% 13)*z^6/jbm^2 + 8/3*(3*alpha^2 + beta^2 - 2*jbk^2 - 1)*z^4/jbm^2 + 8*z^2/jbm^2;

end
