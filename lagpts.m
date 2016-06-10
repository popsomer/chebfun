function [x, w, v] = lagpts(n, int, meth)
%LAGPTS  Laguerre points and Gauss-Laguerre Quadrature Weights.
%   LAGPTS(N) returns N Laguerre points X in (0,inf).
%
%   [X, W] = LAGPTS(N) returns also a row vector W of weights for Gauss-Laguerre
%   quadrature. [X, W, V] = LAGPTS(N) returns in addition a column vector V
%   of the barycentric weights corresponding to X.
%
%   LAGPTS(N, D) scales the nodes and weights for the semi-infinite domain D.
%   D can be either a domain object or a vector with two components.
%
%   [X, W] = LAGPTS(N, METHOD) allows the user to select which method to use.
%   METHOD = 'GW' will use the traditional Golub-Welsch eigenvalue method,
%   which is best for when N is small. METHOD = 'FAST' will use the
%   Glaser-Liu-Rokhlin fast algorithm, which is much faster for large N.
%   By default LAGPTS uses 'GW' when N < 128. METHOD = 'RH' will use asymptotics
%   of Laguerre polynomials.
%
% References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23:221-230, 1969,
%   [2] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the calculation
%       of the roots of special functions", SIAM Journal on Scientific 
%       Computing", 29(4):1420-1438, 2007.
%   [3] P. Opsomer, (in preparation).
%   [4] M. Vanlessen, "Strong asymptotics of Laguerre-Type orthogonal
%       polynomials and applications in Random Matrix Theory", Constr. Approx.,
%       25:125-175, 2007.
%
% See also CHEBPTS, LEGPTS, HERMPTS, and JACPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%
% 'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'FAST' by Nick Hale, March 2010 - algorithm adapted from [2].
% 'RH' by Peter Opsomer, June 2016 - algorithm adapted from [3], based on [4].

% Defaults:
method = 'default';
interval = [0, inf];

if ( n < 0 )
    error('CHEBFUN:lagpts:n', ...
        'First input should be a positive integer.');
end

% Return empty vector if n = 0.
if ( n == 0 )
    [x, w, v] = deal([]);
    return
end

% Check the inputs
if ( nargin > 1 )
    if ( nargin == 3 )
        interval = int;
        method = meth;
    elseif ( nargin == 2 )
        if ( ischar(int) )
            method = int;
        else
            interval = int;
        end
    end
    if ( ~any(strcmpi(method, {'default', 'GW', 'fast', 'RH'})) )
        error('CHEBFUN:lagpts:inputs', 'Unrecognised input string %s.', method);
    end
    if ( numel(interval) > 2 )
        warning('CHEBFUN:lagpts:domain',...
            'Piecewise intervals not supported and will be ignored.');
        interval = interval([1, end]);
    end
end

if ( sum(isinf(interval)) ~= 1 )
    error('CHEBFUN:lagpts:inf', 'LAGPTS only supports semi-infinite domains.');
end

% decide to use GW or FAST
if ( strcmpi(method,'GW') || ( ( n < 128 ) && strcmpi(method,'default') ) )
    % GW, see [1]
    
    alpha = 2*(1:n)-1;  beta = 1:n-1;     % 3-term recurrence coeffs
    T = diag(beta,1) + diag(alpha) + diag(beta,-1);  % Jacobi matrix
    [V, D] = eig(T);                      % eigenvalue decomposition
    [x, indx] = sort(diag(D));            % Laguerre points
    w = V(1,indx).^2;                     % Quadrature weights
    v = sqrt(x).*abs(V(1,indx)).';        % Barycentric weights
    v = v./max(v); 
    v(2:2:n) = -v(2:2:n);
    
elseif ( strcmpi(method,'fast') || strcmpi(method,'default') )
    % Fast, see [2]
    [x, ders] = alg0_Lag(n);              % Nodes and L_n'(x)
    w = exp(-x)./(x.*ders.^2); w = w';    % Quadrature weights
    v = exp(-x/2)./ders;                  % Barycentric weights
    v = -v./max(abs(v));
    
else
    % RH, see [3] and [4]
    [x, w] = alg_rh(n);                   % Nodes and quadrature weights
    v = sqrt(w.*x);                       % Barycentric weights
    v = -v./max(abs(v));
    
end
w = (1/sum(w))*w;                         % Normalise so that sum(w) = 1

% Nonstandard interval
if ( ~all(interval == [0, inf]) )
    a = interval(1); 
    b = interval(2);
    if ( isinf(b) )
        x = x + a;
        w = w*exp(-a);
    else
        x = -x + b;
        w = w*exp(b);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% Routines for FAST algorithm %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, ders] = alg0_Lag(n)
ders = zeros(n, 1);
xs = 1/(2*n+1);
n1 = 20;
n1 = min(n1, n);
x = zeros(n, 1);
for k = 1:n1
    [xs, ders(k)] = alg3_Lag(n, xs);
    x(k) = xs;
    xs = 1.1*xs;
end
[x, ders] = alg1_Lag(x, ders, n, n1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roots, ders] = alg1_Lag(roots, ders, n, n1)

% number of terms in Taylor expansion
m = 30;

% initialise
hh1 = ones(m+1, 1); 
zz = zeros(m, 1); 
u = zeros(1, m+1); 
up = zeros(1, m+1);

x = roots(n1);
for j = n1:(n - 1)
    
    % initial approx
    h = rk2_Lag(pi/2, -pi/2, x, n) - x;
    
    % scaling:
    M = 1/h; 
    M2 = M^2; 
    M3 = M^3; 
    M4 = M^4;
    
    % recurrence relation for Laguerre polynomials
    r = x*(n + .5 - .25*x);  
    p = x^2;
    u(1:2) = [0 ; ders(j)/M];
    u(3) = -.5*u(2)/(M*x) - (n + .5 - .25*x)*u(1)/(x*M^2);
    u(4) = -u(3)/(M*x) + ( -(1+r)*u(2)/6/M^2 - (n+.5-.5*x)*u(1)/M^3 ) / p;
    up(1:3) = [u(2) ; 2*u(3)*M ; 3*u(4)*M];
    
    for k = 2:(m - 2)
        u(k+3) = ( -x*(2*k+1)*(k+1)*u(k+2)/M - (k*k+r)*u(k+1)/M2 - ...
                   (n+.5-.5*x)*u(k)/M3 + .25*u(k-1)/M4 ) / (p*(k+2)*(k+1));
        up(k+2) = (k+2)*u(k+3)*M;
    end
    up(m+1) = 0;
    
    % Flip for more accuracy in inner product calculation.
    u = u(m+1:-1:1);  
    up = up(m+1:-1:1);
    
    % Newton iteration
    hh = hh1; 
    hh(end) = M;    
    step = inf;  
    l = 0;
    if ( M == 1 )
        Mhzz = (M*h) + zz;
        hh = [M ; cumprod(Mhzz)];
        hh = hh(end:-1:1);
    end
    while ( (abs(step) > eps) && (l < 10) )
        l = l + 1;
        step = (u*hh)/(up*hh);
        h = h - step;
        Mhzz = (M*h) + zz;
        % Powers of h (This is the fastest way!)
        hh = [M ; cumprod(Mhzz)];     
        % Flip for more accuracy in inner product
        hh = hh(end:-1:1);          
    end
    
    % Update
    x = x + h;
    roots(j+1) = x;
    ders(j+1) = up*hh;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1, d1] = alg3_Lag(n, xs)
[u, up] = eval_Lag(xs, n);
theta = atan(sqrt(xs/(n + .5 - .25*xs))*up/u);
x1 = rk2_Lag(theta, -pi/2, xs, n);

% Newton iteration
step = inf;  
l = 0;
while ( (abs(step) > eps || abs(u) > eps) && (l < 200) )
    l = l + 1;
    [u, up] = eval_Lag(x1, n);
    step = u/up;
    x1 = x1 - step;
end

[ignored, d1] = eval_Lag(x1, n);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate Laguerre polynomial via recurrence
function [L, Lp] = eval_Lag(x, n)
Lm2 = 0; 
Lm1 = exp(-x/2); 
Lpm2 = 0; 
Lpm1 = 0;
for k = 0:n-1
    L = ( (2*k+1-x).*Lm1 - k*Lm2 ) / (k + 1);
    Lp = ( (2*k+1-x).*Lpm1 - Lm1 - k*Lpm2 ) / (k + 1);
    Lm2 = Lm1; 
    Lm1 = L;
    Lpm2 = Lpm1; 
    Lpm1 = Lp;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runge-Kutta for Laguerre Equation
function x = rk2_Lag(t, tn, x, n)
m = 10; 
h = (tn - t)/m;
for j = 1:m
    f1 = (n + .5 - .25*x);
    k1 = -h/( sqrt(f1/x) + .25*(1/x-.25/f1)*sin(2*t) );
    t = t + h;  
    x = x + k1;   
    f1 = (n + .5 - .25*x);
    k2 = -h/( sqrt(f1/x) + .25*(1/x-.25/f1)*sin(2*t) );
    x = x + .5*(k2 - k1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% Routines for RH algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w] = alg_rh(n)

x = [ besselroots(0, 3).^2/(4*n + 2); zeros(n-3, 1) ];
w = zeros(n, 1);

% This is a heuristic for the number of terms in the expansions.
T = ceil(25/log(n) );
UQ0 = getUQalpha0(T-1);
UQ1 = getUQalpha1(T-1);

% factor = -(1-1/(n+1))^(n+1)*(1+1/n)^(-1/2)/(4*n)/exp(-1-2*log(2))...
%     /sqrt(1-4i*sum((Ul0(1,2,1:(T-1),1) + Ur0(1,2,1:(T-1),1))./(n+1) ...
%     .^reshape(1:(T-1),[1,1,T-1]) ) )*sqrt(1-4i*sum((Ul0(1,2,1:(T-1),1) + ...
%     Ur0(1,2,1:(T-1),1))./n.^reshape(1:(T-1),[1,1,T-1]) ) );

factor = -(1-1/(n+1))^(n+1)*(1+1/n)^(-1/2)/(4*n)/exp(-1-2*log(2))...
    /sqrt(1 - 4i*sum((UQ0(2,2,1:T-1,1) + UQ0(1,2,1:T-1,1))./(n + 1) ...
    .^reshape(1:T-1,[1,1,T-1]) ) )*sqrt(1 - 4i*sum((UQ0(2,2,1:T-1,1) + ...
    UQ0(1,2,1:T-1,1))./n.^reshape(1:T-1,[1,1,T-1]) ) );
% dPoly = @(y) poly(n, x, 0 , UQ0) - poly(n, x, 1, UQ1)*sqrt(n + 1);
% poly = @(np, y, alpha, Ur, Ul) pl(np, y, alpha, Ur, Ul);
poly = @pl;
% dPoly = @(y) poly(n, y, 0 , UQ0) - poly(n, y, 1, UQ1)*sqrt(n + 1);
% [FIXME] Ensure analytic continuation of all functions to avoid adding eps*1i.
dPoly = @(y) pl(n-1,(1+eps*1i)*y, 1, UQ1)*sqrt(n);
% The expansion near zero is employed for x/(4n) < 0.2 as a heuristic, otherwise
% the expansion near 4n. This leads to the given bound for the index by using an
% % approximation of the nodes and of the Bessel zeros.
% The expansion in the lens is cheaper, but was less accurate.
% % for k = 1:n/2
ls = zeros(n,1);
for k = 1:n
    if ( k > 3 )
        x(k) = 3*x(k-1) - 3*x(k-2) + x(k-3);
    end
%     if ( x(k) > 0.2*n )
    if ( x(k) >= 0.8*n )
        poly = @pr;
        % [FIXME] Enable n-1 near y=4n by analytic continuation of all functions.
        dPoly = @(y) pr(n, y, 0 , UQ0) - pr(n, y, 1, UQ1)*sqrt(n + 1);
    end
    step = x(k);
    l = 0;
    ov = inf;
    ox = x(k);
    % [FIXME] Accuracy of the expansions up to machine precision would lower this bound.
    % [FIXME] Accuracy of the expansions up to machine precision lowers this bound.
    while ( ( abs(step) > eps*400*x(k) ) && ( l < 20))%5 ) )
        l = l + 1;
        pe = poly(n, x(k), 0, UQ0);
        % poly' = (p*exp(-Q/2) )' = exp(-Q/2)*(p' -p/2) with orthonormal p
        step = pe/(dPoly( x(k) ) - pe/2);
        if (abs(pe) > abs(ov) + 1e-12)
%             error('Function values increase in Newton method.')
            x(k) = ox;
            break
        end
%         counter = counter + 1;
        ox = x(k);
        x(k) = x(k) -step;
        ov = pe;
    end
    ls(k) = l;
%     if ( l == 5 )
%         warning('Newton method did not convergence as expected.');
%     end
%     w(k) = factor*dPoly( x(k) )*poly(n+1, x(k), 0, UQ0)/exp( x(k) );
    w(k) = factor/dPoly( x(k) )/poly(n+1, x(k), 0, UQ0)/exp( x(k) );
%     [k,counter, x(k), w(k)]
    if ( w(k) == 0 ) && ( k > 1 ) && ( w(k-1) > 0 )
        warning( ['Weights are below realmin*eps from k = ' num2str(k) '.'] );
    end
end
figure; plot(ls);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the expansion of the orthonormal polynomial near zero without e^(x/2)
function p = pl(np, y, alpha, UQ)
z = (1+eps*1i)*y/4/np;
pb = -1i*(pi/2 + sqrt(z).*sqrt(1 - z) - acos(sqrt(z) ) ); % = sqrtphitn

% RL = [1, 0] + [1, 0]*sum(RLkf(z,1:T-1,n)./n.^repmat(reshape(1:T-1,[1,1,T-1]),[2,2,1]),3);
RL = [1, 0];
% if ( ~isempty(UQ) )
%     mx = size(UQ,3) - 1;
%     RL = RL + [1, 0]*sum(RLkf(z, np, UQ)./np.^repmat(reshape(1:mx, [1,1,mx]), [2,2,1]), 3);
% end
% for k = 1:size(UQ,3)
%     Rlt = [0, 0];
%     if abs(z)^(-ti/2) > 1/sqrt(eps)
mk = size(UQ,3); % mk may be zero
if abs(z)^(-mk/2) > 1/sqrt(eps)
	% z is too close to the endpoint so use series expansion.
    for k = 1:mk
        for n = 0:8-k
            RL = RL + UQ(4,:,k,n+1)*z^n/np^k;
        end
    end
else % z is far enough to use formulas with pole cancellation.
    Rko = zeros(1,2,mk);
    sL = zeros(2,2,mk);
    for k = 1:mk
        for m = 1:ceil(3*k/2)
            Rko(:,:,k) = Rko(:,:,k) + UQ(2,:,k,m)./z.^m + UQ(1,:,k,m)./(z-1).^m;
        end
        phi = 2*z - 1 + 2*sqrt(z)*sqrt(z - 1);
        brac = prod(4*alpha^2 - (2*(1:k-1) - 1).^2 )/2^(2*k-2)...
            /factorial(k - 1)/4^k/pb^k;
%         sL(:, :, k) = [1/2^alpha, 0; 0, 2^alpha]*[sqrt(2*z - 1 + ...
%             2*sqrt(z)*sqrt(z - 1) ), 1i/sqrt(2*z - 1 + 2*sqrt(z)*sqrt(z - 1) );
%             -1i/sqrt(2*z - 1 + 2*sqrt(z)*sqrt(z - 1) ), sqrt(2*z - 1 + ...
%             2*sqrt(z)*sqrt(z - 1) ) ]/2/z^(1/4)/(z-1)^(1/4)* ...
%             [(1 - 2*z - 2*sqrt(z)*sqrt(z - 1) )^(alpha/2), 0; 0, ...
%             (1 - 2*z - 2*sqrt(z)*sqrt(z - 1) )^(-alpha/2)];
        sL(:, :, k) = [1/2^alpha, 0; 0, 2^alpha]*[sqrt(phi), 1i/sqrt(phi);
            -1i/sqrt(phi), sqrt(phi) ]/2/z^(1/4)/(z-1)^(1/4)* ...
            [(-phi)^(alpha/2), 0; 0, (-phi)^(-alpha/2)];
        sL(:, :, k) = brac*sL(:, :, k)*[ ( (-1)^k)/k*(alpha^2 + ...
            k/2 - 1/4), (k - 1/2)*1i; -( (-1)^k)*(k - 1/2)*1i, (alpha^2 + ...
            k/2 - 1/4)/k]/sL(:, :, k);
%         sL(:, :, k) = sL(:, :, k)*prod(4*alpha^2-(2*(1:k-1)-1).^2 )/(2^(2*k-2)*factorial(k-1))/4^k/pb^k
        sL(:, :, k) = sL(:, :, k) - mod(k + 1,2)*brac*(4*alpha^2 + 2*k - 1)*...
            eye(2)/2/k;
        
        RL = RL + (Rko(1,:,k) - sL(1,:,k) )/np^k;
        for m = 1:k-1
            RL = RL - Rko(1,:,k-m)*sL(:,:,m)/np^k;
        end
    end
end
RL = RL*[2^(-alpha), 0; 0, 2^alpha]*[sin( (alpha + 1)/2*acos(2*z - 1) - ...
    pi*alpha/2), cos( (alpha + 1)/2*acos(2*z - 1) - pi*alpha/2) ; ...
    -1i*sin( (alpha - 1)/2*acos(2*z - 1) - pi*alpha/2), ...
    -1i*cos( (alpha - 1)/2*acos(2*z - 1) - pi*alpha/2)];
% p = (4*np)^(-1/2 - alpha/2)*exp(np*(1 + 2*log(2)))*sqrt(2)*2^alpha...
p = real( (4*np)^(-1/2 - alpha/2)*sqrt(2)*2^alpha...
    /sqrt(1 - 4i*4^alpha*sum((UQ(2,2,1:mk,1) + UQ(1,2,1:mk,1))./np ...
    .^reshape(1:mk,[1,1,mk]) ) )*(-1)^np*sqrt(1i*np*pb)/z^(1/4)/ ...
    (1 - z)^(1/4)*z^(-alpha/2)*RL*[besselj(alpha,2i*np*pb); ...
    (besselj(alpha-1,2i*np*pb) - alpha/(2i*np*pb)*besselj(alpha, 2i*np*pb) )] );

% p = real( (4*np)^(-1/2 - alpha/2)*sqrt(2)*2^alpha...
%     /sqrt(1 - 4i*sum((UQ(2,2,1:mk,1) + UQ(1,2,1:mk,1))./np ...
%     .^reshape(1:mk,[1,1,mk]) ) )*(-1)^np*sqrt(1i*np*pb)/z^(1/4)/ ...
%     (1 - z)^(1/4)*z^(-alpha/2)*RL*[sin( (alpha + 1)/2*acos(2*z - 1) - ...
%     pi*alpha/2), cos( (alpha + 1)/2*acos(2*z - 1) - pi*alpha/2) ; ...
%     -1i*sin( (alpha - 1)/2*acos(2*z - 1) - pi*alpha/2), ...
%     -1i*cos( (alpha - 1)/2*acos(2*z - 1) - pi*alpha/2)]*...
%     [2^(-alpha)*besselj(alpha,2i*np*pb); 2^alpha*(besselj(alpha-1,2i*np*pb) - ...
%     alpha/(2i*np*pb)*besselj(alpha, 2i*np*pb) ) ] );

% p = (4*np)^(-1/2 - alpha/2)*exp(np*(1 + 2*log(2)))*sqrt(2)*2^alpha...
%     /sqrt(1 - 4i*sum((UQ(2,2,1:mk,1) + UQ(1,2,1:mk,1))./np ...
%     .^reshape(1:mk,[1,1,mk]) ) )*RL;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the expansion of the orthonormal polynomial near 4n without e^(x/2)
function p = pr(np, y, alpha, UQ)
z = (1+eps*1i)*y/4/np;
mxi = 2i*( sqrt(z).*sqrt(1 - z) - acos(sqrt(z) ) ); % = -xin

RR = [1, 0];
mk = size(UQ,3); % mk may be zero
if abs(z)^(-mk/2) > 1/sqrt(eps)
	% z is too close to the endpoint so use series expansion.
    for k = 1:mk
        for n = 0:8-k
            RR = RR + UQ(3,:,k,n+1)*(z-1)^n/np^k;
        end
    end
else % z is far enough to use formulas with pole cancellation.
    Rko = zeros(1,2,mk);
    sR = zeros(2,2,mk);
    for k = 1:mk
        for m = 1:ceil(3*k/2)
            Rko(:,:,k) = Rko(:,:,k) + UQ(2,:,k,m)./z.^m + UQ(1,:,k,m)./(z-1).^m;
        end
        phi = 2*z - 1 + 2*sqrt(z)*sqrt(z - 1);
        mu = 3*gamma(3*k - 1/2)*2^k/27^k/sqrt(pi)/gamma(k*2);
        sR(:, :, k) = [1/2^alpha, 0; 0, 2^alpha]*[sqrt(phi), 1i/sqrt(phi);
            -1i/sqrt(phi), sqrt(phi) ]/2/z^(1/4)/(z-1)^(1/4)* ...
            [phi^(alpha/2), 0; 0, phi^(-alpha/2)];
        sR(:, :, k) = 1/2/mxi^k*sR(:, :, k)*[-(-1)^k*mu/6/k, 1i*mu; ...
            -1i*(-1)^k*mu, -mu/6/k]/sR(:, :, k);
        sR(:, :, k) = sR(:, :, k) + mod(k + 1,2)*mu/6/k/mxi^k*eye(2);
        
        RR = RR + (Rko(1,:,k) - sR(1,:,k) )/np^k;
        for m = 1:k-1
            RR = RR - Rko(1,:,k-m)*sR(:,:,m)/np^k;
        end
    end
end

fn = (np*3/2*mxi)^(2/3);

RR = RR*[2^(-alpha), 0; 0, 2^alpha]*[cos( (alpha + 1)/2*acos(2*z - 1) ), ...
    -1i*sin( (alpha + 1)/2*acos(2*z - 1) ) ; -1i*cos( (alpha - 1)/2* ...
    acos(2*z - 1) ), -sin( (alpha - 1)/2*acos(2*z - 1) )];

p = real( (4*np)^(-1/2 - alpha/2)*sqrt(2)*2^alpha...
    /sqrt(1 - 4i*4^alpha*sum((UQ(2,2,1:mk,1) + UQ(1,2,1:mk,1))./np ...
    .^reshape(1:mk,[1,1,mk]) ) )/z^(1/4)/(z - 1)^(1/4)*z^(-alpha/2)* ...
    RR*[fn^(1/4)*airy(0,fn); fn^(-1/4)*airy(1,fn) ] );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the expansion matrices of R for alpha = 0, with the Laurent series first.
% UQ(:,:,:,:,3:4) gives the Taylor series of R near z = 0 and 1 resp. and we
% take more terms for lower T as they will be multiplied with n^(-T). This is 
% hard-coded for speed and code length, but can be made to get arbitrary orders.
function UQ = getUQalpha0(mk)
% UQ = nan*zeros(2, 2, mk, 8, 4);
% UQ = nan*zeros(2, mk, 8, 4);
UQ = zeros(4, 2, mk, 8);
if ( mk == 0 )
    return
end
UQ(1,:, 1, 1:2) = cat(4, ...
[0.015625, 0.05729166666666667i], ...
[-0.02604166666666667, 0.02604166666666667i]);
UQ(2, :, 1, 1:1) = cat(4, ...
[-0.015625, -0.015625i]);
UQ(3, :, 1, 1:8) = cat(4, ...
[-0.009523809523809525, -0.007142857142857144i], ...
[0.01063492063492063, 0.009365079365079364i], ...
[-0.01129622758194187, -0.01048155019583591i], ...
[0.01174829614829615, 0.01117002997002997i], ...
[-0.01208266968103703, -0.01164522908686174i], ...
[0.01234301148871137, 0.0119973041110328i], ...
[-0.01255318649195314, -0.01227111943654102i], ...
[0.01272750016968636, 0.01249169963638477i]);
UQ(4, :, 1, 1:8) = cat(4, ...
[-3.469446951953614e-18, -0.08333333333333333i], ...
[-0.03888888888888889, -0.0388888888888889i], ...
[-0.07116402116402117, -0.0044973544973545i], ...
[-0.1008289241622575, 0.02615520282186949i], ...
[-0.129241088129977, 0.05509753620864732i], ...
[-0.1569653787325745, 0.08313723470337227i], ...
[-0.1842728591546934, 0.1106475197021934i], ...
[-0.2113083089153498, 0.1378214216323702i]);
if ( mk == 1 )
    return
end
UQ(1,:, 2, 1:3) = cat(4, ...
[-0.0001627604166666668, -0.001085069444444443i], ...
[-0.0007052951388888891, 0.01323784722222223i], ...
[-0.001898871527777778, 0.02278645833333334i]);
UQ(2, :, 2, 1:1) = cat(4, ...
[0.0001627604166666672, 0.004557291666666669i]);
UQ(3, :, 2, 1:7) = cat(4, ...
[1.984126984127046e-05, 0.004007936507936511i], ...
[-3.439153439153486e-05, -0.004221019721019722i], ...
[4.5413316841889e-05, 0.004329315657887089i], ...
[-5.407654074320793e-05, -0.004392654646940364i], ...
[6.109774830863335e-05, 0.0044332342004771i], ...
[-6.692939666830364e-05, -0.004460948839100688i], ...
[7.18695260864394e-05, 0.004480794599931599i]);
UQ(4, :, 2, 1:7) = cat(4, ...
[-0.003472222222222222, -0.006944444444444446i], ...
[-0.003240740740740739, -0.03657407407407406i], ...
[-0.0003747795414462069, -0.09005731922398588i], ...
[0.004772192827748389, -0.1667621987066432i], ...
[0.01205177881103808, -0.2664715497122905i], ...
[0.02138984808385162, -0.3890960516718894i], ...
[0.03274424783742833, -0.5345928880817182i]);
if ( mk == 2 )
    return
end
UQ(1,:, 3, 1:5) = cat(4, ...
[-0.0002585517035590275, -0.000218577443817516i], ...
[-1.465597270447586e-05, -8.356541763117039e-05i], ...
[0.0003698466736593355, 0.006447629575376164i], ...
[-0.002262821903935187, 0.02017682864342207i], ...
[-0.008014160909770449, 0.008014160909770449i]);
UQ(2, :, 3, 1:2) = cat(4, ...
[0.0002585517035590279, -0.0009774102105034725i], ...
[0.0005722045898437501, 0.0005722045898437501i]);
UQ(3, :, 3, 1:6) = cat(4, ...
[0.0005110818194151536, -0.0007267403892403877i], ...
[-0.001045363705958944, 0.000206871642466881i], ...
[0.001590456473343028, 0.0003328986988216683i], ...
[-0.002141460647635705, -0.0008821041836667568i], ...
[0.002696059409107524, 0.001436399596567729i], ...
[-0.003253020840537458, -0.001993718633874486i]);
UQ(4, :, 3, 1:6) = cat(4, ...
[0.005362654320987655, 0.008043981481481482i], ...
[0.02883322310405644, 0.02605544532627866i], ...
[0.0915386169900059, 0.05208755878894767i], ...
[0.2227888285411434, 0.07472004547236026i], ...
[0.4599905388515398, 0.07442783509439153i], ...
[0.8486152325952079, 0.02362629629145646i]);
if ( mk == 3 )
    return
end
UQ(1,:, 4, 1:6) = cat(4, ...
[0.0001123610837959949, 3.220671979488066e-05i], ...
[-4.159894009185921e-05, 1.006970189726202e-05i], ...
[0.0001014470072930732, -0.000150605189947433i], ...
[-0.0001985441019505633, 0.004846322487411191i], ...
[-0.0008181122595390668, 0.02270678924434961i], ...
[-0.000793068006696034, 0.01903363216070482i]);
UQ(2, :, 4, 1:2) = cat(4, ...
[-0.0001123610837959949, -0.0002556506498360341i], ...
[-1.490116119384768e-05, -0.000452995300292969i]);
UQ(3, :, 4, 1:5) = cat(4, ...
[-2.026867991649663e-06, -0.000666664945504233i], ...
[6.757898749300978e-06, 0.001116765488720164i], ...
[-1.295611489560738e-05, -0.001569093104321315i], ...
[2.011295619003682e-05, 0.002022143985881445i], ...
[-2.796022042244435e-05, -0.002475382540803074i]);
UQ(4, :, 4, 1:5) = cat(4, ...
[0.0006884162808641975, 0.0009178883744855927i], ...
[0.002518294998040369, 0.01827557013031548i], ...
[0.005174830777647138, 0.1020198851574237i], ...
[0.006473841918350666, 0.3550440568543063i], ...
[0.001501174033247706, 0.9524842232974917i]);
if ( mk == 4 )
    return
end
UQ(1,:, 5, 1:8) = cat(4, ...
[6.169113234728014e-05, 9.749909808489811e-05i], ...
[-2.735243584977629e-05, -7.023189392787169e-05i], ...
[2.342076600182224e-05, 6.224353508098193e-05i], ...
[-7.31386404171321e-05, -0.0001481373301066077i], ...
[0.0001987962105882525, 0.004849931367753502i], ...
[-0.001539873713001479, 0.03068677518409481i], ...
[-0.01282539667078748, 0.04227630754444762i], ...
[-0.01377542605380878, 0.01377542605380878i]);
UQ(2, :, 5, 1:3) = cat(4, ...
[-6.169113234727972e-05, 0.0002803040936650563i], ...
[-0.0001287894944349925, 0.000105773409207662i], ...
[-0.0001108925789594651, -0.0001108925789594651i]);
UQ(3, :, 5, 1:4) = cat(4, ...
[-0.0002131630878045043, 0.0004025677940905445i], ...
[0.00053567814019176, -0.0003181332313542115i], ...
[-0.0009676950522129104, 0.0001238977582304063i], ...
[0.001509301781329904, 0.0001802687708021319i]);
UQ(4, :, 5, 1:4) = cat(4, ...
[-0.003136156886880289, -0.003920196108600354i], ...
[-0.03190107388066272, -0.03048894498589488i], ...
[-0.1730340465711127, -0.1312375428294408i], ...
[-0.6637971470977652, -0.3913221799067465i]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the expansion matrices of R for alpha = 1.
function UQ = getUQalpha1(mk)
UQ = zeros(4, 2, mk, 8);
if ( mk == 0 )
    return
end
UQ(1,:, 1, 1:2) = cat(4, ...
[-0.046875, 0.06119791666666667i], ...
[-0.02604166666666667, 0.006510416666666667i]);
UQ(2, :, 1, 1:1) = cat(4, ...
[0.046875, 0.01171875i]);
UQ(3, :, 1, 1:8) = cat(4, ...
[0.04047619047619047, -0.03928571428571431i], ...
[-0.04222222222222222, -0.008015873015873025i], ...
[0.04308472479901051, 0.009451041022469598i], ...
[-0.04361026909598338, -0.01006952095523524i], ...
[0.04396981497716192, 0.0104078944789149i], ...
[-0.04423440188249792, -0.01062009860203378i], ...
[0.04443907573247043, 0.01076543684375636i], ...
[-0.04460324252123175, -0.01087128029713879i]);
UQ(4, :, 1, 1:8) = cat(4, ...
[0.08333333333333333, -0.125i], ...
[0.03333333333333333, -0.03333333333333333i], ...
[-0.003174603174603177, -0.03690476190476191i], ...
[-0.03481481481481481, -0.03396825396825398i], ...
[-0.06428731762065096, -0.02897306397306398i], ...
[-0.09264242400750337, -0.0231811821335631i], ...
[-0.1203555135830268, -0.01703582895646388i], ...
[-0.147667867682724, -0.01071718871208446i]);
if ( mk == 1 )
    return
end
UQ(1,:, 2, 1:3) = cat(4, ...
[0.0004882812499999989, 0.01323784722222223i], ...
[-0.01177300347222222, 0.02886284722222222i], ...
[-0.02468532986111112, 0.01139322916666667i]);
UQ(2, :, 2, 1:1) = cat(4, ...
[-0.00048828125, -0.001953125000000002i]);
UQ(3, :, 2, 1:7) = cat(4, ...
[-0.0004166666666666676, -0.0009722222222222267i], ...
[0.0005622895622895622, 0.001434463684463688i], ...
[-0.0006256188256188247, -0.001612467162467158i], ...
[0.0006562001666763589, 0.001701906189049048i], ...
[-0.0006715318212236979, -0.00175433095900883i], ...
[0.0006790103374176761, 0.001788299650036054i], ...
[-0.0006821631748549443, -0.001811910501552051i]);
UQ(4, :, 2, 1:7) = cat(4, ...
[0.003472222222222224, -0.01041666666666667i], ...
[0.03611111111111111, 0.009722222222222219i], ...
[0.09497354497354499, 0.01021825396825397i], ...
[0.1792151675485009, -0.003174603174603194i], ...
[0.2885276040831597, -0.02893518518518523i], ...
[0.422774248874778, -0.06655856841571131i], ...
[0.5818842080253015, -0.1158436958397275i]);
if ( mk == 2 )
    return
end
UQ(1,:, 3, 1:5) = cat(4, ...
[-0.0005976359049479167, -0.000400510246371044i], ...
[0.0008296636887538567, 0.006527672284915132i], ...
[-0.007273111225646224, 0.02336585433394821i], ...
[-0.01923398618344908, 0.01777258037049094i], ...
[-0.008014160909770449, 0.002003540227442612i]);
UQ(2, :, 3, 1:2) = cat(4, ...
[0.0005976359049479174, 0.0009695688883463542i], ...
[-0.0008010864257812502, -0.0002002716064453126i]);
UQ(3, :, 3, 1:6) = cat(4, ...
[-0.000734979989146656, 0.0005610623173123191i], ...
[0.001516718193622956, -0.0004153882471144348i], ...
[-0.002315841844741004, 0.0002340840897612615i], ...
[0.003122487606361571, -4.037681508993723e-05i], ...
[-0.003932667625282309, -0.0001584559939677844i], ...
[0.004744560105086139, 0.0003596418409404624i]);
UQ(4, :, 3, 1:6) = cat(4, ...
[-0.002681327160493828, -0.0004340277777777698i], ...
[-0.01297949735449737, 0.001124338624338644i], ...
[-0.0269951499118166, 0.01863839285714289i], ...
[-0.03246736545347658, 0.06273458794292135i], ...
[-0.008994495608252214, 0.1415801037403221i], ...
[0.07190356645934227, 0.2611881760337451i]);
if ( mk == 3 )
    return
end
UQ(1,:, 4, 1:6) = cat(4, ...
[4.915484675654708e-05, -1.047197192785386e-05i], ...
[-6.56338876166932e-05, -0.0001160778626492966i], ...
[0.000152441503579723, 0.004768120801007313i], ...
[-0.004777727892369407, 0.02510245091630599i], ...
[-0.0228570547614078, 0.0302200650972594i], ...
[-0.01982670016740085, 0.009516816080352408i]);
UQ(2, :, 4, 1:2) = cat(4, ...
[-4.915484675654808e-05, -0.0003443859241626885i], ...
[0.0003713369369506837, 0.0002336502075195314i]);
UQ(3, :, 4, 1:5) = cat(4, ...
[0.0003729230795475838, -0.000150200751078725i], ...
[-0.0007500605012473177, -7.445206877896294e-05i], ...
[0.001127942655094661, 0.0003047918362998525i], ...
[-0.001505607402838352, -0.00053700578303835i], ...
[0.001882808665107142, 0.0007698736891847038i]);
UQ(4, :, 4, 1:5) = cat(4, ...
[-0.0002294720936213962, 0.0003351658950617355i], ...
[-0.01716940402704293, -0.005803066211052317i], ...
[-0.1025708529296493, -0.02748385600382132i], ...
[-0.3638886994207212, -0.06729307069816329i], ...
[-0.9842194497752994, -0.1111429035402386i]);
if ( mk == 4 )
    return
end
UQ(1,:, 5, 1:8) = cat(4, ...
[0.0001335134292826214, 4.24281383919103e-05i], ...
[-0.0001154726233560362, -1.947892172393933e-05i], ...
[0.0001007355757637324, -7.562938905696472e-05i], ...
[1.598066688109512e-05, 0.004724954527204368i], ...
[-0.004761903800468913, 0.03318970787565053i], ...
[-0.03071486300933196, 0.05726074925221438i], ...
[-0.04132627816142633, 0.03194473800409104i], ...
[-0.01377542605380878, 0.003443856513452194i]);
UQ(2, :, 5, 1:3) = cat(4, ...
[-0.0001335134292826253, -3.838322964715599e-06i], ...
[-5.663062135378524e-05, -0.0001810832569996517i], ...
[0.0001355353742837906, 3.388384357094766e-05i]);
UQ(3, :, 5, 1:4) = cat(4, ...
[9.757071356383969e-05, -0.0001007646299157378i], ...
[-0.0003330216038310101, 0.0002083032501474955i], ...
[0.0007064505509263123, -0.000281880394320355i], ...
[-0.001217631725961165, 0.0003210534227464252i]);
UQ(4, :, 5, 1:4) = cat(4, ...
[0.0007840392217200875, 2.86840117026876e-05i], ...
[0.01519177406366366, 0.001661329349769767i], ...
[0.07916956633223091, -0.008926499105987652i], ...
[0.2500575428553918, -0.095433899880876i]);

end
