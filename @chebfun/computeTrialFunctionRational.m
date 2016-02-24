function [p, q, rh, pqh, h, interpSuccess] = computeTrialFunctionRational(f, xk, m, n)

% Vector of alternating signs.
N = m + n;
sigma = ones(N + 2, 1);
sigma(2:2:end) = -1;

% Orthogonal matrix with respect to <,>_{xk}.
[C, ignored] = qr(fliplr(vander(xk)));

fk = feval(f, xk);
dom = f.domain([1, end]);
% Left and right rational interpolation matrices.
ZL = C(:,m+2:N+2).'*diag(fk)*C(:,1:n+1);
ZR = C(:,m+2:N+2).'*diag(sigma)*C(:,1:n+1);

% Solve generalized eigenvalue problem.
[v, d] = eig(ZL, ZR);

% Compute all possible qk and and look for ones with unchanged sign.
qk_all = C(:,1:n+1)*v;
pos =  find(abs(sum(sign(qk_all))) == N + 2);  % Sign changes of each qk.
interpSuccess = 1;

% This shouldn't happen, theoretically
if ( (length(pos) > 1) )
    error('CHEBFUN:CHEBFUN:remez:eigensolver', ...
        'More than one vector doesn''t change sign');
    interpSuccess = 0;
end

if ( isempty(pos) )
    disp('poles on the approximation domain');
    %[~, pos] = max(abs(sum(sign(qk_all))));
    %plusSign = sum(qk_all(:, pos) > 0);
    %minusSign = sum(qk_all(:, pos) < 0);
    %if ( plusSign > minusSign )
    %    qk_all(:, pos) = abs(qk_all(:, pos));
    %else
    %    qk_all(:, pos) = -abs(qk_all(:, pos));
    %end
    interpSuccess = 0;
end

if (interpSuccess == 1)
    
    qk = qk_all(:,pos);       % Keep qk with unchanged sign.
    h = d(pos, pos);          % Levelled reference error.
    disp(h);
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

    rk = feval(p, xk)./feval(q, xk);
    %absInterpError = norm( (fk - rk) - h*sigma, inf);
 
    nn = round(length(xk)/2);
    fvals = fk - h*sigma;
    xx = xk; xx(nn) = [];
    fx = fvals; fx(nn) = [];
    A = berrut(xx, fx, m, n);
    v = null(A);
    rh = @(t) bary(t, fx, xx, v);
    pqh = @(x) p(x)./q(x);
    
else
    % we won't use these values, since the eigensolver failed to provide a
    % valid, pole free, interpolant; hence, we'll go ahead and perturb the
    % previous, valid, reference
    p = 0;
    q = 0;
    rh = 0;
    pqh = 0;
    h = 0;
end

%r = chebfun(fh, dom, 'splitting', 'on');

%if (  absInterpError > 1e-7 )
%    str = sprintf( 'abs interp error: %.16g', absInterpError);
%    disp(str)
%    interpSuccess = 0;
%else
%    interpSuccess = 1;
%end


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
