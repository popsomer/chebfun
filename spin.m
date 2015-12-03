function u = spin(varargin)
%SPIN  Solve a (t,x)-PDE with periodicity in space, using a Fourier spectral 
%method for x and an exponential time-diffenrencing scheme for t.
%
%   U = SPIN(PDECHAR) solves the PDE specified by the string PDECHAR, and plots
%   a movie of the solution as it computes it. The space and time intervals, and 
%   the initial condition are chosen to produce beautiful movies. Six default 
%   PDEs are available: 'KS' for Kuramoto-Sivashinsky equation, 'AC' for 
%   Allen-Cahn equation, 'KdV' for Korteweg-de Vries equation, 'Burg' for 
%   viscous Burgers equation, 'CH' for Cahn-Hilliard equation, and 'NLS' for 
%   nonlinear Schrodinger equation. The output U is a CHEBFUN at the final time.
%   See Remark 1 and Examples 1-6.
%
%   U = SPIN(PDECHAR, 'SLOW') does the same but produces a slower movie.
%
%   U = SPIN(PDECHAR, TSPAN, U0) solves the PDE specified by the string 
%   PDECHAR on TSPAN x U0.DOMAIN, with initial condition a CHEBFUN U0, and plots
%   a movie of the solution as it computes it. TSPAN is a vector of time chunks. 
%   If it is of length 2, i.e., TSPAN = [0 TF] for some scalar TF>0, the output 
%   U is a CHEBFUN at the final time TF. To obtain solutions at specific times 
%   0,T1,T2,..,TF (increasing and evenly spaced), use TSPAN = [0 T1 T2 ... TF]. 
%   See Example 7.  
%
%   U = SPIN(PDECHAR, TSPAN, U0, PREF) allows one to use the preferences 
%   specified by the SPINPREF object PREF. See HELP/SPINPREF.
%
%   U = SPIN(PDEFUNLIN, PDEFUNNONLIN, TSPAN, U0) solves the PDE specified by the 
%   function handles PDEFUNLIN an PDEFUNNONLIN, correponding to the linear and 
%   nonlinear parts of the PDE, and plots a movie of the solution as it computes 
%   it. See Remark 2 and Example 8.
%
%   U = SPIN(PDEFUNLIN, PDEFUNNONLIN, TSPAN, U0, PREF) allows one to use the 
%   preferences specified by the SPINPREF object PREF. 
%
% Remark 1: Available strings PDECHAR are
%
%    - 'KS' for Kuramoto-Sivashinsky equation 
%
%           u_t = -u*u_x - u_xx - u_xxxx,
% 
%    - 'AC' for Allen-Cahn equation 
%
%           u_t = 5e-3*u_xx - u^3 + u,
%
%    - 'KdV' for Korteweg-de Vries equation 
%
%           u_t = -u*u_x - u_xxx,
%
%    - 'Burg' for viscous Burgers equation 
%
%           u_t = 1e-3*u_xx - u*u_x,
%
%    - 'CH' for Cahn-Hilliard equation 
%
%           u_t = 1e-2*(-u_xx - 1e-3*u_xxxx + (u^3)_xx),
%
%    - and 'NLS' for the focusing nonlinear Schrodinger equation 
%
%           i*u_t = -u_xx - |u|^2*u.
%
% Remark 2: If PDEFUNNONLIN involves a spatial derivative, it has to be of the 
% form @(u) param*diff(f(u),m), where param is a scalar, f is a nonlinear 
% function of u that does not involve any derivative, and m is any derivative 
% order. It m is 0, it can be of any form. The following syntaxes are allowed:
%
%    pdefunnonlin = @(u) .5*diff(u.^2);
%    pdefunnonlin = @(u) diff(u.^2 + u.^3,2);
%    pdefunnonlin = @(u) exp(u) + cos(sin(2*u));
%
% The following syntax are not:
%
%    pdefunnonlin = @(u) u.*diff(u); 
%    pdefunnonlin = @(u) diff(u.^2,2) + diff(u.^3,2);
%    pdefunnonlin = @(u) diff(u.^2 + diff(u),2);
%
% Example 1: Kuramoto-Sivashinsky (chaotic attractor)
%
%    u = spin('KS');
%
%    solves the Kuramoto-Sivashinsky equation on [0 32*pi] from t=0 to t=300, 
%    with intial condition cos(x/16).*(1 + sin(x/16)). This is equvalent to 
%    typing
%
%    dom = [0 32*pi]; tspan = [0 300];
%    u0 = chebfun('cos(x/16).*(1 + sin(x/16))', dom, 'trig');
%    u = spin('KS', tspan, u0);
%  
% Example 2: Allen-Cahn equation (mestable patterns)
%
%    u = spin('AC');
%
%    solves the Allen-Cahn equation on [-pi pi] from t=0 to t=300, with initial
%    condition tanh(2*sin(x))+3*exp(-27.*(x-4.2).^2)-3*exp(-23.5.*(x-pi/2).^2) 
%    +3*exp(-38.*(x-5.4).^2). This is equivalent to typing
%
%    dom = [-pi pi]; tspan = [0 300];
%    u0 = chebfun(@(x) tanh(2*sin(x)) + 3*exp(-27.*(x-4.2).^2) ...
%       - 3*exp(-23.5.*(x-pi/2).^2) + 3*exp(-38.*(x-5.4).^2), [0 2*pi], 'trig');
%    u = spin('AC', tspan, u0);
%
% Example 3: Korteweg-de Vries equation (two-soliton solution)
%
%    spin('KdV');
%
%    solves the Korteweg-de Vries equation on [-pi pi] from t=0 to t=2*pi*3/A,
%    with initial condition 3*A^2*sech(.5*A*(x+2)).^2+3*B^2*sech(.5*B*(x+1)).^2,
%    and with A = 25 and B = 16. This is equivalent to typing
%
%    dom = [-pi pi]; tspan = [0 2*pi*3/A]; A = 25^2; B = 16^2; 
%    u0 = @(x) 3*A*sech(.5*sqrt(A)*(x+2)).^2 + 3*B*sech(.5*sqrt(B)*(x+1)).^2;
%    u0 = chebfun(u0, dom, 'trig');
%    u = spin('KdV', tspan, u0);
%
% Example 4: Viscous Burgers equation (shock formation and dissipation)
%
%    spin('Burg');
%
%    solves the viscous Burgers equation on [-1 1] from t=0 to t=20, with 
%    initial condition (1-x.^2).*exp(-30.*(x+1/2).^2. This is equivalent to 
%    typing
%
%    dom = [-1 1]; tspan = [0 20];
%    u0 = chebfun('(1-x.^2).*exp(-30.*(x+1/2).^2)', dom, 'trig');
%    u = spin('Burg', tspan, u0);
%
% Example 5: Cahn-Hilliard equation (metastable patterns)
%
%    spin('CH');
%
%    solves the Cahn-Hilliard equation on [-1 1] from t=0 to t=100, with 
%    initial condition (sin(3*pi*x)).^5-sin(pi*x). This is equivalent to 
%    typing
%
%    dom = [-1 1]; tspan = [0 100];
%    u0 = chebfun('(sin(3*pi*x)).^5-sin(pi*x)', dom, 'trig');
%    u = spin('CH', tspan, u0);
%
% Example 6: Nonlinear Schrodinger equation (breather solution)
%
%    spin('NLS');
%
%    solves the Nonlinear Schrodinger equation on [-pi pi] from t=0 to t=1e5,  
%    with initial condition (sin(3*pi*x)).^5-sin(pi*x), and plots the absolute 
%    value of the solution. This is equivalent to typing
%
%    dom = [-pi pi]; tspan = [0 20]; A = 2; B = 1;
%    u0 = @(x) (2*B^2./(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*x)) - 1)*A;
%    u0 = chebfun(u0, dom, 'trig');
%    u = spin('NLS', tspan, u0);
%   
% Example 7: Kuramoto-Sivashinsky with output at different times
%
%    dom = [0 32*pi]; tspan = 0:10:200;
%    u0 = chebfun('cos(x/16).*(1 + sin(x/16))', dom, 'trig');
%    u = spin('KS', tspan, u0);
%
% Example 8: Kuramoto-Sivashinsky with function handles
%
%    dom = [0 32*pi]; tspan = [0 200];
%    u0 = chebfun('cos(x/16).*(1 + sin(x/16))', dom, 'trig');
%    L = @(u) -diff(u,2) - diff(u,4);
%    N = @(u) -.5*diff(u.^2);
%    u = spin(L, N, tspan, u0);

% Author: Hadrien Montanelli.

%% Inputs:

% Parse Inputs:
pref = [];
count = 0;
iterplot = 20;
for j = 1:nargin
    item =  varargin{j};
    if ( isa(item, 'char') && j == 1 )
        pdefunLin = item;
        pdefunNonlin = item;
        % In that case, we are in DEMO mode:
        if ( nargin <= 2 )
            [tspan, u0, pref] = parseInputs(item);
            % 'SLOW' option passed:
            if ( nargin == 2 )
                if ( strcmpi(varargin{2}, 'SLOW') == 1 )
                    iterplot = 4;
                else
                    error('SPIN:spin', 'Second argument should be ''SLOW''.') 
                end
            end
        end
    elseif ( isa(item, 'function_handle') && count < 1 )
        pdefunLin = item;
        count = 1;
    elseif ( isa(item, 'function_handle') )
        pdefunNonlin = item;
    elseif ( isa(item, 'double') ) 
        tspan = item;
    elseif ( isa(item, 'chebfun') )
        u0 = item;
    elseif ( isa(item, 'spinpref') )
        pref = item;
    end
end

% Space interval and final time TF:
dom = u0.domain; 
TF = tspan(end);

%% Preferences:

% Create a SPINPREF object if none was passed:
if ( isempty(pref) )
   pref = spinpref; 
end

% Dealiasing (0=NO, 1=YES):
dealias = pref.dealias;

% Error tolerance:
errTol = pref.errTol;

% Points for complex means:
M = pref.M;

% Minimum and maximum grid points:
if ( isempty(pref.N) == 1 )
    % Adaptivity in space: 
    u0 = chebfun(u0, dom, 'trig','eps', errTol);
    Nmin = min(max(256, 2^(1+floor(log2(length(u0))))), pref.Nmax);
    Nmax = pref.Nmax;
else
    % Non-adpativity in space, i.e., use the N given by the user:
    Nmin = pref.N;
    Nmax = pref.N;
end

% Plotting options ('movie', or []):
plottingstyle = pref.plotting;

% Create a time-stepping scheme:
scheme = pref.scheme;
scheme = scheme();

% Get the number of steps of the method:
q = scheme.steps;

% Limits of y-axis if specified:
if ( isempty(pref.Ylim) == 0 )
    ylim1 = pref.Ylim(1);
    ylim2 = pref.Ylim(2);
end

%% Pre-processing:

% Start with N=Nmin:
N = Nmin;

% Indexes for dealiasing:
ind = false(N, 1);
ind(floor(N/2)+1-ceil(N/6):floor(N/2)+ceil(N/6)) = 1;
    
% Operators (linear part L, nonlinear part Nv and Nc):
L = getLinearPart(pdefunLin, N, dom);
Nop = getNonlinearPart(pdefunNonlin, N, dom);
Nv = Nop{1};
Nc = Nop{2};

% Time-step dt:
dtmin = 5e-6;
dtmax = 1e-1;
if ( isreal(L) == 1 )
    dtstab = round(600/(max(max(abs(L)))));
else
    dtstab = round(20/(max(max(abs(L)))));
end
if ( isempty(pref.dt) == 1 )
    % Adaptive in time:
    dtstab = min(max(dtstab, dtmin), dtmax);
    dt = TF - tspan(end-1);
    while ( dt > dtstab )
        dt = dt/2;
    end
else
    % Not adpative in time, i.e., use the dt given by the user:
    dt = pref.dt;
end

% Compute linear part for complex means:
LR = computeLR(L, N, M, dt);

% Compute coefficients for the time-stepping scheme:
[A, B, U, V, E] = computeCoeffs(scheme, L, LR, dt);

% Set-up spatial grid, and initial condition:
xx = trigpts(N, dom);
v = u0(xx);
c = fft(v);
vscale = max(abs(v));

% Get enough initial data when using a multistep scheme:
if ( q > 1 )
    c = spinscheme.startMultistep(c, Nc, Nv, L, LR, dt, q);
end

% Plot initial condition if using MOVIE option:
if ( strcmpi(plottingstyle, 'movie') == 1 )
    figure, p = plot(xx, v, 'linewidth', 3);
    set(gca, 'FontSize', 16)
    if ( isempty(pref.Ylim) == 1 )
        ylim1 = min(v) - .1*vscale;
        ylim2 = max(v) + .1*vscale;
    end
    axis([xx(1), xx(end), ylim1, ylim2])
    xlabel('x'), ylabel('u(t,x)'), grid on
    title(sprintf('t = %.3f', 0)), drawnow
    disp('Type <enter> when ready.'), pause
end

% Values to output:
vout{1} = v;

% Values to plot if using WATERFALL option:
if ( strcmpi(plottingstyle, 'waterfall') == 1 )
    vwater{1} = v;
    twater = 0;
end

%% Time-stepping loop:

iter = 0;
tnew = q*dt;
while ( tnew <= TF  + dt/2 )

    % One step in time:
    cnew = spinscheme.oneStep(c, Nc, Nv, A, B, U, V, E);
    
    % Dealiasing procedure:
    if ( dealias == 1 )
        cnew(ind,1) = 0;
    end
        
    % Check happiness:
    vnew = ifft(cnew(:,1));
    u = trigtech({vnew, trigtech.vals2coeffs(vnew)});
    options = trigtech.techPref();
    options.eps = errTol;
    try
        ishappy = happinessCheck(u,[],[],[],options);
    catch ME
        if ( strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:standardCheck:nanEval') )
            error('SPIN:spin', ['The solution blew up. Try to double ' ...
                'the number of points N or/and half the time-step dt.']);
        else
            rethrow(ME)
        end
    end
    
    % If happy, or if we're using the maximum number of points, go on:
    if ( ishappy == 1 || ( ishappy == 0 && N >= Nmax ) )
        iter = iter + 1;
        c = cnew;
        % Plot every ITERPLOT iterations if using MOVIE option:
        if ( strcmpi(plottingstyle, 'movie') == 1 && mod(iter,iterplot) == 0 )
            if ( isempty(pref.Ylim) == 0 )
                plotmovie(xx, tnew, vnew, dt, N, p, ylim1, ylim2);
            else
                [ylim1, ylim2] = plotmovie(xx, tnew, vnew, dt, N, p, ...
                    ylim1, ylim2);
            end
        % Store the values every ITERPLOT iterations if using WATERFALL:
        elseif ( strcmpi(plottingstyle, 'waterfall') == 1 && ...
            mod(iter,iterplot) == 0 )
            vwater{iter/iterplot + 1} = vnew;
            twater = [twater, tnew];
        end
        % Output the solution if TNEW correponds to one of the entries of TSPAN:
        if ( any(abs(tspan-tnew) < 1e-10) )
            [~, pos] = max(abs(tspan-tnew) < 1e-10);
            % Use a cell because the VNEW's can have different sizes:
            vout{pos} = vnew; 
        end
        % Update the time:
        tnew = tnew + dt;
        continue
   
    % If not happy, redo the last time step:
    elseif ( ishappy == 0 && N < Nmax )
        
        % Double the number of points N, if adaptive in space:
        if ( isempty(pref.N) == 1 )
            % Create a new vector of coefficients of size 2*N:
            N = 2*N;
            xx = trigpts(N, dom);
            v = ifft(c);
            u = trigtech({v, trigtech.vals2coeffs(v)});
            v = feval(u, trigpts(N));
            c = fft(v);
            % Update quantities which depend on N:
            L = getLinearPart(pdefunLin, N, dom);
            Nop = getNonlinearPart(pdefunNonlin, N, dom);
            Nv = Nop{1};
            Nc = Nop{2};
            LR = computeLR(L, N, M, dt);
            [A, B, U, V, E] = computeCoeffs(scheme, L, LR, dt);
            ind = false(N, 1);
            ind(floor(N/2)+1-ceil(N/6):floor(N/2)+ceil(N/6)) = 1;
        end
        
        % Half the time-step dt, if adaptive in time:
        if ( isempty(pref.dt) == 1 )
            % Half the time-step:
            dt = max(dt/2, dtmin);
            % Update quantities which depend on dt:
            LR = computeLR(L, N, M, dt);
            [A, B, U, V, E] = computeCoeffs(scheme, L, LR, dt);
        end
        
    end
  
end

%% Post-processing:

% Be sure that the solution at t=TF has been plotted if using MOVIE option:
if ( strcmpi(plottingstyle, 'movie') == 1 )
    plotmovie(xx, TF, vnew, dt, N, p, ylim1, ylim2);
end

% Use WATERFALL if using WATERFALL option:
if ( strcmpi(plottingstyle, 'waterfall') == 1 )
    uwater = []; 
    for k = 1:size(vwater, 2)
        uwater = [ uwater, chebfun(vwater{k}, dom, 'trig') ];  %#ok<*AGROW>
    end
    figure, waterfall(uwater, twater), axis([xx(1), xx(end), 0, TF])
    set(gca, 'FontSize', 16), box on
    xlabel('x'), ylabel('t'), zlabel('u(t,x)')
    view([10 70])
end

% Output a CHEBFUN from values VOUT:
if ( length(tspan) == 2 )
    u = chebfun(vout{end}, dom, 'trig');
else
    u = []; 
    for k = 1:size(vout, 2)
        u = [ u, chebfun(vout{k}, dom, 'trig') ]; 
    end
end

end

%% Parse Inputs for DEMO mode:

function [tspan, u0, pref] = parseInputs(pdeobj)

pref = spinpref();
if ( strcmpi(pdeobj, 'KS') == 1 )
    u0 = chebfun('cos(x/16).*(1 + sin((x-1)/16))', [0 32*pi], 'trig');
    tspan = [0 300];
    pref.dt = 5e-2;
elseif ( strcmpi(pdeobj, 'AC') == 1 )
    u0 = chebfun(@(x) tanh(2*sin(x)) + 3*exp(-27.*(x-4.2).^2) ...
        - 3*exp(-23.5.*(x-pi/2).^2) + 3*exp(-38.*(x-5.4).^2), [0 2*pi], 'trig');
    pref.dt = 5e-2;
    tspan = [0 300];
elseif ( strcmpi(pdeobj, 'KdV') == 1 )
    A = 25^2; B = 16^2;
    u0 = @(x) 3*A*sech(.5*sqrt(A)*x).^2 + 3*B*sech(.5*sqrt(B)*(x-1)).^2;
    u0 = chebfun(u0, [-pi pi], 'trig');
    tspan = [0 2*pi*3/A];
    pref.dt = 7e-6;
elseif ( strcmpi(pdeobj, 'Burg') == 1 )
    u0 = chebfun('(1-x.^2).*exp(-30.*(x+1/2).^2)', [-1 1], 'trig');
    pref.dt = 1e-2;
    tspan = [0 30];
elseif ( strcmpi(pdeobj, 'CH') == 1 )
    u0 = chebfun('(sin(4*pi*x)).^5-sin(pi*x)', [-1 1], 'trig');
    pref.dt = 2e-2;
    tspan = [0 70];
elseif ( strcmpi(pdeobj, 'NLS') == 1 )
    A = 2; B = 1;
    u0 = @(x) (2*B^2./(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*x)) - 1)*A;
    u0 = chebfun(u0, [-pi pi], 'trig');
    pref.dt = 2e-3;
    tspan = [0 60];
else
    error('SPIN:getLinearPart', 'Unrecognized PDE.')
end

end