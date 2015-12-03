function Nop = getNonlinearPart(pdeobj,N,dom)
%GETNONLINEARPART   Get the nonlinear part of the PDE.
%   Nop = GETNONLINEARPART(PDEOBJ,N,DOM), where PDEOBJ is a string or a function 
%   handle, and  N is the number of points to discretize the space interval DOM,
%   outputs a cell-array NOP that reprents the nonlinear part of the PDE 
%   specified by PDEOBJ. NOP{1} reprents the nonlinear part in value space and 
%   is a function handle, while NOP{2} represents the nonlinear part in
%   coefficient space and is a vector.
%
% Notes: When PDEOBJ is a function handle and involves a spatial derivative, 
% it has to be of the form @(u) param*diff(f(u),m), where param is a double,
% f is a nonlinear function of u that does not involve any derivative, and m 
% is any derivative order. It it doesn't have any derivative, it can be of any 
% form. The following syntax are allowed:
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

% Author: Hadrien Montanelli.

% Predefined cases:
if ( isa(pdeobj, 'char') )
    if ( mod(N, 2) == 0 )
        k = [ 0:N/2-1, 0, -N/2+1:-1 ]'/(dom(2) - dom(1))*(2*pi);
    else
        k = [ 0:(N+1)/2-1, -(N+1)/2+1:-1 ]'/(dom(2) - dom(1))*(2*pi);
    end
    if ( strcmpi(pdeobj, 'KS') == 1 )
        Nv = @(z) z.^2;
        Nc = @(k) -0.5i*k;
    elseif ( strcmpi(pdeobj, 'AC') == 1 )
        Nv = @(z) z - z.^3;
        Nc = @(k) 1;
    elseif ( strcmpi(pdeobj, 'KdV') == 1 )
        Nv = @(z) z.^2;
        Nc = @(k) -0.5i*k;
    elseif ( strcmpi(pdeobj, 'Burg') == 1 )
        Nv = @(z) z.^2;
        Nc = @(k) -0.5i*k;
    elseif ( strcmpi(pdeobj, 'CH') == 1 )
        % In that case, since it involves an even derivative, don't zero the 
        % N/2 coefficient when N is even (see TRIGSPEC/DIFFMAT):
        if ( mod(N, 2) == 0 )
            k = [ 0:N/2-1, -N/2:-1 ]'/(dom(2) - dom(1))*(2*pi);
        end
        D = 0.01;
        Nv = @(z) z.^3;
        Nc = @(k) -D*k.^2;
    elseif ( strcmpi(pdeobj, 'NLS') == 1 )
        lambda = 1;
        Nv = @(z) lambda*1i*abs(z).^2.*z;
        Nc = @(k) 1;
    else
        error('SPIN:getNonlinearPart', 'Unrecognized PDE.')
    end
    Nc = Nc(k);
    
% Use REGEXP:
elseif ( isa(pdeobj, 'function_handle') )
    func = functions(pdeobj);
    wrk = func.workspace{1};
    names = fieldnames(wrk);
    if ( ~isempty(names) )
        SIZE = size(names,1);
        for k = 1:SIZE 
          eval(sprintf('%s = wrk.(names{k});',names{k}));
        end
    end
    Nstr = func2str(pdeobj);
    Nv = '';
    Nc = '';
    LENGTH = length(Nstr);
    END = regexp(Nstr,'diff(','end');
    if ( isempty(END) )
        Nv = pdeobj;
        Nc = 1;
    else
        COMMA = regexp(Nstr,',','start');
        PARENTHESIS = 1;
        COUNTER = 0;
        for k = END+1:LENGTH
            if ( strcmp(Nstr(k),'(') )
                PARENTHESIS = PARENTHESIS + 1;
                Nv = [Nv,Nstr(k)]; %#ok<*AGROW>
            elseif ( strcmp(Nstr(k),')') )
                COUNTER = COUNTER + 1;
                if ( COUNTER < PARENTHESIS )
                    Nv = [Nv,Nstr(k)];
                end
            elseif ( strcmp(Nstr(k),',') )
                break
            else
                Nv = [Nv,Nstr(k)];
            end
        end
        Nv = ['@(u)',Nv];
        Nv = eval(Nv);
        if ( ~isempty(COMMA) )
            Nc = [Nc,Nstr(1:END),'u',Nstr(COMMA:LENGTH)];
        else
            Nc = [Nc,Nstr(1:END),'u)'];
        end
        Nc = eval(Nc);
        % Since Nc is always linear, call GETLINEARPART:
        Nc = getLinearPart(Nc,N,dom);
    end
end

Nop = {Nv, Nc};

end