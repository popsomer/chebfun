function L = getLinearPart(pdeobj,N,dom)
%GETLINEARPART   Get the linear part of the PDE.
%   L = GETLINEARPART(PDEOBJ,N,DOM), where PDEOBJ is a string or a function 
%   handle, and  N is the number of points to discretize the space interval DOM,
%   outputs a vector L that reprents the linear part of the PDE specified by 
%   PDEOBJ. 

% Author: Hadrien Montanelli.

% Predefined cases:
if ( isa(pdeobj, 'char') )
    if ( mod(N, 2) == 0 )
        k = [ 0:N/2-1, -N/2:-1 ]'/(dom(2) - dom(1))*(2*pi);
    else
        k = [ 0:(N+1)/2-1, -(N+1)/2+1:-1 ]'/(dom(2) - dom(1))*(2*pi);
    end
    if ( strcmpi(pdeobj, 'KS') == 1 )
        L = @(k) k.^2 - k.^4;
        L = L(k);
    elseif ( strcmpi(pdeobj, 'AC') == 1 )
        epsilon = 5e-3;
        L = @(k) -epsilon*k.^2;
        L = L(k); 
    elseif ( strcmpi(pdeobj, 'KdV') == 1 )
        % In that case, since it involves an odd derivative, zero the N/2
        % coefficients when N is even (see TRIGSPEC/DIFFMAT):
        if ( mod(N, 2) == 0 )
            k = [ 0:N/2-1, 0, -N/2+1:-1 ]'/(dom(2) - dom(1))*(2*pi);
        end
        L = @(k) 1i*k.^3;
        L = L(k);
    elseif ( strcmpi(pdeobj, 'Burg') == 1 )
        epsilon = 1e-3;
        L = @(k) -epsilon*k.^2;
        L = L(k);
    elseif ( strcmpi(pdeobj, 'CH') == 1 )
        D = 0.01;
        gamma = 0.001;
        L = @(k) D*(k.^2 - gamma*k.^4);
        L = L(k);
    elseif ( strcmpi(pdeobj, 'NLS') == 1 )
        L = @(k) -1i*k.^2;
        L = L(k);
    else
        error('SPIN:getLinearPart', 'Unrecognized PDE.')
    end
    
% Use TRIGSPEC to get L in (diagonal) matrix form:
elseif ( isa(pdeobj, 'function_handle') )
    L = chebop(pdeobj,dom);
    L = linop(L);
    prefs = cheboppref();
    prefs.discretization = @trigspec;
    L = full(diag(matrix(L,N,prefs)));
    if ( mod(N,2) == 0 )   
        L = fftshift(L);
    else
        L = ifftshift(L);
    end
end

end