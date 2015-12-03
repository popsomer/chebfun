classdef spinscheme
%SPINSCHEME   Abstract class for representing timestepping schemes.
%
% See also EGLM433, ETDRK4, EXPRK5S8, KROGSTAD.

% Author: Hadrien Montanelli.
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        order           % order of the method
        internalStages  % number of internal stages
        steps           % number of previous time steps used (1 if one-step
                        % method, > 1 if multistep method)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = true )
        
        % Compute coefficients for the time-stepping (step 1, different for
        % each time-stepping schemme):
        [A, B, U, V, E] = computeCoeffs(scheme, L, LR, dt)
        
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Timestepping:
        un = oneStep(un, Nc, Nv, A, B, U, V, E)
        
        %  Get enough initial data when using a multistep scheme:
        c = startMultistep(c, Nc, Nv, L, LR, dt, q)
        
        % Precompute coefficients for the time-stepping (step 2, same for each
        % time-stepping scheme):
        [A, B, U, V, E] = computeMissingCoeffs(A, B, C, U, V, L, phi1, phit, dt)
        
    end
    
end