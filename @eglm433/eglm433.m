classdef eglm433 < spinscheme
%EGLM433   Implements the EGLM433 scheme.
%
% See also SPINSCHEME, ETDRK4, EXPRK5S8, KROGSTAD, PECEC433.

% Author: Hadrien Montanelli.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function scheme = eglm433()
            scheme.order = 4;
            scheme.internalStages = 3;
            scheme.steps = 3;
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public )

        % Compute coefficients for time-stepping:
        [A, B, U, V, E] = computeCoeffs(scheme, L, LR, dt)
        
    end
    
end