classdef pecec433 < spinscheme
%PECEC433  Implements the PECEC433 scheme.
%
% See also SPINSCHEME, EGLM433, ETDRK4, EXPRK5S8, KROGSTAD.

% Author: Hadrien Montanelli.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function scheme = pecec433()
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