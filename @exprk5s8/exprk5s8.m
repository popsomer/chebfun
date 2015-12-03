classdef exprk5s8 < spinscheme
%ETDRK4   Implements the EXPRK5S8 scheme.
%
% See also SPINSCHEME, EGLM433, ETDRK4, KROGSTAD, PECEC433.

% Author: Hadrien Montanelli.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function scheme = exprk5s8()
            scheme.order = 5;
            scheme.internalStages = 8;
            scheme.steps = 1;
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