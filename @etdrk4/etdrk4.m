classdef etdrk4 < spinscheme
%ETDRK4   Implements the ETDRK4 scheme.
%
% See also SPINSCHEME, EGLM433, EXPRK5S8, KROGSTAD, PECEC433.

% Author: Hadrien Montanelli.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function scheme = etdrk4()
            scheme.order = 4;
            scheme.internalStages = 4;
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