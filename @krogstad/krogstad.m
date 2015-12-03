classdef krogstad < spinscheme
%KROGSTAD  Implements the KROGSTAD scheme.
%
% See also SPINSCHEME, EGLM433, ETDRK4, EXPRK5S8, PECEC433.

% Author: Hadrien Montanelli.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function scheme = krogstad()
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