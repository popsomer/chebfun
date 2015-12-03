classdef spinpref
%SPINPREF   Class for managing SPIN preferences.
%
% Available preferences ([] = defaults):
%
%   dealias                   * If 1, use the 2/3-rule to zero high wavenumbers.
%     [0]                       No dealiasing by default.
% 
%   dt                        * Time-step for time discretization. Default is 
%     []                        empty, i.e., automatically chosen by the code to 
%                               achieve stability. 
%
%   errTol                    * Desired accuracy on the solution.
%     [eps]
%
%   M                         * Number of points for complex means to evaluate
%     [64]                      the phi-functions.
%
%   N                         * Number points for spatial discretization when 
%     []                        using a fixed grid. Default is empty, i.e., 
%                               adaptively chosen by the code to achieve errTol.
%
%   Nmax                      * Maximum number of points for spatial 
%     [4096]                    discretization when using an adaptive grid.
%              
%   plotting                  * Plotting options: 'movie' for plotting a 
%     ['movie']                 movie of the solution, or 'waterfall' to use 
%      waterfall                CHEBFUN WATERFALL command.
%
%   scheme                    * Time-stepping scheme.
%     [@etdrk4]
%      @exprk5s8
%      @krogstad
%      @eglm433
%
%   Ylim                      * Limit of the y-axis when 'plotting' is 'movie'.
%     []                        Default is empty, i.e., automatically chosen by
%                               the code. 
%
% See also SPIN.

% Author: Hadrien Montanelli.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        dealias      % To use dealiasing with 2/3-rule
        dt           % Time-step for time discretizaion
        errTol       % Desired accuracy on the solution
        M            % Number of points for complex means
        N            % Number of points for spatial discretization when using
                     % a fixed gird
        Nmax         % Maximum number of points for spatial discretization 
                     % when using an adaptive grid
        plotting     % Plotting options
        scheme       % Time-stepping scheme
        Ylim         % Limit of the y-axis
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinpref(varargin) 
            if ( nargin == 0 )
                pref.dealias = 0;
                pref.dt = [];
                pref.errTol = eps;
                pref.M = 64;
                pref.N = [];
                pref.Nmax = 4096;
                pref.plotting = 'movie';
                pref.scheme = @etdrk4;
                pref.Ylim = [];
            else
                pref = spinpref();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end