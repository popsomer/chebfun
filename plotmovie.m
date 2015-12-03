function [ylim1, ylim2] = plotmovie(xx, t, vnew, dt, N, p, ylim1, ylim2)
%PLOTMOVIE   Plot a movie of the solution.
%   [YLIM1, YLIM2] = PLOTMOVIE(XX, T, VNEW, DT, N, P, YLIM1, YLIM2) plots the 
%   solution VNEW at points XX and time T. N is the number of points in space, 
%   DT is the timestep, P is handle to the figure, and YLIM1 and YLIM2 are 
%   limits for the y-axis.

% Author: Hadrien Montanelli.

% Update the data:
set(p, 'xdata', xx)
if ( isreal(vnew) == 1 )
    set(p, 'ydata', vnew)
else
    vnew = real(vnew);
    set(p, 'ydata', vnew)
end

% Change axes if necessary:
if ( nargout == 2 )
    minvnew = min(vnew);
    maxvnew = max(vnew);
    if ( maxvnew > ylim2 )
        vscalenew = max(abs(minvnew), maxvnew);
        ylim2 = maxvnew + .1*vscalenew;
    end
    if ( minvnew < ylim1 )
        vscalenew = max(abs(minvnew), maxvnew);
        ylim1 = minvnew - .1*vscalenew;
    end
end

% Update title and plot:
title(sprintf('N = %i, dt = %1.0e, t = %.4f', N, dt, t))
axis([xx(1), xx(end), ylim1, ylim2]), drawnow

end