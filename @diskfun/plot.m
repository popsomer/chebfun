function varargout = plot( f, varargin )
%PLOT  Surface plot of a DISKFUN.
%
%   PLOT(F) if F is a real-valued DISKFUN then this is the surface plot and is
%   the same as surf(F). If F is a complex valued then this returns a domain
%   colouring plot of F.
%
%   PLOT(F) if F is a complex-valued SPHEREFUN then we do Wegert's phase portrait
%   plots.
%
%   PLOT(F, S) Plotting with option string plots the column and row slices, and
%   pivot locations used in the construction of F.
%
%   When the first argument in options is a string giving details about
%   linestyle, markerstyle or colour then pivot locations are plotted. Various
%   line types, plot symbols and colors may be obtained with plot(F,S) where S
%   is a character string made from one element from any or all the following 3
%   columns, similar as in the usual plot command:
%
%           b     blue          .     point              -     solid
%           g     green         o     circle             :     dotted
%           r     red           x     x-mark             --    dashed
%           c     cyan          +     plus               -.    dashdot
%           m     magenta       *     star             (none)  no line
%           y     yellow        s     square
%           k     black         d     diamond
%                               v     triangle (down)
%                               ^     triangle (up)
%                               <     triangle (left)
%                               >     triangle (right)
%                               p     pentagram
%                               h     hexagram
%
%   For phase portraits see: E. Wegert, Visual Complex Functions: An
%   introduction with Phase Portraits, Springer Basel, 2012, or for MATLAB code
%   to produce many different styles of phase portraits go to:
%   http://www.visual.wegert.com
% 
% See also SURF, MESH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Make a user option?
plot_full_grid = false;

if ( ~isempty(varargin) )
    
    % Plot the pivots on the surface of the sphere.    
    if ( length(varargin{1}) < 5 )
        dom = f.domain;
        holdState = ishold;
        if ~holdState
            hold on;
        end
        N = 100; %used to generate solid disk as well as plot slicing lines
        th = trigpts(N, dom); th=[th; dom(2)];
        % If the plot is not being added to another then plot a solid 
        % sphere so the lines are more easily discernable.
        if ~holdState
            %
            % Generate a unit disk
            N = 100; 
            th = trigpts(N, dom); th=[th; dom(2)];
            scl = 0.99; %plot disk slightly smaller so unit lines show up more clearly
            r = scl*exp(1i*th);
            
            clr = [255 255 204]/255;
            %plot(r, 'color',clr)
            fill(real(r),imag(r), clr, 'Edgecolor', 'None'); %is there a way to eliminate border?
            
            %
            %[XX,YY,ZZ] = sphere(101);
            % Color of the sphere will be yellow:
            
            % Plot the sphere, make it slightly smaller than unit so lines
            % show up more clearly.
            %scl = 0.99;
            %surf(scl*XX,scl*YY,scl*ZZ,1+0*XX,'EdgeColor','None','FaceColor',clr);
        end

        %% Column, row, pivot plot

        % Only option with <=3 letters is a colour, marker, line
        ll = regexp( varargin{1}, '[-:.]+', 'match' );
        cc = regexp( varargin{1}, '[bgrcmykw]', 'match' );       % color
        mm = regexp( varargin{1}, '[.ox+*sdv^<>ph]', 'match' );  % marker
        
        if ( ~isempty(ll) )
            if ( strcmpi(ll{1},'.') )
                % so we have . first. Don't plot a line.
                ll = {};
            elseif ( strcmpi(ll{1},'.-') )
                % so we have . first. Don't plot a line.
                ll{1} = '-';
            end
        end
        plotline = ~isempty(ll);  % plot row and col pivot lines?
        if ( isempty(mm) )
            mm{1}= '.';
        end
        if ( isempty(cc) ) 
            cc{1}= 'b';
        end
        % Plot the crosses and pivots on domain.

        % Full grid on the disk
        m = length(f.cols);
        n = length(f.rows);
        [TT, RR] = meshgrid([trigpts(n,dom(1:2)); dom(2)],linspace(dom(3),dom(4),m/2+1));

        % Plot pivots:
        
        % Convert pivots to Cartesian coordinates
        pivots = f.pivotLocations;
        pivotsCart = zeros(size(pivots,1),2);
            pivotsCart(:,1) = (pivots(:,2)).*cos(pivots(:,1)); %x=rcosth
            pivotsCart(:,2) = (pivots(:,2)).*sin(pivots(:,1)); %y=rsinth
            XX = RR.*cos(TT);
            YY = RR.*sin(TT);
      
        % Also plot points marking the pivots shifted by pi in longitude
        % since these are technically also included in the GE algorithm.
        pivotsCart = [pivotsCart;[-pivotsCart(:,1:2)]];
            
        defaultopts = { 'MarkerSize', 7 };
        extraopts = { 'Marker', mm{:}, 'LineStyle', 'none', 'Color', cc{:} };
        if ( length(varargin) < 2 )
            h = plot( pivotsCart(:,1), pivotsCart(:,2), ...
                        extraopts{:}, defaultopts{:} );
        else
            h = plot( pivotsCart(:,1), pivotsCart(:,2), ...
                       extraopts{:}, defaultopts{:}, ...
                       varargin{2:end} );
        end
        if ( plotline )
            % Use parametrization for great circles passing through each
            % pivot location and the poles for the column pivot lines and
            % horizontal cirlces at the right latitude for the row pivot
            % lines.
            
                
            colslices = [];
            rowCircs = [];
            
            
            for k=1:size(pivots,1) %note on pole: if r=0 the result will be zero point 
               % if abs(pivots(k,2))<100*eps  %faster to just assign it
                    %colslices=[colslices; 0 0; nan nan];
                   % rowCircs=[rowCircs; 0 0; nan nan];
                 if abs(pivotsCart(k,1))>100*eps  %special case if x=0 
                    colslices = [colslices; [-cos(pivots(k,1)); cos(pivots(k,1))]...
                        pivotsCart(k,2)/pivotsCart(k,1)*[-cos(pivots(k,1)); cos(pivots(k,1))]];
                    colslices = [colslices ; nan nan];
                    %rowcircs use theta pts set up earlier
                    rowCircs = [rowCircs; pivots(k,2)*exp(1i*th)];
                    rowCircs = [rowCircs;nan];
               else   %case where x=0 so slope is undef
                    colslices = [colslices; zeros(2,1) [-1 ; 1]];
                    colslices = [colslices; nan nan];
                    rowCircs = [rowCircs; pivots(k,2)*exp(1i*th)];
                    rowCircs = [rowCircs;nan];
                    end 
                end
            
            
            
            opts = {}; 
            if ( ~isempty(ll) )
                opts = { opts{:}, 'LineStyle', ll{:}, 'color', cc{:} };
            end
            % Plot column lines:
            plot(colslices(:,1),colslices(:,2), opts{:} );
            % Plot row lines:
            plot(rowCircs, opts{:} );            
        end
        if plot_full_grid
            % Plot grayed out grid
            clrGrid = [192 192 192]/256;
            plot(XX,YY,'-','Color',clrGrid);
            plot(XX',YY','-','Color',clrGrid);
        end
        if ~holdState
            hold off;
        end
    else
        %% Standard surface plot 
        h = surf(f, varargin{:});
    end
else
    h = plot@separableApprox( f );
end

if ( nargout > 0 )
    varargout = { h }; 
end

end