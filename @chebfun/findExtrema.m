function rts = findExtrema(f, p, q, rh, h, xk)
% finds all the local maxima and minima
% of f-p./q
% xk is the present reference
% rh is a handle to p/q
% h is the current leveled error

err_handle = @(x) feval(f, x) - rh(x);
rts = [];
%doms = unique(sort([f.domain(:); xk]));
doms = unique([f.domain(1); xk; f.domain(end)]);
%doms = sort([doms; 0]);

% Initial trial
if ( isempty(xk) )
    ek = chebfun(@(x) err_handle(x), doms.', 'splitting', 'on');
    rts = roots(diff(ek), 'nobreaks');
    rts = unique([f.domain(1); rts; f.domain(end)]);
end
   
if ( ~isempty(xk) )
    for k = 1:length(doms)-1
        ek = chebfun(@(x) err_handle(x), [doms(k), doms(k+1)], 200, 'eps', 1e-9);     
        ek = simplify(ek);
        %length(ek)
    %         if (length(ek) > 4000 )
    %             %plot(ek);
    %             %disp( 'reconstructing' )
    %             ek = chebfun(@(x) err_handle(x), [doms(k), doms(k+1)], 'eps', 1e-9 ); 
    %             %length(ek)    
    %             %plot(ek)
    %             %drawnow
    %             %pause()
    %         end        

        rts = [rts; roots(diff(ek), 'nobreaks')];  %#ok<AGROW>
    end    
end

% Append end points of the domain:
rts = [f.domain(1) ; rts; f.domain(end)];

end

