function f = computePhi(l)
%COMPUTEPHI   Compute the phi functions.
%   f = COMPUTEPHI(L) returns a function handle to the phi function of index L.

% Author: Hadrien Montanelli.

% Trivial case:
if ( l == 0 )
   f = @(z) exp(z); 
   
% Compute them recursively using the recurrence formula 
% phi_{j}(z) = 1/z*(phi_{j-1}(z) - 1/(j-1)!).
else
   ff = computePhi(l-1);
   f = @(z) (ff(z) - 1/factorial(l-1))./z; 
end

end