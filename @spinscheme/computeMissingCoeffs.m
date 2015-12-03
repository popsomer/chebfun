function [A, B, U, V, E] = computeMissingCoeffs(A, B, C, U, V, L, phi1, ...
    phit, dt)
%COMPUTEMISSINGCOEFFS   Compute the missing oefficients using the summation 
% properties of the coefficients.
%
%   [A, B, U, V, E] = COMPUTEMISSINGCOEFFS(A, B, C, U, V, L, PHI1, PHIT, DT) 
%   uses the row summing properties to computes the AI1 and B1 coefficients, 
%   also computes the E quantities, using the linear part L of the opeartor, 
%   the PHI1 and PHIT functions, and the timestep DT.

% Author: Hadrien Montanelli.

% Number of internal stages S and number of steps used Q:
s = size(A,1); 
q = size(U,2) + 1;

% Precompute the coefficients Ai1 using the row summing property.
for i = 2:s
    A{i,1} = phit{1,i};
    for j = 2:i-1
        if ( ~isempty(A{i,j}) )
            A{i,1} = A{i,1} - A{i,j};
        end
    end
    for j = 1:q-1
        if ( ~isempty(U{i,j}) )
            A{i,1} = A{i,1} - U{i,j};
        end
    end
end

% Precompute the coefficient B1 using the row summing property.
B{1} = phi1;
for i = 2:s
    if ( ~isempty(B{i}) )
        B{1} = B{1} - B{i};
    end
end
for i = 1:q-1
    if ( ~isempty(V{i}) )
        B{1} = B{1} - V{i};
    end
end

% Precompute the E quantities.
E = cell(s+1, 1);
for i = 1:s
   E{i} = exp(C(i)*dt*L);
end
E{s+1} = exp(dt*L);

% Multiply by time-step:
A = cellfun(@(A) dt*A, A, 'UniformOutput', 0);
B = cellfun(@(B) dt*B, B, 'UniformOutput', 0);
U = cellfun(@(U) dt*U, U, 'UniformOutput', 0);
V = cellfun(@(V) dt*V, V, 'UniformOutput', 0);

end