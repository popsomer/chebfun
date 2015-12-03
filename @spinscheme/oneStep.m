function usol = oneStep(usol, Nc, Nv, A, B, U, V, E)
%ONESTEP   Compute solution at tn+1 from solution at previous times.
%   USOL = ONESTEP(USOL, NC, NV, A, B, U, V, E) updates the solution USOL using 
%   the coefficients A, B, U, V, and E, and the nonlinear part of the operator 
%   in coefficient space NC and in value space NV.

% Author: Hadrien Montanelli.

s = size(A, 1); 
q = size(U, 2) + 1;
N = cell(s, 1);
c = cell(s, 1);

% Compute internal stages:
c{1} = usol(:,1);
N{1} = Nc.*fft(Nv(ifft(c{1})));
for i = 2:s
    c{i} = E{i}.*c{1};
    for j = 1:i-1
        if ( ~isempty(A{i,j}) )
            c{i} = c{i} + A{i,j}.*N{j};
        end
    end
    for j = 1:q-1
        if ( ~isempty(U{i,j}) )
            c{i} = c{i} + U{i,j}.*(Nc.*fft(Nv(ifft(usol(:,j+1)))));
        end
    end
    N{i} = Nc.*fft(Nv(ifft(c{i})));
end

% Compute solution at tn:
c{1} = E{s+1}.*c{1};
for i = 1:s
    if ( ~isempty(B{i}) )
        c{1} = c{1} + B{i}.*N{i};
    end
end
for i = 1:q-1
    if ( ~isempty(V{i}) )
        c{1} = c{1} + V{i}.*(Nc.*fft(Nv(ifft(usol(:,i+1)))));
    end
end

% Update:
if ( q == 1 )
    usol = c{1};
else
    usol = circshift(usol, 1, 2);
    usol(:,1) = c{1};
end
 
end