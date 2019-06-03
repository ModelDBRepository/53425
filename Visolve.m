function Vin = Visolve(Vext,rm,ra,dx);

% solve linear system of equations for intracelular potential

N = length(Vext);

A = sparse([],[],[],N,N,3*N);
A = A + spdiags(repmat(1 + 2*(rm/dx)/(ra*dx),[N 1]),0,N,N) - spdiags(repmat((rm/dx)/(ra*dx),[N 1]),1,N,N) - spdiags(repmat((rm/dx)/(ra*dx),[N-1 1]),-1,N,N);

Vin = [A\Vext']';                        % intracellular potential (mV)

return