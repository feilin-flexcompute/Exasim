
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu(1)*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
L = mu(2);
s = (2*pi*pi)*sin(pi*x1/L)*sin(pi*x2/L)/(L*L);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0); 
end

function u0 = initu(x, mu, eta)
% x1 = x(1);
% x2 = x(2);
% L = mu(2);
%u0 = sin(pi*x1/L)*sin(pi*x2/L);
u0 = sym(0.0);
end

