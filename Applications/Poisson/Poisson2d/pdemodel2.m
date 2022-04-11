function pde = pdemodel2
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.uhat = @uhat;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
s = (2*pi*pi)*sin(pi*x1)*sin(pi*x2);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(uhat, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0); 
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end

function uhat = uhat(u, q, w, v, x, t, mu, eta, uhat, n, tau, u2, q2, w2, v2)
f1 = mu*(q(1)*n(1) + q(2)*n(2));
f2 = mu*(q2(1)*n(1) + q2(2)*n(2));
tau1 = 1/tau;
uhat = 0.5*(u+u2) + 0.5*tau1*(f1-f2);
end