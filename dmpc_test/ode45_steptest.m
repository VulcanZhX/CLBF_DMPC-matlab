clear;
u_discrete = [2 0.9 2 1.1 2.8 1];

syms t;
uf(t) = 0*t;
for i = 1:length(u_discrete)
    uf(t) = uf(t) + u_discrete(i)*(heaviside(t-(i-1)*0.1) - heaviside(t-i*0.1));
end

uf2 = matlabFunction(uf);
f = @(x, u) -x + u;
[t, x_sol]=ode45(@(t, x)f(x, uf2(t)) , [0 length(u_discrete)*0.1], 2);
plot(t, x_sol);
figure;
t = 0:0.001:length(u_discrete)*0.1;
plot(t, uf2(t));