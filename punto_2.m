

f_Rosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

gradf_Rosenbrock = @(x) [400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2 ; 200*( x(2) - x(1)^2 )];

Hessianf_Rosenbrock = @(x) [ 1200*x(1)^2 - 400*x(2) + 2 , -400*x(1); -400*x(1) , 200];

rho = 0.5;

c1 = 1e-4;

x_0_1 = [-1.2 ; 1];

kmax = 100;

btmax = 150;

tolgrad = 1e-5;

[xk, fk, gradfk_norm, k, xseq ,btseq] = Modified_Newton_method(x_0_1, f_Rosenbrock, gradf_Rosenbrock, Hessianf_Rosenbrock, kmax, tolgrad ,c1, rho, btmax)