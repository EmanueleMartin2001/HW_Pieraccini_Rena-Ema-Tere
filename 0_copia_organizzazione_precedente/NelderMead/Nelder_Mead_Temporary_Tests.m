%% test 1
% f = @(x) x(1)+x(2)^2
% X = [2,3; 4,15; 5,10]
% tol = 10e-10
% kmax = 100
%it goes towards -infinity as we expect

%% test 2 Rosenbrock function

% starting simplex: X = [1.2,1.2; -1.2,1; 0,1]
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
X = [1.2,1.2; -1.2,1; 0,1];
% kmax = 200;
% tol = 10^(-6);
% rho = 1;
% chi = 2;
% gamma = 0.5;
% sigma = 0.5;
[x_opt, fx_opt, k, xseq, singular] = Nelder_Mead(X, f)


% % starting simplex: X = [8,14; -100,1; 0,16] (higher distance from x_opt)
% f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
% X = [8,14; -100,1; 0,16];
% kmax = 500;
% tol = 10^(-12);
% rho = 1;
% chi = 2;
% gamma = 0.5;
% sigma = 0.5;
% [x_opt, fx_opt, k, x_seq, singular] = Nelder_Mead(X, f, kmax, tol, rho, chi, gamma, sigma)