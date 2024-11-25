



%%%%%%%% FIRST POINT %%%%%%%%

seed = min([341965, 343316, 284817]);

rng(seed);

%%%%%%%% END FIRST POIN %%%%%%%% 





%%%%%%%% SECOND POINT %%%%%%%%

d = 3;    % alternative: 3,4,5

n = 10^d

%% for the first function we choose the Chainedd Rosenbrock function

% function construction:

f1 = @(x) sum(100 * (x(1:end-1).^2 -x(2:end)).^2 + (x(1:end-1)-1).^2);

% gradient construction:

f1_grad_comp = cell(n,1);

f1_grad_comp{1} = @(x) 400*x(1)^3 -400*x(1)*x(2) + 2*x(1) - 2;
for i = 2:1:(n-1)
    f1_grad_comp{i} = @(x) 400*x(i)^3 -200*x(i-1)^2 - 400*x(i+1)*x(i) + 202*x(i) -2;
end
f1_grad_comp{n} = @(x) -200*(x(n-1)^2 - x(n));

gradf1 = @(x) cell2mat(cellfun(@(f1_grad_comp) f1_grad_comp(x), f1_grad_comp, 'UniformOutput', false));

% computing the Hessian in sparse mode

A1 = cell(1, n);
B1 = cell(1, n);

A1{1} = @(x) 1200*x(1)^2 -400*x(2) + 2;
B1{1} = @(x) -400*x(1);
for i = 2:1:(n-1) 
    A1{i} = @(x) 1200*x(i)^2 - 400*x(i+1) + 202;
    B1{i} = @(x) -400*x(i);
end
A1{n} = @(x) 200;
B1{n} = @(x) 0; %forced by spdiags

Hessf1 = @(x) spdiags([cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false)), ...
                           cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)), ...
                           cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false))], ...
                          [-1, 0, 1], n, n);

%%%%%%%% END SECOND POINT





%%%%%%%% THIRD POINT %%%%%%%%

% construction of the test point x0 for f1:

x_f1 = zeros(n,1);
for i= 1:1:n
    if mod(i,2) == 1
        x_f1(i) = -1.2;
    else
        x_f1(i) = 1.0;
    end
end

% construction of the 10 points 

X_f1 = (x_f1*(ones(n,1))')';       % matrix that copy for each column the vector x_f1
X_f1 = X_f1(:,1:1:10);             % rescale of the matrix

error = rand(n,10);     % matrix of random variable to add to the starting point 
X_f1 = X_f1 + error;    % each column of this vector represent a starting point

%%%%%%%% END THIRD POINT %%%%%%%%

rho = 0.5;
c = 1e-4;

kmax = 10;
tolgrad = 1e-5;
btmax = 1;

% calling the method:

[xk, fk, gradfk_norm, k, xseq, btseq] = ...
    Modified_Newton_method(X_f1(:,1), f1, gradf1, Hessf1, ...
    kmax, tolgrad, c, rho, btmax)

[xk1, fk1, gradfk_norm1, k1, xseq1, btseq1] = ...
    newton_bcktrck(X_f1(:,1), f1, gradf1, Hessf1, ...
    kmax, tolgrad, c, rho, btmax)

norm(xk-xk1)