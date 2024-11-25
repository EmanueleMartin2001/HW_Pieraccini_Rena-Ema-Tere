



%%%%%%%% FIRST POINT %%%%%%%%

seed = min([341965, 343316, 284817]);

rng(seed);

%%%%%%%% END FIRST POIN %%%%%%%% 



%%%%%%%% SECOND POINT %%%%%%%%

d = 1;    % alternative: 3,4,5

n = 10^d;

%% for the first function we choose the Chainedd Rosenbrock function

% function construction:

f1 = @(x) sum(100 * (x(1:end-1).^2 -x(2:end)).^2 + (x(1:end-1)-1).^2);

% gradient construction:

f1_grad_comp = cell(1,n);

f1_grad_comp{1} = @(x) 400*x(1)^3 -400*x(1)*x(2) + 2*x(1) - 1;
for i = 2:1:(n-1)
    f1_grad_comp{i} = @(x) 400*x(i)^3 -200*x(i-1)^2 - 400*x(i+1)*x(i) + 202*x(i);
end
f1_grad_comp{n} = @(x) -200*(x(n-1)^2 - x(n));

gradf1 = @(x) cell2mat(cellfun(@(f1_grad_comp) f1_grad_comp(x), f1_grad_comp, 'UniformOutput', false)); 

% hessian construction

% from calcolus i can see that the hessian is almost diagonal 

A1 = cell(1, n);
B1 = cell(1, n-1);

A1{1} = @(x) 1200*x(1)^2 -400*x(2) + 2;
B1{1} = @(x) -400*x(1);
for i = 2:1:(n-1) 
    A1{i} = @(x) 1200*x(i)^2 - 400*x(i+1) + 202;
    B1{i} = @(x) -400*x(i);
end
A1{n} = @(x) 200;

Hessian_f1 = @(x) diag(cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false))) + diag(cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false)),-1) + diag(cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false)),+1); 



%%%%%%%% END SECOND POINT





%%%%%%%% THIRD POINT %%%%%%%%

