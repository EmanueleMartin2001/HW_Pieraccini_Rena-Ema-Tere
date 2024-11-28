



%%%%%%%% FIRST POINT %%%%%%%%

seed = min([341965, 343316, 284817]);

rng(seed);

%%%%%%%% END FIRST POIN %%%%%%%% 





%%%%%%%% SECOND POINT %%%%%%%%

d = 3;    % alternative: 3,4,5

n = 10^d;

%% for the first function we choose the Chainedd Rosenbrock function

[f1,gradf1,Hessf1] = first_function(n);

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

kmax = 1000;
tolgrad = 1e-6;
btmax = 40;


% calling the method:


[xk, fk, gradfk_norm, k, xseq, btseq] = ...
    Modified_Newton_method(X_f1(:,2), f1, gradf1, Hessf1, ...
    kmax, tolgrad, c, rho, btmax);

[xk1, fk1, gradfk_norm1, k1, xseq1, btseq1] = ...
    newton_bcktrck(X_f1(:,2), f1, gradf1, Hessf1, ...
    kmax, tolgrad, c, rho, btmax);

norm(xk-xk1)

