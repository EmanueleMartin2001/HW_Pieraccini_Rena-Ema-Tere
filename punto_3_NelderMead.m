%%%%%%%% FIRST POINT %%%%%%%%

seed = min([341965, 343316, 284817]);

rng(seed);

%%%%%%%% END FIRST POINT %%%%%%%% 





%%%%%%%% SECOND POINT %%%%%%%%

n = 25; %alternatives 10,25,50

f1 = first_function(n); % Problem 1

f2 = second_function(n); % Problem 27

f3 = third_function(n,10); % Problem 64

%%%%%%%% END SECOND POINT





%%%%%%%% THIRD POINT %%%%%%%%

% creating the starting symplex for f1

x_f1 = zeros(n,1);
for i= 1:1:n
    if mod(i,2) == 1
        x_f1(i) = -1.2;
    else
        x_f1(i) = 1.0;
    end
end

% Here, starting from a given point x0,
% we construct a matrix with n+1 rows:
% the first row corresponds to x0 itself,
% while the remaining n rows represent perturbations of x0
% along its n different coordinates.
% Thus we have the starting symplex.

alpha_1 = 1; %explore term (decides how big should the starting symplex be)
X1 = zeros(n+1,n);
X1(1,:) = x_f1;
for i = 2:n+1
    X1(i,:) = x_f1;
    X1(i,i-1) = X1(i,i-1) + alpha_1; % perturbation of 1 coordinate by alpha_1
end



% creating the starting symplex for f2

x_f2 = zeros(n,1);
for l = 1:1:n
    x_f2(l) = l;
end

alpha_2 = 1; 
X2 = zeros(n+1,n);
X2(1,:) = x_f2;
for i = 2:n+1
    X2(i,:) = x_f2;
    X2(i,i-1) = X2(i,i-1) + alpha_2; % perturbation of 1 coordinate by alpha_2
end

%%%%%%%% END THIRD POINT %%%%%%%%



%%%%%%% START FOURTH POINT %%%%%%%%

%% NELDER MEAD METHOD

kmax = 100000;

tol = 1e-3;

rho = 1;

chi = 2;

gamma = 0.5;

sigma = 0.5;

% the contruction of the 10 starting points is embedeed in the following
% for cycle

for i = 1:1:10

    x_f1 = x_f1 + 2*rand(n,1)-1; % x_f1 perturbed in the hypercube

    alpha_1 = 1; %explore term (decides how big should the starting symplex be)
    X1 = zeros(n+1,n);
    X1(1,:) = x_f1;
    for j = 2:n+1
        X1(j,:) = x_f1;
        X1(j,j-1) = X1(j,j-1) + alpha_1; % perturbation of 1 coordinate by alpha_1
    end

    disp(['**** NELDER MEAD METHOD FOR THE FIRST FUNCTION, STARTING SYMPLEX ', num2str(i), ': STARTED *****']);
    tic;
    [x1_opt, fx1_opt, k1, x1seq, f1_seq, singular1] = Nelder_Mead(X1, f1, kmax, tol, rho, chi, gamma, sigma);
    t = toc;

    disp(['**** NELDER MEAD METHOD FOR THE FIRST FUNCTION, STARTING SYMPLEX ', num2str(i), ': FINISHED *****']);

    disp(['Time: ', num2str(t), ' seconds']);

    disp('**** NELDER MEAD METHOD : RESULTS *****')
    disp('************************************')
    disp(['rho: ',num2str(rho), '  chi: ', num2str(chi), '  gamma: ',num2str(gamma),'  sigma: ', num2str(sigma)])
    disp(['fx1_opt: ', num2str(fx1_opt)])
    disp(['N. of Iterations: ', num2str(k1),'/',num2str(kmax), ';'])
    disp('************************************')

    if k1 == kmax || singular1 == 1
        result_first_function(i) = 0;
        disp('FAIL')
        disp('************************************')
    else
        result_first_function(i) = 1;
        disp('SUCCESS')
        disp('************************************')
    end
    disp(' ')

end

figure; 
semilogy(1:k1, f1_seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
grid on;
xlabel('Iterations (k)');
ylabel('Values of the Rosenbrock function'); 
result_first_function