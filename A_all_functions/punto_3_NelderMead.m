%%%%%%%% FIRST POINT %%%%%%%%

seed = min([341965, 343316, 284817]);

rng(seed);

%%%%%%%% END FIRST POINT %%%%%%%% 





%%%%%%%% SECOND POINT %%%%%%%%

n = 50; %alternatives 10,25,50

% f1 = first_function(n); % Problem 1

%f2 = second_function_75(n) % Problem 75

f3 = third_function(n); % Problem 16

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

x_f2 = -ones(n,1)*1.2;
x_f2(end) = -1;

alpha_2 = 1; 
X2 = zeros(n+1,n);
X2(1,:) = x_f2;
for i = 2:n+1
    X2(i,:) = x_f2;
    X2(i,i-1) = X2(i,i-1) + alpha_2; % perturbation of 1 coordinate by alpha_2
end

x_f3 = ones(n,1);


%%%%%%%% END THIRD POINT %%%%%%%%



%%%%%%% START FOURTH POINT %%%%%%%%

%% NELDER MEAD METHOD

kmax = 100000;

tol = 1e-4;

rho = 1;

chi = 2;

gamma = 0.5;

sigma = 0.5;

%%%  FIRST FUNCTION  %%%
% the construction of the 10 starting points is embedeed in the following
% for cycle
% for i = 1:1:10
% 
%     x_f1 = x_f1 + 2*rand(n,1)-1; % x_f1 perturbed in the hypercube
% 
%     alpha_1 = 1; %explore term (decides how big should the starting symplex be)
%     X1 = zeros(n+1,n);
%     X1(1,:) = x_f1;
%     for j = 2:n+1
%         X1(j,:) = x_f1;
%         X1(j,j-1) = X1(j,j-1) + alpha_1; % perturbation of 1 coordinate by alpha_1
%     end
% 
%     disp(['**** NELDER MEAD METHOD FOR THE FIRST FUNCTION, STARTING SYMPLEX ', num2str(i), ': STARTED *****']);
%     tic;
%     [x1_opt, fx1_opt, k1, x1seq, f1_seq, singular1] = Nelder_Mead(X1, f1, kmax, tol, rho, chi, gamma, sigma);
%     t = toc;
% 
%     disp(['**** NELDER MEAD METHOD FOR THE FIRST FUNCTION, STARTING SYMPLEX ', num2str(i), ': FINISHED *****']);
% 
%     disp(['Time: ', num2str(t), ' seconds']);
% 
%     disp('**** NELDER MEAD METHOD : RESULTS *****')
%     disp('************************************')
%     disp(['rho: ',num2str(rho), '  chi: ', num2str(chi), '  gamma: ',num2str(gamma),'  sigma: ', num2str(sigma)])
%     disp(['fx1_opt: ', num2str(fx1_opt)])
%     disp(['N. of Iterations: ', num2str(k1),'/',num2str(kmax), ';'])
%     disp('************************************')
% 
%     if k1 == kmax || singular1 == 1 || abs(fx1_opt) > tol
%         result_first_function(i) = 0;
%         disp('FAIL')
%         disp('************************************')
%     else
%         result_first_function(i) = 1;
%         disp('SUCCESS')
%         disp('************************************')
%     end
%     disp(' ')
% 
% end
% 
% figure; 
% semilogy(1:k1, f1_seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
% grid on;
% xlabel('Iterations (k)');
% ylabel('Values of the Rosenbrock function'); 
% result_first_function = result_first_function';
% result_first_function


%%%  SECOND FUNCTION  %%%
% % the construction of the 10 starting points is embedeed in the following
% % for cycle
% for i = 1:1:10
%     x_f2 = x_f2 + 2*rand(n,1)-1; % x_f2 perturbed in the hypercube
%     alpha_2 = 1; %explore term (decides how big should the starting symplex be)
%     X2 = zeros(n+1,n);
%     X2(1,:) = x_f2;
%     for j = 2:n+1
%         X2(j,:) = x_f2;
%         X2(j,j-1) = X2(j,j-1) + alpha_2; % perturbation of 1 coordinate by alpha_2
%     end
% 
%     disp(['**** NELDER MEAD METHOD FOR THE SECOND FUNCTION, STARTING SYMPLEX ', num2str(i), ': STARTED *****']);
%     tic;
%     [x2_opt, fx2_opt, k2, x2seq, f2_seq, singular2] = Nelder_Mead(X2, f2, kmax, tol, rho, chi, gamma, sigma);
%     t = toc;
% 
%     disp(['**** NELDER MEAD METHOD FOR THE SECOND FUNCTION, STARTING SYMPLEX ', num2str(i), ': FINISHED *****']);
% 
%     disp(['Time: ', num2str(t), ' seconds']);
% 
%     disp('**** NELDER MEAD METHOD : RESULTS *****')
%     disp('************************************')
%     disp(['rho: ',num2str(rho), '  chi: ', num2str(chi), '  gamma: ',num2str(gamma),'  sigma: ', num2str(sigma)])
%     disp(['fx2_opt: ', num2str(fx2_opt)])
%     disp(['N. of Iterations: ', num2str(k2),'/',num2str(kmax), ';'])
%     disp('************************************')
% 
%     if k2 == kmax || singular2 == 1 || abs(fx2_opt) > tol
%         result_second_function(i) = 0;
%         disp('FAIL')
%         disp('************************************')
%     else
%         result_second_function(i) = 1;
%         disp('SUCCESS')
%         disp('************************************')
%     end
%     disp(' ')
% 
% end
% 
% figure; 
% semilogy(1:k2, f2_seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
% hold on; 
% yline(1, '-', 'Color', [1, 0, 0], 'LineWidth', 2); 
% grid on;
% xlabel('Iterations (k)');
% ylabel('Values of the function 75'); 
% hold off;
% result_second_function'


%%%  THIRD FUNCTION  %%%
% the construction of the 10 starting points is embedeed in the following
% for cycle
for i = 1:1:10

    x_f3 = x_f3 + 2*rand(n,1)-1; % x_f3 perturbed in the hypercube

    alpha_3 = 1; %explore term (decides how big should the starting symplex be)
    X3 = zeros(n+1,n);
    X3(1,:) = x_f3;
    for j = 2:n+1
        X3(j,:) = x_f3;
        X3(j,j-1) = X3(j,j-1) + alpha_3; % perturbation of 1 coordinate by alpha_3
    end

    disp(['**** NELDER MEAD METHOD FOR THE THIRD FUNCTION, STARTING SYMPLEX ', num2str(i), ': STARTED *****']);
    tic;
    [x3_opt, fx3_opt, k3, x3seq, f3_seq, singular3] = Nelder_Mead(X3, f3, kmax, tol, rho, chi, gamma, sigma);
    t = toc;

    disp(['**** NELDER MEAD METHOD FOR THE THIRD FUNCTION, STARTING SYMPLEX ', num2str(i), ': FINISHED *****']);

    disp(['Time: ', num2str(t), ' seconds']);

    disp('**** NELDER MEAD METHOD : RESULTS *****')
    disp('************************************')
    disp(['rho: ',num2str(rho), '  chi: ', num2str(chi), '  gamma: ',num2str(gamma),'  sigma: ', num2str(sigma)])
    disp(['fx3_opt: ', num2str(fx3_opt)])
    disp(['N. of Iterations: ', num2str(k3),'/',num2str(kmax), ';'])
    disp('************************************')
    % n = 10 => f(x3_opt) = -8.0514
    % n = 25 => f(x3_opt) = -16.1379
    % n = 50 => f(x3_opt) = -27.8949
    if n == 50
        if k3 == kmax || singular3 == 1 || (fx3_opt + 8.0514) > tol*10
            result_third_function(i) = 0;
            disp('FAIL')
            disp('************************************')
        else
            result_third_function(i) = 1;
            disp('SUCCESS')
            disp('************************************')
        end
        disp(' ')
    end
    if n == 25
        if k3 == kmax || singular3 == 1 || (fx3_opt + 16.1379) > tol*150
            result_third_function(i) = 0;
            disp('FAIL')
            disp('************************************')
        else
            result_third_function(i) = 1;
            disp('SUCCESS')
            disp('************************************')
        end
        disp(' ')
    end
    if n == 50
        if k3 == kmax || singular3 == 1 || (fx3_opt + 27.8949) > tol*300
            result_third_function(i) = 0;
            disp('FAIL')
            disp('************************************')
        else
            result_third_function(i) = 1;
            disp('SUCCESS')
            disp('************************************')
        end
        disp(' ')
    end
end

figure; 
semilogy(1:k3, f3_seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
grid on;
xlabel('Iterations (k)');
ylabel('Values of the Banded trigonometric problem'); 
result_third_function'

