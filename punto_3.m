
clc;
clear;
close all;



%%%%%%%% FIRST POINT %%%%%%%%

seed = min([341965, 343316, 284817]);

rng(seed);

%%%%%%%% END FIRST POINT %%%%%%%% 





%%%%%%%% SECOND POINT %%%%%%%%

d = 4;    % alternative: 3,4,5

n = 10^d;

 %[f1,gradf1,Hessf1] = first_function(n); % Problem 1

 [f2,gradf2, Hessf2] = second_function(n); % Problem 31

 %[f3, gradf3, Hessf3] = third_function(n); % Problem 16

%%%%%%%% END SECOND POINT





%%%%%%%% THIRD POINT %%%%%%%%

% construction of the test point x0 for f1:

% x_f1 = zeros(n,1);
% for i= 1:1:n
%     if mod(i,2) == 1
%         x_f1(i) = -1.2;
%     else
%         x_f1(i) = 1.0;
%     end
% end

% construction of the test point for f2

x_f2 = -ones(n,1);

% construction of the test point for f3

% x_f3 = ones(n,1);


% construction of the 10 points for f1

% X_f1 = repmat(x_f1, 1, 10);       % matrix that copy for each column the vector x_f1
% X_f1 = X_f1(:,1:1:10);             % rescale of the matrix
% 
% error = rand(n,10);     % matrix of random variable to add to the starting point 
% X_f1 = X_f1 + error;    % each column of this vector represent a starting point

% construction of the 10 points for f2 

X_f2 = repmat(x_f2, 1, 10);
X_f2 = X_f2(:,1:1:10);

% error = rand(n,10);
% X_f2 = X_f2 + error;

% construction of the 10 points for f3 

% X_f3 = repmat(x_f3, 1, 10);
% X_f3 = X_f3(:,1:1:10);
% 
% error = rand(n,10);
% X_f3 = X_f3 + error;


%%%%%%%% END THIRD POINT %%%%%%%%



%%%%%%% START FOURTH POINT %%%%%%%%

%% MODIFIED NEWTON METHOD

rho = 0.5;
c = 1e-4;

kmax = 1000;
tolgrad = 1e-4;
btmax = 40;
%type_tao = 'Gershgorin';
%type_tao = 'Eigen';
type_tao = 'Cholesky';


% calling the method:

result_first_function = 1000*ones(10,1);
result_second_function = 1000*ones(10,1);
result_third_function = 1000*ones(10,1);

time_1 = zeros(10,1);
time_2 = zeros(10,1);
time_3 = zeros(10,1);

iteration_1 = zeros(10,1);
iteration_2 = zeros(10,1);
iteration_3 = zeros(10,1);

number_tao_1 = zeros(10,1);
number_tao_2 = zeros(10,1);
number_tao_3 = zeros(10,1);

conv_rate_1 = zeros(10,1);
conv_rate_2 = zeros(10,1);
conv_rate_3 = zeros(10,1);

soltol1 = 10^-5;

% for i = 1:1:10
% 
%     disp(['**** MODIFIED NEWTON METHOD FOR THE FIRST FUNCTION, POINT ', num2str(i), ': STARTED *****']);
%     tic;
%     [x1k, f1k, gradf1k_norm, k1, x1seq,f1seq, bt1seq, taoseq1] = ...
%         Modified_Newton_method(X_f1(:,i), f1, gradf1, Hessf1, ...
%         kmax, tolgrad, c, rho, btmax, type_tao);
%     t = toc;
% 
%     disp(['**** MODIFIED NEWTON METHOD FOR THE FIRST FUNCTION, POINT ', num2str(i), ': FINISHED *****']);
% 
%     disp(['Time: ', num2str(t), ' seconds']);
% 
%     disp('**** MODIFIED NEWTON METHOD : RESULTS *****');
%     disp('************************************');
%     disp(['N. tao used: ', num2str(nnz(taoseq1))]);
%     disp(['f(xk): ', num2str(f1k)]);
%     disp(['gradfk_norm: ', num2str(gradf1k_norm)]);
% 
%     disp(['N. of Iterations: ', num2str(k1),'/',num2str(kmax), ';']);
%     disp(['Rate of convergence: ', num2str(convergence_rate(x1seq)), ';']);
%     disp('************************************');
% 
%     if k1 == kmax
%         result_first_function(i) = 0;
%         disp('FAIL')
%         disp('************************************')
%     else
%         if (norm(x1k-ones(n,1)) < soltol1)
%             result_first_function(i) = 1;
%             disp('SUCCESS')
%             disp('************************************')
%             success_k1 = k1;
%             success_f1seq = f1seq;
%             success_taoseq1 = taoseq1;
%         else
%             result_first_function(i) = 0;
%             disp('FAIL')
%             disp('************************************')
%         end
%     end
%     disp(' ')
%    time_1(i) = t;
%    iteration_1(i) = k1;
%    number_tao_1(i) = nnz(taoseq1);
%    conv_rate_1(i) = convergence_rate(x1seq);
% 
% end
% 
% figure; 
% semilogy(1:success_k1, success_f1seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
% grid on;
% xlabel('Iterations (k)');
% ylabel('Values of the Rosenbrock function'); 
% 
% figure;
% bar(1:success_k1, success_taoseq1, 'FaceColor', 'blue', 'EdgeColor', 'black')
% grid on;
% xlabel('Iterations (k)');
% ylabel('Tao values for the Rosenbrock function'); 
% 
% disp(' ')
% disp(' ')
% disp(' ')
% disp('******************************************')
% 
% disp('**** RESULTS FOR THE FIRST FUNCTION *****')
% disp(['N. of success: ', num2str(sum(result_first_function))])
% disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_first_function.*iteration_1)/sum(result_first_function))),'/',num2str(kmax), ';'])
% disp(['Mean N. tao used (in case of success): ', num2str(round(sum(result_first_function.*number_tao_1)/sum(result_first_function)))])
% disp(['Mean convergence rate (in case of success): ', num2str(sum(result_first_function.*conv_rate_1)/sum(result_first_function))])
% disp('******************************************')


for i = 1:1:10

    disp(['**** MODIFIED NEWTON METHOD FOR THE SECOND FUNCTION, POINT ', num2str(i), ': STARTED *****']);
    tic;
    [x2k, f2k, gradf2k_norm, k2, x2seq, f2seq, b2tseq, taoseq2] = ...
        Modified_Newton_method(X_f2(:,i), f2, gradf2, Hessf2, ...
        kmax, tolgrad, c, rho, btmax, type_tao);
    t = toc;

    disp(['**** MODIFIED NEWTON METHOD FOR THE SECOND FUNCTION, POINT ', num2str(i), ': FINISHED *****']);

    disp(['Time: ', num2str(t), ' seconds']);

    disp('**** MODIFIED NEWTON METHOD : RESULTS *****')
    disp('************************************')
    disp(['N. tao used: ', num2str(nnz(taoseq2))])
    disp(['f(xk): ', num2str(f2k)])
    disp(['gradfk_norm: ', num2str(gradf2k_norm)])

    disp(['N. of Iterations: ', num2str(k2),'/',num2str(kmax), ';'])
    disp(['Rate of convergence: ', num2str(convergence_rate(x2seq)), ';'])
    disp('************************************')

    if (k2 == kmax || f2k > 10^-1)
        result_second_function(i) = 0;
        disp('FAIL')
        disp('************************************')
    else
        result_second_function(i) = 1;
        disp('SUCCESS')
        disp('************************************')
        success_k2 = k2;
        success_f2seq = f2seq;
        success_taoseq2 = taoseq2;
    end
    disp(' ')
   time_2(i) = t;
   iteration_2(i) = k2;
   number_tao_2(i) = nnz(taoseq2);
   conv_rate_2(i) = convergence_rate(x2seq);

end

figure; 
semilogy(1:success_k2, success_f2seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
grid on;
xlabel('Iterations (k)');
ylabel('Values of the Broyden tridiagonal function'); 

figure;
bar(1:success_k2, success_taoseq2, 'FaceColor', 'blue', 'EdgeColor', 'black')
grid on;
xlabel('Iterations (k)');
ylabel('Tao values for the Broyden tridiagonal function'); 

disp(' ')
disp(' ')
disp(' ')
disp('******************************************')

disp('**** RESULTS FOR THE SECOND FUNCTION *****')
disp(['N. of success: ', num2str(sum(result_second_function))])
disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_second_function.*iteration_2)/sum(result_second_function))),'/',num2str(kmax), ';'])
disp(['Mean N. tao used (in case of success): ', num2str(round(sum(result_second_function.*number_tao_2)/sum(result_second_function)))])
disp(['Mean convergence rate (in case of success): ', num2str(sum(result_second_function.*conv_rate_2)/sum(result_second_function))])
disp('******************************************')

% for i = 1:1:10
% 
%     disp(['**** MODIFIED NEWTON METHOD FOR THE THIRD FUNCTION, POINT ', num2str(i), ': STARTED *****']);
%     tic;
%     [x3k, f3k, gradf3k_norm, k3, x3seq, f3seq, b3tseq, taoseq3] = ...
%         Modified_Newton_method(X_f3(:,i), f3, gradf3, Hessf3, ...
%         kmax, tolgrad, c, rho, btmax, type_tao);
%     t = toc;
% 
%     disp(['**** MODIFIED NEWTON METHOD FOR THE THIRD FUNCTION, POINT ', num2str(i), ': FINISHED *****']);
% 
%     disp(['Time: ', num2str(t), ' seconds']);
% 
%     disp('**** MODIFIED NEWTON METHOD : RESULTS *****')
%     disp('************************************')
%     disp(['N. tao used: ', num2str(nnz(taoseq3))])
%     disp(['f(xk): ', num2str(f3k)])
%     disp(['gradfk_norm: ', num2str(gradf3k_norm)])
%     disp(['N. of Iterations: ', num2str(k3),'/',num2str(kmax), ';'])
%     disp(['Rate of convergence: ', num2str(convergence_rate(x3seq)), ';'])
%     disp('************************************')
% 
%     if k3 == kmax
%         result_third_function(i) = 0;
%         disp('FAIL')
%         disp('************************************')
%     else
%         result_third_function(i) = 1;
%         disp('SUCCESS')
%         disp('************************************')
%     end
%     disp(' ')
% 
%     time_3(i) = t;
%     iteration_3(i) = k3;
%     number_tao_3(i) = nnz(taoseq3);
%     conv_rate_3(i) = convergence_rate(x3seq);
% 
% end
% 
% figure; 
% plot(1:k3, f3seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
% grid on;
% xlabel('Iterations (k)');
% ylabel('Values for the Banded trigonometric problem'); 
% 
% figure;
% hold on;
% bar(1:k3, taoseq3 .* (taoseq3 >= 0), 'FaceColor', 'blue', 'EdgeColor', 'black'); % positive value
% bar(1:k3, taoseq3 .* (taoseq3 < 0), 'FaceColor', 'red', 'EdgeColor', 'black'); % negative
% hold off
% grid on;
% xlabel('Iterations (k)');
% ylabel('Tao values for the Banded trigonometric problem'); 
% 
% disp(' ')
% disp(' ')
% disp(' ')
% disp('******************************************')
% 
% disp('**** RESULTS FOR THE THIRD FUNCTION *****')
% disp(['N. of success: ', num2str(sum(result_third_function))])
% disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_third_function.*iteration_3)/sum(result_third_function))),'/',num2str(kmax), ';'])
% disp(['Mean N. tao used (in case of success): ', num2str(round(sum(result_third_function.*number_tao_3)/sum(result_third_function)))])
% disp(['Mean convergence rate (in case of success): ', num2str(sum(result_third_function.*conv_rate_3)/sum(result_third_function))])
% disp('******************************************')
