% clc;
% clear;
% close all;

%% HIGH LEVEL PARAMETERS INITIALIZATION

seed = min([341965, 343316, 284817]);

rng(seed);


h = 10^(-12); % alternative 2,4,6,8,10,12 

esponenti = 2:2:12;
H = 10.^(-esponenti); %vettore con gli esponenti da 2 a 12 pari

dimensioni = [10^3,10^4,10^5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS SETTINGS 

rho = 0.6; %backtracking parameters
c = 1e-4;
btmax = 60;

kmax = 1000; %stopping conditions
tolgrad = 1e-3; 
delta_step = 1e-5; 

% type_tao = 'Gershgorin'; % methods used to make the hessian def pos
% type_tao = 'Eigen';      
type_tao = 'Cholesky'; %usually works better


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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FIRST FUNCTION (1)
% 
% %statistics for the table plot
% 
% avg_f_best = [];
% avg_norm_grad = [];
% avg_n_iter = [];
% avg_conv_time = [];
% n_success = [];
% avg_roc = [];
% avg_tao_used = [];
% 
% %the simulation
% for l = 1:1:3
% 
%     n = dimensioni(l);
% 
%     [f1,gradf1,Hessf1] = first_function_1(n, "simplified_forward", h, false); % Problem 1
% 
%     % construction of the test point x0 for f1:
% 
%     x_f1 = zeros(n,1);
%     for i= 1:1:n
%         if mod(i,2) == 1
%             x_f1(i) = -1.2;
%         else
%             x_f1(i) = 1.0;
%         end
%     end
% 
%     % construction of the other 10 points for f1
% 
%     X_f1 = repmat(x_f1, 1, 10);       % matrix that copy for each column the vector x_f1
%     X_f1 = X_f1(:,1:1:10);             % rescale of the matrix
% 
%     error = rand(n,10);     % matrix of random variable to add to the starting point 
%     X_f1 = X_f1 + error;    % each column of this vector represent a starting point
%     X_f1 = [x_f1, X_f1];    %adding the point without noise
% 
%     for i = 1:1:11
% 
%         disp(['**** MODIFIED NEWTON METHOD FOR THE FIRST FUNCTION, POINT ', num2str(i), ': STARTED *****']);
%         tic;
%         [x1k, f1k, gradf1k_norm, k1, x1seq,f1seq, bt1seq, taoseq1, gradf_k1, cos_grad1, fail1] = ...
%             Modified_Newton_method(X_f1(:,i), f1, gradf1, Hessf1, ...
%             kmax, tolgrad, delta_step , c, rho, btmax, type_tao);
%         t = toc;
%         disp(['**** MODIFIED NEWTON METHOD FOR THE FIRST FUNCTION, POINT ', num2str(i), ': FINISHED *****']);
% 
%         disp(['Time: ', num2str(t), ' seconds']);
% 
%         disp('**** MODIFIED NEWTON METHOD : RESULTS *****');
%         disp('************************************');
%         disp(['N. tao used: ', num2str(nnz(taoseq1))]);
%         disp(['f(xk): ', num2str(f1k)]);
%         disp(['gradfk_norm: ', num2str(gradf1k_norm)]);
%         disp(['N. of Iterations: ', num2str(k1),'/',num2str(kmax), ';']);
%         disp(['Rate of convergence: ', num2str(convergence_rate(x1seq)), ';']);
%         disp('************************************');
% 
%         if fail1 == "kmax"
%             result_first_function(i) = 0;
%             disp('FAIL: maximum number of iteration reached')
%             disp('************************************')
%         elseif fail1 == "cosine"
%             result_first_function(i) = 0;
%             disp('FAIL: consecutive directions almost identical -> no improvment')
%             disp('************************************')
%         elseif fail1 == "btmax"
%             result_first_function(i) = 0;
%             disp('FAIL: btmax reached -> could not satisfy the Armijo condition')
%             disp('************************************')
%         else %fail1 == "success"
%             result_first_function(i) = 1;
%             disp('SUCCESS')
%             disp('************************************')
%         end
%         disp(' ')
% 
%         time_1(i) = t;
%         iteration_1(i) = k1;
%         number_tao_1(i) = nnz(taoseq1);
%         conv_rate_1(i) = convergence_rate(x1seq);
%         fbest_iter(i) = f1k; %to save it in the table
%         conv_time_iter(i) = t; %to save it in the table
%         norm_grad_iter(i) = gradf1k_norm; %to save it in the table
%     end
%     % if sum(result_first_function) > 0
%     %     figure; 
%     %     semilogy(1:success_k1, success_f1seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
%     %     grid on;
%     %     xlabel('Iterations (k)');
%     %     ylabel('Values of the Rosenbrock function'); 
%     % 
%     %     figure;
%     %     bar(1:success_k1, success_taoseq1, 'FaceColor', 'blue', 'EdgeColor', 'black')
%     %     grid on;
%     %     xlabel('Iterations (k)');
%     %     ylabel('Tao values for the Rosenbrock function'); 
%     % 
%     % end
%     disp(' ')
%     disp(' ')
%     disp(' ')
%     disp('******************************************')
% 
%     disp('**** RESULTS FOR THE FIRST FUNCTION *****')
%     disp(['N. of success: ', num2str(sum(result_first_function)),'/11'])
% 
%     valid_idx = ~isnan(conv_rate_1); %ignore when the algorithm does not converge
%     avg_roc(l) = sum(result_first_function(valid_idx) .* conv_rate_1(valid_idx)) / sum(result_first_function(valid_idx));
%     disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_first_function.*iteration_1)/sum(result_first_function))),'/',num2str(kmax), ';'])
%     disp(['Mean N. tao used (in case of success): ', num2str(round(sum(result_first_function.*number_tao_1)/sum(result_first_function)))])
%     disp(['Mean convergence rate (in case of success): ', num2str(avg_roc(l))])
%     disp('******************************************')
% 
%     avg_f_best(l) = mean(fbest_iter);
%     avg_norm_grad(l) = mean(norm_grad_iter);
%     avg_n_iter(l) = round(mean(iteration_1));
%     avg_conv_time(l) = mean(conv_time_iter);
%     n_success(l) = sum(result_first_function);
%     avg_roc(l) = sum(result_first_function.*conv_rate_1)/sum(result_first_function);
%     avg_tao_used(l) = round(sum(result_first_function.*number_tao_1)/sum(result_first_function));
% 
% end
% 
% %table plot
% clmn_labels = ["avg f(best)", "avg norm grad", "avg n iter", "avg time (s)", "n of success/11", "avg roc (in case of success)", "n tao used (in case of success)"];
% row_names = ["1000", "10000", "100000"];
% T1 = table(avg_f_best', avg_norm_grad', avg_n_iter', avg_conv_time', n_success', ...
%     avg_roc', avg_tao_used', 'VariableNames', clmn_labels, 'RowNames',row_names);
% disp(T1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SECOND FUNCTION (75)
% 
% %statistics for the table plot
% 
% avg_f_best = [];
% avg_norm_grad = [];
% avg_n_iter = [];
% avg_conv_time = [];
% n_success = [];
% avg_roc = [];
% avg_tao_used = [];
% 
% %the simulation
% for l = 1:1:3
% 
%     n = dimensioni(l);
% 
%     [f2, gradf2, Hessf2] = second_function_75(n,"simplified_centered", h, false); % Problem 75
% 
%     % construction of the test point for f2
% 
%     x_f2 = -ones(n,1);
% 
%     % construction of the test point for f2_75
% 
%     x_f2 = -ones(n,1)*1.2;
%     x_f2(end) = -1;
% 
%     % construction of the 10 points for f2 
% 
%     X_f2 = repmat(x_f2, 1, 10);
%     X_f2 = X_f2(:,1:1:10);
% 
%     error = rand(n,10);
%     X_f2 = X_f2 + error;
%     X_f2 = [x_f2, X_f2];
% 
%     fbest_iter = []; %to save it in the table
%     conv_time_iter = []; %to save it in the table
%     norm_grad_iter = []; %to save it in the table
%     for i = 1:1:11 %function 75
% 
%         disp(['**** MODIFIED NEWTON METHOD FOR THE SECOND FUNCTION, POINT ', num2str(i), ': STARTED *****']);
%         tic;
%         [x2k, f2k, gradf2k_norm, k2, x2seq, f2seq, b2tseq, taoseq2, gradf_k2, cos_grad2, fail2] = ...
%             Modified_Newton_method(X_f2(:,i), f2, gradf2, Hessf2, ...
%             kmax, tolgrad, delta_step, c, rho, btmax, type_tao);
%         t = toc;
% 
%         disp(['**** MODIFIED NEWTON METHOD FOR THE SECOND FUNCTION, POINT ', num2str(i), ': FINISHED *****']);
% 
%         disp(['Time: ', num2str(t), ' seconds']);
% 
%         disp('**** MODIFIED NEWTON METHOD : RESULTS *****')
%         disp('************************************')
%         disp(['N. tao used: ', num2str(nnz(taoseq2))])
%         disp(['f(xk): ', num2str(f2k)])
%         disp(['gradfk_norm: ', num2str(gradf2k_norm)])
%         disp(['N. of Iterations: ', num2str(k2),'/',num2str(kmax), ';'])
%         disp(['Rate of convergence: ', num2str(convergence_rate(x2seq)), ';'])
%         disp('************************************')
% 
%         if fail2 == "kmax"
%             result_second_function(i) = 0;
%             disp('FAIL: maximum number of iteration reached')
%             disp('************************************')
%         elseif fail2 == "cosine"
%             result_second_function(i) = 0;
%             disp('FAIL: consecutive directions almost identical -> no improvment')
%             disp('************************************')
%         elseif fail2 == "btmax"
%             result_second_function(i) = 0;
%             disp('FAIL: btmax reached -> could not satisfy the Armijo condition')
%             disp('************************************')
%         else %fail2 == "success"
%             result_second_function(i) = 1;
%             disp('SUCCESS')
%             disp('************************************')
%         end
%         disp(' ')
% 
%         time_2(i) = t;
%         iteration_2(i) = k2;
%         number_tao_2(i) = nnz(taoseq2);
%         conv_rate_2(i) = convergence_rate(x2seq);
%         fbest_iter(i) = f2k; %to save it in the table
%         conv_time_iter(i) = t; %to save it in the table
%         norm_grad_iter(i) = gradf2k_norm; %to save it in the table  
%     end
% 
%     % figure; 
%     % semilogy(1:success_k2, success_f2seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
%     % grid on;
%     % xlabel('Iterations (k)');
%     % ylabel('Values of the function 75'); 
%     % 
%     % figure;
%     % bar(1:success_k2, success_taoseq2, 'FaceColor', 'blue', 'EdgeColor', 'black')
%     % grid on;
%     % xlabel('Iterations (k)');
%     % ylabel('Tao values for the function 75'); 
% 
%     disp(' ')
%     disp(' ')
%     disp(' ')
%     disp('******************************************')
% 
%     valid_idx = ~isnan(conv_rate_2); %ignore when the algorithm does not converge
%     avg_roc(l) = sum(result_second_function(valid_idx) .* conv_rate_2(valid_idx)) / sum(result_second_function(valid_idx));
%     disp('**** RESULTS FOR THE SECOND FUNCTION *****')
%     disp(['N. of success: ', num2str(sum(result_second_function))])
%     disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_second_function.*iteration_2)/sum(result_second_function))),'/',num2str(kmax), ';'])
%     disp(['Mean N. tao used (in case of success): ', num2str(round(sum(result_second_function.*number_tao_2)/sum(result_second_function)))])
%     disp(['Mean convergence rate (in case of success): ', num2str(avg_roc(l))])
%     disp('******************************************')
% 
%     avg_f_best(l) = mean(fbest_iter);
%     avg_norm_grad(l) = mean(norm_grad_iter);
%     avg_n_iter(l) = round(mean(iteration_2));
%     avg_conv_time(l) = mean(conv_time_iter);
%     n_success(l) = sum(result_second_function);
%     avg_roc(l) = sum(result_second_function.*conv_rate_2)/sum(result_second_function);
%     avg_tao_used(l) = round(sum(result_second_function.*number_tao_2)/sum(result_second_function));
% end
% 
% %table plot
% clmn_labels = ["avg f(best)", "avg norm grad", "avg n iter", "avg time (s)", "n of success/11", "avg roc (in case of success)", "n tao used (in case of success)"];
% row_names = ["1000", "10000", "100000"];
% T2 = table(avg_f_best', avg_norm_grad', avg_n_iter', avg_conv_time', n_success', ...
%     avg_roc', avg_tao_used', 'VariableNames', clmn_labels, 'RowNames',row_names);
% disp(T2)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIRD FUNCTION (16)

%statistics for the table plot

avg_f_best = [];
avg_norm_grad = [];
avg_n_iter = [];
avg_conv_time = [];
n_success = [];
avg_roc = [];
avg_tao_used = [];

%the simulation
for l = 1:1:3

    n = dimensioni(l);

    [f3, gradf3, Hessf3] = third_function_3(n, "simplified_forward", h, true); % Problem 16

    % construction of the test point for f3

    x_f3 = ones(n,1);

    % construction of the other 10 points for f3 

    X_f3 = repmat(x_f3, 1, 10);
    X_f3 = X_f3(:,1:1:10);

    error = rand(n,10);
    X_f3 = X_f3 + error;
    X_f3 = [x_f3, X_f3];


    fbest_iter = []; %to save it in the table
    conv_time_iter = []; %to save it in the table
    norm_grad_iter = []; %to save it in the table
    for i = 1:1:11

        disp(['**** MODIFIED NEWTON METHOD FOR THE THIRD FUNCTION, POINT ', num2str(i), ': STARTED *****']);
        tic;
        [x3k, f3k, gradf3k_norm, k3, x3seq, f3seq, b3tseq, taoseq3, gradf_k3, cos_grad3, fail3] = ...
            Modified_Newton_method(X_f3(:,i), f3, gradf3, Hessf3, ...
            kmax, tolgrad, delta_step, c, rho, btmax, type_tao);
        t = toc;

        disp(['**** MODIFIED NEWTON METHOD FOR THE THIRD FUNCTION, POINT ', num2str(i), ': FINISHED *****']);

        disp(['Time: ', num2str(t), ' seconds']);

        disp('**** MODIFIED NEWTON METHOD : RESULTS *****')
        disp('************************************')
        disp(['N. tao used: ', num2str(nnz(taoseq3))])
        disp(['f(xk): ', num2str(f3k)])
        disp(['gradfk_norm: ', num2str(gradf3k_norm)])
        disp(['N. of Iterations: ', num2str(k3),'/',num2str(kmax), ';'])
        disp(['Rate of convergence: ', num2str(convergence_rate(x3seq)), ';'])
        disp('************************************')

        if fail3 == "kmax"
            result_third_function(i) = 0;
            disp('FAIL: maximum number of iteration reached')
            disp('************************************')
        elseif fail3 == "cosine"
            result_third_function(i) = 0;
            disp('FAIL: consecutive directions almost identical -> no improvment')
            disp('************************************')
        elseif fail3 == "btmax"
            result_third_function(i) = 0;
            disp('FAIL: btmax reached -> could not satisfy the Armijo condition')
            disp('************************************')
        else %fail3 == "success"
            result_third_function(i) = 1;
            disp('SUCCESS')
            disp('************************************')
        end
        disp(' ')

        time_3(i) = t;
        iteration_3(i) = k3;
        number_tao_3(i) = nnz(taoseq3);
        conv_rate_3(i) = convergence_rate(x3seq);
        fbest_iter(i) = f3k; %to save it in the table
        conv_time_iter(i) = t; %to save it in the table
        norm_grad_iter(i) = gradf3k_norm; %to save it in the table
    end

    %figure; 
    %plot(1:k3, f3seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
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

    disp(' ')
    disp(' ')
    disp(' ')
    disp('******************************************')

    disp('**** RESULTS FOR THE THIRD FUNCTION *****')
    disp(['N. of success: ', num2str(sum(result_third_function))])
    valid_idx = ~isnan(conv_rate_3); %ignore when the algorithm does not converge
    avg_roc(l) = sum(result_third_function(valid_idx) .* conv_rate_3(valid_idx)) / sum(result_third_function(valid_idx));
    disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_third_function.*iteration_3)/sum(result_third_function))),'/',num2str(kmax), ';'])
    disp(['Mean N. tao used (in case of success): ', num2str(round(sum(result_third_function.*number_tao_3)/sum(result_third_function)))])
    disp(['Mean convergence rate (in case of success): ', num2str(avg_roc(l))])
    disp('******************************************')

    avg_f_best(l) = mean(fbest_iter);
    avg_norm_grad(l) = mean(norm_grad_iter);
    avg_n_iter(l) = round(mean(iteration_3));
    avg_conv_time(l) = mean(conv_time_iter);
    n_success(l) = sum(result_third_function);
    valid_idx_tao = ~isnan(number_tao_3);
    avg_tao_used(l) = round(sum(result_third_function(valid_idx_tao).*number_tao_3(valid_idx_tao))/sum(result_third_function(valid_idx_tao)));
end

%table plot
clmn_labels = ["avg f(best)", "avg norm grad", "avg n iter", "avg time (s)", "n of success/11", "avg roc (in case of success)", "n tao used (in case of success)"];
row_names = ["1000", "10000", "100000"];
T3 = table(avg_f_best', avg_norm_grad', avg_n_iter', avg_conv_time', n_success', ...
    avg_roc', avg_tao_used', 'VariableNames', clmn_labels, 'RowNames',row_names);
disp(T3)



rho
c
btmax
delta_step
tolgrad

beep