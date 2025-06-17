clc;
clear;
close all;

%% HIGH LEVEL PARAMETERS INITIALIZATION

seed = min([341965, 343316, 284817]);

rng(seed);

dimension = [10,25,50];

% Starting from a given point x0,
% we construct a matrix with n+1 rows:
% the first row corresponds to x0 itself,
% while the remaining n rows represent perturbations of x0
% along its n different coordinates.
% Thus we have the starting symplex.

%% PARAMETERS SETTINGS 
kmax = 100000;

tol = 1e-4;

rho = 1;

chi = 2;

gamma = 0.5;

sigma = 0.5;

%% TESTS

%%%  FIRST FUNCTION  %%%
% the construction of the 10 starting points is embedeed in the following
% for cycle
avg_f_best = [];
avg_n_iter = [];
avg_conv_time = [];
n_success = [];
for l = 1:1:3

    n = dimension(l);

    x_f1 = zeros(n,1);
    for i= 1:1:n
       if mod(i,2) == 1
          x_f1(i) = -1.2;
       else
          x_f1(i) = 1.0;
       end
    end

    f1 = function_1(n); % Problem 1

    for i = 1:1:11
        if i > 1 %if i == 1 the starting point is unchanged
            x_f1 = x_f1 + 2*rand(n,1)-1; % x_f1 perturbed in the hypercube
        end
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
        iteration(i) = k1;
        fbest_iter(i) = fx1_opt; %to save it in the table
        conv_time_iter(i) = t; %to save it in the table
    end    
    figure; 
    f1_seq = f1_seq + abs(min(f1_seq)) + 1e-3;
    plot(1:k1, f1_seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
    grid on;
    xlabel('Iterations (k)');
    ylabel('Values of the Chained Rosenbrock function');
    set(gcf, 'Color', 'w');  % Imposta lo sfondo bianco
    text(0.5, 1.025, ['Dimension: n = ', num2str(n)], 'Units', 'normalized', ... %dimension on top of the plot
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
    disp(' ')
    disp(' ')
    disp(' ')
    disp('******************************************')

    disp('**** RESULTS FOR THE FIRST FUNCTION *****')
    disp(['N. of success: ', num2str(sum(result_first_function)),'/11'])
    disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_first_function.*iteration)/sum(result_first_function))),'/',num2str(kmax), ';'])
    disp('******************************************')

    avg_f_best(l) = mean(fbest_iter);
    avg_n_iter(l) = round(mean(iteration));
    avg_conv_time(l) = mean(conv_time_iter);
    n_success(l) = sum(result_first_function);
end

%table plot
clmn_labels = ["avg f(best)", "avg n iter", "avg time (s)", "n of success/11"];
row_names = ["10", "25", "50"];
T1 = table(avg_f_best', avg_n_iter', avg_conv_time', n_success', ...
     'VariableNames', clmn_labels, 'RowNames',row_names);
disp(T1)

% %%%  SECOND FUNCTION  %%%
% % the construction of the 10 starting points is embedeed in the following
% % for cycle
% 
% avg_f_best = [];
% avg_n_iter = [];
% avg_conv_time = [];
% n_success = [];
% for l = 1:1:3
% 
%     n = dimension(l);
% 
%     x_f2 = -ones(n,1)*1.2;  %starting point
%     x_f2(end) = -1;
% 
%     f2 = function_75(n); % Problem 75
% 
%     for i = 1:1:11
%         if i > 1 % i ==1 leave the starting point unchanged
%             x_f2 = x_f2 + 2*rand(n,1)-1; % x_f2 perturbed in the hypercube
%         end
%         alpha_2 = 1; %explore term (decides how big should the starting symplex be)
%         X2 = zeros(n+1,n);
%         X2(1,:) = x_f2;
%         for j = 2:n+1
%             X2(j,:) = x_f2;
%             X2(j,j-1) = X2(j,j-1) + alpha_2; % perturbation of 1 coordinate by alpha_2
%         end
% 
%         disp(['**** NELDER MEAD METHOD FOR THE SECOND FUNCTION, STARTING SYMPLEX ', num2str(i), ': STARTED *****']);
%         tic;
%         [x2_opt, fx2_opt, k2, x2seq, f2_seq, singular2] = Nelder_Mead(X2, f2, kmax, tol, rho, chi, gamma, sigma);
%         t = toc;
% 
%         disp(['**** NELDER MEAD METHOD FOR THE SECOND FUNCTION, STARTING SYMPLEX ', num2str(i), ': FINISHED *****']);
% 
%         disp(['Time: ', num2str(t), ' seconds']);
% 
%         disp('**** NELDER MEAD METHOD : RESULTS *****')
%         disp('************************************')
%         disp(['rho: ',num2str(rho), '  chi: ', num2str(chi), '  gamma: ',num2str(gamma),'  sigma: ', num2str(sigma)])
%         disp(['fx2_opt: ', num2str(fx2_opt)])
%         disp(['N. of Iterations: ', num2str(k2),'/',num2str(kmax), ';'])
%         disp('************************************')
% 
%         if k2 == kmax || singular2 == 1 
%             result_second_function(i) = 0;
%             disp('FAIL')
%             disp('************************************')
%         else
%             result_second_function(i) = 1;
%             disp('SUCCESS')
%             disp('************************************')
%         end
%         disp(' ')
%         iteration(i) = k2;
%         fbest_iter(i) = fx2_opt; %to save it in the table
%         conv_time_iter(i) = t; %to save it in the table
%     end
%     figure; 
%     f2_seq = f2_seq + abs(min(f2_seq)) + 1e-3;
%     semilogy(1:k2, f2_seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
%     grid on;
%     xlabel('Iterations (k)');
%     ylabel('Values of the Chained Rosenbrock function');
%     set(gcf, 'Color', 'w');  % Imposta lo sfondo bianco
%     text(0.5, 1.025, ['Dimension: n = ', num2str(n)], 'Units', 'normalized', ... %dimension on top of the plot
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
%     disp(' ')
%     disp(' ')
%     disp(' ')
%     disp('******************************************')
% 
%     disp('**** RESULTS FOR THE SECOND FUNCTION *****')
%     disp(['N. of success: ', num2str(sum(result_second_function)),'/11'])
%     disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_second_function.*iteration)/sum(result_second_function))),'/',num2str(kmax), ';'])
%     disp('******************************************')
% 
%     avg_f_best(l) = mean(fbest_iter);
%     avg_n_iter(l) = round(mean(iteration));
%     avg_conv_time(l) = mean(conv_time_iter);
%     n_success(l) = sum(result_second_function);
% end
% 
% %table plot
% clmn_labels = ["avg f(best)", "avg n iter", "avg time (s)", "n of success/11"];
% row_names = ["10", "25", "50"];
% T2 = table(avg_f_best', avg_n_iter', avg_conv_time', n_success', ...
%      'VariableNames', clmn_labels, 'RowNames',row_names);
% disp(T2)


% %%%  THIRD FUNCTION  %%%
% % the construction of the 10 starting points is embedeed in the following
% % for cycle
% 
% avg_f_best = [];
% avg_n_iter = [];
% avg_conv_time = [];
% n_success = [];
% for l = 1:1:3
% 
%     n = dimension(l);
% 
%     x_f3 = ones(n,1);
% 
%     f3 = function_16(n); % Problem 16
% 
%     for i = 1:1:11
% 
%         if i > 1 %for i == 1 the starting point is untouched
%             x_f3 = x_f3 + 2*rand(n,1)-1; % x_f3 perturbed in the hypercube
%         end
%         alpha_3 = 1; %explore term (decides how big should the starting symplex be)
%         X3 = zeros(n+1,n);
%         X3(1,:) = x_f3;
%         for j = 2:n+1
%             X3(j,:) = x_f3;
%             X3(j,j-1) = X3(j,j-1) + alpha_3; % perturbation of 1 coordinate by alpha_3
%         end
% 
%         disp(['**** NELDER MEAD METHOD FOR THE THIRD FUNCTION, STARTING SYMPLEX ', num2str(i), ': STARTED *****']);
%         tic;
%         [x3_opt, fx3_opt, k3, x3seq, f3_seq, singular3] = Nelder_Mead(X3, f3, kmax, tol, rho, chi, gamma, sigma);
%         t = toc;
% 
%         disp(['**** NELDER MEAD METHOD FOR THE THIRD FUNCTION, STARTING SYMPLEX ', num2str(i), ': FINISHED *****']);
% 
%         disp(['Time: ', num2str(t), ' seconds']);
% 
%         disp('**** NELDER MEAD METHOD : RESULTS *****')
%         disp('************************************')
%         disp(['rho: ',num2str(rho), '  chi: ', num2str(chi), '  gamma: ',num2str(gamma),'  sigma: ', num2str(sigma)])
%         disp(['fx3_opt: ', num2str(fx3_opt)])
%         disp(['N. of Iterations: ', num2str(k3),'/',num2str(kmax), ';'])
%         disp('************************************')
%         % n = 10 => f(x3_opt) = -8.0514
%         % n = 25 => f(x3_opt) = -16.1379
%         % n = 50 => f(x3_opt) = -27.8949
%         if k3 == kmax || singular3 == 1
%             result_third_function(i) = 0;
%             disp('FAIL')
%             disp('************************************')
%         else
%             result_third_function(i) = 1;
%             disp('SUCCESS')
%             disp('************************************')
%         end
%         disp(' ')
%         iteration(i) = k3;
%         fbest_iter(i) = fx3_opt; %to save it in the table
%         conv_time_iter(i) = t; %to save it in the table
%     end
%     figure; 
%     f3_seq = f3_seq + abs(min(f3_seq)) + 1e-3;
%     semilogy(1:k3, f3_seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
%     grid on;
%     xlabel('Iterations (k)');
%     ylabel('Values of the Banded trigonometric problem');
%     set(gcf, 'Color', 'w');  % Imposta lo sfondo bianco
%     text(0.5, 1.025, ['Dimension: n = ', num2str(n)], 'Units', 'normalized', ... %dimension on top of the plot
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
%     disp(' ')
%     disp(' ')
%     disp(' ')
%     disp('******************************************')
% 
%     disp('**** RESULTS FOR THE THIRD FUNCTION *****')
%     disp(['N. of success: ', num2str(sum(result_third_function)),'/11'])
%     disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_third_function.*iteration)/sum(result_third_function))),'/',num2str(kmax), ';'])
%     disp('******************************************')
% 
%     avg_f_best(l) = mean(fbest_iter);
%     avg_n_iter(l) = round(mean(iteration));
%     avg_conv_time(l) = mean(conv_time_iter);
%     n_success(l) = sum(result_third_function);
% end
% 
% %table plot
% clmn_labels = ["avg f(best)", "avg n iter", "avg time (s)", "n of success/11"];
% row_names = ["10", "25", "50"];
% T3 = table(avg_f_best', avg_n_iter', avg_conv_time', n_success', ...
%      'VariableNames', clmn_labels, 'RowNames',row_names);
% 
% disp(T3)
