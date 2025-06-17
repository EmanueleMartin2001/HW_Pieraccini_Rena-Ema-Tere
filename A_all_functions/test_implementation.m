clc;
clear;
close all;

%% HIGH LEVEL PARAMETERS INITIALIZATION

seed = min([341965, 343316, 284817]);

rng(seed);

functions= {@function_1,@function_75,@function_16};
esponenti = 2:2:12;
H = 10.^(-esponenti); %vettore con gli esponenti da 2 a 12 pari
dimensioni = [10^3,10^4,10^5];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS SETTINGS 

test_function_index=2; 

method= "simplified_centered";

h = 10^(-2); % possible values 2,4,6,8,10,12 

adaptive = true; 

rho = 0.65; %backtracking parameters
c = 0.01;
btmax = 60;

kmax = 1000; %stopping conditions
tolgrad = 1e-3; 
delta_step = 1e-5; 
  
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING FUNCTION 

%statistics for the table plot (total)
avg_f_best = [];
avg_norm_grad = [];
avg_n_iter = [];
avg_conv_time = [];
n_success = [];
avg_roc = [];
avg_tao_used = [];

%statistics for the table plot (in case of convergence)
conv_avg_f_best = [];
conv_avg_norm_grad = [];
conv_avg_n_iter = [];
conv_avg_conv_time = [];
conv_n_success = [];
conv_avg_roc = [];
conv_avg_tao_used = [];



%the simulation
for l = 1:1:3

    n = dimensioni(l);

    [f1,gradf1,Hessf1] = functions{test_function_index}(n, method, h, adaptive); 
    % construction of the test point x0 for f1:

    x_f1 = ones(n,1);
    for i= 1:1:n
        if mod(i,2) == 1
            x_f1(i) = -1.2;
        else
            x_f1(i) = 1.0;
        end
    end

    % construction of the other 10 points for f1

    X_f1 = repmat(x_f1, 1, 10);       % matrix that copy for each column the vector x_f1
    X_f1 = X_f1(:,1:1:10);             % rescale of the matrix

    error = rand(n,10);     % matrix of random variable to add to the starting point 
    X_f1 = X_f1 + error;    % each column of this vector represent a starting point
    X_f1 = [x_f1, X_f1];    %adding the point without noise

    for i = 1:1:11

        disp(['**** MODIFIED NEWTON METHOD FOR THE ',num2str(test_function_index),' FUNCTION, POINT ', num2str(i), ': STARTED *****']);
        tic;
        [x1k, f1k, gradf1k_norm, k1, x1seq,f1seq, bt1seq, taoseq1, gradf_k1, fail1] = ...
            Modified_Newton_method(X_f1(:,i), f1, gradf1, Hessf1, ...
            kmax, tolgrad, delta_step , c, rho, btmax, type_tao);
        t = toc;
        disp(['**** MODIFIED NEWTON METHOD FOR THE ', num2str(test_function_index),' FUNCTION, POINT ', num2str(i), ': FINISHED *****']);

        disp(['Time: ', num2str(t), ' seconds']);

        disp('**** MODIFIED NEWTON METHOD : RESULTS *****');
        disp('************************************');
        disp(['N. tao used: ', num2str(nnz(taoseq1))]);
        disp(['f(xk): ', num2str(f1k)]);
        disp(['gradfk_norm: ', num2str(gradf1k_norm)]);
        disp(['N. of Iterations: ', num2str(k1),'/',num2str(kmax), ';']);
        disp(['Rate of convergence: ', num2str(convergence_rate(x1seq)), ';']);
        disp('************************************');

        if fail1 == "kmax"
            result_first_function(i) = 0;
            disp('FAIL: maximum number of iteration reached')
            disp('************************************')
        elseif fail1 == "cosine"
            result_first_function(i) = 0;
            disp('FAIL: consecutive directions almost identical -> no improvment')
            disp('************************************')
        elseif fail1 == "btmax"
            result_first_function(i) = 0;
            disp('FAIL: btmax reached -> could not satisfy the Armijo condition')
            disp('************************************')
        else %fail1 == "success"
            result_first_function(i) = 1;
            disp('SUCCESS')
            disp('************************************')

        end
        disp(' ')

        time_1(i) = t;
        iteration_1(i) = k1;
        number_tao_1(i) = nnz(taoseq1);
        conv_rate_1(i) = convergence_rate(x1seq);
        fbest_iter(i) = f1k; %to save it in the table
        conv_time_iter(i) = t; %to save it in the table
        norm_grad_iter(i) = gradf1k_norm; %to save it in the table
    end
    % if sum(result_first_function) > 0
    %     figure; 
    %     semilogy(1:success_k1, success_f1seq, 'LineWidth', 2, 'Color', [0.6, 0.2, 0.8]);
    %     grid on;
    %     xlabel('Iterations (k)');
    %     ylabel('Values of the Rosenbrock function'); 
    % 
    %     figure;
    %     bar(1:success_k1, success_taoseq1, 'FaceColor', 'blue', 'EdgeColor', 'black')
    %     grid on;
    %     xlabel('Iterations (k)');
    %     ylabel('Tao values for the Rosenbrock function'); 
    % 
    % end
    disp(' ')
    disp(' ')
    disp(' ')
    disp('******************************************')

    disp('**** RESULTS FOR THE FIRST FUNCTION *****')
    disp(['N. of success: ', num2str(sum(result_first_function)),'/11'])

    valid_idx = ~isnan(conv_rate_1); %ignore when the algorithm does not converge
    avg_roc(l) = sum(result_first_function(valid_idx) .* conv_rate_1(valid_idx)) / sum(result_first_function(valid_idx));
    disp(['Mean N. of Iterations (in case of success): ', num2str(round(sum(result_first_function.*iteration_1)/sum(result_first_function))),'/',num2str(kmax), ';'])
    disp(['Mean N. tao used (in case of success): ', num2str(round(sum(result_first_function.*number_tao_1)/sum(result_first_function)))])
    disp(['Mean convergence rate (in case of success): ', num2str(avg_roc(l))])
    disp('******************************************')

    %total stats
    avg_f_best(l) = mean(fbest_iter);
    avg_norm_grad(l) = mean(norm_grad_iter);
    avg_n_iter(l) = round(mean(iteration_1));
    avg_conv_time(l) = mean(conv_time_iter);
    n_success(l) = sum(result_first_function);
    avg_roc(l) = sum(result_first_function.*conv_rate_1)/sum(result_first_function);
    avg_tao_used(l) = round(sum(result_first_function.*number_tao_1)/sum(result_first_function));
    
    %stats in case of convergence
   
    conv_avg_f_best(l) = dot(fbest_iter,result_first_function)/n_success(l);
    conv_avg_norm_grad(l) = dot(norm_grad_iter,result_first_function)/n_success(l);
    conv_avg_n_iter(l) = round(dot(iteration_1,result_first_function)/n_success(l));
    conv_avg_conv_time(l) = dot(conv_time_iter,result_first_function)/n_success(l); 
    conv_n_success(l) = n_success(l);
    conv_avg_roc(l) = sum(result_first_function.*conv_rate_1)/sum(result_first_function);
    conv_avg_tao_used(l) = round(sum(result_first_function.*number_tao_1)/sum(result_first_function));

end

%table plot
disp("***********TOTAL***********")
clmn_labels = ["avg f(best)", "avg norm grad", "avg n iter", "avg time (s)", "n of success/11", "avg roc (in case of success)", "n tao used (in case of success)"];
row_names = ["1000" , "10000", "100000"];
T1 = table(avg_f_best', avg_norm_grad', avg_n_iter', avg_conv_time', n_success', ...
    avg_roc', avg_tao_used', 'VariableNames', clmn_labels, 'RowNames',row_names);
disp(T1)

%table plot in case of convergence
disp("***********IN CASE OF CONVERGENCE***********")
disp("")
T2 = table(conv_avg_f_best', conv_avg_norm_grad', conv_avg_n_iter', conv_avg_conv_time', conv_n_success', ...
    conv_avg_roc', conv_avg_tao_used', 'VariableNames', clmn_labels, 'RowNames',row_names);
disp(T2)

h
rho
c
btmax
delta_step
tolgrad

beep