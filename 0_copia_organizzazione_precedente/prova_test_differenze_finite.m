
clc;
clear;
close all;



%%%%%%%% FIRST POINT %%%%%%%%

seed = min([341965, 343316, 284817]);
rng(seed);

%%%%%%%% END FIRST POINT %%%%%%%% 





%%%%%%%% SECOND POINT %%%%%%%%

d =3;    % alternative: 3,4,5

n = 10^d;

H = 10.^(-2*(1:6)); %vector with all the increments to be tested

methods= ["forward","centered","backward"];

for a = 1:length(H) %for every increment
    h=H(a); 

    for j = 1:length(methods) %for every method
        method=methods(j);

        [f1,gradf1,Hessf1] = third_function_3(n,method,h); % Problem 1

        
        % construction of the test point x0 for f1:
        X_f1 = zeros(n,1);
        for i= 1:1:n
            if mod(i,2) == 1
                X_f1(i) = -1.2;
            else
                X_f1(i) = 1.0;
            end
        end
    
        
        % MODIFIED NEWTON METHOD
        rho = 0.5;
        c = 1e-4;
        
        kmax = 1000;
        tolgrad = 1e-4;
        btmax = 40;
        % type_tao = 'Gershgorin';
        % type_tao = 'Eigen';
        type_tao = 'Cholesky';
                
        % calling the method:
        
        %vectors ( ,1) containing the result for every possibile
        %combination of h and methods
        result_first_function = 1000*ones(length(H)*length(methods),1);
        time_1 = zeros(length(H)*length(methods),1);
        iteration_1 = zeros(length(H)*length(methods),1);

        soltol1 = 10^-5;


        disp(['**** MODIFIED NEWTON METHOD FOR THE FIRST FUNCTION, INCREMENT ', num2str(h),' METHOD ',method ,': STARTED *****']);
        tic;
        [x1k, f1k, gradf1k_norm, k1, x1seq,f1seq, bt1seq, taoseq1] = ...
            Modified_Newton_method(X_f1, f1, gradf1, Hessf1, ...
            kmax, tolgrad, c, rho, btmax, type_tao);
        t = toc;
    
        disp(['**** MODIFIED NEWTON METHOD FOR THE FIRST FUNCTION, INCREMENT ', num2str(h),' METHOD ',method , ': FINISHED *****']);


        disp(['Time: ', num2str(t), ' seconds']);

        disp('**** MODIFIED NEWTON METHOD : RESULTS *****');
        disp('************************************');
        disp(['N. tao used: ', num2str(nnz(taoseq1))]);
        disp(['f(xk): ', num2str(f1k)]);
        disp(['gradfk_norm: ', num2str(gradf1k_norm)]);
    
        disp(['N. of Iterations: ', num2str(k1),'/',num2str(kmax), ';']);
        disp(['Rate of convergence: ', num2str(convergence_rate(x1seq)), ';']);
        disp('************************************');

        if k1 == kmax
        result_first_function((a-1)*3+j) = 0;
        disp('FAIL')
        disp('************************************')
        else
            if (norm(x1k-ones(n,1)) < soltol1)
                result_first_function((a-1)*3+j) = 1;
                disp('SUCCESS')
                disp('************************************')
                success_k1 = k1;
                success_f1seq = f1seq;
                success_taoseq1 = taoseq1;
            else
                result_first_function((a-1)*3+j) = 0;
                disp('FAIL')
                disp('************************************')
            end
        end

        disp(' ')
       time_1((a-1)*3+j) = t;
       iteration_1((a-1)*3+j) = k1;
       number_tao_1((a-1)*3+j) = nnz(taoseq1);
       conv_rate_1((a-1)*3+j) = convergence_rate(x1seq);


    
    end
end



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
