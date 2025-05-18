clc
clear all

%% MODIFIED NEWTON METHOD

f_Rosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

gradf_Rosenbrock = @(x) [400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2 ; 200*( x(2) - x(1)^2 )];

Hessianf_Rosenbrock = @(x) [ 1200*x(1)^2 - 400*x(2) + 2 , -400*x(1); -400*x(1) , 200];

rho = 0.5;

c1 = 1e-4;

type_tao = 'Cholesky';

x_0_1 = [1.2 ; 1.2];
x_0_2 = [-1.2 ; 1];

kmax = 100;

btmax = 150;

tolgrad = 1e-5;
delta_step = 1e-5;

%% RUN MODIFIED NEWTON METHOD

disp('**** MODIFIED NEWTON METHOD FOR POINT 1 : START *****')
tic;
[xk_mnm, fk_mnm, gradfk_norm, k_mnm, xseq_mnm , f_seq1, btseq] = Modified_Newton_method(x_0_1, f_Rosenbrock, gradf_Rosenbrock, Hessianf_Rosenbrock, kmax, tolgrad, delta_step, c1, rho, btmax, type_tao);
t = toc;
k_plot = k_mnm;

disp('**** MODIFIED NEWTON METHOD FOR POINT 1: FINISHED *****')

disp(['Time: ', num2str(t), ' seconds']);

disp('**** MODIFIED NEWTON METHOD : RESULTS *****')
disp('************************************')
disp(['rho: ', mat2str(rho), '  c: ', mat2str(c1)])
disp(['x_opt: ', mat2str(xk_mnm)])
disp(['f(x_opt): ', num2str(fk_mnm)])
disp(['N. of Iterations: ', num2str(k_mnm),'/',num2str(kmax), ';'])
disp('************************************')

% MNM for the second initial point

disp('**** MODIFIED NEWTON METHOD FOR POINT 2: START *****')
tic;
[xk_mnm, fk_mnm, gradfk_norm, k_mnm, xseq_mnm, f_seq2, btseq] = Modified_Newton_method(x_0_2, f_Rosenbrock, gradf_Rosenbrock, Hessianf_Rosenbrock, kmax, tolgrad, delta_step, c1, rho, btmax, type_tao);
t = toc;
if k_mnm > k_plot
    f_seq1 = [f_seq1,zeros(1, k_mnm-k_plot)];
    k_plot = k_mnm;
else
     f_seq2 = [f_seq2,zeros(1,k_plot-k_mnm)];
end

disp('**** MODIFIED NEWTON METHOD FOR POINT 2: FINISHED *****')

disp(['Time: ', num2str(t), ' seconds']);
format short
disp('**** MODIFIED NEWTON METHOD : RESULTS *****')
disp(['rho: ', mat2str(rho), '  c: ', mat2str(c1)])
disp(['x_opt: ', mat2str(xk_mnm)])
disp(['f(x_opt): ', num2str(fk_mnm)])
disp(['N. of Iterations: ', num2str(k_mnm),'/',num2str(kmax), ';'])
disp('************************************')
%plot
set(gcf, 'Color', 'w');  % Imposta lo sfondo bianco
figure; 
semilogy(1:k_plot, f_seq1, 'LineWidth', 2, 'Color', [0.6, 0.0, 0.2]);
hold on
set(gcf, 'Color', 'w');  % Imposta lo sfondo bianco
grid on;
xlabel('Iterations (k)');
ylabel('Values of the Rosenbrock function');
semilogy(1:k_plot, f_seq2, 'LineWidth', 2, 'Color', [0.1, 0.5, 0.5]);
format short

%% PLOTS (BACKTRACK)

f_meshgrid = @(X, Y) arrayfun(@(i) f_Rosenbrock([X(i); Y(i)]), 1:numel(X));

% Creation of the meshgrid for the contour-plot
[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
% Computation of the values of f for each point of the mesh
Z = reshape(f_meshgrid(X, Y), size(X));

fk_mnm_seq = arrayfun(@(i) f_Rosenbrock(xseq_mnm(:, i)), 1:size(xseq_mnm, 2));

% Plots

% Simple Plot
fig1 = figure();
% Contour plot with curve levels for each point in xseq
[C1, ~] = contour(X, Y, Z);
hold on
% plot of the points in xseq
plot(xseq_mnm(1,:), xseq_mnm(2,:), '--*')

hold off

% More interesting Plot
fig2 = figure();
% Contour plot with curve levels for each point in xseq
% ATTENTION: actually, the mesh [X, Y] is too coarse for plotting the last
% level curves corresponding to the last point in xseq (check it by zooming
% the image).
[C2, ~] = contour(X, Y, Z, fk_mnm_seq);
hold on
% plot of the points in xseq
plot(xseq_mnm(1, :), xseq_mnm(2, :), '--*')
hold off

% Much more interesting plot
fig4 = figure();
surf(X, Y, Z,'EdgeColor','none')
hold on
plot3(xseq_mnm(1, :), xseq_mnm(2, :), fk_mnm_seq, 'r--*')
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% NELDER MEAD METHOD

f_Rosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

kmax = 1000;

tol = 1e-5;

rho = 1;

chi = 2;

gamma = 0.5;

sigma = 0.5;



%% RUN NELDER-MEAD METHOD

k_plot = 0; %used for plotting the results
% starting symplex 1
X1 = [1.2,1.2; 2.2,1.2; 1.2,2.2];

disp('**** NELDER MEAD METHOD FOR SYMPLEX 1: STARTED *****')

tic;
[x_opt, fx_opt, k, x_seq, f_seq1, singular] = Nelder_Mead_method(X1, f_Rosenbrock, kmax, tol, rho, chi, gamma, sigma);
t = toc;

disp('**** NELDER MEAD METHOD FOR SYMPLEX 1: FINISHED *****')

disp(['Time: ', num2str(t), ' seconds']);

disp('**** NELDER MEAD METHOD FOR SYMPLEX 1: RESULTS *****')
disp('************************************')
disp(['rho: ',num2str(rho), '  chi: ', num2str(chi), '  gamma: ',num2str(gamma),'  sigma: ', num2str(sigma)])
disp(['x_opt: ', mat2str(x_opt)])
disp(['f(x_opt): ', num2str(fx_opt)])
disp(['N. of Iterations: ', num2str(k),'/',num2str(kmax), ';'])
disp('************************************')

k_plot = k;

% starting symplex 2
X2 = [-1.2,1; -0.2,1; -1.2, 2];

disp('**** NELDER MEAD METHOD FOR SYMPLEX 2: STARTED *****')

tic;
[x_opt, fx_opt, k, x_seq, f_seq2, singular] = Nelder_Mead_method(X2, f_Rosenbrock, kmax, tol, rho, chi, gamma, sigma);

if k > k_plot
    f_seq1 = [f_seq1,zeros(1,k-k_plot)];
    k_plot = k;
else
    f_seq2 = [f_seq2,zeros(1,k_plot-k)];
end

t = toc;

disp('**** NELDER MEAD METHOD FOR SYMPLEX 2: FINISHED *****')

disp(['Time: ', num2str(t), ' seconds']);

disp('**** NELDER MEAD METHOD FOR SYMPLEX 2: RESULTS *****')
disp('************************************')
disp(['rho: ',num2str(rho), '  chi: ', num2str(chi), '  gamma: ',num2str(gamma),'  sigma: ', num2str(sigma)])
disp(['x_opt: ', mat2str(x_opt)])
disp(['f(x_opt): ', num2str(fx_opt)])
disp(['N. of Iterations: ', num2str(k),'/',num2str(kmax), ';'])
disp('************************************')

%plot
set(gcf, 'Color', 'w');  % Imposta lo sfondo bianco
figure; 
semilogy(1:k_plot, f_seq1, 'LineWidth', 2, 'Color', [0.6, 0.0, 0.2]);
hold on
set(gcf, 'Color', 'w');  % Imposta lo sfondo bianco
grid on;
xlabel('Iterations (k)');
ylabel('Values of the Rosenbrock function');
semilogy(1:k_plot, f_seq2, 'LineWidth', 2, 'Color', [0.1, 0.5, 0.5]);

