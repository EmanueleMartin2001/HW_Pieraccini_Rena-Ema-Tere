clc
clear

%NB!
% The function used in this section, are slight modifications of the
% original Nelder Mead function, in order to show graphs and plots, to
% represent visually how the solution is searched by the algorithm.

% These changes are not already implemented in the original function
% because, in order to make the plots, it's necessary to save all the
% values found during the process, wich makes the resolution much slower
% and is opposite of the aim of our code.


% starting simplex: X = [1.2,1.2; -1.2,1; 0,1]
%(This symplex is chosen randomly just for exposition purposes)
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
X = [-5,5; -5.2,-5; 3,1];


%function that generates a video of the evolution of the symplex in the 2D 
%space w.r.t the iterations. This is helpful to show how this method
%searches the optimum value
[x_opt, fx_opt, k, xseq, singular] = Funz_video_nelder(X, f);

%function that prints the value of the best node of the symplex at each
%iterations. Due to an incredibly fast decrease of such value, the print is
%divided into blocks of 20 iterations each
[x_opt, fx_opt, k, x_seq, singular] = Nelder_Mead_stampa_valori(X, f);
