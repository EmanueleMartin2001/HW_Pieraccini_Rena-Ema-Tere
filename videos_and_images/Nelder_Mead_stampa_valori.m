%function that prints the value of the best node of the symplex at each
%iterations. Due to an incredibly fast decrease of such value, the print is
%divided into blocks of 20 iterations each


%NB!
% The function used in this section, are slight modifications of the
% original Nelder Mead function, in order to show graphs and plots, to
% represent visually how the solution is searched by the algorithm.

% These changes are not already implemented in the original function
% because, in order to make the plots, it's necessary to save all the
% values found during the process, wich makes the resolution much slower
% and is opposite of the aim of our code.

function [x_opt, fx_opt, k, x_seq, singular] = Nelder_Mead_stampa_valori(X, f, varargin) 


%% input parsing
p = inputParser;

% adding optional parameters with default values

addOptional(p, 'kmax', 200);  
addOptional(p, 'tol', 10^(-6));
addOptional(p, 'rho', 1);
addOptional(p, 'chi', 2);
addOptional(p, 'gamma', 0.5);
addOptional(p, 'sigma', 0.5);

parse(p, varargin{:}); % add those variables as input
kmax= p.Results.kmax;
tol = p.Results.tol;
rho = p.Results.rho;
chi= p.Results.chi;
gamma = p.Results.gamma;
sigma = p.Results.sigma;


% inizialization for default output values
x_opt = [];
fx_opt = [];
k = 0;
x_seq = [];
singular = false;

%% 0 singularity check
%at first we have to check that the simplex is not singular
n = size(X,2);
simplex_vectors = ones(n);
for i = 2:n+1
    simplex_vectors(i-1,:) = X(i,:)-X(1,:);
end
if (det(simplex_vectors) == 0) %checks the non singularity
    singular = true;
    return
end


%% 1 initialization
%first we have to reorder: we can save the values of the function evalutations
f_evaluation = ones(n+1,1);
for i = 1:n+1
    f_evaluation(i) = f(X(i,:)); %evaluates the function in the point x_i
end

[f_evaluation,new_indices] = sort(f_evaluation); %in ascending order

values_indices(:,1) = f_evaluation;
values_indices(:,2) = new_indices; %this matrix contains 2 columns;


%% the algorithm
while k < kmax &&... % check ix x_{1} and x_{n+1} are close or f(x_{1}) is close to f(x_{n+1})
        (abs((values_indices(1,1)-values_indices(n+1,1))/values_indices(1,1)) > tol...
        || norm(X(values_indices(n+1,2),:)-X(values_indices(1,2),:),1)/norm(X(values_indices(1,2),:),1) > tol)
    x_baricenter = mean(X([1:values_indices(n+1,2)-1, values_indices(n+1,2)+1:end],:),1);
    %centroid of the simplex (not taking in account x_{k+1}, the worst point)



%--------------------------------------------------------------------
%-----------------------only differences-------------------------

    valori_stampa(k+1)= values_indices(1,1);

%--------------------------------------------------------------------
%--------------------------------------------------------------------





    xr = x_baricenter + rho*(x_baricenter - X(values_indices(n+1,2),:));
    f_xr = f(xr); %point reflected evaluated
    k = k + 1;
    if values_indices(1,1) <= f_xr && f_xr < values_indices(n,1) %%% REFLECTION %%%
        X(values_indices(n+1,2),:) = xr; %x{n+1} --> xr
        values_indices(n+1,1) = f_xr;    %f_x{n+1} --> f_{xr}  
        values_indices = sortrows(values_indices,1); %points reordered
        x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
        
        continue
    elseif f_xr < values_indices(1,1) %%% EXPANTION %%%
        xe = x_baricenter + chi*(xr-x_baricenter);
        f_xe = f(xe); %point expanded evaluated
        if f_xe < f_xr
            X(values_indices(n+1,2),:) = xe; %x_{n+1} --> xe
            values_indices(n+1,1) = f_xe;    %f_x{n+1} --> f_{xe}  
            values_indices = sortrows(values_indices,1); %points reordered
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
            continue
        else
            X(values_indices(n+1,2),:) = xr; %x_{n+1} --> xr
            values_indices(n+1,1) = f_xr;    %f_x{n+1} --> f_{xr}  
            values_indices = sortrows(values_indices,1); %points reordered
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
            continue
        end
    else %%% CONTRACTION %%%  (f_xr >= values_indices(n,1))
        if values_indices(n+1,1) < f_xr %x_{n+1} is better then xr
            xc = x_baricenter - gamma*(x_baricenter - X(values_indices(n+1,2),:));
        else %xr is better then x_{n+1}
            xc = x_baricenter - gamma*(x_baricenter-xr);
        end
        f_xc = f(xc);
        if f_xc < values_indices(n+1,1)
            X(values_indices(n+1,2),:) = xc; %x_{n+1} --> xc
            values_indices(n+1,1) = f_xc;    %f_x{n+1} --> f_{xc}  
            values_indices = sortrows(values_indices,1); %points reordered
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
            continue
        else %%% SHRINKAGE %%%
            for i = 2:n+1
                X(values_indices(i,2),:) = X(values_indices(1,2),:) ...
                + sigma*(X(values_indices(i,2),:) - X(values_indices(1,2),:));
                values_indices(i,1) = f(X(values_indices(i,2),:));
            end
            values_indices = sortrows(values_indices,1);
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
        end
    end    
end
x_opt = X(values_indices(1,2),:); % best point
fx_opt = f(x_opt);



%--------------------------------------------------------------------
%-----------------------only differences-------------------------


%plot complessivo
figure;
plot(linspace(0,k-1,k), valori_stampa,'-b', 'LineWidth', 1.5,'MarkerSize', 5, 'MarkerFaceColor', 'g');
%plot togliendo i primi due valori
figure;
plot(linspace(2,k-1,k-2), valori_stampa(3:end),'-b', 'LineWidth', 1.5,'MarkerSize', 5, 'MarkerFaceColor', 'b');

%plot a blocchi
dim_blocco=20;
blocchi=floor(k/dim_blocco);
for i=1:blocchi    
    figure;
    plot(linspace((i-1)*dim_blocco,(i-1)*dim_blocco+dim_blocco-1,dim_blocco),...
        valori_stampa((i-1)*dim_blocco+1:(i-1)*dim_blocco+dim_blocco),...
        '-ob', 'LineWidth', 1.5,'MarkerSize', 5, 'MarkerFaceColor', 'b');


%--------------------------------------------------------------------
%--------------------------------------------------------------------



end

end