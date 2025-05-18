function [x_opt, fx_opt, k, x_seq, f_seq, singular] = Nelder_Mead_method(X, f, varargin) 

%% INPUTS  
% X = matrix which rows are the starting points of the simplex
% f = is the function handle that we want to work with

% If the following parameters are not provided, they are set to default values
% kmax = maximum number of iterations allowed           (defualt = 200)
% tol = relative tolerance value for stopping criterion (default = 10^(-6))
% rho = reflection's coefficient                        (default = 1)
% chi = expantion's coefficient                         (default = 2)
% gamma = contraction's coefficient                     (default = 0.5)   
% sigma = shrinkage's coefficient                       (default = 0.5)

%% OUTPUTS
% x_opt = best point found
% fx_opt = f(x_opt) (i.e. function evaluated in the best point)
% x_seq = sequence of best points found during the algorithm
% f_seq = values of x_seq (i.e. f(x_seq))
% k = number of iterations performed
% singular = boolean value which tells if the starting simplex is singular (1=true)

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
f_seq = [];
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
%new_indices(k) is the index of the k-th best point of the X matrix
%so new_indices(n+1) is the worst point (x_{n+1} in the slides) 
%f_evaluation(k) is f(x_{k})
values_indices(:,1) = f_evaluation;
values_indices(:,2) = new_indices; %this matrix contains 2 columns;
%a row of this matrix contains f(x) and his corresponding index in the
%matrix X

%% the algorithm

f_bar = mean(values_indices(:,1)); % VARIANCE STOP CRIT
variance = sum((values_indices(:,1)-f_bar).^2)/n; %VARIANCE STOP CRIT
% while k < kmax &&... % check ix x_{1} and x_{n+1} are close or f(x_{1}) is close to f(x_{n+1})
        (abs((values_indices(1,1)-values_indices(n+1,1))/values_indices(1,1)) > tol...
        || norm(X(values_indices(n+1,2),:)-X(values_indices(1,2),:),1)/norm(X(values_indices(1,2),:),1) > tol)
while k < kmax && variance > tol  % VARIANCE STOP CRIT
    f_bar = mean(values_indices(:,1));
    variance = sum((values_indices(:,1)-f_bar).^2)/n;
    
    x_baricenter = mean(X([1:values_indices(n+1,2)-1, values_indices(n+1,2)+1:end],:),1);
    %centroid of the simplex (not taking in account x_{k+1}, the worst point)

    xr = x_baricenter + rho*(x_baricenter - X(values_indices(n+1,2),:));
    f_xr = f(xr); %point reflected evaluated
    k = k + 1;
    if values_indices(1,1) <= f_xr && f_xr < values_indices(n,1) %%% REFLECTION %%%
        X(values_indices(n+1,2),:) = xr; %x{n+1} --> xr
        values_indices(n+1,1) = f_xr;    %f_x{n+1} --> f_{xr}  
        values_indices = sortrows(values_indices,1); %points reordered
        x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
        f_seq(k) = f(x_seq(k,:));
        continue
    elseif f_xr < values_indices(1,1) %%% EXPANTION %%%
        xe = x_baricenter + chi*(xr-x_baricenter);
        f_xe = f(xe); %point expanded evaluated
        if f_xe < f_xr
            X(values_indices(n+1,2),:) = xe; %x_{n+1} --> xe
            values_indices(n+1,1) = f_xe;    %f_x{n+1} --> f_{xe}  
            values_indices = sortrows(values_indices,1); %points reordered
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
            f_seq(k) = f(x_seq(k,:));

            continue
        else
            X(values_indices(n+1,2),:) = xr; %x_{n+1} --> xr
            values_indices(n+1,1) = f_xr;    %f_x{n+1} --> f_{xr}  
            values_indices = sortrows(values_indices,1); %points reordered
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
            f_seq(k) = f(x_seq(k,:));

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
            f_seq(k) = f(x_seq(k,:));
            continue
        else %%% SHRINKAGE %%%
            for i = 2:n+1
                X(values_indices(i,2),:) = X(values_indices(1,2),:) ...
                + sigma*(X(values_indices(i,2),:) - X(values_indices(1,2),:));
                values_indices(i,1) = f(X(values_indices(i,2),:));
            end
            values_indices = sortrows(values_indices,1);
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
            f_seq(k) = f(x_seq(k,:));
        end
    end    
end
x_opt = X(values_indices(1,2),:); % best point
fx_opt = f(x_opt);

end




