function [x_opt,singular] = Nelder_Meada(f, X, tol, kmax, rho, chi, gamma, sigma) 
%[a,b] = Nelder_Meada(f,X,tol,kmax, rho = 1, chi = 2, gamma = 0.5, sigma = 0.5)

n = size(X,2);
singular = "false";
x_opt = 0;

%% 0 singularity check
%at first we have to check that the symplex is not singular
symplex_vectors = ones(n);
for i = 2:n+1
    symplex_vectors(i-1,:) = X(i,:)-X(1,:)
end
if (det(symplex_vectors) == 0) %checks the non singularity
    x_opt = 0;
    singular = "true";
    return
end


%% 1 initialization
%first we have to reorder: we can save the values of the function evalutations
f_evaluation = ones(n+1,1),
for i = 1:n+1
    f_evaluation(i) = f(X(i,:)); %evaluates the function in the point x_i
end

[f_evaluation,new_indices] = sort(f_evaluation) %in ascending order
%new_indices(k) is the index of the k-th best point of the X matrix
%so new_indices(n+1) is the worst point (x_{n+1} in the slides) 
%f_evaluation(k) is f(x_{k})
values_indices(:,1) = f_evaluation
values_indices(:,2) = new_indices %this matrix contains 2 columns;
values_indices %just to check
%a row of this matrix contains f(x) and his corresponding index in the
%matrix X



k = 0
x_baricenter = mean(X([1:values_indices(n+1,2)-1, values_indices(n+1,2)+1:end],:),1)

%% the algorithm
while k < kmax  % and tolerance condition
    X

% X = [2,3; 4,15; 5,10]
% mean(X([1,3],:),1)
    x_baricenter = mean(X([1:values_indices(n+1,2)-1, values_indices(n+1,2)+1:end],:),1)
    %centroid of the simplex  %%%%%%%%% WRONG
    xr = x_baricenter + rho*(x_baricenter - X(values_indices(n+1,2),:))
    f_xr = f(xr); %point reflected evaluated
    k = k + 1
    if values_indices(1,1) <= f_xr && f_xr < values_indices(n,1) %%% REFLECTION %%%
        X(values_indices(n+1,2),:) = xr %x{n+1} --> xr
        values_indices(n+1,1) = f_xr    %f_x{n+1} --> f_{xr}  
        values_indices = sortrows(values_indices,1) %points reordered
        continue
    elseif f_xr < values_indices(1,1) %%% EXPANTION %%%
        xe = x_baricenter + chi*(xr-x_baricenter)
        f_xe = f(xe); %point expanded evaluated
        if f_xe < f_xr
            X(values_indices(n+1,2),:) = xe %x_{n+1} --> xe
            values_indices(n+1,1) = f_xe    %f_x{n+1} --> f_{xe}  
            values_indices = sortrows(values_indices,1) %points reordered
            continue
        else
            X(values_indices(n+1,2),:) = xr %x_{n+1} --> xr
            values_indices(n+1,1) = f_xr    %f_x{n+1} --> f_{xr}  
            values_indices = sortrows(values_indices,1) %points reordered
            continue
        end
    else %%% CONTRACTION %%%  (f_xr >= values_indices(n,1))
        if values_indices(n+1,1) < f_xr %x_{n+1} is better then xr
            xc = x_baricenter - gamma*(x_baricenter - X(values_indices(n+1,2),:))
        else %xr is better then x_{n+1}
            xc = x_baricenter - gamma*(x_baricenter-xr)
        end
        f_xc = f(xc)
        if f_xc < values_indices(n+1,1)
            X(values_indices(n+1,2),:) = xc %x_{n+1} --> xc
            values_indices(n+1,1) = f_xc    %f_x{n+1} --> f_{xc}  
            values_indices = sortrows(values_indices,1) %points reordered
            continue
        else %%% SHRINKAGE %%%
            for i = 2:n+1
                X(values_indices(i,2),:) = X(values_indices(1,2),:) ...
                + sigma*(X(values_indices(i,2),:) - X(values_indices(1,2),:));
                values_indices(i,1) = f(X(values_indices(i,2),:))
            end
            values_indices = sortrows(values_indices,1);
        end
    end    
end
x_opt = X(values_indices(1,2),:) %best point found so far



end




%% test 1
% f = @(x) x(1)+x(2)^2
% X = [2,3; 4,15; 5,10]
% tol = 10e-10
% kmax = 100
%it goes towards -infinity as we expect

%% test 2 Rosenbrock function
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
X = [1.2,1.2; 1,0; 0,1]
kmax = 1000;
tol = 10^(-14);
rho = 1;
chi = 2;
gamma = 0.5;
sigma = 0.5;
[a,b] = Nelder_Meada(f,X,tol,kmax, rho, chi, gamma, sigma)