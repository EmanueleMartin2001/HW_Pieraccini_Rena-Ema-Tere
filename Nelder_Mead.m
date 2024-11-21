function [x_opt,singular] = Nelder_Meada(f, X, tol, kmax, rho, chi, gamma, sigma) 

n = size(X,2);
singular = false;
x_opt = 0

%% 0 singularity check
%at first we have to check that the symplex is not singular
symplex_vectors = ones(n);
for i = 2:n+1
    symplex_vectors(i-1,:) = X(i,:)-X(1,:)
end
if (det(symplex_vectors) == 0) %checks the non singularity
    x_opt = 0;
    singular = true;
    return
end


% no debbuged part
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


k = 0
%% the algorithm
while k < kmax : % and tolerance condition
    x_baricenter = mean(X,1) %centroid of the simplex
    x_r = x_baricenter + rho*(x_baricenter - X(new_indices(n+1),:))
    % 
    if f_evaluation(1) <= f(x_r) && f(x_r) < f_evaluation(n) :
        %substitute x_{r} with x_{n+1}
        X(new_indices(n+1),:) = x_r
        continue;
    end

end



end





f = @(x) x(1)+x(2)^2
X = [2,3; 4,15; 5,10]
tol = 10e-10
kmax = 100
[a,b] = Nelder_Meada(f,X,tol,kmax)