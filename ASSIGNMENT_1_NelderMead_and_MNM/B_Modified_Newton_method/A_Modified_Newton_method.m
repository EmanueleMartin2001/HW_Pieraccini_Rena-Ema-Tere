function [xk, fk, gradfk_norm, k, xseq, fseq, btseq, taoseq, gradfk, cos_grad, fail] = ...
    Modified_Newton_method(x0, f, gradf, Hessf, ...
    kmax, tolgrad, delta_step, c1, rho, btmax, type_tao)
%
% [xk, fk, gradfk_norm, k, xseq, btseq] = ...
    % newton_bcktrck(x0, f, gradf, Hessf, ...
    % kmax, tolgrad, c1, rho, btmax)
%
% Function that performs the Newton optimization method, using
% backtracking strategy for the step-length selection.
% 
% There is also the check in case of an Hessian not definite positive
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% Hessf = function handle that describes the Hessian of f;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy.
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the elements xk of the 
% sequence
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.
% gradfk = gradient of the last iteration
% cos_grad = cosine between two consecutives gradients (used as stopping
% criterion)
% fail = string that indicates the kind of insuccess has occured:
% "success", "btmax", "kmax", "cosine".
%

% Function handle for the armijo condition
farmijo = @(fk, alpha, c1_gradfk_pk) ...
    fk + alpha * c1_gradfk_pk;

% Initializations
n = length(x0);
xseq = zeros(n, kmax);
fseq = zeros(1,kmax);
btseq = zeros(1, kmax);
taoseq = zeros(1,kmax);

xk = x0;
fk = f(xk);
gradfk = gradf(xk);
k = 0;
gradfk_norm = norm(gradfk);
delta = sqrt(eps);
cos_grad = 0;
step_norm = delta_step + 1; %so that the first while condition is always satisfied
fail = "success";
while k < kmax && (gradfk_norm >= sqrt(n)*tolgrad || step_norm > sqrt(n)*delta_step)
    %stopping condition given by the norm of the gradient and the norm of
    %the step x_{k+1}-x_{k}
    Hk = Hessf(xk);   % compute the Hessian
    switch type_tao 
        case 'Cholesky'

            [tao,R] = CholeskyAddIdentity(Hk); % it returns the cholescky decomp. of Hk+taoI that we can use for solving the L.S.
            pk = R\(R'\(-gradfk));
            
        case 'Gershgorin'

            tao = Gershgorin_approx(Hk, delta);
            Bk = Hk + tao*diag(ones(length(x0),1)); % computation of Bk as a positive definite matrix
            [pk, ~, ~, iterk, ~] = pcg(Bk, -gradfk);

        case 'Eigen'

            tao = Eigen_tao(Hk,delta);
            Bk = Hk + tao*diag(ones(length(x0),1)); % computation of Bk as a positive definite matrix
            [pk, ~, ~, iterk, ~] = pcg(Bk, -gradfk);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Reset the value of alpha
    alpha = 1;
    
    % Compute the candidate new xk
    xnew = xk + alpha * pk;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    c1_gradfk_pk = c1 * gradfk' * pk;
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
    end
    if bt == btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        fail = "btmax";
        break
    end
    
    % Update xk, fk, gradfk_norm
    xk = xnew;
    fk = fnew;
    if mod(k,10) == mod(5,10) % to not check it at any iterations
        cos_grad = abs(gradfk'*gradf(xk))./(norm(gradfk)*norm(gradf(xk)));
        if cos_grad < 10^(-3) %it means that the new direction is still bad (similar to the previous one)
            fail = "cosine";
            break
        end
    end
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);
    
    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    if k > 1
        step_norm = norm(xseq(:,k)-xseq(:,k-1),2);
    end
    % Store current fk in fseq
    fseq(k) = fk;
    % Store bt iterations in btseq
    btseq(k) = bt;
    % Store current tao in taoseq
    taoseq(k) = tao;

end

% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
fseq = fseq(1:k);
taoseq = taoseq(1:k);

% "Add" x0 at the beginning of xseq (otherwise the first el. is x1)
xseq = [x0, xseq];
if k == kmax
    fail = "kmax";
end

end