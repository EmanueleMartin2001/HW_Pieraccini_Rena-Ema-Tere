function [tao,R] = CholeskyAddIdentity(Hessian)

    Beta = 10^-3;
    kmax = 1000;

    n = length(Hessian(1,:));
    
    % extract the min on the diagonal

    a = min(diag(Hessian));

    % compare

    if a > 0
        tao = 0;
    else
        tao = -a + Beta;
    end
        
    % computing A

    k = 1;
    flag = 1;       %initial pivoting 
    while(k < kmax && flag > 0 )    

        [R,flag] = chol(Hessian + tao *speye(n));   % flag is 0 if the matrix is positive definite

        k = k + 1;

        tao = max([2*tao;Beta]);
    end


end

