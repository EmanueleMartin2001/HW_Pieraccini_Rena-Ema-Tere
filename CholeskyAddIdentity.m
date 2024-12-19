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
    while(k > 0)
        try
            R = chol(Hessian + tao *diag(ones(n,1)));
            k = -1;
        catch
            tao = max([2*tao;Beta]);
        end

        if k > kmax
            k = -1;
        end
    end


end

