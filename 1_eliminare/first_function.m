function [f1, gradf1, Hessf1] = first_function(n)

    f1 = @(x) sum(100 * (x(1:end-1).^2 -x(2:end)).^2 + (x(1:end-1)-1).^2);

    % gradient construction:
    
    f1_grad_comp = cell(n,1);
    
    f1_grad_comp{1} = @(x) 400*x(1)^3 -400*x(1)*x(2) + 2*x(1) - 2;

    for i = 2:1:(n-1)
        f1_grad_comp{i} = @(x) 400*x(i)^3 -200*x(i-1)^2 - 400*x(i+1)*x(i) + 202*x(i) -2;
    end
    f1_grad_comp{n} = @(x) -200*(x(n-1)^2 - x(n));
    
    gradf1 = @(x) cell2mat(cellfun(@(f1_grad_comp) f1_grad_comp(x), f1_grad_comp, 'UniformOutput', false));


    
    
    % computing the Hessian in sparse mode
    
    A1 = cell(n,1);
    B1 = cell(n,1);
    C1 = cell(n,1);
    
    A1{1} = @(x) 1200*x(1)^2 -400*x(2) + 2;
    B1{1} = @(x) -400*x(1);
    C1{1} = @(x) 0; %forced by spdiags
    for i = 2:1:(n-1) 
        A1{i} = @(x) 1200*x(i)^2 - 400*x(i+1) + 202;
        B1{i} = @(x) -400*x(i);
        C1{i} = @(x) -400*x(i-1);
    end
    A1{n} = @(x) 200;
    B1{n} = @(x) 0; %forced by spdiags
    C1{n} = @(x) -400*x(n-1);
    
    Hessf1 = @(x) spdiags (cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false)),-1, n,n)+ ...
                    spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n)+...
                    spdiags (cell2mat(cellfun(@(C1) C1(x), C1, 'UniformOutput', false)),+1, n,n);
end


