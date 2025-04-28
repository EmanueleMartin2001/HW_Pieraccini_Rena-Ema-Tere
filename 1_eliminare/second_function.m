function [f2, gradf2, Hessf2] = second_function(n)

fk = cell(n,1);
fk{1} = @(x) (3-2*x(1))*x(1) -2*x(2) + 1;
for k = 2:1:n-1
    fk{k} = @(x) (3-2*x(k))*x(k) -x(k-1) -2*x(k+1) + 1;
end
fk{n} = @(x) (3-2*x(n))*x(n) -x(n-1) + 1;

f2 = @(x) 0.5*sum(cell2mat(cellfun(@(fk) fk(x), fk, 'UniformOutput', false)).^2);

deriv_j = cell(n,1);
deriv_j{1} = @(x) fk{1}(x) * (3-4*x(1)) + fk{2}(x) * (-2);
for k = 2:1:n-1
    deriv_j{k} = @(x)  fk{k-1}(x) * (-1) + fk{k}(x) * (3-4*x(k)) + fk{k+1}(x) * (-2);
end
deriv_j{n} = @(x)  fk{n-1}(x) * (-1) + fk{n}(x) * (3-4*x(n));

gradf2 = @(x) cell2mat(cellfun(@(deriv_j) deriv_j(x), deriv_j, 'UniformOutput', false));

% computing the Hessian in sparse mode
    
    A1 = cell(n,1);
    B1 = cell(n,1);
    C1 = cell(n,1);
    
    A1{1} = @(x) (3-4*x(1))^2 + (-2)*(-2) + fk{1}(x)*(-4);
    B1{1} = @(x) 0;
    C1{1} = @(x) (3-4*x(1))*(-1) + (-2)*(3-4*x(2)); %forced by spdiags
    for i = 2:1:(n-1) 
        A1{i} = @(x) (-1)*(-1) + (3-4*x(k))^2 + (-2)*(-2) + fk{k}(x)*(-4);
        B1{i} = @(x) (3-4*x(k-1))*(-1) + (-2)*(3-4*x(k));
        C1{i} = @(x) (3-4*x(k))*(-1) + (-2)*(3-4*x(k+1));
    end
    A1{n} = @(x) (-1)*(-1) + (3-4*x(n))^2 + fk{n}(x)*(-4);
    B1{n} = @(x) (3-4*x(n-1))*(-1) + (-2)*(3-4*x(n)); %forced by spdiags
    C1{n} = @(x) 0;
    
    Hessf2 = @(x) spdiags (cell2mat(cellfun(@(B1) B1(x), C1, 'UniformOutput', false)),-1, n,n)+ ...
                    spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n)+...
                    spdiags (cell2mat(cellfun(@(C1) C1(x), B1, 'UniformOutput', false)),+1, n,n);

end

