function [f2, gradf2, Hessf2] = second_function_75(n)


fk = cell(n,1);
fk{1} = @(x) x(1) - 1;
for k = 2:1:n
    fk{k} = @(x) 10*(k-1)*(x(k)-x(k-1));
end

f2 = @(x) 0.5*sum(cell2mat(cellfun(@(fk) fk(x), fk, 'UniformOutput', false)).^2);

deriv_j = cell(n,1);
deriv_j{1} = @(x) x(1)-1-200*(x(2)-x(1))^3;
for k = 2:1:n-1
    deriv_j{k} = @(x) 200*(k-1)^2*(x(k)-x(k-1))^3 - 200*k^2*(x(k+1)-x(k))^3;
end
deriv_j{n} = @(x) 200*(n-1)^2*(x(n)-x(n-1))^3;


gradf2 = @(x) cell2mat(cellfun(@(deriv_j) deriv_j(x), deriv_j, 'UniformOutput', false));

% computing the Hessian in sparse mode
    
    A1 = cell(n,1);
    B1 = cell(n,1);
    C1 = cell(n,1);
    
    A1{1} = @(x) 1+600*(x(2)-x(1))^2; %diag
    B1{1} = @(x) -600*(x(2)-x(1))^2; %under
    C1{1} = @(x) 0; %forced by spdiags %over

    for i = 2:1:(n-1) 
        A1{i} = @(x) 600*(i-1)^2*(x(i)-x(i-1))^2 + 600*i^2*(x(i+1)-x(i))^2; 
        B1{i} = @(x) -600*i^2*(x(i+1)-x(i))^2;
        C1{i} = @(x) -600*(i-1)^2*(x(i)-x(i-1))^2;
    end
    A1{n} = @(x) 600*(n-1)^2*(x(n)-x(n-1))^2;
    B1{n} = @(x) 0; %forced by spdiags
    C1{n} = @(x) -600*(n-1)^2*(x(n)-x(n-1))^2; 
    
    
    Hessf2 = @(x) spdiags (cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false)),-1, n,n)+ ...
                    spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n)+...
                    spdiags (cell2mat(cellfun(@(C1) C1(x), C1, 'UniformOutput', false)),+1, n,n);

end