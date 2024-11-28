function [f2,gradf2,Hessf2] = second_function(n)

fk = cell(n+1,1);
for k = 1:1:n
    fk{k} = @(x) 1/sqrt(100000)*(x(k)-1);
end
fk{n+1} = @(x) sum(x.^2)-1/4;

f2 = @(x) 0.5*sum(cell2mat(cellfun(@(fk) fk(x), fk, 'UniformOutput', false)).^2);


grad_com = cell(n,1);
for k = 1:1:n
    grad_com{k} = @(x) 1/sqrt(100000) + 2*x(k);
end

gradf2 = @(x) cell2mat(cellfun(@(grad_com) grad_com(x), grad_com, 'UniformOutput', false));

A = 2*ones(n,1);

Hessf2 = @(x) spdiags (A,0, n,n);




end

