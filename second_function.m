function [f2,gradf2,Hessf2] = second_function(n)

fk = cell(n+1,1);
sqrt100000 = sqrt(100000);
for k = 1:1:n
    fk{k} = @(x) 1/sqrt100000*(x(k)-1);
end
fk{n+1} = @(x) sum(x.^2)-1/4;

f2 = @(x) 0.5*sum(cell2mat(cellfun(@(fk) fk(x), fk, 'UniformOutput', false)).^2);

gradf2 = @(x) 1/10000 * (x-1) + 2*x*(sum(x.^2)-0.25);

Hessf2 = @(x) 1/10000*eye(n) + 2*eye(n)*(sum(x.^2)-0.25) + 4*x*x';

end

