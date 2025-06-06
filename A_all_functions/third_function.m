function [f3, gradf3, Hessf3] = third_function(n)

    fk = cell(n,1);
    fk{1} = @(x) 1*((1-cos(x(1)) + 0 - sin(x(2))));
    for k = 2:1:n-1
        fk{k} = @(x) k* ( (1-cos(x(k))) + sin(x(k-1)) - sin(x(k+1)) )   ;
    end
    fk{n} = @(x) n* ((1-cos(x(n))) + sin(x(n-1)));

    f3 = @(x) sum(cell2mat(cellfun(@(fk) fk(x), fk, 'UniformOutput', false)));

    deriv_j = cell(n,1);
    deriv_j{1} = @(x) 1*sin(x(1)) + 2*cos(x(1))   ;
    for k = 2:1:n-1
        deriv_j{k} = @(x) k*sin(x(k)) + 2*cos(x(k))   ;
    end
    deriv_j{n} = @(x) n*sin(x(n)) - (n-1)*cos(x(n))   ;

    gradf3 = @(x) cell2mat(cellfun(@(deriv_j) deriv_j(x), deriv_j, 'UniformOutput', false));



    A2 = cell(n,1);   % diagonal
    A2{1} = @(x) cos(x(1)) - 2* sin(x(1));
    for k = 1:1:n
        A2{k} = @(x) k*cos(x(k)) - 2*sin(x(k));
    end
    A2{n} = @(x) n*cos(x(n)) + (n-1)*sin(x(n));

    Hessf3 = @(x) spdiags (cell2mat(cellfun(@(A2) A2(x), A2, 'UniformOutput', false)),0, n,n);


end

