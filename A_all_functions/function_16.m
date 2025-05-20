function [f3, gradf3, Hessf3] = function_16(n,method, h, adaptive)

    %method can be exact or simplified_forward
    if nargin == 1
        method = "exact";
    end
    if nargin == 2
        h = 0.0001;
        adaptive = false
    end
    if nargin == 3
        adaptive = false
    end

    %building the function
    fk = cell(n,1);
    fk{1} = @(x) 1*((1-cos(x(1)) + 0 - sin(x(2))));
    for k = 2:1:n-1
        fk{k} = @(x) k* ( (1-cos(x(k))) + sin(x(k-1)) - sin(x(k+1)) )   ;
    end
    fk{n} = @(x) n* ((1-cos(x(n))) + sin(x(n-1)));

    f3 = @(x) sum(cell2mat(cellfun(@(fk) fk(x), fk, 'UniformOutput', false)));
    % gradient construction:

    switch method
        case "exact"
            
            deriv_j = cell(n,1);
            deriv_j{1} = @(x) 1*sin(x(1)) + 2*cos(x(1))   ;
            for k = 2:1:n-1
                deriv_j{k} = @(x) k*sin(x(k)) + 2*cos(x(k))   ;
            end
            deriv_j{n} = @(x) n*sin(x(n)) - (n-1)*cos(x(n))   ;
        
            gradf3 = @(x) cell2mat(cellfun(@(deriv_j) deriv_j(x), deriv_j, 'UniformOutput', false));
        
            %hessian 

            A2 = cell(n,1);   % diagonal
            A2{1} = @(x) cos(x(1)) - 2* sin(x(1));
            for k = 2:1:(n-1) 
                A2{k} = @(x) k*cos(x(k)) - 2*sin(x(k));
            end
            A2{n} = @(x) n*cos(x(n)) + (n-1)*sin(x(n));
        
            Hessf3 = @(x) spdiags (cell2mat(cellfun(@(A2) A2(x), A2, 'UniformOutput', false)),0, n,n);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case "simplified_forward"
            if adaptive == false
                hi = @(x)h;
            else
                hi = @(x)h * abs(x);
            end 

            deriv_j = cell(n,1);
            for i = 1:1:n-1
                deriv_j{i} = @(x) (i*(cos(x(i)) - cos(x(i)+hi(x(i)) )) + 2*(sin(x(i)+hi(x(i)) ) - sin(x(i))))/hi(x(i));
            end
            deriv_j{n} = @(x) ((n-1)*(-sin(x(n)+hi(x(n)) ) + sin(x(n))) + n*(-cos(x(n)+hi(x(n)) ) + cos(x(n))))/hi(x(n));
        
            gradf3 = @(x) cell2mat(cellfun(@(deriv_j) deriv_j(x), deriv_j, 'UniformOutput', false));
        
            %hessian 
            
            if adaptive == false
                hi = @(x) h; %for the hessian instead of h we use sqrt(h)
            else
                hi = @(x) h * abs(x);
            end 
            A2 = cell(n,1);   % diagonal
            for i = 1:1:n-1
                A2{i} = @(x) (-i*(cos(x(i)+2*hi(x(i))) - 2*cos(x(i)+hi(x(i))) + cos(x(i)))...
                             +2*(sin(x(i)+2*hi(x(i))) - 2*sin(x(i)+hi(x(i))) + sin(x(i))))/(hi(x(i))^2);
            end
            A2{n} = @(x) (n*(-cos(x(n)+2*hi(x(n))) + 2*cos(x(n)+hi(x(n))) - cos(x(n)))...
                         +(n-1)*(-sin(x(n)+2*hi(x(n))) + 2*sin(x(n)+hi(x(n))) - sin(x(n))))/(hi(x(n))^2);
        
            Hessf3 = @(x) spdiags (cell2mat(cellfun(@(A2) A2(x), A2, 'UniformOutput', false)),0, n,n);


        case "simplified_centered"

            if adaptive == false
                hi = @(x)h;
            else
                hi = @(x)h * abs(x);
            end 
            deriv_j = cell(n,1);
            deriv_j{1} = @(x) 1*sin(x(1))*(1+hi(x(1))^2) + 2*cos(x(1))*(1+hi(x(1))^2);
            for k = 2:1:n-1
                deriv_j{k} = @(x) (k*sin(x(k)) + 2*cos(x(k)))*(1+hi(x(k))^2)  ;
            end
            deriv_j{n} = @(x) (n*sin(x(n)) - (n-1)*cos(x(n)))*(1+hi(x(n))^2)   ;
        
            gradf3 = @(x) cell2mat(cellfun(@(deriv_j) deriv_j(x), deriv_j, 'UniformOutput', false));
        
            %hessian 

            A2 = cell(n,1);   % diagonal
            A2{1} = @(x) (cos(x(1)) - 2* sin(x(1)))*(1+hi(x(1))^2)^2;
            for k = 2:1:(n-1) 
                A2{k} = @(x) (k*cos(x(k)) - 2*sin(x(k)))*(1+hi(x(k))^2)^2;
            end
            A2{n} = @(x) (n*cos(x(n)) + (n-1)*sin(x(n)))*(1+hi(x(n))^2)^2;
        
            Hessf3 = @(x) spdiags (cell2mat(cellfun(@(A2) A2(x), A2, 'UniformOutput', false)),0, n,n);

        otherwise
            error('Metodo non valido! Usa "exact", o "simplified_forward.');
    end
    
    
end
