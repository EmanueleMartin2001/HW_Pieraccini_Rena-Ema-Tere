function [f1, gradf1, Hessf1] = third_function_3(n,method,h)

    % By default, this function calculates the gradient and the Hessian
    % using the exact method.
    % Otherwise, if it is specified by the user in the input to use a finite
    % difference approximation method, it will be used the corresponding
    % one (between 'centered','forward','backward'). 
    % By default, incrementation h will be assigned to 0.0001

    if nargin == 1
        method = "exact"; % Metodo di default
    end

    if nargin == 2
        h = 0.0001; % Incremento di default
    end
    
    %building the function
    fk = cell(n,1);
    fk{1} = @(x) 1*((1-cos(x(1)) + 0 - sin(x(2))));
    for k = 2:1:n-1
        fk{k} = @(x) k* ( (1-cos(x(k))) + sin(x(k-1)) - sin(x(k+1)) )   ;
    end
    fk{n} = @(x) n* ((1-cos(x(n))) + sin(x(n-1)));

    f1 = @(x) sum(cell2mat(cellfun(@(fk) fk(x), fk, 'UniformOutput', false)));
    % gradient construction:

    switch method
        case "exact"
            
            f3=f1;
            
            deriv_j = cell(n,1);
            deriv_j{1} = @(x) 1*sin(x(1)) + 2*cos(x(1))   ;
            for k = 2:1:n-1
                deriv_j{k} = @(x) k*sin(x(k)) + 2*cos(x(k))   ;
            end
            deriv_j{n} = @(x) n*sin(x(n)) - (n-1)*cos(x(n))   ;
        
            gradf1 = @(x) cell2mat(cellfun(@(deriv_j) deriv_j(x), deriv_j, 'UniformOutput', false));
        
        
        
            A2 = cell(n,1);   % diagonal
            A2{1} = @(x) cos(x(1)) - 2* sin(x(1));
            for k = 1:1:n
                A2{k} = @(x) k*cos(x(k)) - 2*sin(x(k));
            end
            A2{n} = @(x) n*cos(x(n)) + (n-1)*sin(x(n));
        
            Hessf1 = @(x) spdiags (cell2mat(cellfun(@(A2) A2(x), A2, 'UniformOutput', false)),0, n,n);


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % The computation of the gradient whit the finite method will use
        % the same steps for every input function given.
        % This is not possible for the computation of the Hessian: due to
        % the sparsity of the Hessian matrix and the its big dimension, it
        % is necessary to leverage on the features of each given function
        % to minimize the computational times and storage memory
        
        
        case "forward"
            gradf1 = forward_difference(n, h, f1);
            
            % calcolo Hessiana 
            
            % computing the Hessian in sparse mode
    
            A1 = cell(n,1);
           
            for i=1:n
                             
                incremento_i=zeros(n,1);
                
                incremento_i(i)=h;
                
                
                A1{i}=@(x) (f1(x+incremento_i+incremento_i)-f1(x+incremento_i)-f1(x+incremento_i)+f1(x))/(h^2);
                
            end
                       
            
            Hessf1 = @(x) spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n);
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case "centered"
            gradf1 = centered_difference(n, h, f1);

            % computing the Hessian in sparse mode
    
            A1 = cell(n,1);
            
            for i=1:n
                             
                incremento_i=zeros(n,1);
                incremento_i(i)=h;
                
                A1{i}=@(x) (f1(x+2*incremento_i)-2*f1(x)+f1(x-2*incremento_i))/(4*(h^2));
                         
            end          
            
            Hessf1 = @(x) spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
        case 'backward'
        gradf1 = centered_difference(n, h, f1);

        % calcolo Hessiana 
            
            % computing the Hessian in sparse mode
    
            A1 = cell(n,1);
            
            for i=1:n
                             
                incremento_i=zeros(n,1);
                
                incremento_i(i)=h;
                
                A1{i}=@(x) (f1(x-incremento_i-incremento_i)-f1(x-incremento_i)-f1(x-incremento_i)+f1(x))/(h^2);
                
            end
            
            Hessf1 = @(x) spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n);
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        
        otherwise
            error('Metodo non valido! Usa "backward", "forward", o "centered".');
    end
    
    
end
