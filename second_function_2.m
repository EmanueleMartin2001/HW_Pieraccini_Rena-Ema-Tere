function [f1, gradf1, Hessf1] = second_function_2(n,method,h)

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
    fk{1} = @(x) (3-2*x(1))*x(1) -2*x(2) + 1;
    for k = 2:1:n-1
        fk{k} = @(x) (3-2*x(k))*x(k) -x(k-1) -2*x(k+1) + 1;
    end
    fk{n} = @(x) (3-2*x(n))*x(n) -x(n-1) + 1;
    
    f1 = @(x) 0.5*sum(cell2mat(cellfun(@(fk) fk(x), fk, 'UniformOutput', false)).^2);

    % gradient construction:

    switch method
        case "exact"
            
            f2=f1;
            
            deriv_j = cell(n,1);
            deriv_j{1} = @(x) fk{1}(x) * (3-4*x(1)) + fk{2}(x) * (-2);
            for k = 2:1:n-1
                deriv_j{k} = @(x)  fk{k-1}(x) * (-1) + fk{k}(x) * (3-4*x(k)) + fk{k+1}(x) * (-2);
            end
            deriv_j{n} = @(x)  fk{n-1}(x) * (-1) + fk{n}(x) * (3-4*x(n));
            
            gradf1 = @(x) cell2mat(cellfun(@(deriv_j) deriv_j(x), deriv_j, 'UniformOutput', false));
            
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
                
                Hessf1 = @(x) spdiags (cell2mat(cellfun(@(B1) B1(x), C1, 'UniformOutput', false)),-1, n,n)+ ...
                                spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n)+...
                                spdiags (cell2mat(cellfun(@(C1) C1(x), B1, 'UniformOutput', false)),+1, n,n);

        
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
            B1 = cell(n,1);
            C1 = cell(n,1);
            
            incremento_primo=zeros(n,1);
            incremento_primo1=zeros(n,1);

            incremento_primo(1)=h;
            incremento_primo1(2)=h;

            
            A1{1}=@(x) (f1(x+incremento_primo+incremento_primo)-f1(x+incremento_primo)-f1(x+incremento_primo)+f1(x))/(h^2);
            B1{1}=@(x) (f1(x+incremento_primo+incremento_primo1)-f1(x+incremento_primo)-f1(x+incremento_primo1)+f1(x))/(h^2);
            C1{1} = @(x) 0; %forced by spdiags

            for i=2:n-1
                             
                incremento_i=zeros(n,1);
                incremento_i0=zeros(n,1);
                incremento_i1=zeros(n,1);

                incremento_i(i)=h;
                incremento_i0(i-1)=h;
                incremento_i1(i+1)=h;

                
                A1{i}=@(x) (f1(x+incremento_i+incremento_i)-f1(x+incremento_i)-f1(x+incremento_i)+f1(x))/(h^2);
                B1{i}=@(x) (f1(x+incremento_i+incremento_i1)-f1(x+incremento_i)-f1(x+incremento_i1)+f1(x))/(h^2);
                C1{i}=@(x) (f1(x+incremento_i+incremento_i0)-f1(x+incremento_i)-f1(x+incremento_i0)+f1(x))/(h^2);

                            
            end
            
            
            incremento_ultimo=zeros(n,1);
            incremento_ultimo1=zeros(n,1);

            incremento_ultimo(n)=h;
            incremento_ultimo1(n-1)=h;
             
            A1{n} = @(x) (f1(x+incremento_ultimo+incremento_ultimo)-f1(x+incremento_ultimo)-f1(x+incremento_ultimo)+f1(x))/(h^2);
            B1{n} = @(x) 0; %forced by spdiags
            C1{n} = @(x) (f1(x+incremento_ultimo+incremento_ultimo1)-f1(x+incremento_ultimo)-f1(x+incremento_ultimo1)+f1(x))/(h^2);
            
            Hessf1 = @(x) spdiags (cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false)),-1, n,n)+ ...
                            spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n)+...
                            spdiags (cell2mat(cellfun(@(C1) C1(x), C1, 'UniformOutput', false)),+1, n,n);
                
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case "centered"
            gradf1 = centered_difference(n, h, f1);

            % computing the Hessian in sparse mode
    
            A1 = cell(n,1);
            B1 = cell(n,1);
            C1 = cell(n,1);
            
            incremento_primo=zeros(n,1);
            incremento_primo1=zeros(n,1);

            incremento_primo(1)=h;
            incremento_primo1(2)=h;

            
            A1{1}=@(x) (f1(x+2*incremento_primo)-2*f1(x)+f1(x-2*incremento_primo))/(4*(h^2));
            B1{1}=@(x) (f1(x+incremento_primo+incremento_primo1)-f1(x+incremento_primo-incremento_primo1)-f1(x-incremento_primo+incremento_primo1)+f1(x-incremento_primo-incremento_primo1))/(4*(h^2));
            C1{1} = @(x) 0; %forced by spdiags

            for i=2:n-1
                             
                incremento_i=zeros(n,1);
                incremento_i0=zeros(n,1);
                incremento_i1=zeros(n,1);

                incremento_i(i)=h;
                incremento_i0(i-1)=h;
                incremento_i1(i+1)=h;

                
                A1{i}=@(x) (f1(x+2*incremento_i)-2*f1(x)+f1(x-2*incremento_i))/(4*(h^2));
                B1{i}=@(x) (f1(x+incremento_i+incremento_i1)+f1(x-incremento_i-incremento_i1)-f1(x-incremento_i+incremento_i1)-f1(x+incremento_i-incremento_i1))/(4*(h^2));
                C1{i}=@(x) (f1(x+incremento_i+incremento_i0)+f1(x-incremento_i-incremento_i0)-f1(x-incremento_i+incremento_i0)-f1(x+incremento_i-incremento_i0))/(4*(h^2));

                            
            end
            
            
            incremento_ultimo=zeros(n,1);
            incremento_ultimo1=zeros(n,1);

            incremento_ultimo(n)=h;
            incremento_ultimo1(n-1)=h;
             
            A1{n} = @(x) (f1(x+2*incremento_ultimo)-2*f1(x)+f1(x-2*incremento_ultimo))/(4*(h^2));
            B1{n} = @(x) 0; %forced by spdiags
            C1{n} = @(x) (f1(x+incremento_ultimo+incremento_ultimo1)+f1(x-incremento_ultimo-incremento_ultimo1)-f1(x-incremento_ultimo+incremento_ultimo1)-f1(x+incremento_ultimo-incremento_ultimo1))/(4*(h^2));
            
            Hessf1 = @(x) spdiags (cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false)),-1, n,n)+ ...
                            spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n)+...
                            spdiags (cell2mat(cellfun(@(C1) C1(x), C1, 'UniformOutput', false)),+1, n,n);
                
            


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
        case 'backward'
        gradf1 = centered_difference(n, h, f1);

        % calcolo Hessiana 
            
            % computing the Hessian in sparse mode
    
            A1 = cell(n,1);
            B1 = cell(n,1);
            C1 = cell(n,1);
            
            incremento_primo=zeros(n,1);
            incremento_primo1=zeros(n,1);

            incremento_primo(1)=h;
            incremento_primo1(2)=h;

            
            A1{1}=@(x) (f1(x-incremento_primo-incremento_primo)-f1(x-incremento_primo)-f1(x-incremento_primo)+f1(x))/(h^2);
            B1{1}=@(x) (f1(x-incremento_primo-incremento_primo1)-f1(x-incremento_primo)-f1(x-incremento_primo1)+f1(x))/(h^2);
            C1{1} = @(x) 0; %forced by spdiags

            for i=2:n-1
                             
                incremento_i=zeros(n,1);
                incremento_i0=zeros(n,1);
                incremento_i1=zeros(n,1);

                incremento_i(i)=h;
                incremento_i0(i-1)=h;
                incremento_i1(i+1)=h;

                
                A1{i}=@(x) (f1(x-incremento_i-incremento_i)-f1(x-incremento_i)-f1(x-incremento_i)+f1(x))/(h^2);
                B1{i}=@(x) (f1(x-incremento_i-incremento_i1)-f1(x-incremento_i)-f1(x-incremento_i1)+f1(x))/(h^2);
                C1{i}=@(x) (f1(x-incremento_i-incremento_i0)-f1(x-incremento_i)-f1(x-incremento_i0)+f1(x))/(h^2);

                            
            end
            
            
            incremento_ultimo=zeros(n,1);
            incremento_ultimo1=zeros(n,1);

            incremento_ultimo(n)=h;
            incremento_ultimo1(n-1)=h;
             
            A1{n} = @(x) (f1(x-incremento_ultimo-incremento_ultimo)-f1(x-incremento_ultimo)-f1(x-incremento_ultimo)+f1(x))/(h^2);
            B1{n} = @(x) 0; %forced by spdiags
            C1{n} = @(x) (f1(x-incremento_ultimo-incremento_ultimo1)-f1(x-incremento_ultimo)-f1(x-incremento_ultimo1)+f1(x))/(h^2);
            
            Hessf1 = @(x) spdiags (cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false)),-1, n,n)+ ...
                            spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n)+...
                            spdiags (cell2mat(cellfun(@(C1) C1(x), C1, 'UniformOutput', false)),+1, n,n);
                
      
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        
        otherwise
            error('Metodo non valido! Usa "backward", "forward", o "centered".');
    end
    
    
end
