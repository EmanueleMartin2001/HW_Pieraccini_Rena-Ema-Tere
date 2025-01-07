function [f1, gradf1, Hessf1] = first_function_1(n,method,h)

    % By default, this function calculates the gradient and the Hessian
    % using the exact method.
    % Otherwise, if it is specified by the user in the input to use a finite
    % difference approximation method, it will be used the corresponding
    % one (between 'centered','forward','backward'). 
    % By default, incrementation h will be assigned to 0.0001

    if nargin == 1
        method = 'exact'; % Metodo di default
    end

    if nargin == 2
        h = 0.0001; % Incremento di default
    end
    
    
    f1 = @(x) sum(100 * (x(1:end-1).^2 -x(2:end)).^2 + (x(1:end-1)-1).^2);

    % gradient construction:

    switch method
        case 'exact'
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
                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % The computation of the gradient whit the finite method will use
        % the same steps for every input function given.
        % This is not possible for the computation of the Hessian: due to
        % the sparsity of the Hessian matrix and the its big dimension, it
        % is necessary to leverage on the features of each given function
        % to minimize the computational times and storage memory
        
        
        case 'forward'
            gradf1 = forward_difference(n, h, f1);
            
            % calcolo Hessiana 
            
            % computing the Hessian in sparse mode
    
            A1 = cell(n,1);
            B1 = cell(n,1);
            C1 = cell(n,1);
            
            incremento_primo=zeros(1, n);
            incremento_primo1=zeros(1, n);

            incremento_primo(1)=h;
            incremento_primo1(2)=h;

            
            A1{1}=@(x) (f1(x+incremento_primo+incremento_primo)-f1(x+incremento_primo)-f1(x+incremento_primo)+f1(x))/(h^2);
            B1{1}=@(x) (f1(x+incremento_primo+incremento_primo1)-f1(x+incremento_primo)-f1(x+incremento_primo1)+f1(x))/(h^2);
            C1{1} = @(x) 0; %forced by spdiags

            for i=2:n-1
                             
                incremento_i=zeros(1, n);
                incremento_i0=zeros(1, n);
                incremento_i1=zeros(1, n);

                incremento_i(i)=h;
                incremento_i0(i-1)=h;
                incremento_i1(i+1)=h;

                
                A1{i}=@(x) (f1(x+incremento_i+incremento_i)-f1(x+incremento_i)-f1(x+incremento_i)+f1(x))/(h^2);
                B1{i}=@(x) (f1(x+incremento_i+incremento_i1)-f1(x+incremento_i)-f1(x+incremento_i1)+f1(x))/(h^2);
                C1{i}=@(x) (f1(x+incremento_i+incremento_i0)-f1(x+incremento_i)-f1(x+incremento_i0)+f1(x))/(h^2);

                            
            end
            
            
            incremento_ultimo=zeros(1, n);
            incremento_ultimo1=zeros(1, n);

            incremento_ultimo(n)=h;
            incremento_ultimo1(n-1)=h;
             
            A1{n} = @(x) (f1(x+incremento_ultimo+incremento_ultimo)-f1(x+incremento_ultimo)-f1(x+incremento_ultimo)+f1(x))/(h^2);
            B1{n} = @(x) 0; %forced by spdiags
            C1{n} = @(x) (f1(x+incremento_ultimo+incremento_ultimo1)-f1(x+incremento_ultimo)-f1(x+incremento_ultimo1)+f1(x))/(h^2);
            
            Hessf1 = @(x) spdiags (cell2mat(cellfun(@(B1) B1(x), B1, 'UniformOutput', false)),-1, n,n)+ ...
                            spdiags (cell2mat(cellfun(@(A1) A1(x), A1, 'UniformOutput', false)),0, n,n)+...
                            spdiags (cell2mat(cellfun(@(C1) C1(x), C1, 'UniformOutput', false)),+1, n,n);
                
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'centered'
            gradf1 = centered_difference(n, h, f1);

            % computing the Hessian in sparse mode
    
            A1 = cell(n,1);
            B1 = cell(n,1);
            C1 = cell(n,1);
            
            incremento_primo=zeros(1, n);
            incremento_primo1=zeros(1, n);

            incremento_primo(1)=h;
            incremento_primo1(2)=h;

            
            A1{1}=@(x) (f1(x-incremento_primo-incremento_primo)-f1(x))/(4*(h^2));
            B1{1}=@(x) (f1(x+incremento_primo+incremento_primo1)-f1(x+incremento_primo-incremento_primo1)-f1(x-incremento_primo+incremento_primo1)+f1(x-incremento_primo-incremento_primo1))/(4*(h^2));
            C1{1} = @(x) 0; %forced by spdiags

            for i=2:n-1
                             
                incremento_i=zeros(1, n);
                incremento_i0=zeros(1, n);
                incremento_i1=zeros(1, n);

                incremento_i(i)=h;
                incremento_i0(i-1)=h;
                incremento_i1(i+1)=h;

                
                A1{i}=@(x) (f1(x-incremento_i-incremento_i)-f1(x))/(4*(h^2));
                B1{i}=@(x) (f1(x+incremento_i+incremento_i1)+f1(x-incremento_i-incremento_i1)-f1(x-incremento_i+incremento_i1)-f1(x+incremento_i-incremento_i1))/(4*(h^2));
                C1{i}=@(x) (f1(x+incremento_i+incremento_i0)+f1(x-incremento_i-incremento_i0)-f1(x-incremento_i+incremento_i0)-f1(x+incremento_i-incremento_i0))/(4*(h^2));

                            
            end
            
            
            incremento_ultimo=zeros(1, n);
            incremento_ultimo1=zeros(1, n);

            incremento_ultimo(n)=h;
            incremento_ultimo1(n-1)=h;
             
            A1{n} = @(x) (f1(x-incremento_ultimo-incremento_ultimo)--f1(x+incremento_ultimo)+f1(x))/(4*(h^2));
            B1{n} = @(x) 0; %forced by spdiags
            C1{n} = @(x) (f1(x+incremento_ultimo+incremento_ultimo1)+f1(x-incremento_ultimo-incremento_ultimo1)-f1(x-incremento_ultimo+incremento_ultimo1)--f1(x+incremento_ultimo-incremento_ultimo1))/(4*(h^2));
            
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
            
            incremento_primo=zeros(1, n);
            incremento_primo1=zeros(1, n);

            incremento_primo(1)=h;
            incremento_primo1(2)=h;

            
            A1{1}=@(x) (f1(x-incremento_primo-incremento_primo)-f1(x-incremento_primo)-f1(x-incremento_primo)+f1(x))/(h^2);
            B1{1}=@(x) (f1(x-incremento_primo-incremento_primo1)-f1(x-incremento_primo)-f1(x-incremento_primo1)+f1(x))/(h^2);
            C1{1} = @(x) 0; %forced by spdiags

            for i=2:n-1
                             
                incremento_i=zeros(1, n);
                incremento_i0=zeros(1, n);
                incremento_i1=zeros(1, n);

                incremento_i(i)=h;
                incremento_i0(i-1)=h;
                incremento_i1(i+1)=h;

                
                A1{i}=@(x) (f1(x-incremento_i-incremento_i)-f1(x-incremento_i)-f1(x-incremento_i)+f1(x))/(h^2);
                B1{i}=@(x) (f1(x-incremento_i-incremento_i1)-f1(x-incremento_i)-f1(x-incremento_i1)+f1(x))/(h^2);
                C1{i}=@(x) (f1(x-incremento_i-incremento_i0)-f1(x-incremento_i)-f1(x-incremento_i0)+f1(x))/(h^2);

                            
            end
            
            
            incremento_ultimo=zeros(1, n);
            incremento_ultimo1=zeros(1, n);

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
