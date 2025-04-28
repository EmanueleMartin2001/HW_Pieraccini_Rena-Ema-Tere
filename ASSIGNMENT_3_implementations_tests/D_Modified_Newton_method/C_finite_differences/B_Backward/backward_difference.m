function [gradf1] = backward_difference(n,h,f)
    
    f_grad_comp = cell(n,1);

    for i=1:n
       
        incremento=zeros(n,1);
        incremento(i)=h;
        f_grad_comp{i} = @(x) (f(x)-f(x-incremento))/h;
    
    end

    gradf1 = @(x) cell2mat(cellfun(@(f_grad_comp) f_grad_comp(x), f_grad_comp, 'UniformOutput', false));


end