function p = convergence_rate(xseq)

    % this function evaluate the rate fo convergence

    if (size(xseq,2)>2)
        
        upper = xseq(:,end) - xseq(:,end-1);
        lower = xseq(:,end-1) - xseq(:,end-2);

        p = log(sum(upper.^2))/log(sum(lower.^2));

        %p = log(norm(upper))/log(norm(lower));

        % we didn't compute the norm becouse computing the sq is costly

    else
        p = 0;
    end

end