f_evaluation = [3,67,1,6,11,43,0]
m = length(f_evaluation)

[f_evaluation,new_indices] = sort(f_evaluation)

values_indices(:,1) = f_evaluation
values_indices(:,2) = new_indices
values_indices(m,1) = 5 %50 è l f_xr
sortrows(values_indices,1)

X = [2,3; 4,15; 5,10]
mean(X([1,3],:),1)