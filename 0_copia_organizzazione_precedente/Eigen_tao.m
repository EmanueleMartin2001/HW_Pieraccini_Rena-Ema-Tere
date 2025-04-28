function tao = Eigen_tao(A,delta)

% This is the case where we use eigen

[eigen_vect, eigen_vals] = eigs(A);

lamda_min = min(diag(eigen_vals));

tao = max(0,delta-lamda_min);

end

