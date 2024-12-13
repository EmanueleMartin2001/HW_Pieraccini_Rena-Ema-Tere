function tao = Gershgorin_approx(A,delta)

% the goal is to find a circle where inside is contained the minimum
% eigenvalue

% delta = sqrt(epslon)

% m represents the center of the circle
[m,i] = min(diag(A));

% find the radius 
r = norm(A(i,:),1)-abs(m);

approx_eigen = m-r;

tao = max(0,delta-approx_eigen);

end

