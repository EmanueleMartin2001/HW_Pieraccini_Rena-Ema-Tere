h = 4
adaptive = true
if adaptive == false
    hi = @(x)h;
else
    hi = @(x)h * abs(x);
end

deriv_j = @(x) x(1)+hi(x(2));

z = deriv_j([3,4])