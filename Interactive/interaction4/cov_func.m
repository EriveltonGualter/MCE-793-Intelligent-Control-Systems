% Covariance function
function k = cov_func(x, xp, p1, p2)
    k = p1*exp( -(x-xp)'*(x-xp) / (2*p2^2) );
end
