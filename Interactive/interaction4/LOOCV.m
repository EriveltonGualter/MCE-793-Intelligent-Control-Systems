% Leave-one-out cross-validation error
function k = LOOCV(fhat, y)
    k = sum((fhat-y).^2);
end
