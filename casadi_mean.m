% CasADi SX variables not compatible with Matlab mean function. This
% implementation fixes the issue
function mean = casadi_mean(numbers)
    import casadi.*
    mean = sum1(numbers)/size(numbers,1);
end