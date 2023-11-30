function [v, lambda, iter] = PowerIteration(A, v0, maxiter, tol)
    v = v0 / norm(v0);
    lambda = v' * A * v;
    r = A * v - lambda * v;
    for iter = 1:maxiter
        if (norm(r)< tol)
            break;
        end
        w = A * v;
        v = w / norm(w);
        lambda = v' * A * v;
        r = A * v - lambda * v;
    end
    v = v / norm(v);
end