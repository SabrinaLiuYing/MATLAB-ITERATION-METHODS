function EigenMethods()
    n = 100;
    maxiter = 10000;
    tol = 1e-4;
    A  = diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1);

    figure();
    v0 = [1; zeros(n-1, 1)];
    [v, lambda, iter] = PowerIteration(A, v0, maxiter, tol);
    plot(v);
    title(sprintf("Power iteration, lambda = %d, iters = %d", lambda, iter));

    figure();
    v0 = ones(n, 1);
    [v, lambda, iter] = RayleighQuotient(A, v0, maxiter, tol);
    plot(v);
    title(sprintf("RQI, lambda = %d, iters = %d", lambda, iter));

    figure();
    [V, Lambda, iter] = QRIteration(A, maxiter, tol);
    plot(diag(Lambda));
    title("Eigenvalues of A");
    
    indices = [20, 40, 60, 80];
    for i = 1:length(indices)
        ind = indices(i);
        figure();
        plot(V(:, ind));
        title(sprintf("QR iteration, lambda= %d, Index= %d, iters= %d", Lambda(ind),ind,iter));
    end
end