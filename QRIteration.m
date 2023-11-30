function [V,Lambda,iter] = QRIteration(A,maxiter,tol)
    n = size(A,1);
    V = eye(n);
    for iter = 1:maxiter
        [Q, R] = qr(A); 
        A = R * Q;
        V = V * Q;
        intol = true;
        for i = 1:n-1
            if abs(A(i+1, i)) > tol
                intol = false;
                break;
            end
        end
        
        if intol
            break;  
        end
    end
    Lambda = diag(A);
end