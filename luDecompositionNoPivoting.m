function [L, U, P] = luDecompositionNoPivoting(A, precision)

    [n, ~] = size(A);
    L = cast(eye(n),precision);
    U = cast(A,precision);
    P = cast(eye(n),precision);

    for k = 1:n-1
        for j = k+1:n
            L(j, k) = U(j, k) / U(k, k);
            U(j, k:n) = U(j, k:n) - L(j, k) * U(k, k:n); 
        end
    end
end
