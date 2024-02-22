function [L, U, P] = luDecompositionWithPivoting(A, precision)
    
    [n, ~] = size(A);
    L = cast(eye(n),precision);
    P = cast(eye(n),precision);
    U = cast(A,precision);
    
    for k = 1:n-1
        [~, m] = max(abs(U(k:n, k)));
        m = m + k - 1;

        % 置换U的行
        if m ~= k
            temp = U(k, :);
            U(k, :) = U(m, :);
            U(m, :) = temp;

            % 置换P的行
            temp = P(k, :);
            P(k, :) = P(m, :);
            P(m, :) = temp;

            % 置换L的行（除第k列之外）
            if k >= 2
                temp = L(k, 1:k-1);
                L(k, 1:k-1) = L(m, 1:k-1);
                L(m, 1:k-1) = temp;
            end
        end

        for j = k+1:n
            L(j, k) = U(j, k) / U(k, k);
            U(j, k:n) = U(j, k:n) - L(j, k) * U(k, k:n);
        end
    end
end
