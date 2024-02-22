function x = solver(L,U,P,b,precision)
    
    L = cast(L,precision);
    U = cast(U,precision);
    P = cast(P,precision);
    b = cast(b,precision);

    % 应用置换矩阵到b
    b = P * b;

    n = length(b);
    y = zeros(n, 1);
    for i = 1:n
        y(i) = (b(i) - L(i, 1:i-1) * y(1:i-1)) / L(i, i);
    end

    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end
