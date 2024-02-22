% [A,b,x] = gravity(100,1,0,4,0.3);
% [A,b,x] = phillips(2500);
[A,b,x] = heat(2500,1);

[n, ~] = size(A);
perturbation = (1e-13)*randn(n,1);
b = b + perturbation;
N = norm(perturbation, 2);
bn = A' * b;
alpha = 10*sqrt(N);

tol = 1e-13;
max_iter = 10;

% r = 2;
% t = 200;
% step = 1;
% M = (r-t)/step + 1;
% R1 = zeros(M);
% R2 = zeros(M);

% count = 1;
% for i = r:step:t
% alpha = i*N;
An = A' * A + alpha*eye(n);

[R1, R2, X1, X2] = deal([], [], [], []);

tic
[L1,U1,P1] = luDecompositionWithPivoting(An,'single');
%[L1,U1,P1] = luDecompositionNoPivoting(An,'single');
x1 = solver(L1,U1,P1,bn,'single');
r1 = bn-An*x1;
X1 = [X1, norm(x1-x)/norm(x)];
R1 = [R1, norm(r1)];
N1 = 0;
while norm(r1)>tol && N1<max_iter
    d1 = solver(L1,U1,P1,r1,'single');
    x1 = x1 + double(d1);
    r1 = bn-An*x1;
    N1 = N1 + 1;
    X1 = [X1, norm(x1-x)/norm(x)];
    R1 = [R1, norm(r1)];
end
toc
% r1 = norm(x1-x, 2);
% R1(count) = r1; 
disp(['混合精度下的迭代次数: ', num2str(N1)]);
disp(['混合精度下的残差范数: ', num2str(norm(r1))]);
disp(['混合精度下的误差范数: ', num2str(norm(x1-x)/norm(x))]);

tic
[L2,U2,P2] = luDecompositionWithPivoting(An,'double');
%[L2,U2,P2] = luDecompositionNoPivoting(An,'double');
x2 = solver(L2,U2,P2,bn,'double');
r2 = bn-An*x2;
X2 = [X2, norm(x2-x)/norm(x)];
R2 = [R2, norm(r2)];
N2 = 0;
while norm(r2)>tol && N2<max_iter
    d2 = solver(L2,U2,P2,r2,'double');
    x2 = x2 + d2;
    r2 = bn-An*x2;
    X2 = [X2, norm(x2-x)/norm(x)];
    R2 = [R2, norm(r2)];
    N2 = N2 + 1;
end
toc
disp(['双精度下的迭代次数: ', num2str(N2)]);
disp(['双精度下的残差范数: ', num2str(norm(r2))]);
disp(['双精度下的误差范数: ', num2str(norm(x2-x)/norm(x))]);

% count = count + 1;
% end

% subplot(1, 2, 1); 
% plot(R1);
% title('混合精度的残差范数');
% subplot(1, 2, 2); 
% plot(X1);
% title('混合精度的相对误差范数');

% subplot(2, 2, 1); 
% plot(R1);
% title('混合精度的残差范数');
% subplot(2, 2, 2); 
% plot(R2);
% title('双精度的残差范数');
% subplot(2, 2, 3); 
% plot(X1);
% title('混合精度的误差范数');
% subplot(2, 2, 4); 
% plot(X2);
% title('双精度的误差范数');

    
