clc,clear,close all
% 设置矩阵维度
n = 1000;
m = n;
% 迭代参数
N = 2000;
tol = 1e-6;

% 设置混合精度参数
single_precision = 'single';
double_precision = 'double';

% A = rand(m, n); % 生成随机矩阵
% x0 = rand(n, 1); % 随机生成精确解
% b = A*x0; % 计算得到右端向量
[A,b,x0] = gravity(n);
alpha = norm(A'*A);
lambda = 0.9/alpha;
epsilon = (1e-4)*randn(m,1);
b1 = b + epsilon;
A2 = cast(A, single_precision);
b2 = cast(b1, single_precision);

%% 双精度
% 初始化
R1 = []; % 储存残差的范数
E1 = []; % 储存相对误差的范数
x1 = ones(n, 1);

tic
for j = 1:N
    r1 = b1 - A*x1;
    %R1(j) = norm(r1);
    R1(j) = log(norm(r1));
    
    d1 = A'*r1;
    x1 = x1 + lambda*d1;
    
    %E1(j) = norm(x1-x0)/norm(x0);
    E1(j) = log(norm(x1-x0)/norm(x0));

    % 检查收敛条件
    if norm(r1) < tol
        break;
    end
end
toc

% 显示结果
subplot(2,3,2);
plot([1:j],R1);
title('双精度下残差图');
subplot(2,3,5);
plot([1:j],E1);
title('双精度下相对误差图');
disp(['双精度下的迭代次数: ', num2str(j)]);
disp(['双精度下的最终残差范数: ', num2str(norm(b-A*x1))]);
disp(['双精度下的最终相对误差范数: ', num2str(norm(x1-x0)/norm(x0))]);

%% 单精度
% 初始化
R2 = []; % 储存残差的范数
E2 = []; % 储存相对误差的范数
x2 = ones(n, 1, single_precision); % 初始解

tic
for p = 1:N
    r2 = b2 - A2*x2;
    %R3(p) = norm(r2);
    R2(p) = log(norm(r2));
    
    d2 = A2'*r2;
    x2 = x2 + lambda*d2;
    
    %E3(p) = norm(x2-x0)/norm(x0);
    E2(p) = log(norm(x2-x0)/norm(x0));

    % 检查收敛条件
    if norm(r2) < tol
        break;
    end
end
toc

% 显示结果
subplot(2,3,3);
plot([1:p],R2);
title('单精度下残差图');
subplot(2,3,6);
plot([1:p],E2);
title('单精度下相对误差图');
disp(['单精度下的迭代次数: ', num2str(p)]);
disp(['单精度下的最终残差范数: ', num2str(norm(b-A*x2))]);
disp(['单精度下的最终相对误差范数: ', num2str(norm(x2-x0)/norm(x0))]);

%% 混合精度
% 初始化
R = []; % 储存残差的范数
E = []; % 储存相对误差的范数
x = ones(n, 1, single_precision); % 初始解

tic
for k = 1:N
    % 使用单精度计算
    Ax = A2*x;
    Ax = double(Ax);
    
    % 使用双精度计算
    r = b1 - Ax;
    %R(k) = norm(r);
    R(k) = log(norm(r));

    % 使用单精度计算
    d = A2'*r;

    % 使用双精度计算
    d = double(d);
    x = double(x);
    x = x + lambda*d;
    %E(k) = norm(x-x0)/norm(x0);
    E(k) = log(norm(x-x0)/norm(x0));

    % 检查收敛条件
    if norm(r) < tol
        break;
    end
end
toc

% 显示结果
subplot(2,3,1);
plot([1:k],R);
title('混合精度下残差图');
subplot(2,3,4);
plot([1:k],E);
title('混合精度下相对误差图');
disp(['混合精度下的迭代次数: ', num2str(k)]);
disp(['混合精度下的最终残差范数: ', num2str(norm(b-A*x))]);
disp(['混合精度下的最终相对误差范数: ', num2str(norm(x-x0)/norm(x0))]);