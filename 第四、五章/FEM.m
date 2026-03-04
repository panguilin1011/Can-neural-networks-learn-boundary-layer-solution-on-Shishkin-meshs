function [u_fem, points] = FEM(epsilon,N)

sigma = min(0.5, 2*epsilon*log(N)); % Shishkin过渡点参数

% Shishkin网格生成（在 x=1 处加密）
x = zeros(N+1,1);
for i = 0:N
    if i <= N/2
        x(i+1) = i*(1 - sigma)/(N/2);
    else
        x(i+1) = 1 - sigma + (i - N/2)*sigma/(N/2);
    end
end

% 初始化矩阵和向量
A = zeros(N+1,N+1);
b = zeros(N+1,1);

% 有限元组装
for i = 1:N
    h = x(i+1) - x(i);
    Ae = epsilon/h * [1 -1; -1 1] + [-1/2 1/2; -1/2 1/2];
    be = h/2 * [1; 1];

    A(i:i+1, i:i+1) = A(i:i+1, i:i+1) + Ae;
    b(i:i+1) = b(i:i+1) + be;
end

% 边界条件
A(1,:) = 0; A(1,1) = 1; b(1) = 0;
A(end,:) = 0; A(end,end) = 1; b(end) = 0;

% 求解线性系统
u_fem = A \ b;
points = x;

%% ===== 计算误差 ===== %%
% 解析解及其导数
u_exact = @(x) (-exp(-1/epsilon) + exp((x - 1)/epsilon)) / (exp(-1/epsilon) - 1) + x;
u_prime = @(x) exp((x - 1)/epsilon) / (epsilon*(exp(-1/epsilon) - 1)) + 1;

% L2误差（梯形法则）
h = diff(x); % 单元长度
e_nodes = u_exact(x) - u_fem; % 节点误差

L2_sq = 0;
for i = 1:N
    L2_sq = L2_sq + (e_nodes(i)^2 + e_nodes(i+1)^2)/2 * h(i);
end
L2_error = sqrt(L2_sq);

% H1误差（两点高斯积分）
H1_semi = 0;
for i = 1:N
    h_i = x(i+1) - x(i);
    u_fem_prime = (u_fem(i+1) - u_fem(i))/h_i; % FEM导数
    
    % 高斯点坐标和权重
    gauss_pts = x(i) + h_i/2*[1 - 1/sqrt(3); 1 + 1/sqrt(3)];
    weights = [h_i/2, h_i/2];
    
    % 计算高斯点处导数误差
    der_error = u_fem_prime - u_prime(gauss_pts);
    H1_semi = H1_semi + sum(weights .* der_error.^2);
end
H1_error = sqrt(L2_error^2 + epsilon*H1_semi);

fprintf('L2误差: %.2e\n', L2_error);
fprintf('H1误差: %.2e\n', H1_error);
end