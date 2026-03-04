function [u, x] = FEM_SUPG(epsilon, N)
format long
% epsilon - 扩散参数
% N - 网格点数（偶数）
%% 构建Shishkin网格
sigma = min(0.5, 2*epsilon*log(N));
x = zeros(N+1,1);
for i=0:N
    if i <= N/2
        x(i+1) = 2*i*(1 - sigma)/N;
    else
        x(i+1) = 1 - sigma + 2*(i - N/2)*sigma/N;
    end
end

%% 初始化全局刚度矩阵和载荷向量
A = zeros(N+1,N+1);
b = zeros(N+1,1);

%% SUPG参数（简化为h/2）
for i=1:N
    h = x(i+1)-x(i);
    tau = h / 2; % 简化稳定参数
    
    %% 元刚度矩阵和载荷向量（线性有限元+SUPG）
    A_local = epsilon/h*[1 -1;-1 1] + 1/2*[-1 1;-1 1] + tau*[1 -1;-1 1]/h;
    b_local = h/2*[1;1] + tau*[-1;1];
    
    %% 组装到全局矩阵
    A(i:i+1,i:i+1) = A(i:i+1,i:i+1) + A_local;
    b(i:i+1) = b(i:i+1) + b_local;
end

%% 施加Dirichlet边界条件
A(1,:) = 0; A(:,1)=0; A(1,1)=1; b(1)=0;
A(end,:) = 0; A(:,end)=0; A(end,end)=1; b(end)=0;

%% 求解线性方程组
u = A\b;

%% ===== 计算误差（按参考代码方法）===== %%
% 解析解及其导数
u_exact = @(x) (-exp(-1/epsilon) + exp((x - 1)/epsilon)) / (exp(-1/epsilon) - 1) + x;
u_prime = @(x) exp((x - 1)/epsilon) / (epsilon*(exp(-1/epsilon) - 1)) + 1;

% L2误差（梯形法则）
h = diff(x); % 单元长度
e_nodes = u_exact(x) - u; % 节点误差

L2_sq = 0;
for i = 1:N
    L2_sq = L2_sq + (e_nodes(i)^2 + e_nodes(i+1)^2)/2 * h(i);
end
L2_error = sqrt(L2_sq);

% H1误差（两点高斯积分）
H1_semi = 0;
for i = 1:N
    h_i = x(i+1) - x(i);
    u_num_prime = (u(i+1) - u(i))/h_i; % 数值解导数
    
    % 高斯点坐标和权重
    gauss_pts = x(i) + h_i/2*[1 - 1/sqrt(3); 1 + 1/sqrt(3)];
    weights = [h_i/2, h_i/2];
    
    % 计算高斯点处导数误差
    der_error = u_num_prime - u_prime(gauss_pts);
    H1_semi = H1_semi + sum(weights .* der_error.^2);
end
H1_error = sqrt(L2_error^2 + epsilon*H1_semi);

fprintf('L2误差: %.4e\n', L2_error);
fprintf('H1误差: %.4e\n', H1_error);

%% 绘图
% figure;
% plot(x,u, x, u_exact);
% xlabel('x'); ylabel('u(x)');
% legend('数值解','解析解','Location', 'northwest');
% title(['SUPG FEM on Shishkin grid, \epsilon=', num2str(epsilon)]);
% grid on;
% 
% %% 误差图
% figure;
% plot(x, abs(u - u_exact), 'LineWidth', 1.5);
% xlabel('x');
% ylabel('|u - u_{exact}|');
% title(['误差图 (SUPG FEM 与解析解误差), \epsilon=', num2str(epsilon)]);
% grid on;
% 
% 

end