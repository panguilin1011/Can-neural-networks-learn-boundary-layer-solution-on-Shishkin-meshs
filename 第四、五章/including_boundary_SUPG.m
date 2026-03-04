clear
format long

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%参数

epsilon = 1e-10;  %扩散系数
N = 20;         %剖分步数
N_iter = 3e5; %迭代次数
eta = 1e-5;    %学习率
beta =0;       %正则化参数
gamma = 0.005;      %罚参数
initial =true;    %是否给初值

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%用于计算的参数

%分界点（过渡点）
tau = Tau(epsilon,N);  

%步长
h = 2*(1-tau)/N;    %(0,1-tau)上的步长
H = 2*tau/N;        %(1-tau,1)上的步长

%节点
x(1:N/2+1) = linspace(0,1-tau,N/2+1);  %(0,1-tau)上的节点坐标
x(N/2+1:N+1) = linspace(1-tau,1,N/2+1);  %(1-tau,1)上的节点坐标

%区间长度
L = diff(x);

%解析解
u_exact = @(x) (-exp(-1/epsilon)+exp((x-1)./epsilon))/(exp(-1/epsilon)-1) + x;
%解析解求导
u_prime = @(x) exp((x-1)./epsilon)/(epsilon*(exp(-1/epsilon)-1)) + 1;
%解析解算的节点值
u_e = u_exact(x);

%大步长区间参数
betah_1 = -(epsilon/h+1);
betah0 = 2*epsilon/h+1;
betah1 = -epsilon/h;

%小步长区间参数
betaH_1 = -(epsilon/H+1);
betaH0 = 2*epsilon/H+1;
betaH1 = -epsilon/H;

%随机数选取（选取为默认种子）
rng('default')

%初始参数选取
F_value = cell(1, N+1);
M = cell(1, N+1);
sigma = cell(1, N+1);
sigma_der = cell(1, N+1);
b2 = rand(3*(N-1), 1);
W2 = rand(3*(N-1), 1);
W3 = rand(1, 3*(N-1));

%是否给初值
if initial == true
    for j = 0:N/2-1
        W2(3*j+1:3*j+3) = [1/h, 2/h, 1/h];
        b2(3*j+1:3*j+3) = [-x(j+1)/h,-x(j+2)*(2/h),-x(j+3)/h];
    end
    W2(3*N/2+1:3*N/2+3) = [1/h,1/h+1/H,1/H];
    b2(3*N/2+1:3*N/2+3) = [-x(N/2+1)/h,-x(N/2+2)*(1/h+1/H),-x(N/2+3)/H];
    for j = N/2+1:N-2
        W2(3*j+1:3*j+3) = [1/H, 2/H, 1/H];
        b2(3*j+1:3*j+3) = [-x(j+1)/H,-x(j+2)*(2/H),-x(j+3)/H];
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%迭代循环
for cycle = 1:N_iter
        
    for j = 1:N+1
        F_value{j} = F(x(j), b2, W2, W3);
    end
    
    for j = 2:N/2
        M{j} = betah_1*F_value{j-1} + ...
        betah0*F_value{j} + ...
        betah1*F_value{j+1} - h;
    end
    M{N/2+1} = betah_1*F_value{N/2} + ...
        (betah0+betaH0)/2*F_value{N/2+1} + ...
        betaH1*F_value{N/2+2} - h;
    for j = N/2+2:N
        M{j} = betaH_1*F_value{j-1} + ...
            betaH0*F_value{j} + ...
            betaH1*F_value{j+1} - H;
    end
    
    for j = 1:N+1
        z = W2*x(j) + b2;
        sigma{j} = relu(z);
        sigma_der{j} = relu_der(z);
    end

%%
%若在不提供任何信息时不对这一部分进行注释
%若在提供足够多的信息时对这一部分进行注释

%     db2 = (F_value{1}*sigma_der{1} + F_value{N+1}*sigma_der{N+1})*gamma;
%     for j = 2:N/2
%         db2 = db2 + M{j}*(...
%             betah_1*sigma_der{j-1} + ...
%             betah0*sigma_der{j} + ...
%             betah1*sigma_der{j+1});
%     end
%     db2 = db2 + M{N/2+1}*(...
%             betah_1*sigma_der{N/2} + ...
%             (betah0+betaH0)/2*sigma_der{N/2+1} + ...
%             betaH1*sigma_der{N/2+2});
%      for j = N/2+2:N
%         db2 = db2 + M{j}*(...
%             betaH_1*sigma_der{j-1} + ...
%             betaH0*sigma_der{j} + ...
%             betaH1*sigma_der{j+1});
%     end
%     db2 = 2*W3'.*db2;
%     db2 = db2 + beta*b2;
% 
%     dw2 = F_value{N+1}*sigma_der{N+1}*x(N+1)*gamma;
%     for j = 2:N/2
%         dw2 = dw2 + M{j}*(...
%             betah_1*sigma_der{j-1}*x(j-1) + ...
%             betah0*sigma_der{j}*x(j) + ...
%             betah1*sigma_der{j+1}*x(j+1));
%     end
%     dw2 = dw2 + M{N/2+1}*(...
%             betah_1*sigma_der{N/2}*x(N/2) + ...
%             (betah0+betaH0)/2*sigma_der{N/2+1}*x(N/2+1) + ...
%             betaH1*sigma_der{N/2+2}*x(N/2+2));
%      for j = N/2+2:N
%         dw2 = dw2 + M{j}*(...
%             betaH_1*sigma_der{j-1}*x(j-1) + ...
%             betaH0*sigma_der{j}*x(j) + ...
%             betaH1*sigma_der{j+1}*x(j+1));
%     end
%     dw2 = 2*W3'.*dw2;
%     dw2 = dw2 + beta*W2;
%%

    dw3 = (F_value{1}*sigma{1} + F_value{N+1}*sigma{N+1})'*gamma;
    for j = 2:N/2
        dw3 = dw3 + (M{j}*(...
            betah_1*sigma{j-1} + ...
            betah0*sigma{j} + ...
            betah1*sigma{j+1}))';
    end
    dw3 = dw3 + (M{N/2+1}*(...
            betah_1*sigma{N/2} + ...
            (betah0+betaH0)/2*sigma{N/2+1} + ...
            betaH1*sigma{N/2+2}))';
     for j = N/2+2:N
        dw3 = dw3 + (M{j}*(...
            betaH_1*sigma{j-1} + ...
            betaH0*sigma{j} + ...
            betaH1*sigma{j+1}))';
    end
    dw3 = 2*dw3;
    dw3 = dw3 + beta*W3;

%     b2 = b2 - eta*db2;
%     W2 = W2 - eta*dw2;
    W3 = W3 - eta*dw3;

    newcost = cost(F_value,M,N,b2,W2,W3,beta,gamma);
    cost_value(cycle) = newcost;
    
    %L2
    e_nodes = u_exact(x) - cell2mat(F_value); % 节点误差

    L2_sq = 0;
    for i = 1:N
      L2_sq = L2_sq + (e_nodes(i)^2 + e_nodes(i+1)^2)/2 * L(i);
    end
    L_2(cycle) = sqrt(L2_sq);

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%数值解

%最后的迭代结果计算的NN节点值
F_values = F(x, b2, W2, W3);

%有限元解节点值
u_fem = FEM_SUPG(epsilon,N);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%绘图
figure(1);
hold on;                     % 先 hold on，下面所有都画在同一张图上

% 1. 精确解（要进图例）
h1 = plot(x, u_e, 'r-', 'LineWidth', 1.5);

% 2. 竖直虚线（不想进图例，关键是 HandleVisibility='off'）
for i = 1:length(x)
    xline(x(i), ':k', 'HandleVisibility', 'off');
end

% 3. NN 曲线（要进图例）
plotting_grid = linspace(0, 1, 10000);
h2 = plot(plotting_grid, F(plotting_grid, b2, W2, W3), ...
          'b--', 'LineWidth', 1.5);

hold off;

% 4. 图例显式绑定两条曲线的句柄
legend([h1 h2], {'精确解', 'NN'}, 'Location', 'northwest');
box on




%绘制cost和L2误差随迭代次数的变化图
figure(2)
yyaxis left;
semilogy(1:N_iter,cost_value,'r-','LineWidth',1.5);
ylabel('cost');
yyaxis right;
semilogy(1:N_iter,L_2,'b-','LineWidth',1.5);
ylabel('L2');
legend('Cost','L^2误差','Location','northeast');

%绘制NN解和解析解的误差图
error = abs(u_e-F_values); %NN解和解析解的绝对值误差
figure(3)
plot(x,error,'LineWidth', 1.5);
legend('error', 'Location', 'northwest');

%绘制NN解与有限元解的误差图
error1 = abs(u_fem'-F_values);
figure(4)
plot(x,error1,'LineWidth', 1.5);
legend('节点绝对值误差', 'Location', 'northwest');

%绘制1-tau到1区间放大的逼近情况
figure(5);
subplot(1, 2, 1);% 左图，0到1-tau区域
u_e_left = u_e(1:N/2+1);  % 解析解在该区域的值
F_left = F_values((1:N/2+1));  % 神经网络近似解在该区域的值
hold on;
plot(x(1:N/2+1), u_e_left, 'r-', 'LineWidth', 1.5); % 解析解，红色实线 
plot(x(1:N/2+1), F_left, 'b--', 'LineWidth', 1.5);  % 神经网络近似解，蓝色虚线
legend( '解析解','NN', 'Location', 'northwest');
xlim([0, 1 - tau]);  % 设置x轴范围
ylim([0,1]);  % 设置y轴范围

subplot(1, 2, 2);% 右图，1-tau到1区域
u_e_right = u_e(N/2+1:N+1);  % 解析解在该区域的值
F_right = F_values(N/2+1:N+1);  % 神经网络近似解在该区域的值
hold on;
plot(x(N/2+1:N+1), u_e_right,'r-', 'LineWidth', 1.5);  % 解析解
plot(x(N/2+1:N+1), F_right, 'b--','LineWidth', 1.5);  % 神经网络近似解
xlim([1 - tau,1]);  % 设置x轴范围
ylim([0,1]);  % 设置y轴范围


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%最终结果输出

fprintf('cost:  %.5e\n', cost_value(end));%Cost值
fprintf('L2误差: %.2e\n', L_2(end)); %L2误差

%最大值误差
max_error = max(abs(F_values-u_e));
fprintf('max误差: %.2e\n',max_error);

%H1误差
H1_semi = 0;
for i = 1:N
    h_i = x(i+1) - x(i);
    u_fem_prime = (F_values(i+1) - F_values(i))/h_i; % FEM导数
    
    % 高斯点坐标和权重
    gauss_pts = x(i) + h_i/2*[1 - 1/sqrt(3); 1 + 1/sqrt(3)];
    weights = [h_i/2, h_i/2];
    
    % 计算高斯点处导数误差
    der_error = u_fem_prime - u_prime(gauss_pts);
    H1_semi = H1_semi + sum(weights .* der_error.^2);
end
error_h = sqrt(L_2(end)^2 + epsilon*H1_semi);
fprintf('H1误差: %.2e\n',error_h);


%%
x_sub = x(x >= 0 & x <= 1-tau);
F_sub = F_values(x >= 0 & x <= 1-tau);
u_e_sub = u_e(x >= 0 & x <= 1-tau);
diff = max(F_sub - u_e_sub) - min(F_sub - u_e_sub);

max_1 = max(abs(F_right-u_e_right));

fprintf('max_1: %.2e\n',max_1);
fprintf('diff: %.2e\n',diff);

%%
u_fem_right = u_fem(N/2+1:N+1)';
u_fem_sub = u_fem(x >= 0 & x <= 1-tau)';
diff_fem = max(u_fem_sub - u_e_sub) - min(u_fem_sub - u_e_sub);
max_fem = max(abs(u_fem_right-u_e_right));
diff_fem_rig = max(u_fem_right - u_e_right) - min(u_fem_right - u_e_right);


fprintf('max_fem: %.2e\n',max_fem);
fprintf('diff_fem: %.2e\n',diff_fem);


%%
diff_rig = max(F_right - u_e_right) - min(F_right - u_e_right);

fprintf('diff_rig: %.2e\n',diff_rig);
fprintf('diff_fem_rig: %.2e\n',diff_fem_rig);

%%
figure(6)
left = linspace(0,1-tau,5001);
right = linspace(1-tau,1,5001);

subplot(1, 2, 1);% 左图，0到1-tau区域
u_e_left = u_exact(left);  % 解析解在该区域的值
F_left = F(left,b2,W2,W3);  % 神经网络近似解在该区域的值
hold on;
plot(left, u_e_left, 'r-', 'LineWidth', 1.5); % 解析解，红色实线 
plot(left, F_left, 'b--', 'LineWidth', 1.5);  % 神经网络近似解，蓝色虚线
% for i = 1: N/2
%     xline(x(i), ':k');
% end
legend( '精确解','NN', 'Location', 'northwest');
title('[0, 1-\tau]区间');
xlim([0, 1 - tau]);  % 设置x轴范围
ylim([0,1]);  % 设置y轴范围

subplot(1, 2, 2);% 右图，1-tau到1区域
u_e_right = u_exact(right);  % 解析解在该区域的值
F_right = F(right,b2,W2,W3);  % 神经网络近似解在该区域的值
hold on;
% for i = 1: N/2
%     xline(x(N/2+i), ':k');
% end
plot(right, u_e_right,'r-', 'LineWidth', 1.5);  % 解析解
plot(right, F_right, 'b--','LineWidth', 1.5);  % 神经网络近似解
title('[1-\tau,1]区间');
xlim([1 - tau,1]);  % 设置x轴范围
ylim([0,1]);  % 设置y轴范围

%%
figure(7)
subplot(1, 2, 1);% 左图，0到1-tau区域
error1_left = error1(1:N/2+1);  % 在该区域误差值
hold on;
plot(x(1:N/2+1), error1_left,  'LineWidth', 1.5); % 解析解，红色实线 
legend( '节点绝对值误差', 'Location', 'northwest');
xlim([0, 1 - tau]);  % 设置x轴范围
ylim([0,1.5e-3]);  % 设置y轴范围
subplot(1, 2, 2);% 右图，1-tau到1区域
error1_right = error1(N/2+1:N+1);  % 在该区域误差值
hold on;
plot(x(N/2+1:N+1), error1_right, 'LineWidth', 1.5);  % 解析解
legend( '节点绝对值误差', 'Location', 'northwest');
xlim([1 - tau,1]);  % 设置x轴范围
ylim([0,1.5e-3]);  % 设置y轴范围