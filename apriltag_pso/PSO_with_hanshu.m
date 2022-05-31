clc
clear
close all

% 导入数据
data1 = xlsread('apriltag.xlsx', 'sheet1', 'A1:D12');

%%
tic;  % 程序运行计时

% 设置相关参数值
E = 0.00000001;  % 允许误差
maxnum = 3000;  % 粒子最大迭代次数
narvs = 12;  % 目标函数的自变量个数，空间维度N，这里是12，代表十二维变量
particlesize = 200;  % 粒子群规模
c1 = 2;  % 每个粒子的个体学习因子，加速度常数1
c2 = 2;  % 每个粒子的社会学习因子，加速度常数2
w = 0.6;  % 惯性因子，一般0.6到0.9
vmax = 0.5;  % 粒子的最大飞行速度，一般取每维变化范围的10%到20%

% 初始化粒子的速度和位置以及目标函数
limit = [-1, 1];
v = vmax * rand(particlesize, narvs);  % 初始化粒子的飞行速度，50行2列的速度矩阵
for i = 1:narvs
    x(:, i) = limit(1) + (limit(2) - limit(1)) .* rand(particlesize, 1);
end

% 初始化粒子群的个体和全局最小值
for i = 1:particlesize  % 完成了对每一个粒子，在各自位置上的适应值
    f(i) = hanshu(data1, (x(i, :))');
end
personalbest_x = x;  % 用于存储对于每一个粒子最佳经历点的自变量x的位置
personalbest_fval = f;  % 同时存储对于每一个粒子的最佳经历点的适应值
[globalbest_fval, i] = min(personalbest_fval);  % min函数返回的第一个是最小值，还有一个就是最小值的下标，这里就是告诉了是在哪个粒子上
globalbest_x = personalbest_x(i, :);  % 这个必定是全局最优点的位置

%%
% 开始迭代运算
k = 1;
while k <= maxnum
    % 更新粒子群的个体和全局最小值
    for i = 1:particlesize
        f(i) = hanshu(data1, (x(i, :))');
        if f(i) < personalbest_fval(i)
            personalbest_fval(i) = f(i);
            personalbest_x(i, :) = x(i, :);
        end
    end
    [globalbest_fval, i] = min(personalbest_fval);
    globalbest_x = personalbest_x(i, :);
    % 更新粒子的速度和位置信息
    for i = 1:particlesize
        v(i, :) = w*v(i, :) + c1*rand*(personalbest_x(i, :) - x(i, :)) + c2*rand*(globalbest_x - x(i, :));
        for j = 1:narvs
            if v(i, j) > vmax
                v(i, j) = vmax;
            elseif v(i, j) < -vmax
                v(i, j) = -vmax;
            end
        end
        x(i, :) = x(i, :) + v(i, :);
    end
    ff(k) = globalbest_fval;  % 为了画适应值曲线图
    % 检查是否符合停止迭代的条件
    if abs(globalbest_fval) < E
        break  % 如果适应值即求取的函数值小于自己预定的一个初值，即跳出循环
    end
    k = k + 1;
end

%%
% 输出结果
figure(1)
plot(1:length(ff), ff)
title("适应值")
disp(['最小值： ', num2str(globalbest_fval)]);
disp(['自变量取值： ', num2str(globalbest_x)]);
% 程序运行结束
toc;