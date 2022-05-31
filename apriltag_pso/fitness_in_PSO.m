clc
clear
close all

% ��������
data1 = xlsread('apriltag.xlsx', 'sheet1', 'A1:D12');

%%
tic;  % �������м�ʱ

% ������ز���ֵ
E = 0.00000001;  % �������
maxnum = 3000;  % ��������������
narvs = 12;  % Ŀ�꺯�����Ա����������ռ�ά��N��������12������ʮ��ά����
particlesize = 200;  % ����Ⱥ��ģ
c1 = 2;  % ÿ�����ӵĸ���ѧϰ���ӣ����ٶȳ���1
c2 = 2;  % ÿ�����ӵ����ѧϰ���ӣ����ٶȳ���2
w = 0.6;  % �������ӣ�һ��0.6��0.9
vmax = 0.5;  % ���ӵ��������ٶȣ�һ��ȡÿά�仯��Χ��10%��20%

% ��ʼ�����ӵ��ٶȺ�λ���Լ�Ŀ�꺯��
limit = [-1, 1];
v = vmax * rand(particlesize, narvs);  % ��ʼ�����ӵķ����ٶȣ�50��2�е��ٶȾ���
for i = 1:narvs
    x(:, i) = limit(1) + (limit(2) - limit(1)) .* rand(particlesize, 1);
end

fitness = @(nx, ny, nz, ox, oy, oz, ax, ay, az, px, py, pz) (data1(1, 1) - nx).^2 + (data1(2, 1) - ny).^2 + (data1(3, 1) - nz).^2 ...
    + (data1(4, 1) - ox).^2 + (data1(5, 1) - oy).^2 + (data1(6, 1) - oz).^2 ...
    + (data1(7, 1) - ax).^2 + (data1(8, 1) - ay).^2 + (data1(9, 1) - az).^2 ...
    + (data1(10, 1) - px).^2 + (data1(11, 1) - py).^2 + (data1(12, 1) - pz).^2 ...
    + (data1(1, 2) - nx).^2 + (data1(2, 2) - ny).^2 + (data1(3, 2) - nz).^2 ...
    + (data1(4, 2) - ox).^2 + (data1(5, 2) - oy).^2 + (data1(6, 2) - oz).^2 ...
    + (data1(7, 2) - ax).^2 + (data1(8, 2) - ay).^2 + (data1(9, 2) - az).^2 ...
    + (data1(10, 2) - px).^2 + (data1(11, 2) - py).^2 + (data1(12, 2) - pz).^2 ...
    + (data1(1, 3) - nx).^2 + (data1(2, 3) - ny).^2 + (data1(3, 3) - nz).^2 ...
    + (data1(4, 3) - ox).^2 + (data1(5, 3) - oy).^2 + (data1(6, 3) - oz).^2 ...
    + (data1(7, 3) - ax).^2 + (data1(8, 3) - ay).^2 + (data1(9, 3) - az).^2 ...
    + (data1(10, 3) - px).^2 + (data1(11, 3) - py).^2 + (data1(12, 3) - pz).^2 ...
    + (data1(1, 4) - nx).^2 + (data1(2, 4) - ny).^2 + (data1(3, 4) - nz).^2 ...
    + (data1(4, 4) - ox).^2 + (data1(5, 4) - oy).^2 + (data1(6, 4) - oz).^2 ...
    + (data1(7, 4) - ax).^2 + (data1(8, 4) - ay).^2 + (data1(9, 4) - az).^2 ...
    + (data1(10, 4) - px).^2 + (data1(11, 4) - py).^2 + (data1(12, 4) - pz).^2;

% ��ʼ������Ⱥ�ĸ����ȫ����Сֵ
for i = 1:particlesize  % ����˶�ÿһ�����ӣ��ڸ���λ���ϵ���Ӧֵ
    f(i) = fitness(x(i, 1), x(i, 2), x(i, 3), x(i, 4), x(i, 5), x(i, 6), x(i, 7), x(i, 8), x(i, 9), x(i, 10), x(i, 11), x(i, 12));
end
personalbest_x = x;  % ���ڴ洢����ÿһ��������Ѿ�������Ա���x��λ��
personalbest_fval = f;  % ͬʱ�洢����ÿһ�����ӵ���Ѿ��������Ӧֵ
[globalbest_fval, i] = min(personalbest_fval);  % min�������صĵ�һ������Сֵ������һ��������Сֵ���±꣬������Ǹ����������ĸ�������
globalbest_x = personalbest_x(i, :);  % ����ض���ȫ�����ŵ��λ��

%%
% ��ʼ��������
k = 1;
while k <= maxnum
    % ��������Ⱥ�ĸ����ȫ����Сֵ
    for i = 1:particlesize
        f(i) = fitness(x(i, 1), x(i, 2), x(i, 3), x(i, 4), x(i, 5), x(i, 6), x(i, 7), x(i, 8), x(i, 9), x(i, 10), x(i, 11), x(i, 12));
        if f(i) < personalbest_fval(i)
            personalbest_fval(i) = f(i);
            personalbest_x(i, :) = x(i, :);
        end
    end
    [globalbest_fval, i] = min(personalbest_fval);
    globalbest_x = personalbest_x(i, :);
    % �������ӵ��ٶȺ�λ����Ϣ
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
    ff(k) = globalbest_fval;  % Ϊ�˻���Ӧֵ����ͼ
    % ����Ƿ����ֹͣ����������
    if abs(globalbest_fval) < E
        break  % �����Ӧֵ����ȡ�ĺ���ֵС���Լ�Ԥ����һ����ֵ��������ѭ��
    end
    k = k + 1;
end

%%
% ������
figure(1)
plot(1:length(ff), ff)
title("��Ӧֵ")
disp(['��Сֵ�� ', num2str(globalbest_fval)]);
disp(['�Ա���ȡֵ�� ', num2str(globalbest_x)]);
% �������н���
toc;