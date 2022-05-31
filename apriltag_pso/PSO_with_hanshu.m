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

% ��ʼ������Ⱥ�ĸ����ȫ����Сֵ
for i = 1:particlesize  % ����˶�ÿһ�����ӣ��ڸ���λ���ϵ���Ӧֵ
    f(i) = hanshu(data1, (x(i, :))');
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
        f(i) = hanshu(data1, (x(i, :))');
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