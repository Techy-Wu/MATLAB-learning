clc %清除命令窗口的内容
clear %清除工作空间的所有变量
clear all %清除工作空间的所有变量，函数，和MEX文件
clf %清除当前的Figure
close %关闭当前的Figure窗口
close all %关闭所有的Figure窗口

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 全局条件 %

% 主变量（ε） % (links 2 P48)

global global_applied_strain global_applied_strain_ratio;
global_applied_strain = 0.002;
global_applied_strain_ratio = 0.02;

% 基础常量 % (links 2 P38)

global global_K global_p0 global_w global_P2 global_Cm global_sigma_p global_delta_t;
global_K = 5.2;
global_p0 = 0.94 * 10 ^ (-4);
global_w = 5 * 10 ^ (-4);
global_P2 = 0.4;
global_Cm = 1.38;
global_sigma_p = 142.5; %44.8;
global_delta_t = 10^(-2);

% 优化后的参数 % (links 2 P38)

global global_lambda global_b global_alpha global_n_alpha global_beta;
global_lambda = 0.42;
global_b = 15.8;
global_alpha = 0.67;
global_n_alpha = 0.33;
global_beta = 0.57;

% 实验测得 % (links 2 P38)

global global_Q global_P1;
global_Q = 274; %180
global_P1 = -19.8; %6.44

% 计算得到 % (links 2 P38)

global global_E global_V global_L global_delta_t;
global_E = 7 * 10 ^ 4;
global_V = 38;
global_L = 32;
global_delta_t = 1 * 10 ^ (-4);

% 计算得到 % (Links 2 P18)

global global_b0 global_kb global_T global_Va global_Ea;
global_b0 = 0.286*10^(-9);
global_kb = 1.38 * 10 ^ (-23);
global_T = 400 + 273; %20 + 273;
global_Va = global_kb * global_T / global_K;
global_Ea = 134;

% 参数表 %
parameter_list = {'global_applied_strain', global_applied_strain;
    'global_applied_strain_ratio', global_applied_strain_ratio;
    'global_K', global_K;
    'global_p0', global_p0;
    'global_w', global_w;
    'global_P2', global_P2;
    'global_Cm', global_Cm;
    'global_sigma_p', global_sigma_p;
    'global_delta_t', global_delta_t;
    'global_lambda', global_lambda;
    'global_b', global_b;
    'global_alpha', global_alpha;
    'global_n_alpha', global_n_alpha;
    'global_beta', global_beta;
    'global_Q', global_Q;
    'global_P1', global_P1;
    'global_E', global_E;
    'global_V', global_V;
    'global_L', global_L;
    'global_delta_t', global_delta_t;
    'global_b0', global_b0;
    'global_kb', global_kb;
    'global_T', global_T;
    'global_Va', global_Va;
    'global_Ea', global_Ea};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 初始条件 %

loopN = 6000;

strain_p = zeros(loopN, 32);
strain_p_ratio = zeros(loopN, 32);
strain_total = zeros(loopN, 1);
ta = zeros(loopN, 32);
ta_ratio = zeros(loopN, 32);
stress = zeros(loopN, 32);
t = zeros(loopN, 1);
stress_total = zeros(loopN, 1);
stress_total_ratio = zeros(loopN, 1);
force = zeros(loopN, 1);
area_total = zeros(loopN, 1);

strain_p_ratio(1, :) = global_applied_strain_ratio;
strain_p(1, :) = global_applied_strain;
ta(1, :) = 0;
area_total(1) = 12;
stress_total(1) = 170;
strain_total(1) = 0.064;

hWaitbar = waitbar(0, '计算中', 'CreateCancelBtn', 'delete(gcbf)');
set(hWaitbar, 'Color', [0.9, 0.9, 0.9]);
compT = 0;
waitbar(compT, hWaitbar, ['完成度', num2str(round(compT, 2) * 100), '%']);

timestr = char(datetime("now", "Format","MM-dd-HHmm"));

ref_para_2 = zeros(loopN, 32);
ref_para_3 = zeros(loopN, 32);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 测试输出功能 %

disp("a = " + 1);
disp("reverse_sinh(a) = " + asinh(1));

%disp(e(1:2,1:16));

save_data(0, 0, timestr, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 计算 %

step = 1;

for i = 1:loopN

    % 求各单元初始应力 %
    fprintf('求各单元初始应力 (循环%d)\n', i);

    parfor index = 1:32
        if step == 1
            stress(step, index) = get_stress(strain_p(step, index), strain_p_ratio(step, index), ta(step, index), parameter_list);
        end
        ta_ratio(step, index) = get_ta_ratio(ta(step, index), strain_p_ratio(step, index), parameter_list);
    end

    % for index = 1:32
    %     disp("stress(" + index + ") = " + stress(step, index));
    % end

    fprintf('完成 (循环%d)\n', i);

    % 叠加应力 %

    fprintf('叠加应力 (循环%d)\n', i);

    parfor index = 1:32
        stress(step, index) = stress(step, index)*(index - 17) * 10 ^ (-4) + stress(step, index);
        % stress(step, index) = stress(step, index) + 10 * sin(10000 * rand * 2 * pi);
        % stress(step, index) = stress(step, index) + sinh(stress(step, index)) / 10;
    end

    fprintf('完成 (循环%d)\n', i);

    % 反代解出各单元应变 %

    fprintf('反代解出各单元应变率 (循环%d)\n', i);

    parfor index = 1:32
        [x, y, z] = get_strain_p_ratio(stress(step, index), strain_p(step, index), ta(step, index), parameter_list);
        strain_p_ratio(step, index) = x;
        ref_para_2(step, index) = y;
        ref_para_3(step, index) = z;
        %strain_p_ratio(step, index) = get_strain_p_ratio(stress(step, index), strain_p(step, index), ta(step, index), parameter_list);
    end

    fprintf('完成 (循环%d)\n', i);

    % 总结计算，并赋值到下一时刻为初始条件 %

    fprintf('总结计算 (循环%d)\n', i);

    step = step + 1;

    % 求时刻 %;
    fprintf('求时刻、应变、有效时间 (循环%d)\n', i);
    t(step) = t(step - 1) + global_delta_t;

    fprintf('求应变、有效时间 (循环%d)\n', i);
    former_strain_p = strain_p(step - 1, :);
    former_ta = ta(step - 1, :);
    parfor index = 1:32
        % 求应变 %
        strain_p(step, index) = former_strain_p(index) + strain_p_ratio(step - 1, index) * get_parameter(parameter_list, 'global_delta_t');
        % 求有效时间 %
        ta(step, index) = former_ta(index) + ta_ratio(step - 1, index) * get_parameter(parameter_list, 'global_delta_t');
    end

    % 求总应力变化率 %
    fprintf('求总应力变化率 (循环%d)\n', i);
    stress_total_ratio(step - 1) = get_stress_total_ratio(strain_p_ratio(step - 1, 1:32), parameter_list); %此处get_stress_total_ratio的step是否需要进1
    % 求总应力 %
    fprintf('求总应力 (循环%d)\n', i);
    stress_total(step) = stress_total(step - 1) + stress_total_ratio(step - 1) * global_delta_t;
    % 求总应变 %
    fprintf('求总应变 (循环%d)\n', i);
    strain_total(step) = get_strain_total(strain_p(step, 1:32));
    % 求横截面积 %
    fprintf('求横截面积 (循环%d)\n', i);
    area_total(step) = get_area_total(area_total(step - 1), strain_total(step - 1));
    % 求平均载荷 %
    fprintf('求平均载荷 (循环%d)\n', i);
    force(step) = get_force(stress_total(step), area_total(step), strain_total(step)); %此处area_total的step是否需要进1
    % 求下一时刻的各单元应力 %
    fprintf('求下一时刻的各单元应力 (循环%d)\n', i);
    former_stress = stress(step - 1, :);
    parfor index = 1:32
        stress(step, index) = get_next_stress(force(step), strain_p(step, index), former_stress(index), get_diff_area(area_total(step), index), parameter_list);
    end

    fprintf('完成 (循环%d)\n', i);
    
    % 更新进度条
    compT = strain_total(step) / 0.2;
    waitbar(compT, hWaitbar, ['完成度', num2str(round(compT, 4) * 100), '% 当前总应变 = ', num2str(strain_total(step))]);
    
    % 保存数据
    save_data(strain_total(step), stress_total(step), timestr, 2);

    % 更新曲线
    subplot(1, 2, 1);
    plot(strain_total(1:step), stress_total(1:step));
    xlabel('Strain');
    ylabel('Stress');
    title('Total SS Curve');
    subplot(1, 2, 2);
    plot(strain_p(1:step, 1), stress(1:step, 1));
    xlabel('Strain');
    ylabel('Stress');
    title('Unit1 SS Curve');

    if strain_total(step) >= 0.2
        fprintf('到达指定应变值\n');
        close all
        break;
    end

end

% 保存参数
fileName = ['parameter2_', timestr, '.txt'];
fileID = fopen(fileName, 'a');
for i = 1:step
    for j = 1:32
        fprintf(fileID, '%-15.5f', ref_para_2(i, j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

fileName = ['parameter3_', timestr, '.txt'];
fileID = fopen(fileName, 'a');
for i = 1:step
    for j = 1:32
        fprintf(fileID, '%-15.5f', ref_para_3(i, j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

fileName = ['strain_p_ratio_', timestr, '.txt'];
fileID = fopen(fileName, 'a');
for i = 1:step
    for j = 1:32
        fprintf(fileID, '%-15.5f', strain_p_ratio(i, j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

fileName = ['strain_p_', timestr, '.txt'];
fileID = fopen(fileName, 'a');
for i = 1:step
    for j = 1:32
        fprintf(fileID, '%-15.5f', strain_p(i, j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

fileName = ['stress_', timestr, '.txt'];
fileID = fopen(fileName, 'a');
for i = 1:step
    for j = 1:32
        fprintf(fileID, '%-15.5f', stress(i, j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

% 显示最终曲线
subplot(1, 2, 1);
plot(strain_total(1:step), stress_total(1:step));
xlabel('Strain');
ylabel('Stress');
title('Total SS Curve');
subplot(1, 2, 2);
plot(strain_p(1:step, 1), stress(1:step, 1));
xlabel('Strain');
ylabel('Stress');
title('Unit1 SS Curve');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 函数定义声明 %

function stress = get_stress(strain_p, strain_p_ratio, ta0, parameter_list)
    global_K = get_parameter(parameter_list, 'global_K');
    global_p0 = get_parameter(parameter_list, 'global_p0');
    global_lambda = get_parameter(parameter_list, 'global_lambda');
    global_beta = get_parameter(parameter_list, 'global_beta');
    global_Q = get_parameter(parameter_list, 'global_Q');
    global_b = get_parameter(parameter_list, 'global_b');
    global_P1 = get_parameter(parameter_list, 'global_P1');
    global_Cm = get_parameter(parameter_list, 'global_Cm');
    global_P2 = get_parameter(parameter_list, 'global_P2');
    global_alpha = get_parameter(parameter_list, 'global_alpha');
    global_n_alpha = get_parameter(parameter_list, 'global_n_alpha');
    global_sigma_p = get_parameter(parameter_list, 'global_sigma_p');
    a = global_K * asinh(strain_p_ratio / global_p0);
    b = global_sigma_p;
    c = (1 + global_lambda * log(strain_p ^ global_beta)) * global_Q * (1 - exp(-global_b * strain_p));  %matlab中，log(x) = ln(x)
    % c = global_Q * (1 - exp(-global_b * strain_p));  %matlab中，log(x) = ln(x)
    d = global_P1 * global_Cm * (1 - exp(-global_P2 * strain_p ^ global_alpha * ta0 ^ global_n_alpha)); %第一次计算ta=0,d=0
    stress = a + b + c + d;
end

function ta = get_ta_ratio(ta, strain_p_ratio, parameter_list)
    global_w = get_parameter(parameter_list, 'global_w');
    ta = 1 - ta / global_w * strain_p_ratio;
end

function [answer, y, z] = get_strain_p_ratio(stress, strain_p, ta0, parameter_list)
    global_K = get_parameter(parameter_list, 'global_K');
    global_p0 = get_parameter(parameter_list, 'global_p0');
    global_lambda = get_parameter(parameter_list, 'global_lambda');
    global_beta = get_parameter(parameter_list, 'global_beta');
    global_Q = get_parameter(parameter_list, 'global_Q');
    global_b = get_parameter(parameter_list, 'global_b');
    global_P1 = get_parameter(parameter_list, 'global_P1');
    global_Cm = get_parameter(parameter_list, 'global_Cm');
    global_P2 = get_parameter(parameter_list, 'global_P2');
    global_alpha = get_parameter(parameter_list, 'global_alpha');
    global_n_alpha = get_parameter(parameter_list, 'global_n_alpha');
    global_sigma_p = get_parameter(parameter_list, 'global_sigma_p');
    %syms strain_p_ratio;
    %a = global_K * asinh(strain_p_ratio / global_p0);
    b = global_sigma_p;
    c = (1 + global_lambda * log(strain_p ^ global_beta)) * global_Q * (1 - exp(-global_b * strain_p));  %matlab中，log(x) = ln(x)
    % c = global_Q * (1 - exp(-global_b * strain_p));  %matlab中，log(x) = ln(x)
    d = global_P1 * global_Cm * (1 - exp(-global_P2 * strain_p ^ global_alpha * ta0 ^ global_n_alpha)); %第一次计算ta=0,d=0
    %result = solve(stress == a + b + c + d, strain_p_ratio);
    %answer = double(result);
    y = b + c;
    z = d;
    %if isreal(answer) ~= true
    %    fprintf("E 复数解");
    %end

    % 弃用牛顿法解方程，由 arcsinh(x) = a <-> sin(a) = x 求解
    a = stress - b - c - d;
    strain_p_ratio = sinh(a / global_K) * global_p0;
    answer = strain_p_ratio;
end

function stress_total_ratio = get_stress_total_ratio(strain_p_ratio, parameter_list)
    %global_V = get_parameter(parameter_list, 'global_V');
    %global_L = get_parameter(parameter_list, 'global_L');
    global_E = get_parameter(parameter_list, 'global_E');
    %a = global_V / global_L;
    global_applied_strain_ratio = get_parameter(parameter_list, 'global_applied_strain_ratio');
    a = global_applied_strain_ratio;
    b = mean(strain_p_ratio(:));
    stress_total_ratio = (a - b) * global_E;
end

function area_total = get_area_total(original_area, strain_p)
    area_total = original_area / (1 + strain_p);
end

function force = get_force(stress_total, area_total, strain_total)
    force = stress_total * area_total / (1 + strain_total);
end

function strain_total = get_strain_total(strain)
    strain_total = sum(strain);
end

function diff_area = get_diff_area(area_total, index) %此处以32份计算
    diff_area = 10 ^ (-4) * area_total * (index - 16) + area_total;
end

function answer = get_next_stress(force, next_strain_p, stress, diff_area, parameter_list)
    % global_E = get_parameter(parameter_list, 'global_E');
    % syms stress;
    % equation = (stress == force * (1 + next_strain_p + stress / global_E) / diff_area);
    % result = solve(equation, stress);
    % answer = double(result);
    global_E = get_parameter(parameter_list, 'global_E');
    result = force * (1 + (next_strain_p + stress / global_E) * 5) / diff_area;
    % result = force * 2 * (next_strain_p + stress / global_E) / diff_area;
    answer = result;
end

function save_data(strain, stress, suffix, title)
    fileName = ['SS-Curve_', suffix, '.txt'];
    fileID = fopen(fileName, 'a');
    if title == 1
        fprintf(fileID, '%-15s%-15s\n', 'Strain', 'Stress');
    else
        fprintf(fileID, '%-15.5f%-15.5f\n', strain, stress);
    end
    fclose(fileID);
end

function result = get_parameter(parameter_list, parameter_name)
    [row, col] = find(strcmp(parameter_name, parameter_list));
    result = parameter_list{row, col + 1};
end