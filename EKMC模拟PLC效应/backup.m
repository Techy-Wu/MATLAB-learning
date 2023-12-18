clc %清除命令窗口的内容
clear %清除工作空间的所有变量
clear all %清除工作空间的所有变量，函数，和MEX文件
clf %清除当前的Figure
close %关闭当前的Figure窗口
close all %关闭所有的Figure窗口

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 全局条件 %

% 主变量（ε） % (links 2 P48)

global global_applied_strain;
global_applied_strain = 2.0 * 10 ^ (-2);

% 基础常量 % (links 2 P38)

global global_K global_p0 global_w global_P2 global_Cm global_sigma_p global_delta_t;
global_K = 5.2;
global_p0 = 0.94 * 10 ^ (-4);
global_w = 5 * 10 ^ (-4);
global_P2 = 0.4;
global_Cm = 1.38;
global_sigma_p = 142.5; %44.8
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
global_Q = 180;
global_P1 = 6.44;

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
global_T = 20 + 273;
global_Va = global_kb * global_T / global_K;
global_Ea = 134;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 初始条件 %

time_array = [10^(-5), 5*10^(-5), 10^(-4), 5*10^(-4), 10^(-3), 5*10^(-3), 10^(-2), 5*10^(-2), 10^(-1)];

for try_time_step = 1:9

    global_delta_t = time_array(try_time_step)

    loopN = 6000;

    strain_p = zeros(loopN, 32); %可能有误
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

    strain_p_ratio(1, :) = 2 * 10 ^ (-2);
    strain_p(1, :) = 0.002;
    ta(1, :) = 0;
    area_total(1) = 12;
    stress_total(1) = 170;
    strain_total(1) = 0.064;

    hWaitbar = waitbar(0, 'Computing，wait...', 'CreateCancelBtn', 'delete(gcbf)');
    set(hWaitbar, 'Color', [0.9, 0.9, 0.9]);
    compT = 0;
    waitbar(compT, hWaitbar, ['Computing...', num2str(round(compT, 2) * 100), '%']);

    timestr = ['dT', num2str(global_delta_t*10^5)];

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

        for index = 1:32
            if step == 1
                stress(step, index) = get_stress(strain_p(step, index), strain_p_ratio(step, index), ta(step, index));
            end
            ta_ratio(step, index) = get_ta_ratio(ta(step, index), strain_p_ratio(step, index));
        end

        % for index = 1:32
        %     disp("stress(" + index + ") = " + stress(step, index));
        % end

        fprintf('完成 (循环%d)\n', i);

        % 叠加应力 %

        fprintf('叠加应力 (循环%d)\n', i);

        for index = 1:32
            stress(step, index) = stress(step, index)*(index - 17) * 10 ^ (-4) + stress(step, index);
        end

        fprintf('完成 (循环%d)\n', i);

        % 反代解出各单元应变 %

        fprintf('反代解出各单元应变率 (循环%d)\n', i);

        for index = 1:32
            strain_p_ratio(step, index) = get_strain_p_ratio(stress(step, index), strain_p(step, index), ta(step, index));
        end

        fprintf('完成 (循环%d)\n', i);

        % 总结计算，并赋值到下一时刻为初始条件 %

        fprintf('总结计算 (循环%d)\n', i);

        step = step + 1;

        fprintf('求时刻、应变、有效时间 (循环%d)\n', i);
        for index = 1:32
            % 求时刻 %;
            t(step) = t(step - 1) + global_delta_t;
            % 求应变 %
            strain_p(step, index) = strain_p(step - 1, index) + strain_p_ratio(step - 1, index) * global_delta_t;
            % 求有效时间 %
            ta(step, index) = ta(step - 1, index) + ta_ratio(step - 1, index) * global_delta_t;
        end

        % 求总应力变化率 %
        fprintf('求总应力变化率 (循环%d)\n', i);
        stress_total_ratio(step - 1) = get_stress_total_ratio(strain_p_ratio(step - 1, 1:32)); %此处get_stress_total_ratio的step是否需要进1
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
        for index = 1:32
            stress(step, index) = get_next_stress(force(step), strain_p(step, index), get_diff_area(area_total(step), index));
        end

        fprintf('完成 (循环%d)\n', i);

        compT = strain_total(step) / 0.2;
        waitbar(compT, hWaitbar, ['Computing...', num2str(round(compT, 4) * 100), '% Strain = ', num2str(strain_total(step))]);

        save_data(stress_total(step), strain_total(step), timestr, 2);

        if strain_total(step) >= 0.2
            fprintf('到达指定应变值\n');
            close all
            break;
        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 函数定义声明 %

function stress = get_stress(strain_p, strain_p_ratio, ta0)
    global global_K global_p0 global_lambda global_beta global_Q global_b global_P1 global_Cm global_P2 global_alpha global_n_alpha global_sigma_p;
    a = global_K * asinh(strain_p_ratio / global_p0);
    b = global_sigma_p;
    c = (1 + global_lambda * log(strain_p ^ global_beta)) * global_Q * (1 - exp(-global_b * strain_p));  %matlab中，log(x) = ln(x)
    d = global_P1 * global_Cm * (1 - exp(-global_P2 * strain_p ^ global_alpha * ta0 ^ global_n_alpha)); %第一次计算ta=0,d=0
    stress = a + b + c + d;
end

function ta = get_ta_ratio(ta, strain_p_ratio)
    global global_w;
    ta = 1 - ta / global_w * strain_p_ratio;
end

function answer = get_strain_p_ratio(stress,strain_p, ta0)
    global global_K global_p0 global_lambda global_beta global_Q global_b global_P1 global_Cm global_P2 global_alpha global_n_alpha global_sigma_p;
    syms strain_p_ratio;
    a = global_K * asinh(strain_p_ratio / global_p0);
    b = global_sigma_p;
    c = (1 + global_lambda * log(strain_p ^ global_beta)) * global_Q * (1 - exp(-global_b * strain_p));  %matlab中，log(x) = ln(x)
    d = global_P1 * global_Cm * (1 - exp(-global_P2 * strain_p ^ global_alpha * ta0 ^ global_n_alpha)); %第一次计算ta=0,d=0
    result = solve(stress == a + b + c + d, strain_p_ratio);
    answer = double(result);
end

function stress_total_ratio = get_stress_total_ratio(strain_p_ratio)
    global global_V global_L global_E;
    %a = global_V / global_L;
    a = 0.02;
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

function answer = get_next_stress(force, next_strain_p, diff_area)
    global global_E;
    syms stress;
    equation = (stress == force * (1 + next_strain_p + stress / global_E) / diff_area);
    result = solve(equation, stress);
    answer = double(result);
end

function ans = arcsinh(x)
    ans = log(x + sqrt(x ^ 2 + 1));
end

function save_data(stress, strain, suffix, title)
    fileName = ['SS_Curve_', suffix, '.txt'];
    fileID = fopen(fileName, 'a');
    if title == 1
        fprintf(fileID, '%-10s%-10s\n', 'Stress', 'Strain');
    else
        fprintf(fileID, '%-10.5f%-10.5f\n', stress, strain);
    end
    fclose(fileID);
end