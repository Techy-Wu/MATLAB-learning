clc %清除命令窗口的内容
clear %清除工作空间的所有变量
clear all %清除工作空间的所有变量，函数，和MEX文件
clf %清除当前的Figure
close %关闭当前的Figure窗口
close all %关闭所有的Figure窗口

% Equations %
%E = Ee + Ep; %应变张量分为弹性和塑性两方面

% 参数表 % (links 2 P48)
% Reference %
applied_strain = [2.0*10^(-2), 1.5*10^(-3), 1.2*10^(-4)];
% Constants %
K = 5.2;
p0 = 0.94*10^(-4);
w = 5*10^(-4);
P2 = 0.4;
Cm = 1.38;
% Optimazed Parameters %
lambda = 0.42;
b = 15.8;
alpha = 0.67;
n_alpha = 0.33;
beta = [0.57, 0.64, 0.67];
% From Experiment%
Q = [180,200,240];
P1 = [6.44, 17.33, -19.80];
% Calculate Conditions %
E = 7*10^4;
V = [38, 2.88, 0.23];
L = 32;
delta_t = 1*10^(-4);

% Substuting %
strain = 