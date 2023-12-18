clc %清除命令窗口的内容
clear %清除工作空间的所有变量
clear all %清除工作空间的所有变量，函数，和MEX文件
clf %清除当前的Figure
close %关闭当前的Figure窗口
close all %关闭所有的Figure窗口

x = -10:0.001:10;
y = power(x,2);

plot(x,y);
xlabel("x");
ylabel("y");
legend("y=x^2");