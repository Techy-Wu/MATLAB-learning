clc
clear

x = 0:0.1:20;
y = sin(x);
z = cos(x);

plot3(x, y, z, 'r');
grid;
box;