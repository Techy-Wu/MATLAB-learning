clc
clear

a = [2,3,4;8,11,9;5,8,13];
%a = [2 2 3; 4 5 6; 7 8 9];

b = det(a);%行列式值
c = inv(a);%矩阵的逆
d = rank(a);%矩阵的秩
e = a';%矩阵的转置

b
c
d
e

a(2,3)
a(1,:)
a(1,1:2)
a(2:3,1:2)
a(:,1)

f = cross(a,c);%矩阵的叉乘
f