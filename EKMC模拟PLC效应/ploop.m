clc
clear

global c d e;
c = 2;
d = 3;
e = 4;

a = [2,3,4,5];
b = [1,1,1,1];

parameters = {'cd', c;
    'd', d;
    'edf', e};

[x,y] = find(strcmp('d',parameters));
[x,y]

disp([x,y]);

parameters{x,y + 1}
disp(parameters{x,y + 1});

ans = parameters{x,y + 1}
disp(ans);

parfor i = 1:4
    a(i) = add(a(i), b(i), parameters);
end

a

function r = add(a, b, parameters)
    %[row, col] = find(parameters == 'd');
    [row, col] = find(strcmp('edf',parameters));
    c = parameters{row, col + 1};
    r = a * b + c;
    %r = a * b;
end