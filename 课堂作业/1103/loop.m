clc
clear

for i = 1:1:10
    fprintf("i = %d\n", i);
end

for i = 1:2:10
    fprintf("i = %d\n", i);
end

for i = 1:1:10
    for j = 1:1:10
        fprintf("i + j = %d\n", i + j);
    end
end

for i = -5:5
    if i >= 0
        fprintf("|i| = %d\n", i);
    else
        fprintf("|i| = %d\n", -i)
    end
end

for i = -5:5
    if i < -3 || i > 3
        fprintf("OutOfRange\n");
    elseif i == 0
        fprintf("i = 0\n");
    elseif i >= 0
        fprintf("|i| = %d\n", i);
    else
        fprintf("|i| = %d\n", -i)
    end
end

for i = 1:20
    if mod(i, 2) == 0
        fprintf("%d is an even.\n", i);
    else
        fprintf("%d is an odd.\n", i);
    end
end