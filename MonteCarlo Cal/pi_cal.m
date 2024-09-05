% Define x /in [-1, +1], y /in [-1, +1]

N = 99999999.0;

x_Nist = randperm(N, N);
y_Nist = randperm(N, N);

x_Nist = x_Nist ./ N;
y_Nist = y_Nist ./ N;

count = 0;

parfor i = 1:N
    if distance(x_Nist(i), y_Nist(i)) <= 1
        count = count + 1;
    end
end

s = count / N * 4.0;

fprintf("pi ≈ %f\n", s)

function result = distance(x, y)
    result = sqrt(x * x + y * y);
end