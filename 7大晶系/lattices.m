clc;
clear;

print_all = true;

% Triclinic Lattice 

clearvars -except print_all;

lattice = build_lattice('Triclinic', 10, 'a', 4.0, 'b', 3.0, 'c', 5.0, 'alpha', 1.2 * pi, 'beta', 0.7 * pi, 'gamma', 0.9 * pi);

try
    if print_all
        subplot(2, 4, 1);
    end
catch
    disp('Program has been ran in sectional mode.')
end
plot3(lattice(:, 1), lattice(:, 2), lattice(:, 3), 'o', 'MarkerFaceColor','g', 'MarkerSize', 10);
title('Triclinic Lattice');
axis square;

%% Monoclinic Lattice

clearvars -except print_all;

lattice = build_lattice('Monoclinic', 10, 'a', 5.0, 'b', 2.0, 'c', 8.0, 'beta', 0.7 * pi);

try
    if print_all
        subplot(2, 4, 2);
    end
catch
    disp('Program has been ran in sectional mode.')
end
plot3(lattice(:, 1), lattice(:, 2), lattice(:, 3), 'o', 'MarkerFaceColor','g', 'MarkerSize', 10);
title('Monoclinic Lattice');
axis square;

%% Trigonal Lattice UNDER CONSTRUCTION

clearvars -except print_all;

lattice = build_lattice('Trigonal', 10, 'a', 5.0, 'alpha', 0.4 * pi);

try
    if print_all
        subplot(2, 4, 3);
    end
catch
    disp('Program has been ran in sectional mode.')
end
plot3(lattice(:, 1), lattice(:, 2), lattice(:, 3), 'o', 'MarkerFaceColor','g', 'MarkerSize', 10);
title('Trigonal Lattice');
axis square;

%% Orthorhombic Lattice

clearvars -except print_all;

lattice = build_lattice('Orthorhombic', 10, 'a', 7.0, 'b', 5.0, 'c', 10.0);

try
    if print_all
        subplot(2, 4, 4);
    end
catch
    disp('Program has been ran in sectional mode.')
end
plot3(lattice(:, 1), lattice(:, 2), lattice(:, 3), 'o', 'MarkerFaceColor','g', 'MarkerSize', 10);
title('Orthorhombic Lattice');
axis square;

%% Tetragonal Lattice

clearvars -except print_all;

lattice = build_lattice('Tetragonal', 10, 'a', 5.0, 'c', 8.0);

try
    if print_all
        subplot(2, 4, 5);
    end
catch
    disp('Program has been ran in sectional mode.')
end
plot3(lattice(:, 1), lattice(:, 2), lattice(:, 3), 'o', 'MarkerFaceColor','g', 'MarkerSize', 10);
title('Tetragonal Lattice');
axis square;

%% Cubic Lattice

clearvars -except print_all;

lattice = build_lattice('Cubic', 10, 'a', 5.0);

try
    if print_all
        subplot(2, 4, 6);
    end
catch
    disp('Program has been ran in sectional mode.')
end
plot3(lattice(:, 1), lattice(:, 2), lattice(:, 3), 'o', 'MarkerFaceColor','g', 'MarkerSize', 10);
title('Cubic Lattice');
axis square;

%% Hexagonal Lattice UNDER CONSTRUCTION

clearvars -except print_all;

lattice = build_lattice('Hexagonal', 10, 'a', 5.0);

try
    if print_all
        subplot(2, 4, 7);
    end
catch
    disp('Program has been ran in sectional mode.')
end
plot3(lattice(:, 1), lattice(:, 2), lattice(:, 3), 'o', 'MarkerFaceColor','g', 'MarkerSize', 10);
title('Hexagonal Lattice');
axis square;

%% Function Defination

function [lattice] = build_lattice(lattice_name, enlarge_limit, varargin)
    if enlarge_limit < 0
        error('Invalid enlarge limitation');
    end
    switch lattice_name
        case 'Triclinic'
            try
                a = get_parameter('a', varargin);
                b = get_parameter('b', varargin);
                c = get_parameter('c', varargin);
                alpha = get_parameter('alpha', varargin);
                beta = get_parameter('beta', varargin);
                gamma = get_parameter('gamma', varargin);
            catch
                error('Parameters invalid or not enough');
            end
            atom1 = [0.0, 0.0, 0.0];
            ux = [a, 0.0, 0.0];
            uy = [b * cos(gamma), b * sin(gamma), 0.0];
            uz = [c * sin(beta), c * sin(alpha), c * cos(beta) * cos(alpha)];
            id = 1;
            for i = -enlarge_limit:enlarge_limit
                for j = -enlarge_limit:enlarge_limit
                    for k = -enlarge_limit:enlarge_limit
                        vector = ux * i + uy * j + uz * k;
                        lattice(id, 1:3) = atom1 + vector;
                        id = id + 1;
                    end
                end
            end

        case 'Monoclinic'
            try
                a = get_parameter('a', varargin);
                b = get_parameter('b', varargin);
                c = get_parameter('c', varargin);
                beta = 0.7*pi;
            catch
                error('Parameters invalid or not enough');
            end
            atom1 = [0.0, 0.0, 0.0];
            ux = [a, 0.0, 0.0];
            uy = [0.0, b, 0.0];
            uz = [c * cos(beta), 0.0, c * sin(beta)];
            id = 1;
            for i = -enlarge_limit:enlarge_limit
                for j = -enlarge_limit:enlarge_limit
                    for k = -enlarge_limit:enlarge_limit
                        vector = ux * i + uy * j + uz * k;
                        lattice(id, 1:3) = atom1 + vector;
                        id = id + 1;
                    end
                end
            end

        case 'Trigonal'
            try
                a = get_parameter('a', varargin);
                alpha = get_parameter('alpha', varargin);
            catch
                error('Parameters invalid or not enough');
            end
            atom1 = [0.0, 0.0, 0.0];
            ux = [a, 0.0, 0.0];
            uy = [a * cos(alpha), a * sin(alpha), 0.0];
            uz = [a * sin(alpha), a * sin(alpha), a * cos(alpha) * cos(alpha)];
            id = 1;
            for i = 0:enlarge_limit
                for j = 0:enlarge_limit
                    for k = 0:enlarge_limit
                        vector = ux * i + uy * j + uz * k;
                        lattice(id, 1:3) = atom1 + vector;
                        id = id + 1;
                    end
                end
            end

        case 'Orthorhombic'
            try
                a = get_parameter('a', varargin);
                b = get_parameter('b', varargin);
                c = get_parameter('c', varargin);
            catch
                error('Parameters invalid or not enough')
            end
            atom1 = [0.0, 0.0, 0.0];
            ux = [a, 0.0, 0.0];
            uy = [0.0, b, 0.0];
            uz = [0.0, 0.0, c];
            id = 1;
            for i = -enlarge_limit:enlarge_limit
                for j = -enlarge_limit:enlarge_limit
                    for k = -enlarge_limit:enlarge_limit
                        vector = ux * i + uy * j + uz * k;
                        lattice(id, 1:3) = atom1 + vector;
                        id = id + 1;
                    end
                end
            end

        case 'Tetragonal'
            try
                a = get_parameter('a', varargin);
                c = get_parameter('c', varargin);
            catch
                error('Parameters invalid or not enough');
            end
            atom1 = [0.0, 0.0, 0.0];
            ux = [a, 0.0, 0.0];
            uy = [0.0, a, 0.0];
            uz = [0.0, 0.0, c];
            id = 1;
            for i = -enlarge_limit:enlarge_limit
                for j = -enlarge_limit:enlarge_limit
                    for k = -enlarge_limit:enlarge_limit
                        vector = ux * i + uy * j + uz * k;
                        lattice(id, 1:3) = atom1 + vector;
                        id = id + 1;
                    end
                end
            end

        case 'Cubic'
            try
                a = get_parameter('a', varargin);
            catch
                error('Parameters invalid or not enough');
            end
            atom1 = [0.0, 0.0, 0.0];
            ux = [a, 0.0, 0.0];
            uy = [0.0, a, 0.0];
            uz = [0.0, 0.0, a];
            id = 1;
            for i = -enlarge_limit:enlarge_limit
                for j = -enlarge_limit:enlarge_limit
                    for k = -enlarge_limit:enlarge_limit
                        vector = ux * i + uy * j + uz * k;
                        lattice(id, 1:3) = atom1 + vector;
                        id = id + 1;
                    end
                end
            end

        case 'Hexagonal'
            try
                a = get_parameter('a', varargin);
                try
                    c = get_parameter('a', varargin);
                catch
                    c = a * 1.633;
                end
            catch
                error('Parameters invalid or not enough');
            end
            atom1 = [0.0, 0.0, 0.0];
            atom2 = [0.5 * a, sqrt(3) / 2 * a, 0.0];
            atom3 = [0.0, 2. / sqrt(3) * a, 0.5 * c];
            atom4 = [-0.5 * a, -a / sqrt(3) / 2, 0.5 * c];
            ux = [a, 0.0, 0.0];
            uy = [0.0, sqrt(3) * a, 0.0];
            uz = [0.0, 0.0, c];
            id = 1;
            for i = -enlarge_limit:enlarge_limit
                for j = -enlarge_limit:enlarge_limit
                    for k = -enlarge_limit:enlarge_limit
                        vector = ux * i + uy * j + uz * k;
                        lattice(id, 1:3) = atom1 + vector;
                        id = id + 1;
                        lattice(id, 1:3) = atom2 + vector;
                        id = id + 1;
                        lattice(id, 1:3) = atom3 + vector;
                        id = id + 1;
                        lattice(id, 1:3) = atom4 + vector;
                        id = id + 1;
                    end
                end
            end

        otherwise
            error('Invalid lattice name');
    end

    function result = get_parameter(parameter_name, parameter_array)
        result = parameter_array{find(strcmp(parameter_name, parameter_array)) + 1};
    end
end