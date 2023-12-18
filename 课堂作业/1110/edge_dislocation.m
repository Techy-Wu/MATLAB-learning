clc
clear

Type = 'SC';
elemment = 'Po';
a0 = 1.0;

atom = [0.0, 0.0, 0.0];

ux = [a0, 0.0, 0.0];
uy = [0.0, a0, 0.0];
uz = [0.0, 0.0, a0];

id = 1;

for i = -10:10
    for j = -10:10
        for k = -10:10
            vector = uz * i + uy * j + uz * k;
            crystal(id, 1:3) = atom + vector;
            id = id + 1;
        end
    end
end

nt = id;

for id = 1:nt
    if abs(crystal(id, 2)) <= 0 & crystal(id, 3) < 0.0
        crystal(id, 1:3) = crystal(nt, 1:3);
        id = id - 1;
        nu = nt - 1;
    end
end

for id = 1:nt
    if crystal(id, 3) < 0.0
        vector = a0 / 20 * ((10 * a0) / crystal(id, 2));
        crystal(id, 2) = crystal(id, 2) - vector;
    end
end

fileID = fopen('SC_lattice_edge_dislocation.pdb', 'w');
formatSpec = '%-6s%-5d %-4s%-1s%-3s  %-4d%-1s   %-8.3f%-8.3f%-8.3f\n';
for i = 1:nt
    fprintf(fileID, formatSpec, 'ATOM', num, element, ' ', ' ', 1, ' ', crystal(id, 1), crystal(id, 2), crystal(id, 3));
end
fclose(fileID);