clc
clear

Type = 'BCC';
element = 'Fe';
a0 = 2.86;
mass = 55.8;
Structure = 1;

atom1 = [0.0, 0.0, 0.0];
atom2 = [0.5 * a0, 0.5 * a0, 0.5 * a0];

ux = a0;
uy = a0;
uz = a0;
uxx = 10;
uyy = 10;
uzz = 10;

id = 0;

formatSpec = '%-6s%-5d %-4s%-1s%-3s  %-4d%-1s   %-8.3f%-8.3f%-8.3f\n';

for i = 0:1
    for j = 0:1
        for k = 0:1
            id = id + 1;
            data(id, 1:3) = atom1 + [ux * i, uy * j, uz * k];
        end
    end
end
id = id + 1;
data(id, 1:3) = atom2;

nt = id;

fileID = fopen('Fe_lattice.pdb', 'a');
num = 0;

for id = 1:nt
    if abs(data(id, 1:3)) <= 20.0
        num = num + 1;
        fprintf(fileID, formatSpec, 'ATOM', num, element, ' ', ' ', 1, ' ', data(id, 1), data(id, 2), data(id, 3));
    end
end

fclose(fileID);