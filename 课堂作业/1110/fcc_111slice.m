clc
clear

Type = 'FCC';
element = 'Cu';
a0 = 3.615;
mass = 63.5460;

atom1 = [0.0 * a0, 0.0 * a0, 0.0 * a0];
atom2 = [0.5 * a0, 0.5 * a0, 0.0 * a0];
atom3 = [0.5 * a0, 0.0 * a0, 0.5 * a0];
atom4 = [0.0 * a0, 0.5 * a0, 0.5 * a0];

ux = [1.0 * a0, 0.0, 0.0];
uy = [0.0, 1.0 * a0, 0.0];
uz = [0.0, 0.0, 1.0 * a0];

id = 0;

formatSpec = '%-6s%-5d %-4s%-1s%-3s  %-4d%-1s   %-8.3f%-8.3f%-8.3f\n';
structure = 2;

%crystal = [1, atom1; 2, atom2; 3, atom3; 4, atom4];

%crystal = [origin;origin*[1,1,2]];

%plot3(origin(:,1), origin(:, 2), origin(:, 3), 'go');
%grid; box;

for i = -10:1:10
    for j = -10:1:10
        for k = -10:1:10
            vector = ux * i + uy * j + uz * k;
            id = id + 1;
            crystal(id, 1:3) = atom1 + vector;
            id = id + 1;
            crystal(id, 1:3) = atom2 + vector;
            %crystal = [crystal; id, atom2 + ux * i + uy * j + uz * k];
            id = id + 1;
            crystal(id, 1:3) = atom3 + vector;
            %crystal = [crystal; id, atom3 + ux * i + uy * j + uz * k];
            id = id + 1;
            crystal(id, 1:3) = atom4 + vector;
            %crystal = [crystal; id, atom4 + ux * i + uy * j + uz * k];
        end
    end
end

nt = id;
num = 0

if structure == 1
    fileID = fopen('FCC_Cu_single_plane.pdb', 'w');
    for id = 1:nt
        if crystal(id, 1:3) >= -10 & crystal(id, 1:3) <= 10
            r = crystal(id, 1) / 1 +  crystal(id, 2) / 1 + crystal(id, 3) / 1;
            if r > 7 & r < 10
                fprintf(fileID, formatSpec, 'ATOM', num, 'Cu', ' ', ' ', 1, ' ', crystal(id, 1), crystal(id, 2), crystal(id, 3));
            end
        end
    end
    fclose(fileID);
elseif structure == 2
    fileID = fopen('FCC_Cu_multi_plane.pdb', 'w');
    for id = 1:nt
        if crystal(id, 1:3) >= -10 & crystal(id, 1:3) <= 10
            r = crystal(id, 1) / 1 +  crystal(id, 2) / 1 + crystal(id, 3) / 1;
            if r > 7 & r < 10
                fprintf(fileID, formatSpec, 'ATOM', num, 'Cu', ' ', ' ', 1, ' ', crystal(id, 1), crystal(id, 2), crystal(id, 3));
            elseif r <= 7
                fprintf(fileID, formatSpec, 'ATOM', num, 'N', ' ', ' ', 1, ' ', crystal(id, 1), crystal(id, 2), crystal(id, 3));
            end
        end
    end
    fclose(fileID);
end
