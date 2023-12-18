clc
clear

Type = 'HCP';
element = 'Mg';
a0 = 3.209;
c0 = a0 * 1.633;
mass = 50.0;

atom1 = [0.0, 0.0, 0.0];
atom2 = [0.5 * a0, sqrt(3) / 2 * a0, 0.0];
atom3 = [0.0, 2.0 / sqrt(3) * a0, 0.5 * c0];
atom4 = [-0.5 * a0, -a0 / sqrt(3) / 2, 0.5 * c0];

origin = [atom1; atom2; atom3; atom4];

ux = [a0, 0.0, 0.0];
uy = [0.0, sqrt(3) * a0, 0.0];
uz = [0.0, 0.0, c0];

formatSpec = '%-6s%-5d %-4s%-1s%-3s  %-4d%-1s   %-8.3f%-8.3f%-8.3f\n';

id = 0;

fileID = fopen('Mg_lattice.pdb', 'w');
num = 0;

for i = -10:10
    for j = -10:10
        for k = -10:10
            vector = ux * i + uy * j + uz * k;
            for m = 1:4
                id = id + 1;
                lattice(id, 1:3) = origin(m) + vector;
                % if m == 1 || m == 2
                %     plot3(lattice(id, 1), lattice(id, 2), lattice(id, 3), 'o', 'MarkerFaceColor', 'r');
                %     hold on;
                % else
                %     plot3(lattice(id, 1), lattice(id, 2), lattice(id, 3), 'o', 'MarkerFaceColor', 'g');
                %     hold on;
                % end

                if lattice(id, 1:3) >= 0.0
                    if lattice(id, 1:3) <= 30.0
                        num = num + 1;
                        if m == 3 || m == 4
                            fprintf(fileID, formatSpec, 'ATOM', num, element, ' ', ' ', 1, ' ', lattice(id, 1), lattice(id, 2), lattice(id, 3));
                        else
                            if mod(k, 2) == 0
                                fprintf(fileID, formatSpec, 'ATOM', num, 'N', ' ', ' ', 1, ' ', lattice(id, 1), lattice(id, 2), lattice(id, 3));
                            else
                                fprintf(fileID, formatSpec, 'ATOM', num, 'H', ' ', ' ', 1, ' ', lattice(id, 1), lattice(id, 2), lattice(id, 3));
                            end
                        end
                    end
                end
            end
        end
    end
end

plot3(lattice(:, 1), lattice(:, 2), lattice(:, 3), 'o', 'MarkerFaceColor', 'g');

nt = id;

% fileID = fopen('Mg_lattice.pdb', 'w');
% num = 0;
% 
% for id = 1:nt
%     if lattice(id, 1:3) >= 0.0
%         if lattice(id, 1:3) <= 30.0
%         num = num + 1;
%             if mod(lattice(id, 3), 1.605) == 0
%                 fprintf(fileID, formatSpec, 'ATOM', num, element, ' ', ' ', 1, ' ', lattice(id, 1), lattice(id, 2), lattice(id, 3));
%             elseif mod(lattice(id, 3) - 1.605, 1.605) == 0
%                 fprintf(fileID, formatSpec, 'ATOM', num, 'N', ' ', ' ', 1, ' ', lattice(id, 1), lattice(id, 2), lattice(id, 3));
%             else
%                 fprintf(fileID, formatSpec, 'ATOM', num, 'C', ' ', ' ', 1, ' ', lattice(id, 1), lattice(id, 2), lattice(id, 3));
%             end
%         end
%     end
% end

fclose(fileID);