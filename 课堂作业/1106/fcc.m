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

plot3(crystal(:, 1), crystal(:, 2), crystal(:, 3), 'o', 'MarkerFaceColor','g', 'MarkerSize', 10);
axis square;
