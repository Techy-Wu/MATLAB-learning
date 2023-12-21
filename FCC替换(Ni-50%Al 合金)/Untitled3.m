%% basic parameters
Type = 'FCC';
element_Ni = 'Ni  ';
element_Al = 'Al  ';
a_Ni = 3.52; % FCC-Ni�ľ�����
a_Al = 4.05; % FCC-Al�ľ�����
a_avg = (a_Ni + a_Al) / 2; % ƽ��������
mass_Ni = 58.69; % Ni��ԭ����
mass_Al = 26.98; % Al��ԭ����

%% ���徧����
ux = a_avg; uy = a_avg; uz = a_avg;

%% ����ԭ������
atom1 = [0.0*a_avg, 0.0*a_avg, 0.0*a_avg];
atom2 = [0.5*a_avg, 0.5*a_avg, 0.0*a_avg];
atom3 = [0.5*a_avg, 0.0*a_avg, 0.5*a_avg];
atom4 = [0.0*a_avg, 0.5*a_avg, 0.5*a_avg];

%% ��������
ID = 0; % ԭ�����
data = zeros(10*10*10*4, 3); % Ԥ����洢�ռ�
if (Type == 'FCC')
    for i = 0:9 % �� x ����ƽ�ƣ�����10������
        for j = 0:9 % �� y ����ƽ�ƣ�����10������
            for k = 0:9 % �� z ����ƽ�ƣ�����10������
                ID = ID + 1; data(ID,1:3) = a_avg*[i, j, k];
                ID = ID + 1; data(ID,1:3) = a_avg*[i+0.5, j+0.5, k];
                ID = ID + 1; data(ID,1:3) = a_avg*[i+0.5, j, k+0.5];
                ID = ID + 1; data(ID,1:3) = a_avg*[i, j+0.5, k+0.5];
            end
        end
    end
end
nt = ID; % ԭ������

%% ����û�50%��Niԭ�ӳ�ΪAlԭ��
%Ni_indices = find(data(:,1) == a_avg*atom1(1) & data(:,2) == a_avg*atom1(2) & data(:,3) == a_avg*atom1(3));
rand_list = randperm(ID);
Ni_indices = rand_list(1:ID / 2);
%Al_indices = Ni_indices(randperm(length(Ni_indices), round(length(Ni_indices)/2)), :);
Al_indices = rand_list(ID / 2 + 1:ID);
%data(Al_indices, :) = repmat([a_avg*atom3(1), a_avg*atom3(2), a_avg*atom3(3)], length(Al_indices), 1);

%% ͳ��Ni��Al��ԭ����
%num_Ni = sum(atoms(:) == 1); % ͳ��Niԭ����
%num_Al = sum(atoms(:) == 2); % ͳ��Alԭ����
num_Ni = ID / 2;
num_Al = ID / 2;
disp(['Niԭ������', num2str(num_Ni)]);
disp(['Alԭ������', num2str(num_Al)]);

%% ���Niԭ�����겢���Ϊ1��Alԭ�����겢���Ϊ2
fid = fopen('Ni_Al_coordinates.txt', 'w');
fprintf(fid, 'Ni atom coordinates:\n');
for i = 1:num_Ni
    fprintf(fid, '%d\t%d\t%d\t%d\n', i, data(Ni_indices(i), 1), data(Ni_indices(i), 2), data(Ni_indices(i), 3));
end
fprintf(fid, 'Al atom coordinates:\n');
for i = 1:num_Al
    fprintf(fid, '%d\t%d\t%d\t%d\n', i, data(Al_indices(i), 1), data(Al_indices(i), 2), data(Al_indices(i), 3));
end
fclose(fid);

%% ���PDB�ļ�
fid = fopen('Ni_Al_structure.pdb', 'w');
fprintf(fid, 'HEADER    FCC Ni/Al Composite\n');
fprintf(fid, 'REMARK    Created by MATLAB\n');
for i = 1:nt
    if ismember(i, Ni_indices)
        element = element_Ni;
    else
        element = element_Al;
    end
    fprintf(fid, '%-6s%5d  %-3s%s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n', 'ATOM', i, element, 'MOL', 1, data(i,1), data(i,2), data(i,3), 1.00, 0.00, element(1:2));
end
fprintf(fid, 'TER\n');
fprintf(fid, 'END\n');
fclose(fid);