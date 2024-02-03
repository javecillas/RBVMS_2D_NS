close all
clear all
clc

path{1} = 'raw_bcs.txt';

data = fileread(path{1});
data = regexp(data,'\n','split');

temp = contains( data , 'PtDisplacement' );
idx = find(temp == 1);

bcs = zeros(length(idx),7);

for i = 1:length(idx)
    % Node ID
    strDel = strcat('PtDisplacement PD',num2str(i),' {Node = ');
    temp2 = erase(char(data(1,idx(i))),strDel);
    bcs(i,1) = str2double(temp2);
    % BC val
    strDel = strcat('                Displacement = ');
    temp2 = erase(char(data(1,idx(i)+1)),strDel);
    bcs(i,2:3) = str2num(char(temp2));
    % Flag
    strDel = strcat('                Flags = ');
    temp2 = erase(char(data(1,idx(i)+2)),strDel);
    temp2 = erase(temp2,'}');
    bcs(i,5:6) = str2num(char(temp2));
end

writematrix(bcs,'bcs.txt','Delimiter','tab');
type 'bcs.txt';

