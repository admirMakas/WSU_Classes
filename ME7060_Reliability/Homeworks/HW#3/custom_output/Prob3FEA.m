clear all
clc
%Area Moment of Inertia
Inrt = 1.33*10^3; % in^4
%Convert length L into inches
L=30*12; % in

%Convert units to lb/in
mux_w = 83.33; % lb/in
sigx_w = 8.33; %lb/in

%Next convert log normal statistics mux and sigx into muy and sigy
sigy_w = log((sigx_w/mux_w)^2 + 1);
muy_w = log(mux_w) - 0.5*sigy_w^2;

%Define statistics for concentrated load
mux_P = 50000; %lbf
sigx_P = 10000; %lbf

%Define statistics for Young's Modulus
mux_E = 29*10^6; %lb/in^2
sigx_E = 1*10^5; %lb/in^2

%Following 2 lines define number of samples(N) and variables(NVar)
N = 10000;
NVar = 3;

%Define random vector U(N, NVar) dimensions
U=rand(N,NVar);

data=[];

%Perform LHS method in loop below
for i = 1:NVar
    index=randperm(N);
    prob = (index'-U(:,i))/N;
    data(:,i)=prob;
end

%Next define distributions for variables P , E, and w
XP = mux_P + sigx_P*icdf('normal', data(:,1), 0, 1);

XE = mux_E + sigx_E*icdf('normal', data(:,2), 0, 1);

Xw = exp(muy_w + sigy_w*icdf('normal', data(:,3), 0, 1));

for j = 1:N
%Generate new input file for ABAQUS run------------------------------------
    fin = fopen('Job-1.inp','r'); %opens .inp file as read only
    fout = fopen('Job-1Test.inp','w'); %creted new file with .inp extention
    idk=0;
while ~feof(fin) %while loop to generate new input file
    idk=idk+1;
    s = fgetl(fin);
    if idk==72
        format = ' %s, 0.29\r\n';
        fprintf(fout, format, num2str(XE(j)));
    elseif idk==95
        format = '_PickedSurf8, P, %s\r\n';
        fprintf(fout, format, num2str(Xw(j)));
    elseif idk==98
        format = 'loadPSet, 2, %s\r\n';
        fprintf(fout, format, num2str(-XP(j)));
    else
        fprintf(fout,'%s\r\n',s); %prints new .inp file 
                                  %with the modified lines
    end
end
fclose(fin);
fclose(fout);
%--------------------------------------------------------------------------
%
%Run ABAQUS with the revised input file
%CANT HAVE SPACES IN FILE PATH!!!!!!!!!
command = ['abaqus job=BeamSol ' ...
'input=Job-1Test.inp ' ...
'interactive'];

[status,~] = dos(command)

%--------------------------------------------------------------------------
%
%Read the results from the .dat file
fid=fopen('BeamSol.dat');
linenum = 359; % Or whichever line you wish to read
fileLine = textscan(fid, '%s', 1, 'delimiter', '\n', ...
'headerlines', linenum-1);
fileLine = char(fileLine{1});
fileLine = fileLine(29:34);
%number = str2double(fileLine(15:end))
fileID = fopen('maxDisp.txt','a');
fprintf(fileID, '%s \r\n', fileLine);
fclose(fileID);
fclose(fid);
%
%Rinse and repeat
end