% fileID = fopen('exp.txt','w');
% A1 = [9.9, 9900];
% A2 = [8.8,  7.7 ; ...
%       8800, 7700];
% %to move to new line typically can use '\n' but for
% %notepad need to use '\r\n' instead
% formatSpec = 'X is %4.2f meters or %8.3f mm\r\n';
% fprintf(fileID, formatSpec , A1, A2);
% fclose(fileID);
%%
%Read value from the .dat file (in this case max displacement)
fid=fopen('Job1.dat');
linenum = 359; % Or whichever line you wish to read
%the function 'textscan' is used to exctract certain stings from a
%predetermined line in a text file. Here we are storing lines of text from 'fid'
%starting at line 308 and below. Next only line 308 is extracted from fileLine by
%using char(fileLine{1})
%then the number is stored under 'number' by converting it from sting to double
%(15:end) will take all the values from those locations and store as number.
fileLine = textscan(fid, '%s', 1, 'delimiter', '\n', ...
'headerlines', linenum-1);
fileLine = char(fileLine{1});
fileLine = fileLine(29:34);
%number = str2double(fileLine(15:end))
fileID = fopen('exp.txt','a');
fprintf(fileID, '%s \n', fileLine);
fclose(fileID);
%%
E=3.9e6;
P=-55000;
w=-95.0;
%Code that reads and modifies the abaqus input file
fin = fopen('Job-1.inp','r'); %opens .inp file as read only
fout = fopen('Job-1Test.inp','w'); %creted new file with .inp extention
idk=0;
while ~feof(fin) %while
    idk=idk+1;
    s = fgetl(fin);
    if idk==72
        format = '%s, 0.29\n'
        fprintf(fout, format, num2str(E))
    elseif idk==95
        format = '_PickedSurf8, P, %s\n'
        fprintf(fout, format, num2str(w))
    elseif idk==98
        format = 'loadPSet, 2, %s\n'
        fprintf(fout, format, num2str(P))
    else
        fprintf(fout,'%s\n',s); %prints new .inp file with the modified lines
    end
    end
fclose(fin);
fclose(fout);