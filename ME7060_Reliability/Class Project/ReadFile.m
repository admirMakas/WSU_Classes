function [ FN ] = ReadFile( fileName )

fid = fopen(fileName,'r');
NodeOutput = '                                              R E A L   E I G E N V A L U E S';
tline = fgetl(fid);
FN = zeros(1,3);

while ischar(tline)
    
    if isempty(strfind(tline,NodeOutput)) == 0
        for i = 1:3
            tline = fgetl(fid);
        end 
        B = textscan(tline,'%f %f %f %f %f %f %f');
        FN(1) = B{4}(1);
        
        B = textscan(fgetl(fid),'%f %f %f %f %f %f %f');
        FN(2) = B{4}(1);
        
        B = textscan(fgetl(fid),'%f %f %f %f %f %f %f');
        FN(3) = B{4}(1);
        
        break
        
    else
        tline = fgetl(fid);
    end
end
fclose(fid);

end