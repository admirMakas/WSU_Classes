function GenerateNewJob(CASE, newFilePath, X,newFileName)

readIn = fopen(CASE,'r');
WriteOut = fopen(newFilePath, 'w');

i = 1;
tline = fgetl(readIn);

while ischar(tline)
    % Double Check the line numbers in windows
    if i == 0 % E
        fprintf(WriteOut,' %s, 0.29 \n', num2str(X{2,2}) );
        
    elseif i == 362 % K  PBUSH         10       K1.000+151.000+151.000+151.000+151.000+151.000+15
        fprintf(WriteOut,'PBUSH*, 10, K,  1.000+15, 1.000+15, %f, 1.000+15, 1.000+15, 1.000+15\n', X(1));
        
    elseif i == 367 %E P   MAT1,2, 75000.00, ,0.330000, 2.7110-6
        fprintf(WriteOut,'MAT1*, 2, %f, , %f, 2.7110-6 \n', X(5), X(4));
        
    else
        fprintf(WriteOut,'%s \n',tline);
    
    end
    i = i+1;
    tline = fgetl(readIn);
end

fclose(readIn);
fclose(WriteOut);

display(sprintf('Finished Creating file %s.', newFileName))

end