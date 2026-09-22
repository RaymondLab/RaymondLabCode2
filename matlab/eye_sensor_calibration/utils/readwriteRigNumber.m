function rigNum = readwriteRigNumber(num, readOnly)
%READWRITERIGNUMBER

filePath = fullfile(userpath, 'rigNum.txt');

if ~isfile(filePath) || ~readOnly
    rignumFile = fopen(filePath, 'w');
    fprintf(rignumFile, num);
    rigNum = num;
else
    rignumFile = fopen(filePath, 'r');
    rigNum = fscanf(rignumFile, '%s');
end

fclose(rignumFile);

end