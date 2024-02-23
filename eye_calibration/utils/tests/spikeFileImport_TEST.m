[~, filenameroot]= fileparts(cd);
fullfilename = fullfile(cd, [filenameroot,'.smr']);

try
    rawMagnetData = importSpike(fullfilename, [4 5 6 10 31]);
catch
    rawMagnetData = importSpike(fullfilename, [2 7 8 11]);
end