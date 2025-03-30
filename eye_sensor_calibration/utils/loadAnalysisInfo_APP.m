%% load analysis information
try
    % don't load from a file if you don't have to
    mag1 = vars.mag1;
    mag2 = vars.mag2;
    vid = vars.vid;
catch
    try
        pathname = cd;
        [~, filenameroot]= fileparts(pathname);
        load(fullfile(cd, [filenameroot '_analysis.mat']))
    catch
        error('Unable to read desaccading data')
    end
end