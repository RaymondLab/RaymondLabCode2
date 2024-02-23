%% Save Analysis Information

vars.mag1 = mag1;
vars.mag2 = mag2;
vars.vid = vid;

pathname = cd;
[~, filenameroot]= fileparts(pathname);
save(fullfile(cd, [filenameroot '_analysis.mat']), 'mag1', 'mag2', 'vid');