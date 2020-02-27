function loaddeps()
    [current, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(current, 'condense'));
    addpath(fullfile(current, 'condense/quality'));
    addpath(fullfile(current, 'centrality'));
end