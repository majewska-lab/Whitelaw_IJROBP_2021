
% Script to analyze coordinates (2D) for longitudinal post-radiation
% analysis

% Choose file and import table
input_dir = uigetdir(); % Select directory with results in .csv
input_list = dir(fullfile(input_dir,'*.csv')); % Grab filenames 
input_names = {input_list.name};
input_names = strtrim(input_names);
input_path = fullfile(input_dir,input_names);

% Create directory to store results:
NN_results_dir = fullfile(input_dir,'NN_results');
mkdir(NN_results_dir);
% Generate file names for NN data to be saved
NN_results_names = strrep(input_names,'.csv','_NN.mat');
NN_results_path = fullfile(NN_results_dir,NN_results_names);

% Input experimental parameters:
time_pts = [-2 0.25 1 7 16]; % pre-define list of timepoints
%time_pts = [-2 2 7 16];
cutoff = 10; % cutoff distance (um) for cell classification
scale = 0.667; % Define scale in um per pixel



Dv1_all = {};
Dvi_all = {};

for i = 1:size(input_names,2)
    % Run NN function to calculate distances (vs 1 or i), and cell classification 
    [x, y, Dv1, Dvi, n_cells_norm, frac_pers, frac_new] = NN_2D_bsw(input_path{i},time_pts,cutoff,scale);
    % Save calculations to new .mat file
    Dv1_all = cat(1,Dv1_all,Dv1);
    Dvi_all = cat(1,Dvi_all,Dvi);
    save(NN_results_path{i},'x','y','Dv1','Dvi','n_cells_norm','frac_pers','frac_new'); 
end 

% Concatenating NN distances for given time points:

% Initilize cell array with empty vectors for distances to go into
Dv1_cat = {};
Dvi_cat = {};

for i=1:size(time_pts,2) 
    Dv1_cat{i} = []; 
    Dvi_cat{i} = [];
end

cells_norm = [];
fracs_pers = [];
fracs_new = [];

for i=1:size(input_names,2) % i is for each sample
    load(NN_results_path{i},'Dv1','Dvi','n_cells_norm','frac_pers','frac_new');
    % Add cell NN distance data for each timepoint into concatenated vector
    for j=1:size(Dv1,2) %j is for each timepoint 
        Dv1_cat{j} = cat(1,Dv1_cat{j},Dv1{j});
        Dvi_cat{j} = cat(1,Dvi_cat{j},Dvi{j});
    end
    % Get mean distance, and fraction cells persisted or new into vectors
    cells_norm = cat(1,cells_norm,n_cells_norm);
    fracs_pers = cat(1,fracs_pers,frac_pers);
    fracs_new = cat(1,fracs_new,frac_new);
end

save(fullfile(input_dir,'NNresults.mat'),'Dv1_all','Dvi_all','Dv1_cat','Dvi_cat','cells_norm','fracs_pers','fracs_new','input_names','time_pts')


function [x, y, Dv1, Dvi, n_cells_norm, frac_pers, frac_new] = NN_2D_bsw(file_path,time_pts,cutoff,scale)

% Function to extract cell location data, calculate and store:
% Coordinates: x,y: cell array with each element corresponding to vector
% for one time point
% Nearest-neighbor distance: 
    % Dv1 (compared to t=1): for old/new classifation
    % Dvi (compared to t=i): for persistent/gone classifcation
% Number of cells normalized to first time point
    % n_cells_norm
% Cell classification results:
    % frac_pers: persistent cells
    % frac_new: new cells


mytable = readtable(file_path);

num_t = size(time_pts,2);
x_coords = mytable.X;
y_coords = mytable.Y;
marker = mytable.Type;

% Make cell array and store time-specific coordinates there
% Also store number of cells
x = {};
y = {};
n_cells = [];
for i = 1:num_t
    ind = marker == i;
    x{i} = scale*x_coords(ind); % Convert from pixels to microns
    y{i} = scale*y_coords(ind);
    n_cells(i) = sum(ind,'all');
end

n_cells_norm = n_cells/n_cells(1);

% Calculate nearest neighbors and mean NN distance:
% v1 = vs. t=1... for every cell at t=i, how close is the nearest in t=1
    % used for old/new classification and also NN cdf plots
Dv1 = {}; 
for i = 1:num_t
    [~, Dv1{i}] = knnsearch([x{1}, y{1}], [x{i}, y{i}]); 
    % suppress index output. not needed.
end

% vi = vs. t=i... for every cell at t=1, how close is the nearest in t=i
    % used for persistent/dropped classification
Dvi = {};
for i =1:num_t
    [~, Dvi{i}] = knnsearch([x{i},y{i}],[x{1},y{1}]);
end
 

% Calculate persistent/dropped cells: cells in t=1 that are <= or > 16px
% (~10 um; or cutoff) from nearest cell in t=i
for i = 1:num_t
    bigD1{i} = Dvi{i} > cutoff; % Build logical index for dropped vs. persistent cells
    smaD1{i} = Dvi{i} <= cutoff; % Big (dropped) vs small (persistance) distance
    x_drop{i} = x{1}(bigD1{i}); % Coordinates of dropped cells
    y_drop{i} = y{1}(bigD1{i});
    x_pers{i} = x{1}(smaD1{i}); % Coordinates of Persistent cells
    y_pers{i} = y{1}(smaD1{i});
end
% Calculate old/new cells: cells in t=i that are <= or > ~10 um from
% nearest cell in t=1
for i = 1:num_t
    bigDi{i} = Dv1{i} > cutoff; % Log Ind for new (big distance) vs old (small distance) cells 
    smaDi{i} = Dv1{i} <= cutoff;
    x_new{i} = x{i}(bigDi{i}); % New cell coordinates
    y_new{i} = y{i}(bigDi{i});
    x_old{i} = x{i}(smaDi{i}); % Old cell coordinates
    y_old{i} = y{i}(smaDi{i});
end


% Calculate the fraction of persistent cells and fraction of new cells 
% at each timepoint
frac_pers = [];
frac_new = [];

for i = 1:num_t
    frac_pers(i) = size(x_pers{i},1)/size(x{1},1);
    frac_new(i) = size(x_new{i},1)/size(x{i},1);
end


end


