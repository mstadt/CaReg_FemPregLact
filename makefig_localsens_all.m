% postprocess output from all the local sensitivity analysis
% output from compute_localsensitivity.m
clear all;

% male
fname = './results_localsens/15-May-2023_localsens_sexORrep-male_notes-malerun.mat';
male_dat = load(fname);

% female
fname = './results_localsens/15-May-2023_localsens_sexORrep-female_notes-femalerun.mat';
female_dat = load(fname);

% preg
fname = './results_localsens/15-May-2023_localsens_sexORrep-preg_notes-pregrun.mat';
preg_dat = load(fname);

% lact
fname = './results_localsens/15-May-2023_localsens_sexORrep-lact_notes-lactrun.mat';
lact_dat = load(fname);

[male_xlabs, male_vals] = getdata(male_dat);
[female_xlabs, female_vals] = getdata(female_dat);
[preg_xlabs, preg_vals] = getdata(preg_dat);
[lact_xlabs, lact_vals] = getdata(lact_dat);

% sanity check: check that xlabels are the same
flag = 0;
if size(male_xlabs) ~= size(female_xlabs)
    flag = 1;
elseif size(male_xlabs) ~= size(preg_xlabs)
    flag = 1;
elseif size(male_xlabs) ~= size(lact_xlabs)
    flag = 1;
end
for ii = 1:size(male_xlabs,2)
    if ~strcmp(male_xlabs{ii}, female_xlabs{ii})
        flag = 1;
    elseif ~strcmp(male_xlabs{ii}, preg_xlabs{ii})
        flag = 1;
    elseif ~strcmp(male_xlabs{ii}, lact_xlabs{ii})
        flag = 1;
    end
end
if flag
    error('xlabels do not match')
end
xlabels = male_xlabs;

% Make figures
temp = [male_vals, female_vals, preg_vals, lact_vals];

clims = [min(temp, [], 'all'), max(temp, [], 'all')];
ylabels = {'[PTH]_p', '[Ca^{2+}]_p', '[1,25(OH)_2D_3]_p'};
fsize = 14;
colmap = parula;
cmissdat = 'w';
labmissdat = '<0.5%';


% figure without removing 
fig_noremove = 0;
if fig_noremove
figure(1)
clf
% male fig
subplot(4,1,1)
h1 = heatmap(xlabels, ylabels, male_vals,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Male')
% female fig
subplot(4,1,2)
h1 = heatmap(xlabels, ylabels, female_vals,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'on');
h1.FontSize = fsize;
ylabel('Female')
% preg fig
subplot(4,1,3)
h1 = heatmap(xlabels, ylabels, preg_vals,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Pregnancy')
% lact fig
subplot(4,1,4)
h1 = heatmap(xlabels, ylabels, lact_vals,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Lactation')
end



% remove Nan values
Nan_male = find_all_Nan(male_vals);
Nan_female = find_all_Nan(female_vals);
temp = intersect(Nan_male, Nan_female);
Nan_preg = find_all_Nan(preg_vals);
temp = intersect(temp, Nan_preg);
Nan_lact = find_all_Nan(lact_vals);
AllNan = intersect(temp, Nan_lact);

[xlabs_rm, male_vals_rm] = removeAllNan(AllNan, xlabels, male_vals);
[~, female_vals_rm] = removeAllNan(AllNan, xlabels, female_vals);
[~, preg_vals_rm] = removeAllNan(AllNan, xlabels, preg_vals);
[~, lact_vals_rm] = removeAllNan(AllNan, xlabels, lact_vals);

% figure with removing 
fig_remove = 1;
if fig_remove
figure(2)
clf
% male fig
subplot(4,1,1)
h1 = heatmap(xlabs_rm, ylabels, male_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Male')
% remove xlabels
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% female fig
subplot(4,1,2)
h1 = heatmap(xlabs_rm, ylabels, female_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Female')
% remove xlabels
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% preg fig
subplot(4,1,3)
h1 = heatmap(xlabs_rm, ylabels, preg_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'on');
h1.FontSize = fsize;
ylabel('Pregnancy')
% remove xlabels
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% lact fig
subplot(4,1,4)
h1 = heatmap(xlabs_rm, ylabels, lact_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Lactation')

sgtitle('Local sensitivity analysis', 'fontsize', 24)
end

%------------------
% functions used
%------------------
function [xlabels, round_data] = getdata(dat)
    frac_sens = dat.frac_change;
    labels = cell(size(dat.param_names));
    for ii = 1:size(dat.param_names, 2)
        labels{ii} = convert_param_name(dat.param_names{ii});
    end
    xlabels = labels;
    round_data = round(frac_sens', 2, 'significant');

    [r,c] = find(abs(round_data) <= 0.5);% r - row values, c - column value
    for ii = 1:length(r)
        round_data(r(ii),c(ii)) = NaN;
    end
end

function AllNan_vals = find_all_Nan(round_data)
    % finds which indices are all Nan values
    PTH_Nan = find(isnan(round_data(1,:)));
    Ca_Nan  = find(isnan(round_data(2,:)));
    temp = intersect(PTH_Nan, Ca_Nan);
    D3_Nan  = find(isnan(round_data(3,:)));
    AllNan_vals = intersect(temp, D3_Nan);
end

function [rmNan_xlabels, rmNan_round_data] = removeAllNan(AllNan_vals, xlabels, round_data)
    % removes columns listed in indices from AllNan_vals
    rmNan_round_data = round_data;
    rmNan_xlabels = xlabels;
    
    rmNan_round_data(:,AllNan_vals) = [];
    for ii = length(AllNan_vals):-1:1
        rmNan_xlabels(AllNan_vals(ii)) = [];
    end
end