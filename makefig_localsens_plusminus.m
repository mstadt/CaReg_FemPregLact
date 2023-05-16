% postprocess output from all the local sensitivity analysis
% output from compute_localsensitivity.m
clear all;

% male
fname = './results_localsens/16-May-2023_localsens_plusminussexORrep-male_notes-plusminus.mat';
male_dat = load(fname);

% female
fname = './results_localsens/16-May-2023_localsens_plusminussexORrep-female_notes-plusminus.mat';
female_dat = load(fname);

% preg
fname = './results_localsens/16-May-2023_localsens_plusminussexORrep-preg_notes-plusminus.mat';
preg_dat = load(fname);

% lact
fname = './results_localsens/16-May-2023_localsens_plusminussexORrep-lact_notes-plusminus.mat';
lact_dat = load(fname);


%%
% increase by 5% sensitivity
[male_xlabs_plus, male_vals_plus] = getdata_plus(male_dat);
[female_xlabs_plus, female_vals_plus] = getdata_plus(female_dat);
[preg_xlabs_plus, preg_vals_plus] = getdata_plus(preg_dat);
[lact_xlabs_plus, lact_vals_plus] = getdata_plus(lact_dat);

% sanity check: check that xlabels are the same
flag = 0;
if size(male_xlabs_plus) ~= size(female_xlabs_plus)
    flag = 1;
elseif size(male_xlabs_plus) ~= size(preg_xlabs_plus)
    flag = 1;
elseif size(male_xlabs_plus) ~= size(lact_xlabs_plus)
    flag = 1;
end
for ii = 1:size(male_xlabs_plus,2)
    if ~strcmp(male_xlabs_plus{ii}, female_xlabs_plus{ii})
        flag = 1;
    elseif ~strcmp(male_xlabs_plus{ii}, preg_xlabs_plus{ii})
        flag = 1;
    elseif ~strcmp(male_xlabs_plus{ii}, lact_xlabs_plus{ii})
        flag = 1;
    end
end
if flag
    error('xlabels do not match')
end
xlabels = male_xlabs_plus;

% Make figures
temp = [male_vals_plus, female_vals_plus, preg_vals_plus, lact_vals_plus];

clims = [min(temp, [], 'all'), max(temp, [], 'all')];
ylabels = {'[PTH]_p', '[Ca^{2+}]_p', '[1,25(OH)_2D_3]_p'};
fsize = 14;
colmap = parula;
cmissdat = 'w';
labmissdat = '<1.0%';


% figure without removing 
% fig_noremove = 0;
% if fig_noremove
% figure(1)
% clf
% % male fig
% subplot(4,1,1)
% h1 = heatmap(xlabels, ylabels, male_vals_plus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'off');
% h1.FontSize = fsize;
% ylabel('Male')
% % female fig
% subplot(4,1,2)
% h1 = heatmap(xlabels, ylabels, female_vals_plus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'on');
% h1.FontSize = fsize;
% ylabel('Female')
% % preg fig
% subplot(4,1,3)
% h1 = heatmap(xlabels, ylabels, preg_vals_plus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'off');
% h1.FontSize = fsize;
% ylabel('Pregnancy')
% % lact fig
% subplot(4,1,4)
% h1 = heatmap(xlabels, ylabels, lact_vals_plus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'off');
% h1.FontSize = fsize;
% ylabel('Lactation')
% end



% remove Nan values
Nan_male = find_all_Nan(male_vals_plus);
Nan_female = find_all_Nan(female_vals_plus);
temp = intersect(Nan_male, Nan_female);
Nan_preg = find_all_Nan(preg_vals_plus);
temp = intersect(temp, Nan_preg);
Nan_lact = find_all_Nan(lact_vals_plus);
AllNan = intersect(temp, Nan_lact);

[xlabs_rm, male_vals_rm] = removeAllNan(AllNan, xlabels, male_vals_plus);
[~, female_vals_rm] = removeAllNan(AllNan, xlabels, female_vals_plus);
[~, preg_vals_rm] = removeAllNan(AllNan, xlabels, preg_vals_plus);
[~, lact_vals_rm] = removeAllNan(AllNan, xlabels, lact_vals_plus);

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
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% female fig
subplot(4,1,2)
h1 = heatmap(xlabs_rm, ylabels, female_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Female')
% remove xlabels
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% preg fig
subplot(4,1,3)
h1 = heatmap(xlabs_rm, ylabels, preg_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'on');
h1.FontSize = fsize;
ylabel('Pregnancy')
% remove xlabels
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% lact fig
subplot(4,1,4)
h1 = heatmap(xlabs_rm, ylabels, lact_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Lactation')

sgtitle('Local sensitivity analysis - 5% increase', 'fontsize', 24)
end

%%
% decrease by 5% sensitivity
[male_xlabs_minus, male_vals_minus] = getdata_minus(male_dat);
[female_xlabs_minus, female_vals_minus] = getdata_minus(female_dat);
[preg_xlabs_minus, preg_vals_minus] = getdata_minus(preg_dat);
[lact_xlabs_minus, lact_vals_minus] = getdata_minus(lact_dat);

% sanity check: check that xlabels are the same
flag = 0;
if size(male_xlabs_minus) ~= size(female_xlabs_minus)
    flag = 1;
elseif size(male_xlabs_minus) ~= size(preg_xlabs_minus)
    flag = 1;
elseif size(male_xlabs_minus) ~= size(lact_xlabs_minus)
    flag = 1;
end
for ii = 1:size(male_xlabs_minus,2)
    if ~strcmp(male_xlabs_minus{ii}, female_xlabs_minus{ii})
        flag = 1;
    elseif ~strcmp(male_xlabs_minus{ii}, preg_xlabs_minus{ii})
        flag = 1;
    elseif ~strcmp(male_xlabs_minus{ii}, lact_xlabs_minus{ii})
        flag = 1;
    end
end
if flag
    error('xlabels do not match')
end
xlabels = male_xlabs_minus;

% Make figures
temp = [male_vals_minus, female_vals_minus, preg_vals_minus, lact_vals_minus];

clims = [min(temp, [], 'all'), max(temp, [], 'all')];
ylabels = {'[PTH]_p', '[Ca^{2+}]_p', '[1,25(OH)_2D_3]_p'};
fsize = 14;
colmap = parula;
cmissdat = 'w';
labmissdat = '<1.0%';


% figure without removing 
% fig_noremove = 0;
% if fig_noremove
% figure(1)
% clf
% % male fig
% subplot(4,1,1)
% h1 = heatmap(xlabels, ylabels, male_vals_minus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'off');
% h1.FontSize = fsize;
% ylabel('Male')
% % female fig
% subplot(4,1,2)
% h1 = heatmap(xlabels, ylabels, female_vals_minus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'on');
% h1.FontSize = fsize;
% ylabel('Female')
% % preg fig
% subplot(4,1,3)
% h1 = heatmap(xlabels, ylabels, preg_vals_minus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'off');
% h1.FontSize = fsize;
% ylabel('Pregnancy')
% % lact fig
% subplot(4,1,4)
% h1 = heatmap(xlabels, ylabels, lact_vals_minus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'off');
% h1.FontSize = fsize;
% ylabel('Lactation')
% end



% remove Nan values
Nan_male = find_all_Nan(male_vals_minus);
Nan_female = find_all_Nan(female_vals_minus);
temp = intersect(Nan_male, Nan_female);
Nan_preg = find_all_Nan(preg_vals_minus);
temp = intersect(temp, Nan_preg);
Nan_lact = find_all_Nan(lact_vals_minus);
AllNan = intersect(temp, Nan_lact);

[xlabs_rm, male_vals_rm] = removeAllNan(AllNan, xlabels, male_vals_minus);
[~, female_vals_rm] = removeAllNan(AllNan, xlabels, female_vals_minus);
[~, preg_vals_rm] = removeAllNan(AllNan, xlabels, preg_vals_minus);
[~, lact_vals_rm] = removeAllNan(AllNan, xlabels, lact_vals_minus);

% figure with removing 
fig_remove = 1;
if fig_remove
figure(3)
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
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% female fig
subplot(4,1,2)
h1 = heatmap(xlabs_rm, ylabels, female_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Female')
% remove xlabels
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% preg fig
subplot(4,1,3)
h1 = heatmap(xlabs_rm, ylabels, preg_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'on');
h1.FontSize = fsize;
ylabel('Pregnancy')
% remove xlabels
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% lact fig
subplot(4,1,4)
h1 = heatmap(xlabs_rm, ylabels, lact_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Lactation')

sgtitle('Local sensitivity analysis - 5% decrease', 'fontsize', 24)
end

%%
% increase and decrease by 5% sensitivity
[male_xlabs_plusminus, male_vals_plusminus] = getdata_plusminus(male_dat);
[female_xlabs_plusminus, female_vals_plusminus] = getdata_plusminus(female_dat);
[preg_xlabs_plusminus, preg_vals_plusminus] = getdata_plusminus(preg_dat);
[lact_xlabs_plusminus, lact_vals_plusminus] = getdata_plusminus(lact_dat);

% sanity check: check that xlabels are the same
flag = 0;
if size(male_xlabs_plusminus) ~= size(female_xlabs_plusminus)
    flag = 1;
elseif size(male_xlabs_plusminus) ~= size(preg_xlabs_plusminus)
    flag = 1;
elseif size(male_xlabs_plusminus) ~= size(lact_xlabs_plusminus)
    flag = 1;
end
for ii = 1:size(male_xlabs_plusminus,2)
    if ~strcmp(male_xlabs_plusminus{ii}, female_xlabs_plusminus{ii})
        flag = 1;
    elseif ~strcmp(male_xlabs_plusminus{ii}, preg_xlabs_plusminus{ii})
        flag = 1;
    elseif ~strcmp(male_xlabs_plusminus{ii}, lact_xlabs_plusminus{ii})
        flag = 1;
    end
end
if flag
    error('xlabels do not match')
end
xlabels = male_xlabs_plusminus;

% Make figures
temp = [male_vals_plusminus, female_vals_plusminus, preg_vals_plusminus, lact_vals_plusminus];

clims = [min(temp, [], 'all'), max(temp, [], 'all')];
ylabels = {'[PTH]_p', '[Ca^{2+}]_p', '[1,25(OH)_2D_3]_p'};
fsize = 14;
colmap = parula;
cmissdat = 'w';
labmissdat = '<1.0%';


% figure without removing 
% fig_noremove = 0;
% if fig_noremove
% figure(1)
% clf
% % male fig
% subplot(4,1,1)
% h1 = heatmap(xlabels, ylabels, male_vals_plusminus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'off');
% h1.FontSize = fsize;
% ylabel('Male')
% % female fig
% subplot(4,1,2)
% h1 = heatmap(xlabels, ylabels, female_vals_plusminus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'on');
% h1.FontSize = fsize;
% ylabel('Female')
% % preg fig
% subplot(4,1,3)
% h1 = heatmap(xlabels, ylabels, preg_vals_plusminus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'off');
% h1.FontSize = fsize;
% ylabel('Pregnancy')
% % lact fig
% subplot(4,1,4)
% h1 = heatmap(xlabels, ylabels, lact_vals_plusminus,...
%                 'colormap', colmap,...
%                 'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
%                 'ColorLimits', clims, 'colorbarvisible', 'off');
% h1.FontSize = fsize;
% ylabel('Lactation')
% end



% remove Nan values
Nan_male = find_all_Nan(male_vals_plusminus);
Nan_female = find_all_Nan(female_vals_plusminus);
temp = intersect(Nan_male, Nan_female);
Nan_preg = find_all_Nan(preg_vals_plusminus);
temp = intersect(temp, Nan_preg);
Nan_lact = find_all_Nan(lact_vals_plusminus);
AllNan = intersect(temp, Nan_lact);

[xlabs_rm, male_vals_rm] = removeAllNan(AllNan, xlabels, male_vals_plusminus);
[~, female_vals_rm] = removeAllNan(AllNan, xlabels, female_vals_plusminus);
[~, preg_vals_rm] = removeAllNan(AllNan, xlabels, preg_vals_plusminus);
[~, lact_vals_rm] = removeAllNan(AllNan, xlabels, lact_vals_plusminus);

% figure with removing 
fig_remove = 1;
if fig_remove
figure(4)
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
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% female fig
subplot(4,1,2)
h1 = heatmap(xlabs_rm, ylabels, female_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Female')
% remove xlabels
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% preg fig
subplot(4,1,3)
h1 = heatmap(xlabs_rm, ylabels, preg_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'on');
h1.FontSize = fsize;
ylabel('Pregnancy')
% remove xlabels
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% lact fig
subplot(4,1,4)
h1 = heatmap(xlabs_rm, ylabels, lact_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
ylabel('Lactation')

sgtitle('Local sensitivity analysis - 5% increase and decrease', 'fontsize', 24)
end

%------------------
% functions used
%------------------
function [xlabels, round_data] = getdata_plus(dat)
    frac_sens = dat.frac_change_plus;
    labels = cell(size(dat.param_names));
    for ii = 1:size(dat.param_names, 2)
        labels{ii} = convert_param_name(dat.param_names{ii});
    end
    xlabels = labels;
    round_data = round(frac_sens', 2, 'significant');

    [r,c] = find(abs(round_data) <= 1.0);% r - row values, c - column value
    for ii = 1:length(r)
        round_data(r(ii),c(ii)) = NaN;
    end
end

function [xlabels, round_data] = getdata_minus(dat)
    frac_sens = dat.frac_change_minus;
    labels = cell(size(dat.param_names));
    for ii = 1:size(dat.param_names, 2)
        labels{ii} = convert_param_name(dat.param_names{ii});
    end
    xlabels = labels;
    round_data = round(frac_sens', 2, 'significant');

    [r,c] = find(abs(round_data) <= 1.0);% r - row values, c - column value
    for ii = 1:length(r)
        round_data(r(ii),c(ii)) = NaN;
    end
end

function [xlabels, round_data] = getdata_plusminus(dat)
    frac_sens = dat.frac_change_plusminus;
    labels = cell(size(dat.param_names));
    for ii = 1:size(dat.param_names, 2)
        labels{ii} = convert_param_name(dat.param_names{ii});
    end
    xlabels = labels;
    round_data = round(frac_sens', 2, 'significant');

    [r,c] = find(abs(round_data) <= 1.0);% r - row values, c - column value
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