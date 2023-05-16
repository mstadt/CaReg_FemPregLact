% This file is used to postprocess output from
% compute_preg2fem_lact2fem_all.m

% Makes Figure 3.5 from manuscript
clear all
% file where results for preg 2 virgin, lact 2 virgin results are saved
fname = './results_preglact_sensitivity/16-May-2023_preg2fem_lact2fem_all_notes-PTHchange.mat';
dat = load(fname);

[preg_xlabs, preg_data] = getdata(dat,'preg');
[lact_xlabs, lact_data] = getdata(dat,'lact');

% sanity check: check that xlabels are the same
flag = 0;
if size(preg_xlabs) ~= size(lact_xlabs)
    flag = 1;
end
for ii = 1:size(preg_xlabs,2)
    if ~strcmp(preg_xlabs, lact_xlabs)
        flag = 1;
    end
end
if flag
    error('xlabels do not match')
end
xlabels = preg_xlabs;

Nan_preg = find_all_Nan(preg_data);
Nan_lact = find_all_Nan(lact_data);
allNan = intersect(Nan_preg, Nan_lact);

[rm_xlabs, rm_pregdat] = removeAllNan(allNan, xlabels, preg_data);
[temp, rm_lactdat] = removeAllNan(allNan, xlabels, lact_data);

if size(rm_xlabs) ~= size(temp)
    error('rm_xlabs do not match')
end

%% make figures
temp = [rm_pregdat, rm_lactdat];
clim = [min(temp,[],'all'), max(temp,[],'all')];
cmissing = 'w';
labmissing = '<1.0%';
fsize = 18;
ylabels = {'[PTH]_p', '[Ca^{2+}]_p', '[1,25(OH)_2D_3]_p'};
cmap = turbo;

figure(1)
clf
subplot(2,1,1)
h = heatmap(rm_xlabs, ylabels, rm_pregdat,...
        'colormap', cmap,...
        'MissingDataColor', cmissing, 'MissingDataLabel', labmissing,...
        'ColorLimits', clim, 'colorbarvisible', 'on');
ylabel('Pregnant to virgin')
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
h.FontSize = fsize;

subplot(2,1,2)
h = heatmap(rm_xlabs, ylabels, rm_lactdat,...
        'colormap', cmap,...
        'MissingDataColor', cmissing, 'MissingDataLabel', labmissing,...
        'ColorLimits', clim, 'colorbarvisible', 'off');
ylabel('Lactating to virgin')
xlabel('Parameters')
h.FontSize = fsize;

sgtitle('Impact of individual maternal adaptations on pregnancy and lactation models', 'fontsize', 24)






%---------------------------
% Functions used
%---------------------------
function [xlabels, round_data] = getdata(dat, pregORlact)
    if strcmp(pregORlact, 'preg')
        frac_sens = dat.preg_frac;
    elseif strcmp(pregORlact, 'lact')
        frac_sens = dat.lact_frac;
    end
    labels = cell(size(dat.param_names));
    for ii = 1:size(dat.param_names,2)
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

