% postprocess results from compute_preg2lact_lact2preg_sensitivity

% makes figure 3.6 from manuscript
clearvars

% Get Data
% filename where results are stored 
fname = './results_preglact_sensitivity/16-May-2023_preg2lactsens_notes-PTHupdate.mat';
dat = load(fname);

[p2l_xlabs, p2l_data] = getdata(dat, 'p2l');
[l2p_xlabs, l2p_data] = getdata(dat, 'l2p');

% sanity check: check that xlabels are the same
flag = 0;
if size(p2l_xlabs) ~= size(l2p_xlabs)
    flag = 1;
end
for ii = 1:size(p2l_xlabs, 2)
    if ~strcmp(p2l_xlabs, l2p_xlabs)
        flag = 1;
    end
end 
if flag
    error('xlabels do not match')
end
xlabels = p2l_xlabs;

Nan_p2l = find_all_Nan(p2l_data);
Nan_l2p = find_all_Nan(l2p_data);
allNan = intersect(Nan_p2l, Nan_l2p);

[rm_xlabs, rm_p2ldat] = removeAllNan(allNan, xlabels, p2l_data);
[~, rm_l2pdat]        = removeAllNan(allNan, xlabels, l2p_data);



%% make figures
temp = [rm_p2ldat, rm_l2pdat];
clim = [min(temp,[],'all'),max(temp,[],'all')];
cmissing = 'w';
labmissing = '<1.0%';
fsize = 18;
ylabels = {'[PTH]_p', '[Ca^{2+}]_p', '[1,25(OH)_2D_3]_p'};
cmap = turbo; %cool;

figure(15)
clf
subplot(2,1,1)
h = heatmap(rm_xlabs, ylabels, rm_p2ldat,...
        'colormap', cmap,...
        'MissingDataColor', cmissing, 'MissingDataLabel', labmissing,...
        'ColorLimits', clim, 'colorbarvisible', 'on');
ylabel('Pregnant to lactating')
%Ax = gca;
%Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
h.FontSize = fsize;


subplot(2,1,2)
h = heatmap(rm_xlabs, ylabels, rm_l2pdat,...
        'colormap', cmap,...
        'MissingDataColor', cmissing, 'MissingDataLabel', labmissing,...
        'ColorLimits', clim, 'colorbarvisible', 'off');
ylabel('Lactating to pregnant')
xlabel('Parameters')
h.FontSize = fsize;

%sgtitle('Differential impact of pregnancy and lactation adaptations', 'fontsize', 24)

%------------------
% functions used
%------------------
function [xlabels, round_data] = getdata(dat, p2lORl2p)
    if strcmp(p2lORl2p, 'p2l')
        frac_sens = dat.preg2lact_frac;
    elseif strcmp(p2lORl2p, 'l2p')
        frac_sens = dat.lact2preg_frac;
    end
    labels = cell(size(dat.diffIDs));
    for ii = 1:length(dat.diffIDs)
        IDval = dat.diffIDs(ii);
        labels{ii} = convert_param_name(dat.param_names{IDval});
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


% xlabels = cell(size(dat.diffIDs));
% for ii = 1:length(dat.diffIDs)
%     IDval = dat.diffIDs(ii);
%     xlabels{ii} = convert_param_name(dat.param_names{IDval});
% end
% ylabels = {'[PTH]_p', '[Ca^{2+}]_p', '[1,25(OH)_2D_3]_p'};
% fsize = 16; % fontsize
% 
% 
% 
% %%% preg2lact figure
% figure(1)
% clf
% round_preg2lact = round(dat.preg2lact_frac', 2, 'significant');
% [r,c] = find(abs(round_preg2lact) < 0.1);
% for ii = 1:length(r)
%     round_preg2lact(r(ii),c(ii)) = NaN;
% end
% h = heatmap(xlabels, ylabels, round_preg2lact,...
%             'colormap', spring,...
%             'MissingDataColor', 'w', 'MissingDataLabel', '<0.1%');
% sortx(h, '[Ca^{2+}]_p', 'ascend',...
%             'MissingPlacement', 'last'); % sort heatmap by calcium
% h.FontSize = fsize;
% title('Impact of individual lactation parameters on pregnancy model plasma concentration')
% 
% 
% %%% lact2preg figure
% figure(2)
% clf
% round_lact2preg = round(dat.lact2preg_frac', 2, 'significant');
% [r,c] = find(abs(round_lact2preg) < 0.1);
% for ii = 1:length(r)
%     round_lact2preg(r(ii),c(ii)) = NaN;
% end
% h = heatmap(xlabels, ylabels, round_lact2preg,...
%             'colormap', spring,...
%             'MissingDataColor', 'w', 'MissingDataLabel', '<0.1%');
% sortx(h, '[Ca^{2+}]_p', 'ascend',...
%             'MissingPlacement', 'last'); % sort heatmap by calcium
% h.FontSize = fsize;
% title('Impact of individual pregnancy parameters on lactation model plasma concentration')
