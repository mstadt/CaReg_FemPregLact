% make figs best on results from driver_hyperparathyroid
clear all;

PTHchange_vals = 1:1:100;

notes_male = 'may16'; %'gammaProdD3';
notes_female = 'may16';
notes_preg = 'PTH2';
notes_lact = 'PTHchange';
in_date_male = '16-May-2023';
in_date_female = '16-May-2023';
in_date = '16-May-2023';

[male_SS, male_vars, male_SS_norm, male_vars_norm]         = get_vals(PTHchange_vals, notes_male, 'male', in_date_male);
[female_SS, female_vars, female_SS_norm, female_vars_norm] = get_vals(PTHchange_vals, notes_female, 'female', in_date_female);
[preg_SS, preg_vars, preg_SS_norm, preg_vars_norm]         = get_vals(PTHchange_vals, notes_preg, 'preg', in_date);
[lact_SS, lact_vars, lact_SS_norm, lact_vars_norm]         = get_vals(PTHchange_vals, notes_lact, 'lact', in_date);

% make figures
lw = 3.5;
ls1 = '-'; ls2 = '--'; ls3 = ':'; ls4 = '-.';
fgca = 18;
fleg = 14;
xlab = 'Normalized PTH synthesis';
graymap = gray(5);
darkgray = graymap(2,:);

% concentrations
% male and female
figure(1)
clf
hold on
inds = [2, 3, 4];
cmap = parula(4);
for ii = 1:length(inds)
    ind = inds(ii);
    plot(PTHchange_vals, male_SS_norm(:,ind), 'linestyle', ls1, 'color', cmap(ii,:), 'linewidth', lw)
    plot(PTHchange_vals, female_SS_norm(:,ind),'linestyle',ls2, 'color',cmap(ii,:), 'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized concentration')
legend('[PTH]_p (male)', '[PTH]_p (female)',...
            '[Ca^{2+}]_p (male)', '[Ca^{2+}]_p (female)',...
            '[1,25(OH)_2D_3]_p (male)', '[1,25(OH)_2D_3]_p (female)',...
                'fontsize', fleg,...
                'location', 'northwest')
title('Impact of hyperparathyroidism on plasma concentrations in male and female')
set(gca, 'fontsize', fgca)
xlim([min(PTHchange_vals), max(PTHchange_vals)])
ylim([1.0, 7])
%yticks(1.0:2:7.0)
grid on

% concentrations
% virgin, preg, lact
figure(2)
clf
hold on
inds = [2, 3, 4];
cmap = parula(4);
for ii = 1:length(inds)
    ind = inds(ii);
    plot(PTHchange_vals, female_SS_norm(:,ind), 'linestyle', ls1, 'color', cmap(ii,:), 'linewidth', lw)
    plot(PTHchange_vals, preg_SS_norm(:,ind),'linestyle',ls2, 'color',cmap(ii,:), 'linewidth',lw)
    plot(PTHchange_vals, lact_SS_norm(:,ind),'linestyle',ls3, 'color',cmap(ii,:), 'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized concentration')
legend('[PTH]_p (virgin)', '[PTH]_p (pregnant)', '[PTH]_p (lactating)',...
            '[Ca^{2+}]_p (virgin)', '[Ca^{2+}]_p (pregnant)', '[Ca^{2+}]_p (lactating)',...
            '[1,25(OH)_2D_3]_p (virgin)', '[1,25(OH)_2D_3]_p (pregnant)', '[1,25(OH)_2D_3]_p (lactating)',...
                'fontsize', fleg,...
                'location', 'northwest')
title('Impact of hyperparathyroidism on plasma concentrations in virgin, pregnant, and lactating rats')
set(gca, 'fontsize', fgca)
xlim([min(PTHchange_vals), max(PTHchange_vals)])
ylim([1.0,7])
grid on


figure(4)
% inward fluxes
clf
subplot(1,2,1)
hold on
inds = [4, 9, 14];
cmap3 = parula(8);
m=2;
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
for jj = 1:length(inds)
    ind = inds(jj);
    plot(PTHchange_vals, female_vars_norm(:,ind), 'linestyle', ls1,'color',cmap3(m*jj,:), 'linewidth',lw)
    plot(PTHchange_vals, preg_vars_norm(:,ind), 'linestyle',ls2, 'color',cmap3(m*jj,:), 'linewidth',lw)
    plot(PTHchange_vals, lact_vars_norm(:,ind), 'linestyle', ls3, 'color', cmap3(m*jj,:), 'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized fluxes')
legend('','Gut absorption (virgin)', 'Gut absorption (pregnancy)', 'Gut absorption (lactation)',...
    'Bone resorption (virgin)', 'Bone resorption (pregnancy)', 'Bone resorption (lactation)',...
    'FastPool to Plasma (virgin)', 'FastPool to Plasma (pregnancy)', 'FastPool to Plasma (lactation)',...
    'fontsize', fleg,...
    'location', 'southwest')
title('Inward calcium fluxes')
set(gca, 'fontsize', fgca)
grid on

%figure(5)
% outward fluxes
% virgin, preg, lact
%clf

subplot(1,2,2)
hold on
inds = [7, 8, 13];
cmap4 =autumn(8);
m=2;
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
for ii = 1:length(inds)
    ind = inds(ii);
    plot(PTHchange_vals, female_vars_norm(:,ind), 'linestyle', ls1,'color',cmap4(m*ii,:), 'linewidth',lw)
    plot(PTHchange_vals, preg_vars_norm(:,ind), 'linestyle',ls2, 'color',cmap4(m*ii,:), 'linewidth',lw)
    plot(PTHchange_vals, lact_vars_norm(:,ind), 'linestyle', ls3, 'color', cmap4(m*ii,:), 'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized fluxes')
legend('','Urine excretion (virgin)', 'Urine excretion (pregnancy)', 'Urine excretion (lactation)',...
    'Bone accretion (virgin)', 'Bone accretion (pregnancy)', 'Bone accretion (lactation)',...
    'Plasma to FastPool (virgin)', 'Plasma to FastPool (pregnancy)', 'Plasma to FastPool (lactation)',...
    'fontsize', fleg,...
    'location', 'southwest')
title('Outward calcium fluxes')
set(gca, 'fontsize', fgca)
grid on

sgtitle('Impact of hyperparathyroid on inward calcium fluxes during pregnancy and lactation', 'fontsize', 20)

% combine bone fluxes together (resorption + Fast2Plasma, Accretion + Plasma2FastPool) 
figure(21)
clf
subplot(1,2,1)
hold on
% bone_in = accretion (8) + plasma2fastpool (13); % REMOVE ACCRETION BC NOT
% INTERACTING WITH PLASMA
bone_in_female = female_vars(:,13); %female_vars(:,8) + female_vars(:,13);
bone_in_female_norm = bone_in_female./female_vars(1,13); %bone_in_female./(female_vars(1,8) + female_vars(1,13));

bone_in_preg =  preg_vars(:,13);
bone_in_preg_norm = bone_in_preg./(preg_vars(1,13));

bone_in_lact = lact_vars(:,13);
bone_in_lact_norm = bone_in_lact./(lact_vars(1,13));

% bone_out = resorption(9) + fast2plasma (14)
bone_out_female = female_vars(:,9) + female_vars(:,14);
bone_out_female_norm = bone_out_female./(female_vars(1,9) + female_vars(1,14));

bone_out_preg = preg_vars(:,9) + preg_vars(:,14);
bone_out_preg_norm = bone_out_preg./(preg_vars(1,9) + preg_vars(1,14));

bone_out_lact = lact_vars(:,9) + lact_vars(:,14);
bone_out_lact_norm = bone_out_lact./(lact_vars(1,9) + lact_vars(1,14));

% inward fluxes
hold on
cmap3 = parula(6);
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
% gut absorption (4)
plot(PTHchange_vals, female_vars_norm(:,4), 'linestyle', ls1,'color',cmap3(1,:), 'linewidth',lw)
plot(PTHchange_vals, preg_vars_norm(:,4), 'linestyle',ls2, 'color',cmap3(1,:), 'linewidth',lw)
plot(PTHchange_vals, lact_vars_norm(:,4), 'linestyle', ls3, 'color', cmap3(1,:), 'linewidth',lw)

% bone out (bone2plasma -- resorption(9) + fast2plasma (14))
plot(PTHchange_vals, bone_out_female_norm, 'linestyle', ls1, 'color', cmap3(4,:), 'linewidth', lw)
plot(PTHchange_vals, bone_out_preg_norm, 'linestyle', ls2, 'color', cmap3(4,:), 'linewidth', lw)
plot(PTHchange_vals, bone_out_lact_norm, 'linestyle', ls3, 'color', cmap3(4,:), 'linewidth', lw)

ylabel('Normalized fluxes')
legend('','Gut absorption (virgin)', 'Gut absorption (pregnancy)', 'Gut absorption (lactation)',...
            'Bone to plasma (virgin)', 'Bone to plasma (pregnancy)', 'Bone to plasma (lactation)',...
            'fontsize', fleg,...
            'location', 'northwest')
title('Inward calcium fluxes')
set(gca, 'fontsize', fgca)
ylim([1.0,1.6])
grid on

subplot(1,2,2)
hold on
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
% outward fluxes
cmap4 = autumn(6);
% urine (7)
plot(PTHchange_vals, female_vars_norm(:,7), 'linestyle', ls1,'color',cmap4(1,:), 'linewidth',lw)
plot(PTHchange_vals, preg_vars_norm(:,7), 'linestyle',ls2, 'color',cmap4(1,:), 'linewidth',lw)
plot(PTHchange_vals, lact_vars_norm(:,7), 'linestyle', ls3, 'color', cmap4(1,:), 'linewidth',lw)

% bone in (plasma to bone -- accretion (8) + plasma2fastpool (13))
plot(PTHchange_vals, bone_in_female_norm, 'linestyle', ls1, 'color', cmap4(5,:), 'linewidth', lw)
plot(PTHchange_vals, bone_in_preg_norm, 'linestyle', ls2, 'color', cmap4(5,:), 'linewidth', lw)
plot(PTHchange_vals, bone_in_lact_norm, 'linestyle', ls3, 'color', cmap4(5,:), 'linewidth', lw)
xlabel(xlab)
ylabel('Normalized fluxes')
grid on

legend('','Urine Excretion (virgin)', 'Urine Excretion (pregnancy)', 'Urine Excretion (lactation)',...
            'Plasma to bone (virgin)', 'Plasma to bone (pregnancy)', 'Plasma to bone (lactation)',...
            'fontsize', fleg,...
            'location', 'northwest')
title('Outward calcium fluxes')
set(gca, 'fontsize', fgca)
ylim([1.0,2.6])
sgtitle('Impact of hyperparathyroid on calcium fluxes during pregnancy and lactation', 'fontsize', 20)

%-----------------------------
% functions used
%-----------------------------

function [SS_vals, SS_vars, SS_vals_norm, SS_vars_norm] = get_vals(PTHchange_vals, notes, sexORrep, in_date)
    SS_vals = zeros(length(PTHchange_vals), 6);
    SS_vars = zeros(length(PTHchange_vals), 14);
    for ii = 1:length(PTHchange_vals)
        PTHchange = PTHchange_vals(ii);
        fname = strcat('./results_hyperPTH/', in_date,'_hyperparathyroid_', 'sexORrep-',sexORrep,...
                                '_PTHchange-', num2str(PTHchange), '_notes-', notes, '.mat');
        dat = load(fname);
        SS_vals(ii,:) = dat.SS;

        vals = dat.valsSS;
        SS_vars(ii,:) = [vals.PTHg_synthesis;              % 1
                                vals.PTHg_exocytosis;       % 2
                                vals.Gut_frac_absorption;   % 3
                                vals.Gut_absorption;        % 4
                                vals.Renal_filtration;      % 5
                                vals.Renal_frac_reab;       % 6 
                                vals.Urine_excretion;       % 7
                                vals.Bone_accretion;        % 8
                                vals.Bone_resorption;       % 9
                                vals.Lambda_PT;             % 10
                                vals.Lambda_TAL;            % 11
                                vals.Lambda_DCT;            % 12
                                vals.Plasma_to_FastPool;    % 13
                                vals.FastPool_to_Plasma]';  % 14  
        fclose('all');
    end

    % make normalized data
    SSbase = SS_vals(1,:);
    SSvars_base = SS_vars(1,:);

    SS_vals_norm = SS_vals./SSbase;
    SS_vars_norm = SS_vars./SSvars_base;
end
