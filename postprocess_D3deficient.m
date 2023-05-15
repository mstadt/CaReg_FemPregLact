% make figs based on results from driver_D3deficiency
clear all;

% user input
D3change_vals = (100:-1:1)./100.0;

fig_D3change = (ones(size(D3change_vals)) - D3change_vals) * 100;
notes_male = 'finalmale';
notes_female = 'finalfemale';
notes_preg = 'finalpreg'; % TO DO: also without Gamma_ac
notes_lact = 'finallact';
in_date_male = '12-May-2023';
in_date_female = '12-May-2023';
in_date_preg = '12-May-2023';
in_date_lact = '12-May-2023';

[male_SS, male_vars, male_SS_norm, male_vars_norm]         = getvals(D3change_vals, notes_male, 'male', in_date_male);
[female_SS, female_vars, female_SS_norm, female_vars_norm] = getvals(D3change_vals, notes_female, 'female', in_date_female);
[preg_SS, preg_vars, preg_SS_norm, preg_vars_norm]         = getvals(D3change_vals, notes_preg, 'preg', in_date_preg);
[lact_SS, lact_vars, lact_SS_norm, lact_vars_norm]         = getvals(D3change_vals, notes_lact, 'lact', in_date_lact);

%%
% make figures
lw = 4.0;
ls1 = '-';
ls2 = '--';
ls3 = ':';
ls4 = '-.';
cmap = parula(12);
graymap = gray(5);
darkgray = graymap(2,:);

f_gca = 18;
fleg = 16;
xlab = 'Inhibition of inactive vitamin D_3 concentration (%)';

% figure(1)
% clf
% hold on
% for ii = 1:6
%     plot(fig_D3change, male_SS_norm(:,ii), 'linestyle', ls1, 'color', cmap(2*ii,:),'linewidth',lw)
%     plot(fig_D3change, female_SS_norm(:,ii), 'linestyle', ls2,'color',cmap(2*ii,:), 'linewidth',lw)
%     plot(fig_D3change, preg_SS_norm(:,ii), 'linestyle',ls3, 'color',cmap(2*ii,:), 'linewidth',lw)
%     plot(fig_D3change, lact_SS_norm(:,ii), 'linestyle', ls4, 'color', cmap(2*ii,:), 'linewidth',lw)
% end
% xlabel(xlab)
% ylabel('Normalized concentrations')
% legend('PTH_g', '', '', '',...
%         'PTH_p', '', '', '',...
%         'Ca_p','','', '',...
%         'D3_p', '', '', '',...
%         'NCa_f', '', '', '',...
%         'NCa_s', '', '', '',...
%         'fontsize',fleg)
% set(gca, 'fontsize', f_gca)
% grid on


% virgin, lactation, pregnancy
figure(2)
clf
hold on
inds = [2, 3, 4];
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
for ii = 1:length(inds)
    %plot(fig_D3change, male_SS_norm(:,ii), 'linestyle', ls1, 'color', cmap(3*ii,:),'linewidth',lw)
    ind = inds(ii);
    plot(fig_D3change, female_SS_norm(:,ind), 'linestyle', ls1,'color',cmap(4*ii-2,:), 'linewidth',lw)
    plot(fig_D3change, preg_SS_norm(:,ind), 'linestyle',ls2, 'color',cmap(4*ii-2,:), 'linewidth',lw)
    plot(fig_D3change, lact_SS_norm(:,ind), 'linestyle', ls3, 'color', cmap(4*ii-2,:), 'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized concentration')
legend('','[PTH]_p (virgin)', '[PTH]_p (pregnancy)', '[PTH]_p (lactation)',...
        '[Ca^{2+}]_p (virgin)','[Ca^{2+}]_p (pregnancy)','[Ca^{2+}]_p (lactation)',...
        '[1,25(OH)_2D_3]_p (virgin)', '[1,25(OH)_2D_3]_p (pregnancy)', '[1,25(OH)_2D_3]_p (lactation)',...
        'fontsize',fleg,...
        'location','northwest')
title('Impact of vitamin D_3 deficiency on plasma concentrations during pregnancy and lactation')
set(gca, 'fontsize', f_gca)
grid on

% virgin, preg, lact fluxes
figure(3)
clf
hold on
inds = [4, 7, 8, 9, 13, 14];
cmap2 = parula(25);
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
for ii = 1:length(inds)
    ind = inds(ii);
    plot(fig_D3change, female_vars_norm(:,ind), 'linestyle', ls1,'color',cmap2(4*ii,:), 'linewidth',lw)
    plot(fig_D3change, preg_vars_norm(:,ind), 'linestyle',ls2, 'color',cmap2(4*ii,:), 'linewidth',lw)
    plot(fig_D3change, lact_vars_norm(:,ind), 'linestyle', ls3, 'color', cmap2(4*ii,:), 'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized fluxes')
legend('','Gut absorption (virgin)', 'Gut absorption (pregnancy)', 'Gut absorption (lactation)',...
    'Urine excretion (virgin)', 'Urine excretion (pregnancy)', 'Urine excretion (lactation)',...
    'Bone accretion (virgin)', 'Bone accretion (pregnancy)', 'Bone accretion (lactation)',...
    'Bone resorption (virgin)', 'Bone resorption (pregnancy)', 'Bone resorption (lactation)',...
    'Plasma to FastPool (virgin)', 'Plasma to FastPool (pregnancy)', 'Plasma to FastPool (lactation)',...
    'FastPool to Plasma (virgin)', 'FastPool to Plasma (pregnancy)', 'FastPool to Plasma (lactation)',...
    'fontsize', fleg,...
    'location', 'southwest')
title('Impact of vitamin D_3 deficiency on calcium fluxes during pregnancy and lactation')
set(gca, 'fontsize', f_gca)
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
    plot(fig_D3change, female_vars_norm(:,ind), 'linestyle', ls1,'color',cmap3(m*jj,:), 'linewidth',lw)
    plot(fig_D3change, preg_vars_norm(:,ind), 'linestyle',ls2, 'color',cmap3(m*jj,:), 'linewidth',lw)
    plot(fig_D3change, lact_vars_norm(:,ind), 'linestyle', ls3, 'color', cmap3(m*jj,:), 'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized fluxes')
legend('','Gut absorption (virgin)', 'Gut absorption (pregnancy)', 'Gut absorption (lactation)',...
    'Bone resorption (virgin)', 'Bone resorption (pregnancy)', 'Bone resorption (lactation)',...
    'FastPool to Plasma (virgin)', 'FastPool to Plasma (pregnancy)', 'FastPool to Plasma (lactation)',...
    'fontsize', fleg,...
    'location', 'southwest')
title('Inward calcium fluxes')
set(gca, 'fontsize', f_gca)
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
    plot(fig_D3change, female_vars_norm(:,ind), 'linestyle', ls1,'color',cmap4(m*ii,:), 'linewidth',lw)
    plot(fig_D3change, preg_vars_norm(:,ind), 'linestyle',ls2, 'color',cmap4(m*ii,:), 'linewidth',lw)
    plot(fig_D3change, lact_vars_norm(:,ind), 'linestyle', ls3, 'color', cmap4(m*ii,:), 'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized fluxes')
legend('','Urine excretion (virgin)', 'Urine excretion (pregnancy)', 'Urine excretion (lactation)',...
    'Bone accretion (virgin)', 'Bone accretion (pregnancy)', 'Bone accretion (lactation)',...
    'Plasma to FastPool (virgin)', 'Plasma to FastPool (pregnancy)', 'Plasma to FastPool (lactation)',...
    'fontsize', fleg,...
    'location', 'southwest')
title('Outward calcium fluxes')
set(gca, 'fontsize', f_gca)
grid on

sgtitle('Impact of vitamin D_3 deficiency on inward calcium fluxes during pregnancy and lactation', 'fontsize', 20)

% sex-specific
figure(10)
clf
hold on
inds = [2, 3, 4];
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
for ii = 1:length(inds)
    ind = inds(ii);
    plot(fig_D3change, male_SS_norm(:,ind), 'linestyle', ls4, 'color', cmap(4*ii-2,:),'linewidth',lw)
    plot(fig_D3change, female_SS_norm(:,ind), 'linestyle', ls1,'color',cmap(4*ii-2,:), 'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized concentration')
legend('','[PTH]_p (male)', '[PTH]_p (female)', ...
        '[Ca^{2+}]_p (male)','[Ca^{2+}]_p (female)',...
        '[1,25(OH)_2D_3]_p (male)', '[1,25(OH)_2D_3]_p (female)',...
        'fontsize',fleg,...
        'location','northwest')
title('Impact of vitamin D_3 deficiency on plasma concentrations in males versus females')
set(gca, 'fontsize', f_gca)
grid on

%%
% combine bone fluxes together (resorption + Fast2Plasma, Accretion + Plasma2FastPool) 
figure(20)
clf
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
subplot(1,2,1)
hold on
cmap3 = parula(6);
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
% gut absorption (4)
plot(fig_D3change, female_vars_norm(:,4), 'linestyle', ls1,'color',cmap3(1,:), 'linewidth',lw)
plot(fig_D3change, preg_vars_norm(:,4), 'linestyle',ls2, 'color',cmap3(1,:), 'linewidth',lw)
plot(fig_D3change, lact_vars_norm(:,4), 'linestyle', ls3, 'color', cmap3(1,:), 'linewidth',lw)

% bone out (bone2plasma -- resorption(9) + fast2plasma (14))
plot(fig_D3change, bone_out_female_norm, 'linestyle', ls1, 'color', cmap3(4,:), 'linewidth', lw)
plot(fig_D3change, bone_out_preg_norm, 'linestyle', ls2, 'color', cmap3(4,:), 'linewidth', lw)
plot(fig_D3change, bone_out_lact_norm, 'linestyle', ls3, 'color', cmap3(4,:), 'linewidth', lw)

xlabel(xlab)
ylabel('Normalized fluxes')
grid on
legend('','Gut absorption (virgin)', 'Gut absorption (pregnancy)', 'Gut absorption (lactation)',...
            'Bone to plasma (virgin)', 'Bone to plasma (pregnancy)', 'Bone to plasma (lactation)',...
            'fontsize', fleg,...
            'location', 'southwest')
title('Inward calcium fluxes')
ylim([0.3,1.1])
set(gca, 'fontsize', f_gca)

% outward fluxes
subplot(1,2,2)
hold on
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
cmap4 = autumn(6);
% urine (7)
plot(fig_D3change, female_vars_norm(:,7), 'linestyle', ls1,'color',cmap4(1,:), 'linewidth',lw)
plot(fig_D3change, preg_vars_norm(:,7), 'linestyle',ls2, 'color',cmap4(1,:), 'linewidth',lw)
plot(fig_D3change, lact_vars_norm(:,7), 'linestyle', ls3, 'color', cmap4(1,:), 'linewidth',lw)

% bone in (plasma to bone -- accretion (8) + plasma2fastpool (13))
plot(fig_D3change, bone_in_female_norm, 'linestyle', ls1, 'color', cmap4(5,:), 'linewidth', lw)
plot(fig_D3change, bone_in_preg_norm, 'linestyle', ls2, 'color', cmap4(5,:), 'linewidth', lw)
plot(fig_D3change, bone_in_lact_norm, 'linestyle', ls3, 'color', cmap4(5,:), 'linewidth', lw)
xlabel(xlab)
ylabel('Normalized fluxes')
ylim([0.3,1.1])
grid on
legend('','Urine Excretion (virgin)', 'Urine Excretion (pregnancy)', 'Urine Excretion (lactation)',...
            'Plasma to bone (virgin)', 'Plasma to bone (pregnancy)', 'Plasma to bone (lactation)',...
            'fontsize', fleg,...
            'location', 'southwest')
title('Outward calcium fluxes')
set(gca, 'fontsize', f_gca)

sgtitle('Impact of vitamin D_3 deficiency on calcium fluxes during pregnancy and lactation', 'fontsize', 20)


% combine bone fluxes together (resorption + Fast2Plasma, Accretion + Plasma2FastPool) 
figure(21)
clf
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
plot(fig_D3change, female_vars_norm(:,4), 'linestyle', ls1,'color',cmap3(1,:), 'linewidth',lw)
plot(fig_D3change, preg_vars_norm(:,4), 'linestyle',ls2, 'color',cmap3(1,:), 'linewidth',lw)
plot(fig_D3change, lact_vars_norm(:,4), 'linestyle', ls3, 'color', cmap3(1,:), 'linewidth',lw)

% bone out (bone2plasma -- resorption(9) + fast2plasma (14))
plot(fig_D3change, bone_out_female_norm, 'linestyle', ls1, 'color', cmap3(4,:), 'linewidth', lw)
plot(fig_D3change, bone_out_preg_norm, 'linestyle', ls2, 'color', cmap3(4,:), 'linewidth', lw)
plot(fig_D3change, bone_out_lact_norm, 'linestyle', ls3, 'color', cmap3(4,:), 'linewidth', lw)


% outward fluxes
cmap4 = autumn(6);
% urine (7)
plot(fig_D3change, female_vars_norm(:,7), 'linestyle', ls1,'color',cmap4(1,:), 'linewidth',lw)
plot(fig_D3change, preg_vars_norm(:,7), 'linestyle',ls2, 'color',cmap4(1,:), 'linewidth',lw)
plot(fig_D3change, lact_vars_norm(:,7), 'linestyle', ls3, 'color', cmap4(1,:), 'linewidth',lw)

% bone in (plasma to bone -- accretion (8) + plasma2fastpool (13))
plot(fig_D3change, bone_in_female_norm, 'linestyle', ls1, 'color', cmap4(5,:), 'linewidth', lw)
plot(fig_D3change, bone_in_preg_norm, 'linestyle', ls2, 'color', cmap4(5,:), 'linewidth', lw)
plot(fig_D3change, bone_in_lact_norm, 'linestyle', ls3, 'color', cmap4(5,:), 'linewidth', lw)
xlabel(xlab)
ylabel('Normalized fluxes')
grid on

xlabel(xlab)
ylabel('Normalized fluxes')
grid on
legend('','Gut absorption (virgin)', 'Gut absorption (pregnancy)', 'Gut absorption (lactation)',...
            'Bone to plasma (virgin)', 'Bone to plasma (pregnancy)', 'Bone to plasma (lactation)',...
            'Urine Excretion (virgin)', 'Urine Excretion (pregnancy)', 'Urine Excretion (lactation)',...
            'Plasma to bone (virgin)', 'Plasma to bone (pregnancy)', 'Plasma to bone (lactation)',...
            'fontsize', fleg,...
            'location', 'southwest')
set(gca, 'fontsize', f_gca)

title('Impact of vitamin D_3 deficiency on calcium fluxes during pregnancy and lactation', 'fontsize', 20)

%%
%------------------------------------------------
% functions used
%------------------------------------------------

function [SSvals, SSvars_vals, SSvals_norm, SSvars_norm] = getvals(D3change_vals, notes, sexORrep, in_date)
% put data into matrices
SSvals = zeros(length(D3change_vals), 6);
SSvars_vals = zeros(length(D3change_vals),14);

for ii = 1:length(D3change_vals)
    D3change = D3change_vals(ii);
    fname = strcat('./results_D3deficient/', in_date,'_D3deficient_', 'sexORrep-',sexORrep,...
                                '_D3change-', num2str(D3change), '_notes-', notes, '.mat');
    dat = load(fname);
    SSvals(ii,:) = dat.SS;
    vals = dat.valsSS;
    SSvars_vals(ii,:) =  [vals.PTHg_synthesis;              % 1
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
end %ii

% make normalized data
SSbase = SSvals(1,:);
SSvarsbase = SSvars_vals(1,:);

SSvals_norm = SSvals./SSbase;
SSvars_norm = SSvars_vals./SSvarsbase;
end % getvarsvals