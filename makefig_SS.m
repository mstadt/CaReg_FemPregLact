% This script creates the baseline steady state figures found in the
% manuscript draft (Fig 3.1 and 3.2)

%-----------------------
% User Input
%-----------------------
compare = 4; % 3 or 4


% load files where SS solutions are stored

% male file
fname1 = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-male_notes-newmale.mat';
SSdat1 = load(fname1);
vals1 = SSdat1.valsSS;
SS1 = SSdat1.SS;
lab1 = 'male';

% female file
fname2 = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-female_notes-newfemale.mat';
SSdat2 = load(fname2);
vals2 = SSdat2.valsSS;
SS2 = SSdat2.SS;
lab2 = 'female';

% pregnancy file
fname3 = './SSbest/16-May-2023_calcium_mod_SS_sexORrep-preg_notes-PTHupdate2.mat';
SSdat3 = load(fname3);
vals3 = SSdat3.valsSS;
SS3 = SSdat3.SS;
lab3 = 'pregnancy';

% lactation file
fname4 = './SSbest/16-May-2023_calcium_mod_SS_sexORrep-lact_notes-PTHupdate2.mat';
SSdat4 = load(fname4);
vals4 = SSdat4.valsSS;
SS4 = SSdat4.SS;
lab4 = 'lactation';

% figure specs
w = 1.0;
cmap = spring(4); %summer(14); %parula(14);
cvals = [cmap(4,:); cmap(3,:); cmap(2,:); cmap(1,:)];%[cmap(14,:); cmap(10,:); cmap(6,:); cmap(3,:)];
%[cmap(1,:); cmap(5,:); cmap(9,:); cmap(13,:)];
ce = 'black'; % error bar color
f_gca = 18;

%---------------------
%---------------------
if compare == 4
    leg_vals = {lab1, lab2, lab3, lab4};
    xvals = [0.5, 2, 3.5, 5.0];
    xnames = leg_vals;
    % get param vals
    for ii = 1:length(SSdat1.param_names)
        assign_name = strcat(SSdat1.param_names{ii}, '_1');
        assignin('base', assign_name, SSdat1.params(ii));
    end
    
    for ii = 1:length(SSdat2.param_names)
        assign_name = strcat(SSdat2.param_names{ii}, '_2');
        assignin('base', assign_name, SSdat2.params(ii));
    end

    for ii = 1:length(SSdat3.param_names)
        assign_name = strcat(SSdat3.param_names{ii}, '_3');
        assignin('base', assign_name, SSdat3.params(ii));
    end

    for ii = 1:length(SSdat4.param_names)
        assign_name = strcat(SSdat4.param_names{ii}, '_4');
        assignin('base', assign_name, SSdat4.params(ii));
    end


elseif compare == 3
    leg_vals = {lab1, lab2, lab3};
    xvals = [0.5, 2, 3.5];
    xnames = leg_vals;
    for ii = 1:length(SSdat1.param_names)
        assign_name = strcat(SSdat1.param_names{ii}, '_1');
        assignin('base', assign_name, SSdat1.params(ii));
    end
    
    for ii = 1:length(SSdat2.param_names)
        assign_name = strcat(SSdat2.param_names{ii}, '_2');
        assignin('base', assign_name, SSdat2.params(ii));
    end

    for ii = 1:length(SSdat3.param_names)
        assign_name = strcat(SSdat3.param_names{ii}, '_3');
        assignin('base', assign_name, SSdat3.params(ii));
    end
end

%% plasma concentrations
figure(1)
clf;
nrows = 1; ncols = 3;
subplot(nrows,ncols,1)
if compare == 3
    convals = [SS1(2)/Vp_1 SS2(2)/Vp_2 SS3(2)/Vp_3];
elseif compare == 4
    convals = [SS1(2)/Vp_1 SS2(2)/Vp_2 SS3(2)/Vp_3 SS4(2)/Vp_4];
end
hold on
for ii = 1:compare
    bar(xvals(ii), convals(ii),w, 'facecolor', cvals(ii,:))
end
ylabel('[PTH]_p (pmol/L)')
%title('Plasma PTH concentration')
xticks(xvals)
xticklabels(xnames)
set(gca, 'fontsize', f_gca)
grid on

subplot(nrows,ncols,2)
if compare == 3
    convals = [SS1(3)/Vp_1 SS2(3)/Vp_2 SS3(3)/Vp_3];
elseif compare == 4
    convals = [SS1(3)/Vp_1 SS2(3)/Vp_2 SS3(3)/Vp_3 SS4(3)/Vp_4];
end
hold on
for ii = 1:compare
    bar(xvals(ii), convals(ii),w, 'facecolor', cvals(ii,:))
end
ylabel('[Ca^{2+}]_p (mmol/L)')
%title('Plasma calcium concentration')
xticks(xvals)
xticklabels(xnames)
set(gca, 'fontsize', f_gca)
ylim([0.0,1.4])
grid on

subplot(nrows,ncols,3)
if compare == 3
    convals = [SS1(4)/Vp_1 SS2(4)/Vp_2 SS3(4)/Vp_3];
elseif compare == 4
    convals = [SS1(4)/Vp_1 SS2(4)/Vp_2 SS3(4)/Vp_3 SS4(4)/Vp_4];
end
hold on
for ii = 1:compare
    bar(xvals(ii), convals(ii),w, 'facecolor', cvals(ii,:))
end
ylabel('[1,25(OH)_2D_3]_p (pmol/L)')
%title('Plasma calcitriol concentration')
xticks(xvals)
xticklabels(xnames)
set(gca, 'fontsize', f_gca)
grid on

AddLetters2Plots(figure(1), {'(a)', '(b)', '(c)'}, 'FontSize', 18)

%% calcium fluxes
% bone together
% Fetus or Milk
figure(2)
mmol2mumol = 1e3;
clf
nrows=2; ncols=3;
%% calcium into plasma
subplot(nrows,ncols,1)
% gut absorption
temp = mmol2mumol*[vals1.Gut_absorption, vals2.Gut_absorption, vals3.Gut_absorption, vals4.Gut_absorption];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{abs} (\mumol/min)')
%title('Intestinal absorption')
xticks(xvals)
xticklabels(xnames)
set(gca,'fontsize',f_gca)
grid on

subplot(nrows,ncols,2)
% bone 2 plasma: resorption + fastpool2plasma
temp1 = mmol2mumol*[vals1.Bone_resorption, vals2.Bone_resorption, vals3.Bone_resorption, vals4.Bone_resorption];
temp2 = mmol2mumol*[vals1.FastPool_to_Plasma, vals2.FastPool_to_Plasma, vals3.FastPool_to_Plasma, vals4.FastPool_to_Plasma];
temp = temp1 + temp2;
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{res} + \Gamma_{f-p} (\mumol/min)')
%title('Bone to Plasma')
xticks(xvals)
xticklabels(xnames)
ylim([0.0 1.2])
yticks(0.0:0.2:1.2)
set(gca,'fontsize',f_gca)
grid on

% subplot(nrows,ncols,3)
% % fast2plasma
% temp = mmol2mumol*[vals1.FastPool_to_Plasma, vals2.FastPool_to_Plasma, vals3.FastPool_to_Plasma, vals4.FastPool_to_Plasma];
% hold on
% for ii = 1:compare
%     bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
% end
% ylabel('\Gamma_{f-p} (\mumol/min)')
% title('Fast Pool to Plasma')
% xticks(xvals)
% xticklabels(xnames)
% set(gca,'fontsize',f_gca)
% grid on


%% calcium OUT of plasma
subplot(nrows,ncols,4)
% urine
temp = mmol2mumol*[vals1.Urine_excretion, vals2.Urine_excretion, vals3.Urine_excretion, vals4.Urine_excretion];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_u (\mumol/min)')
%title('Urine excretion')
xticks(xvals)
xticklabels(xnames)
ylim([0.0,0.06])
yticks(0.0:0.01:0.06)
set(gca,'fontsize',f_gca)
grid on

subplot(nrows,ncols,5)
% plasma 2 fastpool + bone accretion
temp1 = mmol2mumol*[vals1.Plasma_to_FastPool, vals2.Plasma_to_FastPool, vals3.Plasma_to_FastPool, vals4.Plasma_to_FastPool];
%temp2 = mmol2mumol*[vals1.Bone_accretion, vals2.Bone_accretion, vals3.Bone_accretion, vals4.Bone_accretion];
temp = temp1; %+ temp2; % accretion doesnt interact with plasma!
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{p-f} (\mumol/min)')
%title('Plasma to Bone')
ylim([0.0, 2.1])
xticks(xvals)
xticklabels(xnames)
ylim([0.0,2.4])
yticks(0.0:0.4:2.4)
set(gca,'fontsize',f_gca)
grid on

% subplot(nrows,ncols,7)
% % plasma 2 fastpool
% temp = mmol2mumol*[vals1.Bone_accretion, vals2.Bone_accretion, vals3.Bone_accretion, vals4.Bone_accretion];
% hold on
% for ii = 1:compare
%     bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
% end
% ylabel('\Gamma_{ac} (\mumol/min)')
% title('Bone accretion')
% xticks(xvals)
% xticklabels(xnames)
% set(gca,'fontsize',f_gca)
% grid on

subplot(nrows,ncols,6)
% fetus or milk
temp = mmol2mumol*[FetusORMilk_1, FetusORMilk_2, FetusORMilk_3, FetusORMilk_4];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{x} (\mumol/min)')
%title('Fetus or milk calcium')
xticks(xvals)
xticklabels(xnames)
ylim([0.0,1.2])
yticks(0.0:0.2:1.2)
set(gca,'fontsize',f_gca)
grid on
AddLetters2Plots(figure(2), {'(a)', '(b)', '(c)', '(d)', '(e)'}, 'FontSize', 18)


% calcium fluxes
figure(3)
mmol2mumol = 1e3;
clf
nrows=2; ncols=4;
%% calcium into plasma
subplot(nrows,ncols,1)
% gut absorption
temp = mmol2mumol*[vals1.Gut_absorption, vals2.Gut_absorption, vals3.Gut_absorption, vals4.Gut_absorption];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{abs} (\mumol/min)')
title('Gut absorption')
xticks(xvals)
xticklabels(xnames)
set(gca,'fontsize',f_gca)
grid on

subplot(nrows,ncols,2)
% resorption
temp = mmol2mumol*[vals1.Bone_resorption, vals2.Bone_resorption, vals3.Bone_resorption, vals4.Bone_resorption];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{res} (\mumol/min)')
title('Bone resorption')
xticks(xvals)
xticklabels(xnames)
set(gca,'fontsize',f_gca)
grid on

subplot(nrows,ncols,3)
% fast2plasma
temp = mmol2mumol*[vals1.FastPool_to_Plasma, vals2.FastPool_to_Plasma, vals3.FastPool_to_Plasma, vals4.FastPool_to_Plasma];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{f-p} (\mumol/min)')
title('Fast Pool to Plasma')
xticks(xvals)
xticklabels(xnames)
set(gca,'fontsize',f_gca)
grid on


%% calcium OUT of plasma
subplot(nrows,ncols,5)
% urine
temp = mmol2mumol*[vals1.Urine_excretion, vals2.Urine_excretion, vals3.Urine_excretion, vals4.Urine_excretion];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_u (\mumol/min)')
title('Urine excretion')
xticks(xvals)
xticklabels(xnames)
set(gca,'fontsize',f_gca)
grid on

subplot(nrows,ncols,6)
% plasma 2 fastpool
temp = mmol2mumol*[vals1.Plasma_to_FastPool, vals2.Plasma_to_FastPool, vals3.Plasma_to_FastPool, vals4.Plasma_to_FastPool];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{p-f} (\mumol/min)')
title('Plasma to Fast Pool')
xticks(xvals)
xticklabels(xnames)
set(gca,'fontsize',f_gca)
grid on

subplot(nrows,ncols,7)
% plasma 2 fastpool
temp = mmol2mumol*[vals1.Bone_accretion, vals2.Bone_accretion, vals3.Bone_accretion, vals4.Bone_accretion];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{ac} (\mumol/min)')
title('Bone accretion')
xticks(xvals)
xticklabels(xnames)
set(gca,'fontsize',f_gca)
grid on

subplot(nrows,ncols,8)
% fetus or milk
temp = mmol2mumol*[FetusORMilk_1, FetusORMilk_2, FetusORMilk_3, FetusORMilk_4];
hold on
for ii = 1:compare
    bar(xvals(ii), temp(ii),w,'facecolor',cvals(ii,:))
end
ylabel('\Gamma_{FetusORMilk} (\mumol/min)')
title('Fetus or milk calcium')
xticks(xvals)
xticklabels(xnames)
set(gca,'fontsize',f_gca)
grid on
