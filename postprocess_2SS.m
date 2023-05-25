clear all;

% which plots
con_vals = 1;
calcium_vals = 1;
kidney_vals = 1;
kidney_impact_vals = 0;
D3_vals = 1;
gut_vals = 1;
PTH_vals = 1;
bone_vals = 1;


% load SS
SS1dat = load('./results_malepreglact/16-May-2023_calcium_mod_SS_sexORrep-malelact_notes-lactmale.mat');
vals1 = SS1dat.valsSS;
SS1 = SS1dat.SS;
lab1 = 'SS1';
SS1dat.sexORrep = 'lact'

SS2dat = load('./SSbest/16-May-2023_calcium_mod_SS_sexORrep-lact_notes-PTHupdate2.mat');
vals2 = SS2dat.valsSS;
SS2 = SS2dat.SS;
lab2 = 'SS2';

% get param vals
for ii = 1:length(SS1dat.param_names)
    assign_name = strcat(SS1dat.param_names{ii}, '_1');
    assignin('base', assign_name, SS1dat.params(ii));
end

for ii = 1:length(SS2dat.param_names)
    assign_name = strcat(SS2dat.param_names{ii}, '_2');
    assignin('base', assign_name, SS2dat.params(ii));
end

% figure specs
w1 = 0.75; w2 = 0.5;
cvals = summer(4);
c1 = cvals(4,:);
c2 = cvals(1,:);
ce = 'black'; % error bar color


if con_vals
    % optimal concentrations for PTHp, Cap, D3p
    % male range
    opt_male = [1.5 7.25 13; % PTHp
                1.1 1.2 1.3; % Cap
                80 165 250]; % D3p

    % female range
    opt_female =[opt_male(1,:)*1.0; % PTHp
                opt_male(2,:)*1.0;  % Cap
                opt_male(3,:)*0.6]; % D3p

    % preg range
    opt_preg  = [opt_female(1,:)*1.8; % PTHp 
                opt_female(2,:)*1.0;  % Cap
                opt_female(3,:)*3.0]; % D3p

    % lact range
    opt_lact  = [opt_female(1,:)*1.8; % PTHp
            opt_female(2,:)*1.0;  % Cap
            opt_female(3,:)*3.5]; % D3p
    
    if strcmp(SS1dat.sexORrep, 'male')
        opt1 = opt_male;
    elseif strcmp(SS1dat.sexORrep, 'female')
        opt1 = opt_female;
    elseif strcmp(SS1dat.sexORrep, 'preg')
        opt1 = opt_preg; 
    elseif strcmp(SS1dat.sexORrep, 'lact')
        opt1 = opt_lact;
    end
    if strcmp(SS2dat.sexORrep, 'male')
        opt2 = opt_male;
    elseif strcmp(SS2dat.sexORrep, 'female')
        opt2 = opt_female;
    elseif strcmp(SS2dat.sexORrep, 'preg')
        opt2 = opt_preg; 
    elseif strcmp(SS2dat.sexORrep, 'lact')
        opt2 = opt_lact;
    end


    figure(1)
    clf
    x = [1, 2];
    nrows = 2;ncols=3;
    subplot(nrows,ncols,1)
    bar(x, [opt1(1,2) opt2(1,2)], w1, 'facecolor',c1)
    hold on
    bar(x, [SS1(2)/Vp_1 SS2(2)/Vp_2],w2,'facecolor',c2)
    errorbar(1, opt1(1,2),opt1(1,2)-opt1(1,1),opt1(1,3)-opt1(1,2),'color',ce)
    errorbar(2, opt2(1,2),opt2(1,2)-opt2(1,1),opt2(1,3)-opt2(1,2),'color',ce)
    xlabel('[PTH]_p')
    ylabel('pmol/L')
    title('PTH plasma concentration')
    xticklabels({lab1, lab2})
    grid on

    subplot(nrows,ncols,2)
    bar(x, [opt1(2,2) opt2(2,2)], w1, 'facecolor',c1)
    hold on
    bar(x, [SS1(3)/Vp_1 SS2(3)/Vp_2],w2,'facecolor',c2)
    errorbar(1, opt1(2,2),opt1(2,2)-opt1(2,1),opt1(2,3)-opt1(2,2),'color',ce)
    errorbar(2, opt2(2,2),opt2(2,2)-opt2(2,1),opt2(2,3)-opt2(2,2),'color',ce)
    xlabel('[Ca^{2+}]_p')
    ylabel('mmol/L')
    title('Ca^{2+} plasma concentration')
    xticklabels({lab1, lab2})
    grid on

    subplot(nrows,ncols,3)
    bar(x, [opt1(3,2) opt2(3,2)], w1, 'facecolor',c1)
    hold on
    bar(x, [SS1(4)/Vp_1 SS2(4)/Vp_2],w2,'facecolor',c2)
    errorbar(1, opt1(3,2),opt1(3,2)-opt1(3,1),opt1(3,3)-opt1(3,2),'color',ce)
    errorbar(2, opt2(3,2),opt2(3,2)-opt2(3,1),opt2(3,3)-opt2(3,2),'color',ce)
    xlabel('[1,25(OH)_2D_3]_p')
    ylabel('pmol/L')
    title('1,25(OH)_2D_3 plasma concentration')
    xticklabels({lab1, lab2})
    grid on

    subplot(nrows,ncols,4)
    bar(x, [SS1(1) SS2(1)],w2,'facecolor',c2)
    xlabel('PTH_g')
    ylabel('pmol')
    title('PTH gland pool')
    xticklabels({lab1, lab2})
    grid on

    subplot(nrows,ncols,5)
    bar(x, [SS1(5) SS2(5)],w2,'facecolor',c2)
    xlabel('NCa_f')
    ylabel('mmol')
    title('Fast bone pool')
    xticklabels({lab1, lab2})
    grid on

    subplot(nrows,ncols,6)
    bar(x, [SS1(6) SS2(6)],w2,'facecolor',c2)
    xlabel('NCa_s')
    ylabel('mmol')
    title('Slow bone pool')
    xticklabels({lab1, lab2})
    grid on

    sgtitle('Concentrations and variable values')

    fprintf('final steady states \n')
    fprintf('           SS1       SS2\n')
    fprintf('PTHg:      %0.3f     %0.3f\n', SS1(1), SS2(1))
    fprintf('PTHp_con:  %0.3f      %0.3f\n', SS1(2)/Vp_1, SS2(2)/Vp_2)
    fprintf('Cap_con:   %0.3f      %0.3f\n', SS1(3)/Vp_1, SS2(3)/Vp_2)
    fprintf('D3p_con:   %0.3f    %0.3f\n', SS1(4)/Vp_1, SS2(4)/Vp_2)
    fprintf('NCaf:      %0.3f      %0.3f\n', SS1(5), SS2(5))
    fprintf('NCas:      %0.3f    %0.3f\n', SS1(6), SS2(6))
end

if calcium_vals
    figure(2)
    clf
    x = [1, 2];
    nrows = 2;ncols=3;
    subplot(nrows,ncols,1)
    bar(x,[vals1.Gut_absorption,vals2.Gut_absorption],w1,'facecolor',c1)
    ylabel('Gut\_absorption')
    title('Gut absorption')
    xticklabels({lab1, lab2})
    grid on
    
    subplot(nrows,ncols,2)
    bar(x,[vals1.Bone_resorption,vals2.Bone_resorption],w1,'facecolor',c1)
    ylabel('Bone\_resorption')
    title('Bone resorption')
    xticklabels({lab1, lab2})
    grid on
    
    subplot(nrows,ncols,3)
    bar(x,[vals1.FastPool_to_Plasma,vals2.FastPool_to_Plasma],w1,'facecolor',c1)
    ylabel('FastPool\_to\_Plasma')
    xticklabels({lab1, lab2})
    title('Fast Pool to Plasma')
    grid on
    
    subplot(nrows,ncols,4)
    bar(x,[vals1.Plasma_to_FastPool,vals2.Plasma_to_FastPool],w1,'facecolor',c1)
    ylabel('Plasma\_to\_FastPool')
    xticklabels({lab1, lab2})
    title('Plasma to Fast Pool')
    grid on
    
    subplot(nrows,ncols,5)
    bar(x,[vals1.Urine_excretion,vals2.Urine_excretion],w1,'facecolor',c1)
    ylabel('Urine\_excretion')
    xticklabels({lab1, lab2})
    title('Urine excretion')
    grid on

    sgtitle('Calcium fluxes')

end

if kidney_vals
    fprintf('Renal calcium transport values \n')
    fprintf('                    SS1                SS2   \n')
    fprintf('Renal filtration:   %0.3d        %0.3d  \n', vals1.Renal_filtration, vals2.Renal_filtration)
    fprintf('PT frac reab:       %0.3f            %0.3f\n',vals1.Lambda_PT,vals2.Lambda_PT)
    fprintf('TAL frac reab:      %0.3f            %0.3f\n', vals1.Lambda_TAL,vals2.Lambda_TAL)
    fprintf('DT frac reab:       %0.3f            %0.3f\n', vals1.Lambda_DCT,vals2.Lambda_DCT)
    fprintf('Renal frac reab:    %0.3f            %0.3f \n', vals1.Renal_frac_reab, vals2.Renal_frac_reab)
    fprintf('Urine excretion:    %0.3d        %0.3d\n', vals1.Urine_excretion, vals2.Urine_excretion)

    GFRmale = 2.0e-3;
    LambdaPT0_male = 0.66; deltaPTmax_male = 0.035;
    LambdaTAL0_male = 0.185; deltaTALmax_male = 0.015;
    LambdaDCT0_male = 0.095; deltaDCTmax_male = 0.015;
    opt_male = [0.97 0.99 1.00; % frac renal reab
                1.2*GFRmale*0.005 1.2*GFRmale*0.0175 1.2*GFRmale*0.03;% urinary excretion ([Ca]_p*GFR*(1 - Frac_reab))
                LambdaPT0_male, LambdaPT0_male + deltaPTmax_male*0.75,LambdaPT0_male + deltaPTmax_male; % PT reabsorption
                LambdaTAL0_male, LambdaTAL0_male + deltaTALmax_male*0.75, LambdaTAL0_male + deltaTALmax_male;% TAL reabsorption
                LambdaDCT0_male, LambdaDCT0_male + deltaDCTmax_male*0.75, LambdaDCT0_male + deltaDCTmax_male]; % DT reabsorpton
                

    opt_female= [opt_male(1,:)*1.0; % frac reab
                opt_male(2,:)*0.75;% urine
                opt_male(3,:)*0.9; % PT frac reab
                opt_male(4,:)*1.175; % TAL frac reab
                opt_male(5,:)*1.3]; % DT frac reab
    
    opt_preg = [opt_female(1,:)*1.0; % frac reab
                opt_female(2,:)*1.7;% urine
                opt_female(3,:)*1.0; % PT frac reab
                opt_female(4,:)*1.0; % TAL frac reab
                opt_female(5,:)*1.0]; % DT frac reab

    opt_lact = [opt_female(1,:)*1.0; % frac reab
                opt_female(2,:)*1.0;% urine
                opt_female(3,:)*1.0; % PT frac reab
                opt_female(4,:)*1.0; % TAL frac reab
                opt_female(5,:)*1.0]; % DT frac reab

    if strcmp(SS1dat.sexORrep, 'male')
        opt1 = opt_male;
    elseif strcmp(SS1dat.sexORrep, 'female')
        opt1 = opt_female;
    elseif strcmp(SS1dat.sexORrep, 'preg')
        opt1 = opt_preg; 
    elseif strcmp(SS1dat.sexORrep, 'lact')
        opt1 = opt_lact;
    end
    if strcmp(SS2dat.sexORrep, 'male')
        opt2 = opt_male;
    elseif strcmp(SS2dat.sexORrep, 'female')
        opt2 = opt_female;
    elseif strcmp(SS2dat.sexORrep, 'preg')
        opt2 = opt_preg; 
    elseif strcmp(SS2dat.sexORrep, 'lact')
        opt2 = opt_lact;
    end

    figure(3)
    clf
    nrows = 2; ncols = 3;
    subplot(nrows,ncols,1)
    x = [1,2];
    bar(x,[vals1.Renal_filtration, vals2.Renal_filtration],w2,'facecolor',c2)
    xticklabels({lab1,lab2})
    ylabel('Renal\_filtration')
    title('Renal filtration')
    grid on

    subplot(nrows,ncols,2)
    bar(x,[opt1(3,2), opt2(3,2)],w1,'facecolor',c1)
    hold on
    bar(x,[vals1.Lambda_PT, vals2.Lambda_PT],w2,'facecolor',c2)
    errorbar(1,opt1(3,2),opt1(3,2)-opt1(3,1), opt1(3,3)-opt1(3,2),'color',ce);
    errorbar(2,opt2(3,2),opt2(3,2)-opt2(3,1), opt2(3,3)-opt2(3,2),'color',ce)
    xticklabels({lab1,lab2})
    ylabel('Lambda\_PT')
    title('Fractional PT reabsorption')
    ylim([0 0.7])
    grid on

    subplot(nrows,ncols,3)
    bar(x,[opt1(4,2), opt2(4,2)],w1,'facecolor',c1)
    hold on
    bar(x,[vals1.Lambda_TAL,vals2.Lambda_TAL],w2,'facecolor',c2)
    errorbar(1,opt1(4,2),opt1(4,2)-opt1(4,1), opt1(4,3)-opt1(4,2),'color',ce)
    errorbar(2,opt2(4,2),opt2(4,2)-opt2(4,1), opt2(4,3)-opt2(4,2),'color',ce)
    xticklabels({lab1,lab2})
    ylabel('Lambda\_TAL')
    title('Fractional TAL reabsorption')
    grid on

    subplot(nrows,ncols,4)
    bar(x,[opt1(5,2), opt2(5,2)],w1,'facecolor',c1)
    hold on
    bar(x,[vals1.Lambda_DCT, vals2.Lambda_DCT],w2,'facecolor',c2)
    errorbar(1,opt1(5,2),opt1(5,2)-opt1(5,1), opt1(5,3)-opt1(5,2),'color',ce)
    errorbar(2,opt2(5,2),opt2(5,2)-opt2(5,1), opt2(5,3)-opt2(5,2),'color',ce)
    xticklabels({lab1,lab2})
    ylabel('Lambda\_DCT')
    title('Fractional DT reabsorption')
    grid on

    subplot(nrows,ncols,5)
    bar(x,[opt1(1,2), opt2(1,2)],w1,'facecolor',c1)
    hold on
    bar(x,[vals1.Renal_frac_reab, vals2.Renal_frac_reab],w2,'facecolor',c2)
    errorbar(1,opt1(1,2),opt1(1,2)-opt1(1,1), opt1(1,3)-opt1(1,2),'color',ce);
    errorbar(2,opt2(1,2),opt2(1,2)-opt2(1,1), opt2(1,3)-opt2(1,2),'color',ce)
    xticklabels({lab1,lab2})
    ylim([0.5 1.01])
    ylabel('Renal\_frac\_reab')
    title('Renal fractional reabsorption')
    grid on

    subplot(nrows,ncols,6)
    bar(x,[opt1(2,2),opt2(2,2)],w1,'facecolor',c1)
    hold on
    bar(x,[vals1.Urine_excretion, vals2.Urine_excretion],w2,'facecolor',c2)
    errorbar(1,opt1(2,2),opt1(2,2)-opt1(2,1), opt1(2,3)-opt1(2,2),'color',ce);
    errorbar(2,opt2(2,2),opt2(2,2)-opt2(2,1), opt2(2,3)-opt2(2,2),'color',ce)
    xticklabels({lab1,lab2})
    ylabel('Urine\_excretion')
    title('Urinary excretion')
    grid on

    sgtitle('Kidney segmental reabsorption')
end

if kidney_impact_vals
    figure(4)
    clf
    nrows = 2; ncols = 3;
    subplot(nrows,ncols,1)
    x = [1,2];
    bar(x,[vals1.delta_PT_PTH, vals2.delta_PT_PTH],w2,'facecolor',c2)
    xticklabels({lab1,lab2})
    ylabel('delta\_PT\_PTH')
    title('PTH impact on PT reab')
    grid on

    subplot(nrows,ncols,2)
    bar(x,[vals1.delta_TAL_Ca, vals2.delta_TAL_Ca],w2,'facecolor',c2)
    xticklabels({lab1,lab2})
    ylabel('delta\_TAL\_Ca')
    title('Ca impact on TAL reab')
    grid on

    subplot(nrows,ncols,3)
    bar(x,[vals1.delta_TAL_PTH,vals2.delta_TAL_PTH],w2,'facecolor',c2)
    xticklabels({lab1,lab2})
    ylabel('delta\_TAL\_PTH')
    title('PTH impact on TAL reab')
    grid on

    subplot(nrows,ncols,4)
    bar(x,[vals1.delta_DCT_PTH, vals2.delta_DCT_PTH],w2,'facecolor',c2)
    xticklabels({lab1,lab2})
    ylabel('delta\_DCT\_PTH')
    title('PTH impact on DT reab')
    grid on

    subplot(nrows,ncols,5)
    bar(x,[vals1.delta_DCT_D3, vals2.delta_DCT_D3],w2,'facecolor',c2)
    xticklabels({lab1,lab2})
    ylabel('delta\_DCT\_D3')
    title('D3 impact on DT reab')
    grid on
    
    sgtitle('Kidney impact vals')

end

if gut_vals
    frac_gut_abs_opt_male = [0.4 0.5 0.6];
    frac_gut_abs_opt_female = frac_gut_abs_opt_male * 1.3;
    frac_gut_abs_opt_preg   = frac_gut_abs_opt_female * 1.3;
    frac_gut_abs_opt_lact   = frac_gut_abs_opt_female * 1.3; % TO DO

    if strcmp(SS1dat.sexORrep, 'male')
        opt1 = frac_gut_abs_opt_male;
    elseif strcmp(SS1dat.sexORrep, 'female')
        opt1 = frac_gut_abs_opt_female;
    elseif strcmp(SS1dat.sexORrep, 'preg')
        opt1 = frac_gut_abs_opt_preg; 
    elseif strcmp(SS1dat.sexORrep, 'lact')
        opt1 = frac_gut_abs_opt_lact;
    end
    if strcmp(SS2dat.sexORrep, 'male')
        opt2 = frac_gut_abs_opt_male;
    elseif strcmp(SS2dat.sexORrep, 'female')
        opt2 = frac_gut_abs_opt_female;
    elseif strcmp(SS2dat.sexORrep, 'preg')
        opt2 = frac_gut_abs_opt_preg; 
    elseif strcmp(SS2dat.sexORrep,'lact')
        opt2 = frac_gut_abs_opt_lact;
    end

    fprintf('\n Gut SS vals \n')
    fprintf('                SS1          SS2 \n')
    fprintf('Frac gut abs:   %0.3f        %0.3f \n', vals1.Gut_frac_absorption, vals2.Gut_frac_absorption)
    fprintf('Net gut abs:    %0.3d    %0.3d \n', vals1.Gut_absorption, vals2.Gut_absorption)

    figure(5)
    clf
    nrows = 1; ncols = 3;
    subplot(nrows,ncols,1)
    x = [1,2];
    bar(x,[vals1.Gut_impact_D3, vals2.Gut_impact_D3],w2,'facecolor',c2)
    xticklabels({lab1,lab2})
    ylabel('Gut\_impact\_D3')
    title('D3 impact on gut absorption')
    grid on

    subplot(nrows,ncols,2)
    bar(x,[opt1(2) opt2(2)], w1, 'facecolor', c1)
    hold on
    bar(x,[vals1.Gut_frac_absorption, vals2.Gut_frac_absorption],w2,'facecolor',c2)
    errorbar(1, opt1(2), opt1(2)-opt1(1), opt1(3)-opt1(2),'color',ce);
    errorbar(2, opt2(2), opt2(2)-opt2(1),opt2(3)-opt2(2),'color',ce)
    xticklabels({lab1,lab2})
    ylabel('Gut\_frac\_absorption')
    title('Fractional gut absorption')
    grid on

    subplot(nrows,ncols,3)
    bar(x,[vals1.Gut_absorption,vals2.Gut_absorption],w2,'facecolor',c2)
    xticklabels({lab1,lab2})
    ylabel('Gut\_absorption')
    title('Net gut absorption')
    grid on
  
    sgtitle('Gut vals')
end

if D3_vals
    figure(6)
    clf
    nrows = 2;ncols=4;
    subplot(nrows,ncols,1)
    bar(x,[vals1.PTH_impact_D3,vals2.PTH_impact_D3],w1,'facecolor',c1)
    ylabel('PTH\_impact\_D3')
    title('PTH impact on D3 synthesis')
    xticklabels({lab1, lab2})
    grid on
    
    subplot(nrows,ncols,2)
    bar(x,[vals1.Ca_impact_D3,vals2.Ca_impact_D3],w1,'facecolor',c1)
    ylabel('Ca\_impact\_D3')
    title('Ca^{2+} impact on D3 synthesis')
    xticklabels({lab1, lab2})
    grid on
    
    subplot(nrows,ncols,3)
    bar(x,[vals1.D3_impact_D3,vals2.D3_impact_D3],w1,'facecolor',c1)
    ylabel('D3\_impact\_D3')
    xticklabels({lab1, lab2})
    title('D3 impact on D3 synthesis')
    grid on
    
    subplot(nrows,ncols,4)
    bar(x,[vals1.Rconv,vals2.Rconv],w1,'facecolor',c1)
    ylabel('Rconv')
    xticklabels({lab1, lab2})
    title('Rconv')
    grid on
    
    subplot(nrows,ncols,5)
    bar(x,[vals1.D3_synthesis,vals2.D3_synthesis],w1,'facecolor',c1)
    ylabel('D3\_synthesis')
    xticklabels({lab1, lab2})
    title('D3 synthesis')
    grid on

    subplot(nrows,ncols,6)
    bar(x,[vals1.PTH_impact_D3_degradation,vals2.PTH_impact_D3_degradation],w1,'facecolor',c1)
    ylabel('PTH\_impact\_D3\_degradation')
    xticklabels({lab1, lab2})
    title('PTH impact on D3 degradation')
    grid on


    subplot(nrows,ncols,7)
    bar(x,[vals1.D3_degradation,vals2.D3_degradation],w1,'facecolor',c1)
    ylabel('D3\_degradation')
    xticklabels({lab1, lab2})
    title('D3 degradation')
    grid on

    sgtitle('D3 values')

end

if PTH_vals
    fprintf('\n PTH SS vals \n')
    fprintf('                        SS1         SS2 \n')
    fprintf('nCanorm:                %0.3f       %0.3f \n',vals1.n_Ca_norm, vals2.n_Ca_norm)
    fprintf('CaSR_PTHg_regulation:   %0.3f       %0.3f \n', vals1.CaSR_PTHg_regulation, vals2.CaSR_PTHg_regulation)
    fprintf('FCa:                    %0.3d   %0.3d \n',vals1.F_Ca,vals2.F_Ca)
    fprintf('PTHg exocytosis:        %0.3f       %0.3f \n', vals1.PTHg_exocytosis, vals2.PTHg_exocytosis)

    figure(7)
    clf
    nrows = 2; ncols = 4;
    subplot(nrows,ncols,1)
    bar(x,[vals1.PTHg_prod_effect_D3, vals2.PTHg_prod_effect_D3], w2,'facecolor',c2)
    ylabel('PTH\_prod\_effect\_D3')
    xticklabels({lab1,lab2})
    title('D3 effect on PTH prod')
    grid on

    subplot(nrows,ncols,2)
    bar(x,[vals1.PTHg_synthesis, vals2.PTHg_synthesis],w2,'facecolor',c2)
    ylabel('PTHg\_synthesis')
    xticklabels({lab1, lab2})
    title('PTHg synthesis')
    grid on

    subplot(nrows,ncols,3)
    bar(x,[vals1.n_Ca_norm, vals2.n_Ca_norm],w2,'facecolor',c2)
    ylabel('nCa\_norm')
    xticklabels({lab1, lab2})
    title('nCa\_norm')
    grid on

    subplot(nrows,ncols,4)
    bar(x,[vals1.CaSR_PTHg_regulation, vals2.CaSR_PTHg_regulation],w2,'facecolor',c2)
    ylabel('CaSR\_PTHg\_regulation')
    xticklabels({lab1, lab2})
    title('CaSR regulation of PTHg')
    grid on


    subplot(nrows,ncols,5)
    bar(x,[vals1.F_Ca, vals2.F_Ca],w2,'facecolor',c2)
    ylabel('F\_Ca')
    xticklabels({lab1, lab2})
    title('F\_Ca')
    grid on

    subplot(nrows,ncols,6)
    bar(x,[vals1.PTHg_exocytosis, vals2.PTHg_exocytosis],w2,'facecolor',c2)
    ylabel('PTHg\_exocytosis')
    xticklabels({lab1, lab2})
    title('PTHg exocytosis')
    grid on

    subplot(nrows,ncols,7)
    bar(x,[vals1.PTHg_degradation, vals2.PTHg_degradation],w2,'facecolor',c2)
    ylabel('PTHg\_degradation')
    xticklabels({lab1, lab2})
    title('PTHg degradation')
    grid on


    sgtitle('PTH values')
end

if bone_vals

    fprintf('\n Bone SS vals \n')
    fprintf('                        SS1         SS2 \n')
    fprintf('Plasma_to_FastPool:     %0.3d     %0.3d \n',vals1.Plasma_to_FastPool, vals2.Plasma_to_FastPool)
    fprintf('FastPool_to_Plasma:     %0.3d     %0.3d \n', vals1.FastPool_to_Plasma, vals2.FastPool_to_Plasma)
    fprintf('Bone_accretion:         %0.3d     %0.3d \n',vals1.Bone_accretion,vals2.Bone_accretion)
    fprintf('Bone_resorption:        %0.3d     %0.3d \n', vals1.Bone_resorption, vals2.Bone_resorption)

    figure(8)
    clf
    nrows = 2; ncols = 2;
    subplot(nrows,ncols,1)
    bar(x,[vals1.Plasma_to_FastPool, vals2.Plasma_to_FastPool], w2,'facecolor',c2)
    ylabel('Plasma\_to\_FastPool')
    xticklabels({lab1,lab2})
    title('Plasma to FastPool')
    grid on

    subplot(nrows,ncols,2)
    bar(x,[vals1.FastPool_to_Plasma, vals2.FastPool_to_Plasma],w2,'facecolor',c2)
    ylabel('FastPool\_to\_Plasma')
    xticklabels({lab1, lab2})
    title('FastPool to Plasma')
    grid on

    subplot(nrows,ncols,3)
    bar(x,[vals1.Bone_accretion, vals2.Bone_accretion],w2,'facecolor',c2)
    ylabel('Bone\_accretion')
    xticklabels({lab1, lab2})
    title('Bone accretion')
    grid on

    subplot(nrows,ncols,4)
    bar(x,[vals1.Bone_resorption, vals2.Bone_resorption],w2,'facecolor',c2)
    ylabel('Bone\_resorption')
    xticklabels({lab1, lab2})
    title('Bone resorption')
    grid on

    sgtitle('Bone values')

end



