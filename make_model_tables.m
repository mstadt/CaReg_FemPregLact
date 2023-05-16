% script to make tables with SS and parameter values
% exports to .xlsx files
clear all;

% filenames
% male file
fn1 = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-male_notes-newmale.mat'; 
dat_male = load(fn1);
% female file
fn2 = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-female_notes-newfemale.mat'; 
dat_female = load(fn2);
% preg file
fn3 = './SSbest/16-May-2023_calcium_mod_SS_sexORrep-preg_notes-PTHupdate2.mat';
dat_preg = load(fn3);
% lact file
fn4 = './SSbest/16-May-2023_calcium_mod_SS_sexORrep-lact_notes-PTHupdate2.mat';
dat_lact = load(fn4);

% check par names match for all simulations
parcheck = dat_male.param_names;
flag = 0;
for ii = 1:size(parcheck,2)
    if ~strcmp(parcheck{ii}, dat_female.param_names{ii})
        flag = 1;
        fprintf('parval %i does not match in female and male \n', ii)
    end
    if ~strcmp(parcheck{ii},dat_preg.param_names{ii})
        flag = 1;
        fprintf('parval %i does not match in preg and male \n', ii)
    end
    if ~strcmp(parcheck{ii},dat_lact.param_names{ii})
        flag = 1;
        fprintf('parval %i does not match in lact and male \n', ii)
    end
end
if flag
    error('parameter names not matching')
end

%% Params table variables
VarNames = {'Male', 'Female', 'Pregnant', 'Lactation',...
                    'Female2Male', 'Preg2Female', 'Lact2Female', 'Lact2Preg'};
RowNames = dat_male.param_names';

ParNames   = dat_male.param_names';
MalePars   = dat_male.params';
FemalePars = dat_female.params';
PregPars   = dat_preg.params';
LactPars   = dat_lact.params';

F2Mratio  = round(FemalePars./MalePars,3,'decimals');
temp = find(F2Mratio == 1.0);
F2Mratio(temp) = NaN;

P2Fratio  = round(PregPars./FemalePars,3,'decimals');
temp = P2Fratio == 1.0;
P2Fratio(temp) = NaN;

L2Fratio  = round(LactPars./FemalePars,3,'decimals');
temp = find(L2Fratio == 1.0);
L2Fratio(temp) = NaN;

L2Pratio  = round(LactPars./PregPars,3,'decimals');
temp = find(L2Pratio == 1.0);
L2Pratio(temp) = NaN;

T_pars = table(MalePars, FemalePars, PregPars, LactPars, ...
                    F2Mratio, P2Fratio, L2Fratio, L2Pratio,...
                    'VariableNames', VarNames, ...
                    'RowNames', RowNames);

% fname = strcat(date, '_model_parameters.xlsx');
% writetable(T_pars, fname,...
%                     'WriteRowNames', true)


%% SS values table variables

VarNames = {'Male', 'Female', 'Pregnant', 'Lactation',...
                    'Female2Male', 'Preg2Female', 'Lact2Female', 'Lact2Preg'};

RowNames = {'[PTHp]', '[Cap]', '[1,25(OH)2D3]',...
                'PTHg', 'PTHp', 'Cap', 'D3p', 'NCaf', 'NCas',...
                'GutFracAbs', 'GutAbs', ...
                'BoneResorption', 'Plasma2FastPool', 'FastPool2Plasma', ...
                'RenalFracReab', 'UrineExcretion', ...
                'BoneAccretion'};

MaleVals   = getvals(dat_male);
FemaleVals = getvals(dat_female);
PregVals   = getvals(dat_preg);
LactVals   = getvals(dat_lact);

F2Mratio  = round(FemaleVals./MaleVals,3,'decimals');
temp = find(F2Mratio == 1.0);
F2Mratio(temp) = NaN;

P2Fratio  = round(PregVals./FemaleVals,3,'decimals');
temp = P2Fratio == 1.0;
P2Fratio(temp) = NaN;

L2Fratio  = round(LactVals./FemaleVals,3,'decimals');
temp = find(L2Fratio == 1.0);
L2Fratio(temp) = NaN;

L2Pratio  = round(LactVals./PregVals,3,'decimals');
temp = find(L2Pratio == 1.0);
L2Pratio(temp) = NaN;

T_vars = table(MaleVals, FemaleVals, PregVals, LactVals, ...
                    F2Mratio, P2Fratio, L2Fratio, L2Pratio,...
                    'VariableNames', VarNames, ...
                    'RowNames', RowNames);
% fname = strcat(date, '_model_variables.xlsx');
% writetable(T_vars, fname,...
%                     'WriteRowNames', true)


fname = strcat(date, '_model_results.xlsx');
writetable(T_pars, fname,'WriteRowNames',true, 'Sheet', 'parameters');
writetable(T_vars, fname, 'WriteRowNames', true, 'Sheet', 'variables');


%--------------------------
% functions used
%--------------------------
function vals = getvals(dat)
    SS = dat.SS;
    params = dat.params;
    param_names = dat.param_names;
    valsSS = dat.valsSS;

    Ind_Vp = 0;
        for ii = 1:size(param_names,2)
            if strcmp(param_names{ii}, 'Vp')
                Ind_Vp = ii;
            end
        end
        if Ind_Vp ~= 18
            fprintf('Ind_VP: %i, expected: 18 \n', Ind_Vp)
            error('Vp index not same')
        end
        Vp = params(Ind_Vp);
        vals = [SS(2)/Vp;
                SS(3)/Vp;
                SS(4)/Vp;
                SS;
                valsSS.Gut_frac_absorption;
                valsSS.Gut_absorption;
                valsSS.Bone_resorption;
                valsSS.Plasma_to_FastPool;
                valsSS.FastPool_to_Plasma;
                valsSS.Renal_frac_reab;
                valsSS.Urine_excretion;
                valsSS.Bone_accretion];
    
end
