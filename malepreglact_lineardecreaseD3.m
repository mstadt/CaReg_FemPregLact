
% Male preg lact analysis with a linear decrease in D3p

clear all;

% get parameter values
sexORrep = 'male';
run('read_in_params');
pars_male = params;
clearvars -except pars_male

sexORrep = 'female';
run('read_in_params');
pars_female = params;
clearvars -except pars_female pars_male

sexORrep = 'preg';
run('read_in_params');
pars_preg = params;
clearvars -except pars_female pars_preg pars_male

sexORrep = 'lact';
run('read_in_params');
pars_lact = params;
clearvars -except pars_female pars_preg pars_lact param_names pars_male


% compute preg 2 female and lact 2 female ratios
preg2vir = pars_preg./pars_female;
preg2vir(37) = -1; % FetusORMilk

lact2vir = pars_lact./pars_female;
lact2vir(37) = -1; % FetusORMilk

% set preg/lact male parameters
par_malepreg = pars_male.*preg2vir;
par_malepreg(37) = pars_preg(37); % set FetusORMilk

par_malelact = pars_male.*lact2vir;
par_malelact(37) = pars_lact(37); % set FetusORMilk

% variables to fix
NCas_fixed = 1;

% figure specs
cmap = parula(4);
lw = 3.0; c1 = cmap(1,:); c2 = cmap(3,:);
%% male pregnancy
% get initial condition
SSfile = './results_malepreglact/16-May-2023_calcium_mod_SS_sexORrep-malepreg_notes-pregmale.mat';
IC = load(SSfile).SS;
Vp_malepreg = par_malepreg(18); % set volume for concentrations

tspan = [0 10000];
options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9);
fprintf('solving male preg 2 female preg model eqns \n')

tspan_decrease = 6000; % how long to decrease D3p con over
[tpreg, ypreg] = ode15s(@(t,y) calcium_mod_fixedD3(t,y,par_malepreg,...
                                        'preg', tspan_decrease,...
                                        'NCas_fixed', NCas_fixed),...
                                        tspan, IC, options);

figure(50)
clf;
nrows = 1; ncols=3;
subplot(nrows,ncols,1)
plot(tpreg, ypreg(:,2)/Vp_malepreg, 'linewidth', lw, 'color', c1)
ylabel('[PTH]_p (pmol/L)')
xlabel('t')
title('Plasma PTH concentration')
grid on

subplot(nrows,ncols,2)
plot(tpreg, ypreg(:,3)/Vp_malepreg, 'linewidth', lw, 'color', c1)
ylabel('[Ca^{2+}]_p (pmol/L)')
xlabel('t')
title('Plasma calcium concentration')
grid on

rep = 'preg';
if strcmp(rep, 'preg')
    startD3p = 481.372; % male preg D3p
    endD3p   = 269.106; % female preg D3p
elseif strcmp(rep, 'lact')
    startD3p = 573.95; % male lact D3p
    endD3p   = 329.17; % female lact D3p
end
m = (startD3p - endD3p)/tspan_decrease; % slope
b = startD3p;

D3p_con_vals = max(endD3p, b - m.*tpreg);
subplot(nrows,ncols,3)
plot(tpreg, D3p_con_vals, 'linewidth',lw,'color',c1)
ylabel('[1,25(OH)_2D_3]_p (pmol/L)')
title('Plasma 1,25(OH)_2D_3 concentration')
grid on

sgtitle('Pregnancy')

%% male lactation
% get initial condition
SSfile = './results_malepreglact/16-May-2023_calcium_mod_SS_sexORrep-malelact_notes-lactmale.mat';
IC = load(SSfile).SS;
Vp_malelact = par_malelact(18); % set volume for concentrations

tspan = [0 10000];
options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9);
fprintf('solving male lact 2 female lact model eqns \n')

tspan_decrease = 6000; % how long to decrease D3p con over
[tlact, ylact] = ode15s(@(t,y) calcium_mod_fixedD3(t,y,par_malelact,...
                                        'lact', tspan_decrease,...
                                        'NCas_fixed', NCas_fixed),...
                                        tspan, IC, options);

figure(51)
clf;
nrows = 1; ncols=3;
lw = 3.0; 
subplot(nrows,ncols,1)
plot(tlact, ylact(:,2)/Vp_malelact, 'linewidth', lw, 'color', c2)
ylabel('[PTH]_p (pmol/L)')
xlabel('t')
title('Plasma PTH concentration')
grid on

subplot(nrows,ncols,2)
plot(tlact, ylact(:,3)/Vp_malelact, 'linewidth', lw, 'color', c2)
ylabel('[Ca^{2+}]_p (pmol/L)')
xlabel('t')
title('Plasma calcium concentration')
grid on

rep = 'lact';
if strcmp(rep, 'preg')
    startD3p = 481.372; % male preg D3p
    endD3p   = 269.106; % female preg D3p
elseif strcmp(rep, 'lact')
    startD3p = 573.95; % male lact D3p
    endD3p   = 329.17; % female lact D3p
end
m = (startD3p - endD3p)/tspan_decrease; % slope
b = startD3p;

D3p_con_vals = max(endD3p, b - m.*tlact);
subplot(nrows,ncols,3)
plot(tlact, D3p_con_vals, 'linewidth',lw,'color',c2)
ylabel('[1,25(OH)_2D_3]_p (pmol/L)')
title('Plasma 1,25(OH)_2D_3 concentration')
grid on

sgtitle('Lactation')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions used

function dydt = calcium_mod_fixedD3(t,y,params, rep, tspan_decrease, varargin)
% calcium homeostasis with fixed D3 for the male preg and male lact
% simulations


% get variable inputs
Ca_fixed     = false; % turn on to have no calcium dynamics
D3_fixed     = true; % D3 is always fixed!
PTHg_fixed   = false; % turn on to have no PTHg dynamics
PTHp_fixed   = false; % turn on to have no PTHp dynamics
NCaf_fixed   = false; % turn on to have no NCaf dynamics
NCas_fixed   = false; % turn on to have no NCas dynamics
EGTA_inject  = false; % turn on if EGTA injection
alt_sim      = false; % alt_sim if need to test equations
for ii = 1:2:length(varargin)
    if strcmp(varargin{ii},'Ca_fixed')
        Ca_fixed = varargin{ii+1};
        if Ca_fixed
            fprintf('NOTE: running model with fixed calcium \n')
        end
    elseif strcmp(varargin{ii},'PTHg_fixed')
        PTHg_fixed = varargin{ii+1};
        if PTHg_fixed
            fprintf('NOTE: running model with fixed PTHg \n')
        end
    elseif strcmp(varargin{ii}, 'PTHp_fixed')
        PTHp_fixed = varargin{ii+1};
        if PTHp_fixed
            fprintf('NOTE: running model with fixed PTHp \n')
        end
    elseif strcmp(varargin{ii}, 'D3_fixed')
        D3_fixed = varargin{ii+1};
        if D3_fixed
            fprintf('NOTE: running model with fixed D3 \n')
        end
    elseif strcmp(varargin{ii}, 'NCaf_fixed')
        NCaf_fixed = varargin{ii+1};
        if NCaf_fixed
            fprintf('NOTE: running model with fixed NCaf \n')
        end
    elseif strcmp(varargin{ii}, 'NCas_fixed')
        NCas_fixed = varargin{ii+1};
%         if NCas_fixed
%             %fprintf('NOTE: running model with fixed NCas \n')
%         end
    elseif strcmp(varargin{ii}, 'EGTA_inject')
       temp = varargin{ii+1};
       EGTA_inject = temp(1);
       k_EGTA_inject = temp(2);
       if k_EGTA_inject == -1
           % function form
           tdose = t-202;
           k_EGTA_inject = 1.6e-4*(exp(0.001*(tdose-10)) + exp(0.0025*(tdose-50)));
       elseif k_EGTA_inject == -2
           % female scaling
           tdose = t-202;
           k_EGTA_inject = 0.65 * (1.6e-4*(exp(0.001*(tdose-10)) + exp(0.0025*(tdose-50))));
       end
    elseif strcmp(varargin{ii},'alt_sim')
        alt_sim = varargin{ii+1};
        if alt_sim
            fprintf('running alt sim \n')
        end
    else
        fprintf('varargin{ii}: %s \n', varargin{ii})
        error('varargin not recognized')
    end
end

% variable names
PTH_g = y(1); % amount of PTH in parathyroid gland pool (pmol)
PTH_p = y(2); % amount of PTH in plasma (pmol)
Ca_p  = y(3); % amount of calcium in plasma (mmol)
D3_p  = y(4); % amount of vitamin D3 in plasma (pmol)
NCa_f = y(5); % amount of calcium in the rapidly excangeable pool (mmol)
NCa_s = y(6); % amount of calcium in the slow pool (mmol)
if EGTA_inject
    EGTAp_con   = y(7);
    CaEGTAp_con = y(8);
end


% parameters -- ordered from sex-specific parameters
k_PTHg_deg = params(1);
rho_exo = params(2);
R = params(3);
k_PTHp_deg = params(4);
Gamma_res_min = params(5);
delta_res_max = params(6);
kappa_b = params(7);
nconv = params(8);
gamma_conv_Ca = params(9);
k_deg_D3 = params(10);
k_pf_Ca = params(11);
k_fp_Ca = params(12);
nPT = params(13);
Cap_ref = params(14);
nTAL = params(15);
k_EGTA_on = params(16);
k_EGTA_off = params(17);
Vp = params(18);
GFR = params(19);
gamma_conv_D3 = params(20);
delta_conv_max = params(21);
k_conv_min = params(22);
D3_inact_p = params(23);
gamma_prod_D3 = params(24);
ICa = params(25);
Gamma_abs0 = params(26);
delta_abs_D3 = params(27);
K_abs_D3 = params(28);
K_D3p_res = params(29);
Lambda_PT0 = params(30);
delta_PT_max = params(31);
Lambda_TAL0 = params(32);
delta_TAL_max = params(33);
delta_DCT_max = params(34);
K_DCT_D3p = params(35);
Lambda_DCT0 = params(36);
FetusORMilk = params(37);
K_Ca_CASR = params(38);
K_conv_PTH = params(39);
k_prod_PTHg = params(40);
K_PTHp_res = params(41);
gamma_deg_PTHp = params(42);
PTHp_ref = params(43);
K_TAL_PTHp = params(44);
K_DCT_PTHp = params(45);
n1_exo = params(46);
n2_exo = params(47);
beta_exo_PTHg = params(48);
gamma_exo_PTHg = params(49);
Gamma_ac = params(50);



% model equations
dydt = zeros(length(y),1);


% change to concentrations
PTHp_con = PTH_p/Vp;
Cap_con = Ca_p/Vp;
%D3p_con  = D3_p/Vp;
% set fixed D3p_con
% male preg D3p = 481.372;
% female preg D3p = 269.106
% male lact D3p = 573.95
% female lact = 329.17
%tspan_decrease = 6000; % minutes, how long to decrease over
if strcmp(rep, 'preg')
    startD3p = 481.372; % male preg D3p
    endD3p   = 269.106; % female preg D3p
elseif strcmp(rep, 'lact')
    startD3p = 573.95; % male lact D3p
    endD3p   = 329.17; % female lact D3p
end
m = (startD3p - endD3p)/tspan_decrease; % slope
b = startD3p;

D3p_con = max(endD3p, b - m*t);
if t>2000
    disp(D3p_con)
end

% Parathyroid gland PTH (PTH_g)
if PTHg_fixed
    dydt(1) = 0;
else
    PTHg_prod_effect_D3 = 1/(1 + gamma_prod_D3 * D3p_con);
    PTHg_basal_synthesis = k_prod_PTHg; % fit for this parameter
    PTHg_synthesis = PTHg_basal_synthesis*PTHg_prod_effect_D3;
    
    PTHg_degradation = k_PTHg_deg*PTH_g;
    n_Ca_norm = n1_exo/(1 + exp(-rho_exo*(R - Cap_con))) + n2_exo;
    CaSR_PTHg_regulation =  Cap_con^n_Ca_norm / (Cap_con^n_Ca_norm + (K_Ca_CASR)^n_Ca_norm);
    F_Ca = beta_exo_PTHg - gamma_exo_PTHg * CaSR_PTHg_regulation;
    PTHg_exocytosis = F_Ca * PTH_g;
    
    dydt(1) = PTHg_synthesis - PTHg_exocytosis - PTHg_degradation;
end

% Plasma PTH amount (PTH_p)
if PTHp_fixed
    dydt(2) = 0;
else
    PTHp_influx = PTHg_exocytosis;
    PTHp_degradation = k_PTHp_deg*PTHp_con;
    
    dydt(2) = PTHp_influx - PTHp_degradation;
end


% EGTA reaction
if EGTA_inject
    %fprintf('**doing EGTA!!!**')
    CaEGTA_form = k_EGTA_on * Cap_con * EGTAp_con;
    CaEGTA_diss = k_EGTA_off * CaEGTAp_con;
    CaEGTA_reaction = CaEGTA_form - CaEGTA_diss;
    dydt(7) = k_EGTA_inject/Vp - CaEGTA_reaction; % no degradataion so just stays high... maybe add this later?
    dydt(8) = CaEGTA_reaction;
end


% Plasma calcium amount (Ca_p)
    % intestinal calcium absorption
    Gut_impact_D3 = D3p_con^2/(K_abs_D3^2 + D3p_con^2);
    Gut_frac_absorption = Gamma_abs0 + (delta_abs_D3)*Gut_impact_D3;
    Gut_absorption = ICa*Gut_frac_absorption;
    
    % bone resorption
    PTHp_res_effect = delta_res_max*0.2*PTHp_con^2/(PTHp_con^2 + K_PTHp_res^2);
    D3_res_effect   = delta_res_max*0.8*D3p_con^2/(D3p_con^2 + K_D3p_res^2);
    Bone_resorption = Gamma_res_min + (PTHp_res_effect + D3_res_effect);

    % renal calcium handling
    Renal_filtration = GFR * Cap_con;
    delta_PT_PTH = delta_PT_max/(1 + (PTHp_con/PTHp_ref)^nPT);
    Lambda_PT = Lambda_PT0 + delta_PT_PTH;
    delta_TAL_Ca = 0.7*delta_TAL_max/(1 + (Cap_con/Cap_ref)^nTAL);
    delta_TAL_PTH = 0.3*delta_TAL_max*PTHp_con/(PTHp_con + K_TAL_PTHp);
    Lambda_TAL = Lambda_TAL0 + delta_TAL_Ca + delta_TAL_PTH;
    delta_DCT_PTH = 0.8*delta_DCT_max*PTHp_con/(PTHp_con + K_DCT_PTHp);
    delta_DCT_D3  = 0.2*delta_DCT_max*D3p_con/(D3p_con + K_DCT_D3p);
    Lambda_DCT = Lambda_DCT0 + delta_DCT_PTH + delta_DCT_D3;
    Renal_frac_reab = min(0.995, Lambda_PT + Lambda_TAL + Lambda_DCT);
    Urine_excretion = (1-Renal_frac_reab)*Renal_filtration;
    
    % bone fast pool
    Plasma_to_FastPool = k_pf_Ca*Cap_con;
    FastPool_to_Plasma = k_fp_Ca*NCa_f;

    if ~EGTA_inject
        CaEGTA_reaction = 0;
    end

if Ca_fixed
    dydt(3) = 0;
else
    dydt(3) = (1 - kappa_b)*(Gut_absorption + Bone_resorption + FastPool_to_Plasma...
            - Plasma_to_FastPool - Urine_excretion - FetusORMilk) - CaEGTA_reaction*Vp;
end


% D3_p
if D3_fixed
    dydt(4) = 0;
else
    PTH_impact_D3 = PTHp_con^nconv/(PTHp_con^nconv + K_conv_PTH^nconv);
    Ca_impact_D3  = 1/(1 + gamma_conv_Ca*Cap_con);
    D3_impact_D3  = 1/(1 + gamma_conv_D3*D3p_con);
    Rconv = k_conv_min + delta_conv_max*PTH_impact_D3*Ca_impact_D3*D3_impact_D3; % eq 3.73 thesis
    D3_synthesis = Rconv * D3_inact_p;

    PTH_impact_D3_degradation = 1./(1 + gamma_deg_PTHp*PTHp_con);
    D3_degradation = k_deg_D3 * D3p_con * PTH_impact_D3_degradation;
    dydt(4) = D3_synthesis - D3_degradation; % eq 3.74 thesis 
end

% NCa_f
if NCaf_fixed
    dydt(5) = 0;
else
    Bone_accretion = Gamma_ac*NCa_f;
    dydt(5) = Plasma_to_FastPool - FastPool_to_Plasma - Bone_accretion;
end

% NCa_s
if NCas_fixed
    dydt(6) = 0;
else
    dydt(6) = Bone_accretion - Bone_resorption;
end

end % calcium_mod