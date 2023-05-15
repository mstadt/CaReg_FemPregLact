function vals = compute_vars(yvals, params)
% Input: yvals -- output from ODE model, params -- parameter vector

% name parameters
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


% variable names
PTHg_vals = yvals(:,1);
PTHp_vals = yvals(:,2); PTHp_con_vals = PTHp_vals./Vp;
Cap_vals = yvals(:,3); Cap_con_vals = Cap_vals./Vp;
D3p_vals = yvals(:,4); D3p_con_vals = D3p_vals./Vp;
NCaf_vals = yvals(:,5);

%PTHg_basal_synthesis = k_prod_PTHg; % fit for this parameter
%PTHg_degradation = k_PTHg_deg*PTH_g;
vals.PTHg_prod_effect_D3 = ones(size(D3p_con_vals))./(ones(size(D3p_con_vals)) + gamma_prod_D3.*D3p_con_vals);
PTHg_basal_synthesis = k_prod_PTHg;
vals.PTHg_synthesis = PTHg_basal_synthesis*vals.PTHg_prod_effect_D3;

vals.n_Ca_norm = n1_exo*ones(size(Cap_con_vals))./(ones(size(Cap_con_vals)) + exp(-rho_exo.*(R*ones(size(Cap_con_vals)) - Cap_con_vals))) + n2_exo;
vals.CaSR_PTHg_regulation =  Cap_con_vals.^vals.n_Ca_norm ./ (Cap_con_vals.^vals.n_Ca_norm + (K_Ca_CASR).^vals.n_Ca_norm);
vals.F_Ca = beta_exo_PTHg - gamma_exo_PTHg .* vals.CaSR_PTHg_regulation;
vals.PTHg_exocytosis = vals.F_Ca .* PTHg_vals;

vals.PTHg_degradation = k_PTHg_deg*PTHg_vals;

%PTHp_degradation
vals.PTHp_degradation = k_PTHp_deg.*PTHp_con_vals;

% Calcium
%gut
vals.Gut_impact_D3 = D3p_con_vals.^2./(K_abs_D3.^2 + D3p_con_vals.^2);
vals.Gut_frac_absorption = Gamma_abs0 + (delta_abs_D3).*vals.Gut_impact_D3;
vals.Gut_absorption = ICa.*vals.Gut_frac_absorption;

% resorption
vals.PTHp_res_effect = 0.2*delta_res_max.*PTHp_con_vals.^2./(K_PTHp_res.^2.*ones(size(PTHp_con_vals)) + PTHp_con_vals.^2);
vals.D3p_res_effect   = 0.8*delta_res_max.*D3p_con_vals.^2./(K_D3p_res.^2.*ones(size(D3p_con_vals)) + D3p_con_vals.^2);
vals.Bone_resorption = Gamma_res_min.*ones(size(D3p_con_vals)) + (vals.PTHp_res_effect + vals.D3p_res_effect);

%vals.Bone_accretion = Gamma_ac*ones(size(PTHp_con_vals));
vals.Plasma_to_FastPool = k_pf_Ca*Cap_con_vals;
vals.FastPool_to_Plasma = k_fp_Ca*NCaf_vals;
% renal
vals.Renal_filtration = GFR .* Cap_con_vals;
vals.delta_PT_PTH = delta_PT_max*ones(size(PTHp_con_vals))./(ones(size(PTHp_con_vals))+(PTHp_con_vals./PTHp_ref).^nPT);
vals.Lambda_PT = Lambda_PT0*ones(size(PTHp_con_vals)) + vals.delta_PT_PTH;
vals.delta_TAL_Ca = 0.7*delta_TAL_max*ones(size(Cap_con_vals))./(ones(size(Cap_con_vals)) + (Cap_con_vals./Cap_ref).^nTAL);
vals.delta_TAL_PTH = 0.3*delta_TAL_max*PTHp_con_vals./(PTHp_con_vals + K_TAL_PTHp);
vals.Lambda_TAL = ones(size(PTHp_con_vals))*Lambda_TAL0 + vals.delta_TAL_PTH + vals.delta_TAL_Ca;
vals.delta_DCT_PTH = 0.8*delta_DCT_max.*PTHp_con_vals./(PTHp_con_vals + K_DCT_PTHp*ones(size(PTHp_con_vals)));
vals.delta_DCT_D3 = 0.2*delta_DCT_max.*D3p_con_vals./(D3p_con_vals + K_DCT_D3p*ones(size(D3p_con_vals)));
vals.Lambda_DCT = ones(size(D3p_con_vals))*Lambda_DCT0 + vals.delta_DCT_PTH + vals.delta_DCT_D3;
vals.Renal_frac_reab = min(0.995, vals.Lambda_PT + vals.Lambda_TAL + vals.Lambda_DCT);
vals.Urine_excretion = (1-vals.Renal_frac_reab).* vals.Renal_filtration;


% D3
vals.PTH_impact_D3 = PTHp_con_vals.^nconv./(PTHp_con_vals.^nconv + K_conv_PTH.^nconv);
vals.Ca_impact_D3  = ones(size(Cap_con_vals))./(ones(size(Cap_con_vals)) + gamma_conv_Ca.*Cap_con_vals);
vals.D3_impact_D3  = ones(size(D3p_con_vals))./(ones(size(D3p_con_vals)) + gamma_conv_D3.*D3p_con_vals);
vals.Rconv = k_conv_min + delta_conv_max.*vals.PTH_impact_D3.*vals.Ca_impact_D3.*vals.D3_impact_D3;
vals.D3_synthesis = vals.Rconv .* D3_inact_p;

vals.PTH_impact_D3_degradation = ones(size(PTHp_con_vals))./(1 + gamma_deg_PTHp.*PTHp_con_vals);
vals.D3_degradation = k_deg_D3 .* D3p_con_vals .* vals.PTH_impact_D3_degradation;

% NCaf
vals.Bone_accretion = Gamma_ac*NCaf_vals;
end