% This is used to see how pregnancy/lactation parameters effect male model

% compute male model SS with preg parameter changes
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

% set preg/lact male
par_malepreg = pars_male.*preg2vir;
par_malepreg(37) = pars_preg(37); % set FetusORMilk

par_malelact = pars_male.*lact2vir;
par_malelact(37) = pars_lact(37); % set FetusORMilk

% compute preg/lact male

% set initial conditions
% use male IC
SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-male_notes-newmale.mat';
SS_IG = load(SSfile).SS;
Vp_malepreg = par_malepreg(18); % set volume for concentrations

tspan = [0 5000];
options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9);
fprintf('solving preg male model eqns \n')
IC = SS_IG;

% variables to fix
D3_fixed   = 0;
Ca_fixed   = 0;
NCas_fixed = 1;
NCaf_fixed = 0;
PTHg_fixed = 0;
PTHp_fixed = 0;
altsim = 0;

[t, y] = ode15s(@(t,y) calcium_mod(t,y,par_malepreg, ...
                    'Ca_fixed', Ca_fixed,...
                    'D3_fixed', D3_fixed,...
                    'PTHg_fixed', PTHg_fixed,...
                    'PTHp_fixed', PTHp_fixed,...
                    'NCaf_fixed', NCaf_fixed,...
                    'NCas_fixed', NCas_fixed, ...
                    'alt_sim', altsim), ...
                    tspan, IC, options);

figure(2)
clf
num_rows = 2; num_cols = 3;
lw = 2; c1 = [0.3010 0.7450 0.9330];
c_gray = uint8([70 78 81]);
subplot(num_rows,num_cols,1)
plot(t, y(:,1),'linewidth', lw, 'color', c1)
ylabel('PTH_g (pmol)')
xlabel('t')
title('Parathyroid gland PTH pool')
yline(0,'color',c_gray,'linewidth',lw)
grid on

subplot(num_rows,num_cols,2)
plot(t, y(:,2)/Vp_malepreg,'linewidth', lw, 'color', c1)
ylabel('[PTH_p] (pmol/L)')
xlabel('t')
title('Plasma PTH concentration')
yline(1.5, 'color', c_gray, 'linewidth', lw) % min of normal range
yline(13, 'color', c_gray,'linewidth', lw) % max of normal range
grid on

subplot(num_rows,num_cols,3)
plot(t, y(:,3)/Vp_malepreg,'linewidth', lw, 'color', c1)
ylabel('[Ca_p] (mmol/L)')
xlabel('t')
title('Plasma calcium concentration')
yline(1.1, 'color', c_gray,'linewidth', lw) % min of normal range
yline(1.3, 'color', c_gray,'linewidth', lw) % max of normal range
grid on

subplot(num_rows,num_cols,4)
plot(t, y(:,4)/Vp_malepreg, 'linewidth', lw, 'color', c1)
ylabel('[D3_p] (pmol/L)')
xlabel('t')
yline(80*0.6, 'color', c_gray,'linewidth', lw) % min of normal range
yline(250*0.6, 'color', c_gray,'linewidth', lw) % max of normal range
title('Plasma vitamin D3 concentration')
grid on

subplot(num_rows,num_cols,5)
plot(t,y(:,5),'linewidth',lw,'color',c1)
ylabel('NCa_f')
xlabel('t')
title('Fast bone pool calcium')
yline(0,'color',c_gray,'linewidth',lw)
grid on

subplot(num_rows,num_cols,6)
plot(t,y(:,6),'linewidth',lw, 'color', c1)
yline(0, 'color',c_gray,'linewidth',lw)
ylabel('NCa_s')
xlabel('Slow bone pool calcium')
grid on

sgtitle('ODE trajectory for preg male model')
SS_IG = y(end,:)';

fprintf('fsolve for preg male SS \n')
options = optimoptions('fsolve', 'Display', 'final', 'MaxFunEvals', 10000, 'MaxIter', 10000);
[SS, residual, exitflag, output] = fsolve(@(y) calcium_mod(0,y,par_malepreg,...
                                                         'Ca_fixed', Ca_fixed,...
                                                         'D3_fixed', D3_fixed,...
                                                         'PTHg_fixed', PTHg_fixed,...
                                                         'PTHp_fixed', PTHp_fixed,...
                                                         'NCaf_fixed', NCaf_fixed,...
                                                         'NCas_fixed', NCas_fixed),...
                                                         SS_IG, options);
% check between ODE and fsolve result
consODE = SS_IG(2:4)/Vp_malepreg;
consfsolve = SS(2:4)/Vp_malepreg;
fracchange = (consODE - consfsolve)./consODE;
if max(abs(fracchange)) > 0.1
    fprintf('maximum ODE to fsolve change: %0.3f \n', max(abs(fracchange)))
    

    fprintf('ODE steady states \n')
    fprintf('           SS_IG\n')
    fprintf('PTHg:      %0.3f\n', SS_IG(1))
    fprintf('PTHp_con:  %0.3f\n', SS_IG(2)/Vp_malepreg)
    fprintf('Cap_con:   %0.3f\n', SS_IG(3)/Vp_malepreg)
    fprintf('D3p_con:   %0.3f\n', SS_IG(4)/Vp_malepreg)
    fprintf('NCaf:      %0.3f\n', SS_IG(5))
    fprintf('NCas:      %0.3f\n', SS_IG(6))

    fprintf('**WARNING: ODE to fsolve change by more than 10 percent *** \n')
end

if exitflag<1
    fprintf('***exitflag indicates error!!*** \n')
end

fprintf('final steady states for preg male \n')
fprintf('           SS\n')
fprintf('PTHg:      %0.3f\n', SS(1))
fprintf('PTHp_con:  %0.3f\n', SS(2)/Vp_malepreg)
fprintf('Cap_con:   %0.3f\n', SS(3)/Vp_malepreg)
fprintf('D3p_con:   %0.3f\n', SS(4)/Vp_malepreg)
fprintf('NCaf:      %0.3f\n', SS(5))
fprintf('NCas:      %0.3f\n', SS(6))

valsSS = compute_vars(SS', par_malepreg);
save_SS = 1;
if save_SS
    notes = input('notes: ');
    fname_save = strcat('./SS/',date,'_calcium_mod_SS', '_sexORrep-malepreg_notes-',notes,'.mat');
    save(fname_save, 'SS', 'valsSS','par_malepreg')

    fprintf('male preg results saved to: \n %s \n', fname_save);
end

%% compute lact
% set initial conditions
% use male IC
SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-male_notes-newmale.mat';
SS_IG = load(SSfile).SS;
Vp_malelact = par_malelact(18); % set volume for concentrations

tspan = [0 5000];
options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9);
fprintf('solving lact male model eqns \n')
IC = SS_IG;

% variables to fix
D3_fixed   = 0;
Ca_fixed   = 0;
NCas_fixed = 1;
NCaf_fixed = 0;
PTHg_fixed = 0;
PTHp_fixed = 0;
altsim = 0;

[t, y] = ode15s(@(t,y) calcium_mod(t,y,par_malelact, ...
                    'Ca_fixed', Ca_fixed,...
                    'D3_fixed', D3_fixed,...
                    'PTHg_fixed', PTHg_fixed,...
                    'PTHp_fixed', PTHp_fixed,...
                    'NCaf_fixed', NCaf_fixed,...
                    'NCas_fixed', NCas_fixed, ...
                    'alt_sim', altsim), ...
                    tspan, IC, options);

figure(3)
clf
num_rows = 2; num_cols = 3;
lw = 2; c1 = [0.3010 0.7450 0.9330];
c_gray = uint8([70 78 81]);
subplot(num_rows,num_cols,1)
plot(t, y(:,1),'linewidth', lw, 'color', c1)
ylabel('PTH_g (pmol)')
xlabel('t')
title('Parathyroid gland PTH pool')
yline(0,'color',c_gray,'linewidth',lw)
grid on

subplot(num_rows,num_cols,2)
plot(t, y(:,2)/Vp_malelact,'linewidth', lw, 'color', c1)
ylabel('[PTH_p] (pmol/L)')
xlabel('t')
title('Plasma PTH concentration')
yline(1.5, 'color', c_gray, 'linewidth', lw) % min of normal range
yline(13, 'color', c_gray,'linewidth', lw) % max of normal range
grid on

subplot(num_rows,num_cols,3)
plot(t, y(:,3)/Vp_malelact,'linewidth', lw, 'color', c1)
ylabel('[Ca_p] (mmol/L)')
xlabel('t')
title('Plasma calcium concentration')
yline(1.1, 'color', c_gray,'linewidth', lw) % min of normal range
yline(1.3, 'color', c_gray,'linewidth', lw) % max of normal range
grid on

subplot(num_rows,num_cols,4)
plot(t, y(:,4)/Vp_malelact, 'linewidth', lw, 'color', c1)
ylabel('[D3_p] (pmol/L)')
xlabel('t')
yline(80*0.6, 'color', c_gray,'linewidth', lw) % min of normal range
yline(250*0.6, 'color', c_gray,'linewidth', lw) % max of normal range
title('Plasma vitamin D3 concentration')
grid on

subplot(num_rows,num_cols,5)
plot(t,y(:,5),'linewidth',lw,'color',c1)
ylabel('NCa_f')
xlabel('t')
title('Fast bone pool calcium')
yline(0,'color',c_gray,'linewidth',lw)
grid on

subplot(num_rows,num_cols,6)
plot(t,y(:,6),'linewidth',lw, 'color', c1)
yline(0, 'color',c_gray,'linewidth',lw)
ylabel('NCa_s')
xlabel('Slow bone pool calcium')
grid on

sgtitle('ODE trajectory for lact male model')
SS_IG = y(end,:)';

fprintf('fsolve for lact male SS \n')
options = optimoptions('fsolve', 'Display', 'final', 'MaxFunEvals', 10000, 'MaxIter', 10000);
[SS, residual, exitflag, output] = fsolve(@(y) calcium_mod(0,y,par_malelact,...
                                                         'Ca_fixed', Ca_fixed,...
                                                         'D3_fixed', D3_fixed,...
                                                         'PTHg_fixed', PTHg_fixed,...
                                                         'PTHp_fixed', PTHp_fixed,...
                                                         'NCaf_fixed', NCaf_fixed,...
                                                         'NCas_fixed', NCas_fixed),...
                                                         SS_IG, options);
% check between ODE and fsolve result
consODE = SS_IG(2:4)/Vp_malelact;
consfsolve = SS(2:4)/Vp_malelact;
fracchange = (consODE - consfsolve)./consODE;
if max(abs(fracchange)) > 0.1
    fprintf('maximum ODE to fsolve change: %0.3f \n', max(abs(fracchange)))
    

    fprintf('ODE steady states \n')
    fprintf('           SS_IG\n')
    fprintf('PTHg:      %0.3f\n', SS_IG(1))
    fprintf('PTHp_con:  %0.3f\n', SS_IG(2)/Vp_malelact)
    fprintf('Cap_con:   %0.3f\n', SS_IG(3)/Vp_malelact)
    fprintf('D3p_con:   %0.3f\n', SS_IG(4)/Vp_malelact)
    fprintf('NCaf:      %0.3f\n', SS_IG(5))
    fprintf('NCas:      %0.3f\n', SS_IG(6))

    fprintf('**WARNING: ODE to fsolve change by more than 10 percent *** \n')
end

if exitflag<1
    fprintf('***exitflag indicates error!!*** \n')
end

fprintf('final steady states for lact male \n')
fprintf('           SS\n')
fprintf('PTHg:      %0.3f\n', SS(1))
fprintf('PTHp_con:  %0.3f\n', SS(2)/Vp_malelact)
fprintf('Cap_con:   %0.3f\n', SS(3)/Vp_malelact)
fprintf('D3p_con:   %0.3f\n', SS(4)/Vp_malelact)
fprintf('NCaf:      %0.3f\n', SS(5))
fprintf('NCas:      %0.3f\n', SS(6))

valsSS = compute_vars(SS', par_malelact);
save_SS = 1;
if save_SS
    notes = input('notes: ');
    fname_save = strcat('./SS/',date,'_calcium_mod_SS', '_sexORrep-malelact_notes-',notes,'.mat');
    save(fname_save, 'SS', 'valsSS','par_malelact')

    fprintf('male lact results saved to: \n %s \n', fname_save);
end


%% lower calcitriol synthesis
par_malepreg()
