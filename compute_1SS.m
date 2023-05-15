clear all;

sexORrep = 'preg'; % options: male, female, preg, lact
altsim = 0;

%----------------
% User input
%----------------
do_ODE = 1; % start with solving ODE model for SS initial guess

dofig = 1; % show ODE fig

% variables to fix
D3_fixed   = 0;
Ca_fixed   = 0;
NCas_fixed = 1;
NCaf_fixed = 0;
PTHg_fixed = 0;
PTHp_fixed = 0;


% initialize parameter values
fprintf('loading %s params \n', sexORrep)
run('read_in_params.m')

% % change parameters here!
%GamFetORMilk = 0.389e-3; params(37) = GamFetORMilk;
%ICa = 1.05*1.84e-3; params(25) = ICa;
%GamAbs0 = 0.55; params(26) = GamAbs0;
%deltaabs = 0.425; params(27) = deltaabs;
%KabsD3 = 265; params(28) = KabsD3;
%--------------------
% End of user input
%--------------------
% %% compute steady state

% set initial conditions
if strcmp(sexORrep,'male')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-male_notes-newmale.mat';
    SS_IG = load(SSfile).SS;
    Vp = Vp_male;
elseif strcmp(sexORrep, 'female')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-female_notes-newfemale.mat';
    SS_IG = load(SSfile).SS;
    Vp = Vp_female;
elseif strcmp(sexORrep, 'preg')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-preg_notes-PTHadjust_preg.mat';
    SS_IG = load(SSfile).SS;
    Vp = Vp_preg;
elseif strcmp(sexORrep, 'lact')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-lact_notes-PTHadjustlact.mat';
    SS_IG = load(SSfile).SS;
    Vp = Vp_lact;
end

%%% find IG for steady state from ODE system
if do_ODE
    tspan = [0 4000];
    
    options = odeset('RelTol',1.0e-6,'AbsTol',1e-9);
    
    fprintf('solving model eqns \n')

    IC = SS_IG;
    
    [t, y] = ode15s(@(t,y) calcium_mod(t,y,params, ...
                        'Ca_fixed', Ca_fixed,...
                        'D3_fixed', D3_fixed,...
                        'PTHg_fixed', PTHg_fixed,...
                        'PTHp_fixed', PTHp_fixed,...
                        'NCaf_fixed', NCaf_fixed,...
                        'NCas_fixed', NCas_fixed, ...
                        'alt_sim', altsim), ...
                       tspan, IC, options);

    if dofig
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
    plot(t, y(:,2)/Vp,'linewidth', lw, 'color', c1)
    ylabel('[PTH_p] (pmol/L)')
    xlabel('t')
    title('Plasma PTH concentration')
    yline(1.5, 'color', c_gray, 'linewidth', lw) % min of normal range
    yline(13, 'color', c_gray,'linewidth', lw) % max of normal range
    grid on
    
    subplot(num_rows,num_cols,3)
    plot(t, y(:,3)/Vp,'linewidth', lw, 'color', c1)
    ylabel('[Ca_p] (mmol/L)')
    xlabel('t')
    title('Plasma calcium concentration')
    yline(1.1, 'color', c_gray,'linewidth', lw) % min of normal range
    yline(1.3, 'color', c_gray,'linewidth', lw) % max of normal range
    grid on
    
    subplot(num_rows,num_cols,4)
    plot(t, y(:,4)/Vp, 'linewidth', lw, 'color', c1)
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
    
    sgtitle('ODE trajectory for model')
    end
    SS_IG = y(end,:)';
end

fprintf('fsolve for %s SS \n', sexORrep)
options = optimoptions('fsolve', 'Display', 'final', 'MaxFunEvals', 10000, 'MaxIter', 10000);
[SS, residual, exitflag, output] = fsolve(@(y) calcium_mod(0,y,params,...
                                                         'Ca_fixed', Ca_fixed,...
                                                         'D3_fixed', D3_fixed,...
                                                         'PTHg_fixed', PTHg_fixed,...
                                                         'PTHp_fixed', PTHp_fixed,...
                                                         'NCaf_fixed', NCaf_fixed,...
                                                         'NCas_fixed', NCas_fixed),...
                                                         SS_IG, options);

% check between ODE and fsolve result
consODE = SS_IG(2:4)/Vp;
consfsolve = SS(2:4)/Vp;
fracchange = (consODE - consfsolve)./consODE;
if max(abs(fracchange)) > 0.1
    fprintf('maximum ODE to fsolve change: %0.3f \n', max(abs(fracchange)))
    

    fprintf('ODE steady states \n')
    fprintf('           SS_IG\n')
    fprintf('PTHg:      %0.3f\n', SS_IG(1))
    fprintf('PTHp_con:  %0.3f\n', SS_IG(2)/Vp)
    fprintf('Cap_con:   %0.3f\n', SS_IG(3)/Vp)
    fprintf('D3p_con:   %0.3f\n', SS_IG(4)/Vp)
    fprintf('NCaf:      %0.3f\n', SS_IG(5))
    fprintf('NCas:      %0.3f\n', SS_IG(6))

    fprintf('**WARNING: ODE to fsolve change by more than 10 percent *** \n')
end

if exitflag<1
    fprintf('***exitflag indicates error!!*** \n')
end

fprintf('final steady states \n')
fprintf('           SS\n')
fprintf('PTHg:      %0.3f\n', SS(1))
fprintf('PTHp_con:  %0.3f\n', SS(2)/Vp)
fprintf('Cap_con:   %0.3f\n', SS(3)/Vp)
fprintf('D3p_con:   %0.3f\n', SS(4)/Vp)
fprintf('NCaf:      %0.3f\n', SS(5))
fprintf('NCas:      %0.3f\n', SS(6))

valsSS = compute_vars(SS', params);
save_SS = input('save SS? (0/1) ');
if save_SS
    notes = input('notes: ');
    fname_save = strcat('./SS/',date,'_calcium_mod_SS', '_sexORrep-', sexORrep, '_notes-',notes,'.mat');
    save(fname_save, 'SS', 'valsSS','params','sexORrep',...
                        'exitflag', 'residual', 'param_names')

    fprintf('results saved to %s \n', fname_save);
end
