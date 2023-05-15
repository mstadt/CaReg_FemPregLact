% predicted impact of hyperparathyroidism
clear all;
sexORrep = 'preg';

% set up baseline parameters
run('read_in_params.m')
base_k_prod_PTHg = params(40);



% PTH synthesis changes
PTHchange_vals = 1:1:100; % multiply
notes = input('notes for fname: ');


% set initial conditions
%%% set ICs
if strcmp(sexORrep, 'male')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-male_notes-newmale.mat';
elseif strcmp(sexORrep, 'female')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-female_notes-newfemale.mat';
elseif strcmp(sexORrep, 'preg')
    SSfile = './SSbest/15-May-2023_calcium_mod_SS_sexORrep-preg_notes-FetORMilk.mat';
elseif strcmp(sexORrep, 'lact')
    SSfile = './SSbest/15-May-2023_calcium_mod_SS_sexORrep-lact_notes-FetORMilkupdate_lact.mat';
end
IG = load(SSfile).SS;

NCas_fixed = true;

options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9);
opts_fsolve = optimoptions('fsolve', 'Display', 'off', ...
                                'MaxFunEvals', 5000, 'MaxIter', 2000);

for ii = 1:length(PTHchange_vals)
    PTHchange = PTHchange_vals(ii);
    k_prod_PTHg = base_k_prod_PTHg * PTHchange;
    params(40) = k_prod_PTHg;
    fprintf('k_prod_PTHg change: %0.1f \n', PTHchange)

    % compute SS vals
    fprintf('solving ODE \n')
    tspan = [0 3000];

    [t, y] = ode15s(@(t,y) calcium_mod(t,y,params,...
                        'NCas_fixed', NCas_fixed),...
                        tspan, IG, options);

    IG = y(end,:)';
    fprintf('solving fsolve \n')
    
    [SS, fval, exitflag, ~] = fsolve(@(y) calcium_mod(0,y,params,...
                                            'NCas_fixed', NCas_fixed),...
                                            IG, opts_fsolve);
    if exitflag < 1
        fprintf('**** %s exitflag indicates error! **** \n', sexORrep)
        Vp = params(18);
        figure(1)
        clf
         num_rows = 2; num_cols = 3;
        lw = 3; c1 = [0.3010 0.7450 0.9330];
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
        title('Fast bone pool clacium')
        yline(0,'color',c_gray,'linewidth',lw)
        grid on
        
        subplot(num_rows,num_cols,6)
        plot(t,y(:,6),'linewidth',lw, 'color', c1)
        yline(0, 'color',c_gray,'linewidth',lw)
        ylabel('NCa_s')
        xlabel('Slow bone pool calcium')
        grid on
        
        temp = strcat('ODE trajectory for ',sexORrep, ' model, D3change = ', D3change);
        sgtitle(temp)
        
        pause(3)
        max(fval)
        cont = input('continue? (0/1) ');
        
        if ~cont
            break
        end
    end % exitflag check

    valsSS = compute_vars(SS', params);

    % export to .mat files for info
    fname_save = strcat('./results_hyperPTH/', date, '_hyperparathyroid_', 'sexORrep-',sexORrep,...
                            '_PTHchange-', num2str(PTHchange), '_notes-', notes, '.mat');
    save(fname_save, 'SS', 'params', 'sexORrep', 'valsSS',...
            'PTHchange', 'k_prod_PTHg')

end