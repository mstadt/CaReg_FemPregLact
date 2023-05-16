% computes local sensitivity analysis of all parameters for given sexORrep
% variable
clear all
sexORrep = 'lact';

if strcmp(sexORrep,'male')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-male_notes-newmale.mat';
elseif strcmp(sexORrep, 'female')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-female_notes-newfemale.mat';
elseif strcmp(sexORrep, 'preg')
    SSfile = './SSbest/16-May-2023_calcium_mod_SS_sexORrep-preg_notes-PTHupdate2.mat';
elseif strcmp(sexORrep, 'lact')
    SSfile = './SSbest/16-May-2023_calcium_mod_SS_sexORrep-lact_notes-PTHupdate2.mat';
end

save_res = 1; %input('save? (0/1) ');
notes = input('notes: ');

dat = load(SSfile);
SS_IG = dat.SS;
pars = dat.params;
param_names = dat.param_names;


SSbase = compute_SS(pars, SS_IG, -1);


sens_plus = zeros(length(pars)-2, 3); % col 1: PTH, col 2: Ca, col 3: calcitriol
sens_minus = zeros(size(sens_plus));
frac_change_plus = zeros(size(sens_plus));
frac_change_minus = zeros(size(sens_minus));
basevals = SSbase(2:4)'; % used to compute fractional changes
frac_change_plusminus = zeros(size(sens_plus));

for ii = 1:50
    disp(ii)
    disp(param_names{ii})
    [SSplus, SSminus] = compute_SS(pars, SS_IG, ii);
    sens_plus(ii, :) = SSplus(2:4);
    sens_minus(ii,:) = SSminus(2:4);
    frac_change_plus(ii,:) = 100.0*(sens_plus(ii,:) - basevals)./basevals;
    frac_change_minus(ii,:) = 100.0*(sens_minus(ii,:) - basevals)./basevals;

    frac_change_plusminus(ii,:) = 100.0*(sens_plus(ii,:) - sens_minus(ii,:))./basevals;
end

%% save results
if save_res
    %notes = input('notes: ');
    fname = strcat('./results_localsens/',date,'_localsens_plusminus',...
                            'sexORrep-', sexORrep, '_notes-', notes, '.mat');
    save(fname, 'pars', 'sexORrep', 'frac_change_plus', 'sens_plus',...
        'sens_minus', 'frac_change_minus', 'basevals', ...
        'frac_change_plusminus',...
        'SSbase', 'param_names');
    fprintf('sensitivity analysis results saved to: \n %s \n', fname)
end


%%%
function [SS_plus, SS_minus] = compute_SS(pars, IC, parchange_ID)
    % pars - par values
    % IG - initial guess
    % parchange_ID - par ID to change, set to -1 for normal SS
    params_plus = pars;

    % compute SS with increase by 5%
    if parchange_ID > -1 
        params_plus(parchange_ID) = 1.05*pars(parchange_ID);
    end
    
    tspan = [0 4000];
    options = odeset('RelTol',1.0e-6,'AbsTol',1e-9);
    [~, y_plus] = ode15s(@(t,y) calcium_mod(t,y,params_plus, ...
                                'NCas_fixed', true), ...
                       tspan, IC, options);
    IG1 = y_plus(end,:)';
    
    options = optimoptions('fsolve', 'Display', 'off');
    [SS_plus, ~, exitflag, ~] = fsolve(@(y) calcium_mod(0,y,params_plus),...
                                                             IG1, options);

    Vp = params_plus(18);
    SS_plus(2:4) = SS_plus(2:4)/Vp; % change to concentration
    if exitflag<1
        fprintf('***exitflag indicates error!!*** \n')
    end


    % compute SS with decrease by 5%
    params_minus = pars;
    if parchange_ID > -1
        params_minus(parchange_ID) = 0.95*pars(parchange_ID);
    end

    tspan = [0 4000];
    options = odeset('RelTol',1.0e-6,'AbsTol',1e-9);
    [~, y_minus] = ode15s(@(t,y) calcium_mod(t,y,params_minus,...
                                        'NCas_fixed', true), ...
                            tspan, IC, options);
    IG2 = y_minus(end,:)';

    options = optimoptions('fsolve', 'Display', 'off');
    [SS_minus, ~, exitflag, ~] = fsolve(@(y) calcium_mod(0,y,params_minus),...
                                        IG2, options);

    Vp = params_minus(18);
    SS_minus(2:4) = SS_minus(2:4)/Vp;
    if exitflag<1
        fprintf('***exitplag indicates error!!*** \n')
    end
    
end
