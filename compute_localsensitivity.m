% computes local sensitivity analysis of all parameters for given sexORrep
% variable
clear all
sexORrep = 'lact';

if strcmp(sexORrep,'male')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-male_notes-newmale.mat';
elseif strcmp(sexORrep, 'female')
    SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-female_notes-newfemale.mat';
elseif strcmp(sexORrep, 'preg')
    SSfile = './SSbest/15-May-2023_calcium_mod_SS_sexORrep-preg_notes-FetORMilk.mat';
elseif strcmp(sexORrep, 'lact')
    SSfile = './SSbest/15-May-2023_calcium_mod_SS_sexORrep-lact_notes-FetORMilkupdate_lact.mat';
end

save_res = 1; %input('save? (0/1) ');
notes = input('notes: ');

dat = load(SSfile);
SS_IG = dat.SS;
pars = dat.params;
param_names = dat.param_names;


SSbase = compute_SS(pars, SS_IG, -1);


sens = zeros(length(pars)-2, 3); % col 1: PTH, col 2: Ca, col 3: calcitriol
frac_change = zeros(size(sens));
basevals = SSbase(2:4)'; % used to compute fractional changes

for ii = 1:50
    disp(ii)
    disp(param_names{ii})
    SSsens = compute_SS(pars, SS_IG, ii);
    sens(ii, :) = SSsens(2:4);
    frac_change(ii,:) = 100.0*(sens(ii,:) - basevals)./basevals;
end

%% save results
if save_res
    %notes = input('notes: ');
    fname = strcat('./results_localsens/',date,'_localsens_',...
                            'sexORrep-', sexORrep, '_notes-', notes, '.mat');
    save(fname, 'pars', 'sexORrep', 'frac_change', 'SSbase', 'sens', 'param_names');
    fprintf('sensitivity analysis results saved to %s \n', fname)
end


%%%
function SS = compute_SS(pars, IC, parchange_ID)
    % pars - par values
    % IG - initial guess
    % parchange_ID - par ID to change, set to -1 for normal SS
    params = pars;
    if parchange_ID > -1 
        params(parchange_ID) = 1.05*pars(parchange_ID);
    %params(parchange_ID) = 1.1*pars(parchange_ID);
    end
    
    tspan = [0 4000];
    options = odeset('RelTol',1.0e-6,'AbsTol',1e-9);
    [~, y] = ode15s(@(t,y) calcium_mod(t,y,params, ...
                                'NCas_fixed', true), ...
                       tspan, IC, options);
    IG = y(end,:)';
    
    options = optimoptions('fsolve', 'Display', 'off');
    [SS, ~, exitflag, ~] = fsolve(@(y) calcium_mod(0,y,params),...
                                                             IG, options);

    Vp = params(18);
    SS(2:4) = SS(2:4)/Vp; % change to concentration
    if exitflag<1
        fprintf('***exitflag indicates error!!*** \n')
    end
end