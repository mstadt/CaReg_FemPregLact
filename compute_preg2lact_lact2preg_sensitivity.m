% driver to compare lactation and pregnancy parameter effects

% get parameter values
sexORrep = 'preg';
run('read_in_params.m');
pars_preg = params;
clearvars -except pars_preg

sexORrep = 'lact';
run('read_in_params.m');
pars_lact = params;
clearvars -except pars_preg pars_lact param_names

%%% compute pregnancy to lactation sensitivity
fprintf('computing pregnancy to lactation sensitivity \n')
SSfile = './SSbest/16-May-2023_calcium_mod_SS_sexORrep-preg_notes-PTHupdate2.mat ';
ICpreg = load(SSfile).SS;
SSfile = './SSbest/16-May-2023_calcium_mod_SS_sexORrep-lact_notes-PTHupdate2.mat';
IClact = load(SSfile).SS;

preglactdiff = abs(pars_preg - pars_lact);
diffIDs = find(preglactdiff > 1e-12);
disp(size(diffIDs));

preg2lact_sens = zeros(length(diffIDs), 3);  % col 1: PTH, col 2: Ca, col 3: calcitriol
preg2lact_frac = zeros(size(preg2lact_sens));
preg_base = compute_SS(ICpreg,pars_lact, pars_preg, -1);

lact2preg_sens = zeros(size(preg2lact_sens));
lact2preg_frac = zeros(size(preg2lact_sens));
lact_base = compute_SS(IClact,pars_preg, pars_lact, -1);
for ii = 1:length(diffIDs)
    disp(ii)
    parID = diffIDs(ii);
    disp(param_names{parID})
    SSsens = compute_SS(ICpreg, pars_lact, pars_preg, parID);
    preg2lact_sens(ii,:) = SSsens(2:4);
    preg2lact_frac(ii,:) = 100.0*(preg2lact_sens(ii,:) - preg_base(2:4)')./preg_base(2:4)';

    SSsens = compute_SS(IClact, pars_preg, pars_lact, parID);
    lact2preg_sens(ii,:) = SSsens(2:4);
    lact2preg_frac(ii,:) = 100.0*(lact2preg_sens(ii,:) - lact_base(2:4)')./lact_base(2:4)';
end

%%% save results
saveres = input('save results? (0/1)');
if saveres
    notes = input('notes: ');
    fname = strcat('./results_preglact_sensitivity/', date, ...
                '_preg2lactsens_', 'notes-', notes, '.mat');
    save(fname, 'param_names','pars_preg', 'pars_lact',...
                'preg2lact_sens', 'preg2lact_frac', 'preg_base',...
                'lact2preg_sens', 'lact2preg_frac', 'lact_base',...
                'diffIDs');

    fprintf('preg/lact sensitivity analysis results saved to %s \n', fname)
end


%-----------------------------------------------------
function SS = compute_SS(IC, pars1, pars2, IDchange)
    % IC -- initial condition for preg or lact model
    % pars1 -- baseline parameters, this set is fixed
    % pars2 -- parameters that will be change
    % IDchange    -- ID of parameter to set to pars1 value, IDchange = -1
    %                       means no change
    if IDchange > -1
        pars2(IDchange) = pars1(IDchange);
    end
    tspan = [0 4000];
    options = odeset('RelTol',1.0e-6,'AbsTol',1e-9);
    [~, y] = ode15s(@(t,y) calcium_mod(t,y,pars2, ...
                                'NCas_fixed', true), ...
                       tspan, IC, options);
    IG = y(end,:)';
    
    options = optimoptions('fsolve', 'Display', 'off');
    [SS, ~, exitflag, ~] = fsolve(@(y) calcium_mod(0,y,pars2),...
                                                             IG, options);

    Vp = pars2(18);
    SS(2:4) = SS(2:4)/Vp;
    if exitflag<1
        fprintf('***exitflag indicates error!!*** \n')
    end
end