% driver to compare impact of changing individual preg/lactation parameters
% to respective virgin values
% NOTE: this index sets it up so all parameters are considered

% get parameter values
sexORrep = 'female';
run('read_in_params');
pars_female = params;
clearvars -except pars_female

sexORrep = 'preg';
run('read_in_params');
pars_preg = params;
clearvars -except pars_female pars_preg

sexORrep = 'lact';
run('read_in_params');
pars_lact = params;
clearvars -except pars_female pars_preg pars_lact param_names

notes = input('notes: ');

%%% compute pregnancy
fprintf('computing virgin to pregnancy sensitivity \n')
% baseline pregnancy SS results
SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-female_notes-newfemale.mat';
IC = load(SSfile).SS;

pregdiff = abs(pars_preg - pars_female);
diffIDs_preg = find(pregdiff > 1e-12);
disp(size(diffIDs_preg))

fem2preg_sens = zeros(length(pars_preg), 3);  % col 1: PTH, col 2: Ca, col 3: calcitriol
fem2preg_frac = zeros(size(fem2preg_sens));% col 1: PTH, col 2: Ca, col 3: calcitriol
female_base = compute_SS(IC,pars_female, pars_preg, -1);
for ii = 1:50
    disp(ii)
    disp(param_names{ii})
    if ismember(ii, diffIDs_preg)
        SSsens = compute_SS(IC, pars_female, pars_preg, ii);
        fem2preg_sens(ii,:) = SSsens(2:4);
        fem2preg_frac(ii,:) = 100.0*(fem2preg_sens(ii,:) - female_base(2:4)')./female_base(2:4)';
    end
end


%%% compute lactation
fprintf('computing virgin to lactation sensitivity \n')
% baseline lactation SS results
SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-female_notes-newfemale.mat';
IC = load(SSfile).SS;

lactdiff = abs(pars_lact - pars_female);
diffIDs_lact = find(lactdiff > 1e-12);  
disp(size(diffIDs_lact))
fem2lact_sens = zeros(length(pars_lact), 3);  % col 1: PTH, col 2: Ca, col 3: calcitriol
fem2lact_frac = zeros(size(fem2lact_sens));
%female_base = compute_SS(IC,pars_female, pars_lact, -1);
for ii = 1:50
    disp(ii)
    disp(param_names{ii})
    if ismember(ii, diffIDs_lact)
        SSsens = compute_SS(IC, pars_female, pars_lact, ii);
        fem2lact_sens(ii,:) = SSsens(2:4);
        fem2lact_frac(ii,:) = 100.0*(fem2lact_sens(ii,:) - female_base(2:4)')./female_base(2:4)';
    end
end


%%% save results
saveres = 1; %input('save results? (0/1)');
if saveres
    fname = strcat('./results_virgin2preglact_sensitivity/', date, ...
                '_vir2preg_vir2lact_all_', 'notes-', notes, '.mat');
    save(fname, 'param_names', 'pars_female', 'pars_preg', 'pars_lact',...
                'fem2preg_sens', 'fem2preg_frac', 'female_base',...
                'fem2lact_sens', 'fem2lact_frac', 'female_base',...
                'diffIDs_preg', 'diffIDs_lact');

    fprintf('virgin to preg/lact sensitivity analysis results saved to %s \n', fname)
end
%-----------------------------------------------------
function SS = compute_SS(IC, pars_female, pars_pregORlact, IDchange)
    % IC -- initial condition for preg or lact model
    % pars_female -- baseline female parameters
    % pars_pregORlact -- baseline preg or lact parameters
    % IDchange    -- ID of parameter to set to female value, IDchange = -1
    %                       means no change
    if IDchange > -1
        pars_female(IDchange) = pars_pregORlact(IDchange);
        %pars_pregORlact(IDchange) = pars_female(IDchange);
    end
    tspan = [0 4000];
    options = odeset('RelTol',1.0e-6,'AbsTol',1e-9);
    [~, y] = ode15s(@(t,y) calcium_mod(t,y,pars_female, ...
                            'NCas_fixed', true), ...
                   tspan, IC, options);
    IG = y(end,:)';
    
    options = optimoptions('fsolve', 'Display', 'off');
    [SS, ~, exitflag, ~] = fsolve(@(y) calcium_mod(0,y,pars_female),...
                                                             IG, options);
    Vp = pars_female(18);
    SS(2:4) = SS(2:4)/Vp;
    if exitflag<1
        fprintf('***exitflag indicates error!!*** \n')
    end
end

