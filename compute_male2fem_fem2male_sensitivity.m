% driver to compare lactation and pregnancy parameter effects

% get parameter values
sexORrep = 'male';
run('read_in_params.m');
pars_male = params;
clearvars -except pars_male

sexORrep = 'female';
run('read_in_params.m');
pars_female = params;
clearvars -except param_names pars_male pars_female

saveres = 1; %input('save results? (0/1)');
notes = input('notes: ');

%%% compute pregnancy to lactation sensitivity
fprintf('computing male to female sensitivity \n')
SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-male_notes-newmale.mat';
ICmale = load(SSfile).SS;
SSfile = './SSbest/12-May-2023_calcium_mod_SS_sexORrep-female_notes-newfemale.mat';
ICfemale = load(SSfile).SS;

mfdiff = abs(pars_male - pars_female);
diffIDs = find(mfdiff > 1e-12);
disp(size(diffIDs));

male2female_sens = zeros(length(diffIDs), 3);  % col 1: PTH, col 2: Ca, col 3: calcitriol
male2female_frac = zeros(size(male2female_sens));
male_base = compute_SS(ICmale,pars_female, pars_male, -1);

female2male_sens = zeros(size(male2female_sens));
female2male_frac = zeros(size(male2female_sens));
female_base = compute_SS(ICfemale,pars_male, pars_female, -1);
for ii = 1:length(diffIDs)
    disp(ii)
    parID = diffIDs(ii);
    disp(param_names{parID})
    SSsens = compute_SS(ICmale, pars_female, pars_male, parID);
    male2female_sens(ii,:) = SSsens(2:4);
    male2female_frac(ii,:) = 100.0*(male2female_sens(ii,:) - male_base(2:4)')./male_base(2:4)';

    SSsens = compute_SS(ICfemale, pars_male, pars_female, parID);
    female2male_sens(ii,:) = SSsens(2:4);
    female2male_frac(ii,:) = 100.0*(female2male_sens(ii,:) - female_base(2:4)')./female_base(2:4)';
end

%%% save results

if saveres
    fname = strcat('./results_male2fem_fem2male/', date, ...
                '_male2femalesens_', 'notes-', notes, '.mat');
    save(fname, 'param_names','pars_male', 'pars_female',...
                'male2female_sens', 'male2female_frac', 'male_base',...
                'female2male_sens', 'female2male_frac', 'female_base',...
                'diffIDs');

    fprintf('male2female & female2male sensitivity analysis results saved to %s \n', fname)
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
    SS(2:4) = SS(2:4)/Vp; % change to concentration
    if exitflag<1
        fprintf('***exitflag indicates error!!*** \n')
    end
end