% need to set sexORpreg variable (options: male, female, preg, lact)

new_pars = 0; % are there new parameters

if strcmp(sexORrep, 'male')
    param_file = 'params_male.txt';
elseif strcmp(sexORrep,'female')
    param_file = 'params_female.txt';
elseif strcmp(sexORrep, 'preg')
    param_file = 'params_pregnancy.txt';
elseif strcmp(sexORrep, 'lact')
    param_file = 'params_lactation.txt';
else
    fprintf('sexORrep: %s \n', sex)
    error('sexORrep is not done')
end
num_check1 = 50;
temp = readcell(param_file);
pars_list1 = temp(:,1);
numpars1 = size(pars_list1,1);
if numpars1~=num_check1
    fprintf('expected: %i, numpars: %i \n', num_check1, numpars1)
    error('sexORrep specific params numpars incorrect')
end

param_names = {};
numpars = numpars1;
params = zeros(1,numpars);
for ii = 1:numpars
    temp = pars_list1{ii};
    parval = temp(find(~isspace(temp)));
    param_names{ii} = extractBefore(parval, '=');
    param_vals{ii}  = extractAfter(parval, '=');
    pname = strcat(param_names{ii},'_', sexORrep);
    assignin('base', pname ,eval(param_vals{ii}))
    params(ii) = eval(param_vals{ii});
end

% optional: easier for writing up model equations and checking params
if new_pars
    fid = fopen('newparamnames.txt','wt');
    for jj = 1:numpars
        fprintf(fid, '%s = params(%d);\n',param_names{jj},jj);
    end
    fprintf('TO DO: copy results from newparamnames.txt into parameters set up in model function \n')
end
fclose('all');
