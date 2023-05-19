clear all
%%% Pregnancy

CaPerDay_mg_perpup = 5

%%%%
npups = 9;
mgCaperday = CaPerDay_mg_perpup*npups; 


gperday = mgCaperday * 1e-3;
molperday = gperday / (2 * 40.078);
molpermin = molperday / (60 * 24);
mmolpermin = molpermin * 1e3;

Preg_GamFetusORMilk = mmolpermin


%% Lactation
%%%%
mgCaperday = 125/9; % 125 mg per day for whole litter
npups = 9;

gperday = npups * mgCaperday * 1e-3;
molperday = gperday / (2 * 40.078);
molpermin = molperday / (60 * 24);
mmolpermin = molpermin * 1e3;

Lact_GamFetusORMilk = mmolpermin






