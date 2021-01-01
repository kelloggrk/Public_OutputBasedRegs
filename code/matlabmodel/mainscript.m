% mainscript.m
% Ryan Kellogg
% Created: 12 April, 2018


%{
Runs models for quantities, emissions, and welfare. Calculates welfare
effects of indexed and non-indexed quantity standards
%}


clear all

% Set directories for repo and dropbox
S = pwd;
test = strcmp(S(end-14:end),'OutputBasedRegs') + strcmp(S(end-14:end),'outputbasedregs');
while test==0
    S = S(1:end-1);
    test = strcmp(S(end-14:end),'OutputBasedRegs') + strcmp(S(end-14:end),'outputbasedregs');
end
clear test
cd(S)
globals         % call path names in globals.m
clear S

% Set up directories for code, data, and output
wdir = strcat(repodir, '/code/matlabmodel');
fdir = strcat(repodir, '/paper/Figures/');
tdir = strcat(repodir, '/paper/Tables/');
rdir = strcat(dropbox, '/output/other/');
gpdir = strcat(dropbox, '/output/fuelprices/');
workspacedir = strcat(dropbox, '/output/model/');
sndir = strcat(repodir, '/paper/SingleNumberTex/');

cd(wdir)        % working directory for model objects


% First step below is to define parameter inputs to all models
% Storing all parameters in structs, with one struct per model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1A: DEFINE INPUTS TO FOOTPRINT MODEL fp IN WHICH
%   Q = VEHICLE FOOTPRINT (square feet)
%   E = VEHICLE FUEL ECONOMY (gallons per 100 miles)
%   eta = SHOCKS TO DEMAND FOR FOOTPRINT
%   Ignore VMT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Direct inputs
infile = strcat(rdir, 'LifetimeDiscMiles.csv');
B_EF = -csvread(infile,0,0) / 100;  % Discounted lifetime VMT (100 miles). BKZ (2013) and AW (2014)
clear infile
slopes = [4.04 3.865556]';  % slopes of abr set in 2012 for trucks and cars, in CO2 g/mile per sqft
gasCO2 = 8.91;          % kg CO2 / gallon
Q0 = 56.8349;       % 2012 sales-weighted average footprint (sq ft)
W0 = 3665.7355;     % 2012 sales-weighted average weight (lbs)
E0 = 24.934446;     % 2012 sales-weighted average mpg
B_QQwIS = -0.99 / 0.43 * 1000 / 100^2;  % B_QQ in units of $ per kg^2. From Ito & Sallee table 3 column (2)
infile = strcat(rdir, 'Footprint_B_EE.csv');
B_EE = -csvread(infile,0,0);    % From 2015 Natl Academies Panel. Units of $ / (gal/100mi)^2
clear infile
B_Qeta = 1;                           % normalization. $ per (sqft * shock)
credit = 36;        % Tesla credit sale in 2012 in $/Mg CO2
conv = 195264;      % EPA & NHTSA conversion for lifetime miles driven / car
phig = 0.42;        % from Kellogg JPubE (2018), externality in $/gal
wtfootelast = 1;    % elast of weight w.r.t footprint (Whitefoot & Skerlos 2012)
mton_ton = 1.10231; % metric ton per ton
infile = strcat(gpdir, 'RealGasPrice2012.csv');
F0 = csvread(infile,0,0);   % 2012 retail gasoline price from EIA, $/gal real 2012
Th = 10;            % regulatory time horizon in years
clear infile

% Derived model parameters
gamma0 = mean(slopes) / (gasCO2*1000) * 100;     % abr slope in gal/100mi per sqft
B_QE = -gamma0 * B_EE;
B_QQw = B_QQwIS / (mton_ton*2)^2;       % convert to $ per lb^2
E0 = 1 / E0 * 100;                      % convert to gal / 100 miles
wtslope = wtfootelast * W0 / Q0;        % change in weight corresponding to chg in footprint
B_QQ = B_QQw * wtslope^2;               % units of $ per ft^4
clear wtslope wtfootelast B_QQw
phigact = credit/1000 * gasCO2;         % actual CAFE shadow value in $/gallon (lambda_g)
phiact = conv/100 * phigact;            % actual CAFE shadow val in $ per (gal/100mi) (lambda)
phi = -B_EF * phig;                     % externality in $ per (gal/100mi)
Eeta = -B_QQ / B_Qeta * Q0;             % expected marginal value of footprint
clear conv

% Directly import fuel price volatility
infile = strcat(gpdir, 'GasolinePriceStdDeviations_Levels_MA.csv');
sigmaFvec = csvread(infile,0,0);
sigmaF = mean(sigmaFvec);       % average vol over 10 year horizon
clear infile

% Import footprint and fuel price data; impute variance of eta. 
% Use implication from model that variation in F does not
% affect footprint choice before 2010 (when fuel econ std was flat)
infile = strcat(rdir, 'AnnualFootprint_PGas.csv');    % columns for pgas and footprint
footpg = csvread(infile,1,1); footvec = footpg(:,2); pgasvec = footpg(:,1);
clear infile footpg
% obtain year-to-year gas price shocks
T = length(footvec);
diff_g = pgasvec(2:T) - pgasvec(1:T-1);
% Run AR1 regression for footprints. Get initial AR1 estimate, residuals,
% and residual variance
[AR1hat,~,diff_f,~,stats] = regress(footvec(2:T), [ones(T-1,1) [1:T-1]' footvec(1:T-1)]);
AR1 = 1 / (T-3) * (T * AR1hat(3) + 1);  % correct AR1 for short T
footvol1 = sqrt(stats(4));  % error std deviation
% Extrapolate volatility to Th years and take average
footvolvec = zeros(Th,1);      % initialize
for t = 1:Th
    footvolvec(t) = footvol1 * (1-AR1^(2*t))^0.5 / (1-AR1^2)^0.5;
end
sigmafoot = mean(footvolvec);  % volatility of footprint
% Footprint volatility assuming random walk
footvol1_rw = std(footvec(2:T) - footvec(1:T-1));
ind = 1:Th;
sigmafoot_rw = footvol1_rw * mean(ind'.^0.5);
% Convert footprint volatility to eta volatility
sigmaeta = sigmafoot / -B_Qeta * B_QQ;          % vol of eta
sigmaeta_rw = sigmafoot_rw / -B_Qeta * B_QQ;    % vol of eta, assuming random walk
sigma = sigmaeta;       % volatility in single index model (same for footprint application)
sigma_rw = sigmaeta_rw; % random walk in single index model
% Correlation between eta and F shocks, with p-val
[R,P] = corrcoef(diff_f,diff_g); R = R(2,1); P = P(2,1);
rho = 0;    % force zero for main analysis
clear footvec footvolvec gasvol1 gasvolvec
clear ind pgasvec T diff* t stats

% Define parameter struct to input into model
paramsfp.BQQ = B_QQ; paramsfp.BQE = B_QE; paramsfp.BEE = B_EE; paramsfp.BEF = B_EF;
paramsfp.BQeta = B_Qeta; paramsfp.phi = phi; paramsfp.phiact = phiact; paramsfp.gamma0 = gamma0;
paramsfp.Q0 = Q0; paramsfp.E0 = E0; paramsfp.Eeta = Eeta; paramsfp.F0 = F0; paramsfp.sigma = sigma;
paramsfp.sigmaeta = sigmaeta; paramsfp.sigmaF = sigmaF; paramsfp.rho = rho;
paramsfp.RP = [R P]; paramsfp.AR1 = AR1; paramsfp.AR1hat = AR1hat;
clear B_QQ B_QE B_EE B_EF B_Qeta phi phiact gamma0 Q0 E0 F0 Eeta sigma...
    sigmaeta sigmaF rho R P AR1 AR1hat

% Model with nonzero correlation between eta and F
paramsfpr = paramsfp;
paramsfpr.rho = paramsfp.RP(1);

% Model with random walk volatility for eta
paramsfprw = paramsfp;
paramsfprw.sigmaeta = sigmaeta_rw;
paramsfprw.sigma = sigma_rw;
clear sigmaeta_rw sigma_rw

% Model with sigmaF derived from monthly rather than MA data
paramsfphf = paramsfp;
infile = strcat(gpdir, 'MeanGasolinePriceStdDeviation_Levels.csv');
paramsfphf.sigmaF = csvread(infile,0,0);
clear infile

% Create single-number tex files for paper
% B_QQw
fid = fopen(strcat(sndir, 'foot_B_QQwIS.tex'),'w');
fprintf(fid,'%8.3f', -B_QQwIS);
fclose(fid); clear B_QQwIS
% B_QQ
fid = fopen(strcat(sndir, 'foot_B_QQ.tex'),'w');
fprintf(fid,'%8.0f', -paramsfp.BQQ);
fclose(fid);
% B_QE
fid = fopen(strcat(sndir, 'foot_B_QE.tex'),'w');
fprintf(fid,'$\\$ %8.1f $', paramsfp.BQE);
fclose(fid);
% B_EE
fid = fopen(strcat(sndir, 'foot_B_EE.tex'),'w');
fprintf(fid,'%8.0f', -paramsfp.BEE);
fclose(fid);
% car slope
fid = fopen(strcat(sndir, 'foot_carslope.tex'),'w');
fprintf(fid,'%8.2f', slopes(2));
fclose(fid);
% truck slope
fid = fopen(strcat(sndir, 'foot_truckslope.tex'),'w');
fprintf(fid,'%8.2f', slopes(1));
fclose(fid); clear slopes
% gasCO2 emissions factor for gasoline
fid = fopen(strcat(sndir, 'foot_gasCO2.tex'),'w');
fprintf(fid,'%8.2f', gasCO2);
fclose(fid); clear gasCO2
% phig
fid = fopen(strcat(sndir, 'foot_phig.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', phig);
fclose(fid);
% Q0
fid = fopen(strcat(sndir, 'foot_Q0.tex'),'w');
fprintf(fid,'%8.1f', paramsfp.Q0);
fclose(fid);
% E0
fid = fopen(strcat(sndir, 'foot_E0.tex'),'w');
fprintf(fid,'%8.2f', paramsfp.E0);
fclose(fid);
% F0
fid = fopen(strcat(sndir, 'foot_F0.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', paramsfp.F0);
fclose(fid);
% W0
fid = fopen(strcat(sndir, 'foot_W0.tex'),'w');
fprintf(fid,'%8.0f', W0);
fclose(fid); clear W0
% credit
fid = fopen(strcat(sndir, 'foot_credit.tex'),'w');
fprintf(fid,'$\\$ %8.0f $', credit);
fclose(fid); clear credit
% sigmaF1
fid = fopen(strcat(sndir, 'foot_sigmaF1.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', sigmaFvec(1));
fclose(fid);
% sigmaF10
fid = fopen(strcat(sndir, 'foot_sigmaF10.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', sigmaFvec(length(sigmaFvec)));
fclose(fid); clear sigmaFvec
% sigmaF
fid = fopen(strcat(sndir, 'foot_sigmaF.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', paramsfp.sigmaF);
fclose(fid);
% high vol sigmaF
fid = fopen(strcat(sndir, 'foot_sigmaFh.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', paramsfphf.sigmaF);
fclose(fid);
% AR1hat
fid = fopen(strcat(sndir, 'foot_AR1hat.tex'),'w');
fprintf(fid,'%8.2f', paramsfp.AR1hat(3));
fclose(fid);
% AR1
fid = fopen(strcat(sndir, 'foot_AR1.tex'),'w');
fprintf(fid,'%8.2f', paramsfp.AR1);
fclose(fid);
% footvol1
fid = fopen(strcat(sndir, 'foot_footvol1.tex'),'w');
fprintf(fid,'%8.2f', footvol1);
fclose(fid); clear footvol1
% footvol1_rw
fid = fopen(strcat(sndir, 'foot_footvol1_rw.tex'),'w');
fprintf(fid,'%8.2f', footvol1_rw);
fclose(fid); clear footvol1_rw
% sigmafoot
fid = fopen(strcat(sndir, 'foot_sigmafoot.tex'),'w');
fprintf(fid,'%8.2f', sigmafoot);
fclose(fid); clear sigmafoot
% sigmafoot_rw
fid = fopen(strcat(sndir, 'foot_sigmafoot_rw.tex'),'w');
fprintf(fid,'%8.2f', sigmafoot_rw);
fclose(fid); clear sigmafoot_rw
% sigmaeta_rw
fid = fopen(strcat(sndir, 'foot_sigmaeta_rw.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', paramsfprw.sigmaeta);
fclose(fid);
% rho
fid = fopen(strcat(sndir, 'foot_rho.tex'),'w');
fprintf(fid,'%8.2f', paramsfp.RP(1));
fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1B: DEFINE INPUTS TO FUEL ECONOMY (VMT) MODEL fe, IN WHICH
%   Q = LIFETIME DISCOUNTED VEHICLE VMT (100 miles)
%   e = VEHICLE FUEL ECONOMY (gallons per 100 miles)
%   E = LIFETIME FUEL USE = Qe (gallons)
%   eta = SHOCKS TO DEMAND FOR VMT
%   Ignore attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct inputs
Q0 = -paramsfp.BEF;     % lifetime VMT (100 miles) equals -B_EF from attribute basing parameters
B_EF = -1;              % effect of $1/gal increase in gas price on marg utility of gallons
B_ee = paramsfp.BEE;    % convert attribute based problem to new notation
e0 = paramsfp.E0;       % 2012 sales-weighted average gal / 100 miles
elastQF = -0.081;       % Gillingham (2018) lit review of studies using odometer readings
B_Qeta = 1;             % normalization
F0 = paramsfp.F0;       % 2012 retail gasoline price from EIA, $/gal real 2012
sigmaF = paramsfp.sigmaF;   % same fuel price volatility as in footprint application

% Derived model parameters
B_EE = B_ee / Q0^2;     % rescale B_ee from attribute basing parameters. $ / gal^2
E0 = Q0 * e0;           % 2012 lifetime discounted gallons
dQdF = elastQF * Q0 / F0;   % convert elasticity to derivative, holds e fixed. Units: 100 miles per ($/gal)
U_QQ = e0 / dQdF;           % Units of $ per (100 miles)^2
MCe = -Q0 * F0 - paramsfp.phiact;    % marginal cost of gal/100mile in 2012 ($ per (gal/100mi))
phi = -B_EF * phig;         % externality in $ per car
phiact = -B_EF * phigact;   % actual reg shadow cost in $ per gal
gamma0 = e0;                % actual intenstity standard in 2012 (gal/100mi)
B_QQ = U_QQ + B_ee * E0^2 / Q0^4 - 2 * MCe * E0 / Q0^3; % $ per (100 miles)^2
B_QE = -B_ee * E0/ Q0^3 + MCe / Q0^2;                   % $ per gallon per 100 miles
dQdeta = -B_Qeta / (B_EE*gamma0^2 + 2*B_QE*gamma0 + B_QQ);  % 100 miles per eta (marg val of vmt)
Eeta = Q0 / dQdeta;         % expected marginal value of 100 miles
clear U_QQ

% Import VMT and fuel price data; impute variance of eta
% First scale VMT so that VMT in 2012 matches Q0. Then subtract off (from
% other years) changes induced by fuel price variation. Then calculate
% residual variance
infile = strcat(rdir, 'AnnualVMT_PGas.csv');    % columns for pgas and vmt
vmtpg = csvread(infile,1,1); vmtvec = vmtpg(:,2); pgasvec = vmtpg(:,1);
infile = strcat(rdir, 'VMT2012.csv');
vmt2012 = csvread(infile,1,0);
clear infile vmtpg
% Scale annual vmt vector and adjust for fuel prices
vmtvec = vmtvec * Q0 / vmt2012;
Fadj = dQdF * (pgasvec - F0);   % VMT change (vs 2012) due to fuel price diff
vmtadj = vmtvec - Fadj;         % adjusted vmt
% Take first differences and compute one-year volatility
T = length(vmtadj);
diff_g = pgasvec(2:T) - pgasvec(1:T-1);
diff = vmtvec(2:T) - vmtvec(1:T-1);
diff_v = vmtadj(2:T) - vmtadj(1:T-1); 
vol1 = std(diff); vmtvol1 = std(diff_v);
% Extrapolate volatility to Th years and take average
% for both raw vmt data and vmt adjusted for fuel price shocks
ind = 1:Th;
volvec = vol1 * ind'.^0.5; vmtvolvec = vmtvol1 * ind'.^0.5;
sigma = mean(volvec); sigmavmt = mean(vmtvolvec);       % volatility of vmt
% Convert vmt volatility to eta volatility
sigma = sigma / dQdeta;             % vol of eta
sigmaeta = sigmavmt / dQdeta;       % vol of eta (adj for fuel prices)
% Correlation between eta and F shocks, with p-val
[R,P] = corrcoef(diff_v,diff_g); R = R(2,1); P = P(2,1);
rho = 0;    % force zero for main analysis
clear vmt2012 Fadj vmtadj T diff vmtvol1 ind vmtvolvec...
    dQdF vol1 volvec dQdeta diff*

% Define parameter struct to input into model
paramsfe.BQQ = B_QQ; paramsfe.BQE = B_QE; paramsfe.BEE = B_EE; paramsfe.BEF = B_EF;
paramsfe.BQeta = B_Qeta; paramsfe.phi = phi; paramsfe.phiact = phiact; paramsfe.gamma0 = gamma0;
paramsfe.Q0 = Q0; paramsfe.E0 = E0; paramsfe.Eeta = Eeta; paramsfe.F0 = F0; paramsfe.sigma = sigma;
paramsfe.sigmaeta = sigmaeta; paramsfe.sigmaF = sigmaF; paramsfe.rho = rho;
paramsfe.RP = [R P]; 

% Model with nonzero correlation between eta and F
paramsfer = paramsfe;
paramsfer.rho = paramsfe.RP(1);



%%%%%%%%%%%
% High elasticity case for vmt
elastQF0 = elastQF;     % save original Gillingham elasticity
sigmavmt0 = sigmavmt;   % save vmt volatility under original Gillingham elasticity
elastQF = -0.141;       % Gillingham (2018) lit review of all studies
dQdF = elastQF * Q0 / F0;   % convert elasticity to derivative, holds e fixed. Units: 100 miles per ($/gal)
U_QQ = e0 / dQdF;           % Units of $ per (100 miles)^2
B_QQ = U_QQ + B_ee * E0^2 / Q0^4 - 2 * MCe * E0 / Q0^3; % $ per (100 miles)^2
dQdeta = -B_Qeta / (B_EE*gamma0^2 + 2*B_QE*gamma0 + B_QQ);  % 100 miles per eta (marg val of vmt)
Eeta = Q0 / dQdeta;         % expected marginal value of 100 miles
clear U_QQ e0 B_ee MCe
Fadj = dQdF * (pgasvec - F0);   % VMT change (vs 2012) due to fuel price diff
vmtadj = vmtvec - Fadj;         % adjusted vmt
% obtain year-to-year gas price shocks
T = length(vmtadj);
diff_g = pgasvec(2:T) - pgasvec(1:T-1);
diff = vmtvec(2:T) - vmtvec(1:T-1);
diff_v = vmtadj(2:T) - vmtadj(1:T-1); 
vol1 = std(diff); vmtvol1 = std(diff_v);
% Extrapolate volatility to Th years and take average
% for both raw vmt data and vmt adjusted for fuel price shocks
ind = 1:Th;
volvec = vol1 * ind'.^0.5; vmtvolvec = vmtvol1 * ind'.^0.5;
sigma = mean(volvec); sigmavmt = mean(vmtvolvec);       % volatility of vmt
% Convert vmt volatility to eta volatility
sigma = sigma / dQdeta;             % vol of eta
sigmaeta = sigmavmt / dQdeta;       % vol of eta (adj for fuel prices)
% Correlation between eta and F shocks, with p-val
[R,P] = corrcoef(diff_v,diff_g); R = R(2,1); P = P(2,1);
rho = 0;    % force zero for main analysis
clear vmtvec pgasvec Th ind Fadj vmtadj vmtvol1 vmtvolvec vol1 volvec
clear dQdeta dQdF diff* sigmavmt T

% Define parameter struct to input into model
paramsfeh.BQQ = B_QQ; paramsfeh.BQE = B_QE; paramsfeh.BEE = B_EE; paramsfeh.BEF = B_EF;
paramsfeh.BQeta = B_Qeta; paramsfeh.phi = phi; paramsfeh.phiact = phiact; paramsfeh.gamma0 = gamma0;
paramsfeh.Q0 = Q0; paramsfeh.E0 = E0; paramsfeh.Eeta = Eeta; paramsfeh.F0 = F0; paramsfeh.sigma = sigma;
paramsfeh.sigmaeta = sigmaeta; paramsfeh.sigmaF = sigmaF; paramsfeh.rho = rho;
paramsfeh.RP = [R P]; 
clear B_QQ B_QE B_EE B_EF B_Qeta phi* gamma0 Q0 E0 F0 Eeta sigma...
    sigmaeta sigmaF rho R P

% Create single-number tex files for paper
% Q0
fid = fopen(strcat(sndir, 'vmt_Q0.tex'),'w');
fprintf(fid,'%8.0f', paramsfe.Q0*100);
fclose(fid);
% VMT elasticity
fid = fopen(strcat(sndir, 'vmt_elast0.tex'),'w');
fprintf(fid,'%8.3f', elastQF0);
fclose(fid); clear elastQF0
% High VMT elasticity
fid = fopen(strcat(sndir, 'vmt_elasth.tex'),'w');
fprintf(fid,'%8.3f', elastQF);
fclose(fid); clear elastQF
% sigmavmt
fid = fopen(strcat(sndir, 'vmt_sigmavmt.tex'),'w');
fprintf(fid,'%8.0f', sigmavmt0*100);
fclose(fid); clear sigmavmt0
% rho
fid = fopen(strcat(sndir, 'vmt_rho.tex'),'w');
fprintf(fid,'%8.2f', paramsfe.RP(1));
fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1C: DEFINE INPUTS TO ELECTRICITY MODEL eg, IN WHICH
%   Q = TOTAL ANNUAL GENERATION (TWh/year)
%   E = TOTAL CO2 EMISSIONS = (millions of metric tons / year (Mmton / year)
%   eta = shock to demand for Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Direct inputs
E0 = 2031.452;      % Baseline Mmton/year. EIA Electric Power Annual for 2015
phi0 = 38;          % Damages per metric ton CO2
phiact = 0;         % Actual regulatory price wedge
gamma0 = 0;         % Actual regulatory slope
B_Qeta = 1;         % Normalization
PctFallE_40perton = 7.9;    % % fall in E for a $40/ton CO2 tax, holding Q fixed. Cullen and Mansur table 2
elastQ = -0.09;     % Ito (2014) elasticity over monthly price lags from 0 to 4; also Deryungina 0-6 month elast.
P0 = 0.1265;        % $/kWh. Average US price to end users in 2015
PCoal0 = 42.58;     % $/short ton. For 2015. EIA Annual Coal Report.
mmBtu_per_short = 19.146;   % mmBtu per short ton for 2015. Monthly Energy Review table A5
Eng = 117.0;        % emission factor for natural gas, lbs CO2/mmBtu (Cullen and Mansur 2017)
Ecoal = 210.86;     % emission factor for coal, lbs CO2/mmBtu (Cullen and Mansur 2017)
Th = 15;            % regulatory time horizon in years

% Import 2015 elec generation and real $ Jan 2016 gas price for 2015
infile = strcat(rdir, 'ElecGen_Pnatgas_2015.csv');    % columns for pgas and vmt
tempin = csvread(infile,1,0);
Q0 = tempin(1) / 1000000;     % baseline elec gen in TWh/year 
F0 = tempin(2);               % baseline natural gas price in $/mmBtu
clear infile tempin

% Derived model parameters
PCoal = PCoal0 / mmBtu_per_short;   % $/mmBtu
phi = phi0 * 1000000;       % Damages per Mmton
dEdphi = -PctFallE_40perton / 100 * E0 / 40;    % Change in E (Mmton/year) per $/ton tax (Q fixed)
dEdphi = dEdphi / mton_ton / 1000000;           % Mmton/year per $/Mmton
B_EE = 1 / dEdphi;          % $ per Mmton^2
clear dEdphi
% dEdQ = dEdQ_GZ / 2000 / mton_ton / 1000000 * 1E9;  % Mmton CO2 per TWh, using Graff Zivin et al
dEdQ = E0/Q0;               % dEdQ assuming CES for overall elec industry
B_QE = -dEdQ * B_EE;        % $ per Mmton per TWh
clear dEdQ_GZ
P0 = P0 * 1E9;              % $ per TWh
dQdP = elastQ * Q0 / P0;    % TWh^2 per $
B_QQ = 1 / dQdP;            % $ per TWh^2
elastQ0 = elastQ;           % store elasticity
clear dQdP elastQ
denom = B_QQ*B_EE - B_QE^2;
dQdeta = -B_EE * B_Qeta / denom;        % Twh per eta (marg val of kWh)
Eeta = Q0 / dQdeta;     % expected marginal value of elec
CR = PCoal / F0;        % coal vs gas cost ratio at baseline with no CO2 price
dPngdphi = (CR*Eng-Ecoal)/mton_ton/2000 / CR;      % gas price change equiv to $1/mton CO2 change
dPngdphi = dPngdphi / 1000000;  % gas price change equiv to $1/Mmton CO2 change. Mmton CO2 / mmBtu
B_EF = 1 / dPngdphi;   % $ per Mmton CO2 per $/mmBtu gas. No negative sign since F = -Png
clear CR dPngdphi PCoal
dQdF = B_QE*B_EF / denom;               % TWh per $/mmBtu
clear denom

% Directly import HH natural gas price volatility
infile = strcat(gpdir, 'MeanHHStdDeviation_Levels.csv');
sigmaF = csvread(infile,0,0);
clear infile

% Import elec gen and natural gas price data; impute variance of eta
% First subtract off changes in gen induced by fuel price variation. Then calculate
% residual variance
infile = strcat(rdir, 'AnnualElecGen_Pnatgas.csv');    % columns for pgas and vmt
genpng = csvread(infile,1,1); pngvec = genpng(:,2); genvec = genpng(:,1)/ 1000000;
clear infile genpng
% Adjust generation for fuel price shocks to isolate demand shocks
pngadj = -dQdF * (pngvec - F0);     % gen change (vs 2015) due to nat gas price diff (F=-Pg)
genvecadj = (genvec - Q0) - pngadj; % isolate year-to-year differences from demand shocks
% Take first differences and compute one-year volatility
T = length(genvecadj); 
diff = genvec(2:T) - genvec(1:T-1);         % generation change
diff_e = genvecadj(2:T) - genvecadj(1:T-1); % gen change due to eta
diff_f = -pngvec(2:T) + pngvec(1:T-1);      % F change (F = -Pg)
vol1 = std(diff); genvol1 = std(diff_e);
% Correlation between eta and F shocks, with p-val
[R,P] = corrcoef(diff_e,diff_f); R = R(2,1); P = P(2,1);
rho = 0;    % force zero for main analysis

% Extrapolate volatility to Th years and take average
ind = 1:Th;
volvec = vol1 * ind'.^0.5; genvolvec = genvol1 * ind'.^0.5;
sigma = mean(volvec) / dQdeta;          % vol of eta
sigmaeta = mean(genvolvec) / dQdeta;    % vol of eta (adj for gas prices)
clear dQdF dQdeta mton_ton pngadj genvecadj diff_e diff
clear genvol1 volvec genvolvec

% Define parameter struct to input into model
paramseg.BQQ = B_QQ; paramseg.BQE = B_QE; paramseg.BEE = B_EE; paramseg.BEF = B_EF;
paramseg.BQeta = B_Qeta; paramseg.phi = phi; paramseg.phiact = phiact; paramseg.gamma0 = gamma0;
paramseg.Q0 = Q0; paramseg.E0 = E0; paramseg.Eeta = Eeta; paramseg.F0 = F0; paramseg.sigma = sigma;
paramseg.sigmaeta = sigmaeta; paramseg.sigmaF = sigmaF; paramseg.rho = rho;
paramseg.RP = [R P];

% Model with nonzero correlation between eta and F
paramsegr = paramseg;
paramsegr.rho = paramseg.RP(1);

% Model with fixed positive correlation between eta and F
paramsegrh = paramseg;
paramsegrh.rho = 0.25;

% Model with half eta volatility
paramseghalf = paramseg;
paramseghalf.sigmaeta = paramseg.sigmaeta / 2;

% Model with double eta volatility
paramsegdoub = paramseg;
paramsegdoub.sigmaeta = paramseg.sigmaeta * 2;


%%%%%%%%%%%
% High elasticity case for electricity
elastQ = -0.27;     % Deryugina et al (2017) ten-year elasticity lower bound
dQdP = elastQ * Q0 / P0;    % TWh^2 per $
B_QQ = 1 / dQdP;            % $ per TWh^2
clear dQdP
denom = B_QQ*B_EE - B_QE^2;
dQdF = B_QE*B_EF / denom;               % TWh per $/mmBtu
dQdeta = -B_EE * B_Qeta / denom;        % Twh per eta (marg val of kWh)
Eeta = Q0 / dQdeta;             % expected marginal value of elec
% Adjust generation for fuel price shocks to isolate demand shocks
pngadj = -dQdF * (pngvec - F0);     % gen change (vs 2015) due to nat gas price diff (F=-Pg)
genvecadj = (genvec - Q0) - pngadj; % isolate year-to-year differences from demand shocks
% Take first differences and compute one-year volatility
diff_e = genvecadj(2:T) - genvecadj(1:T-1); % gen change due to eta
genvol1 = std(diff_e);
% Correlation between eta and F shocks, with p-val
[R,P] = corrcoef(diff_e,diff_f); R = R(2,1); P = P(2,1);
rho = 0;    % force zero for main analysis

% Extrapolate volatility to Th years and take average
ind = 1:Th;
volvec = vol1 * ind'.^0.5; genvolvec = genvol1 * ind'.^0.5;
sigma = mean(volvec) / dQdeta;          % vol of eta
sigmaeta = mean(genvolvec) / dQdeta;    % vol of eta (adj for gas prices)
clear denom dQdF dQdeta pngvec genvec pngadj genvecadj diff_e diff_f 
clear vol1 genvol1 ind T Th genvolvec volvec

% Define parameter struct to input into model
paramsegh.BQQ = B_QQ; paramsegh.BQE = B_QE; paramsegh.BEE = B_EE; paramsegh.BEF = B_EF;
paramsegh.BQeta = B_Qeta; paramsegh.phi = phi; paramsegh.phiact = phiact; paramsegh.gamma0 = gamma0;
paramsegh.Q0 = Q0; paramsegh.E0 = E0; paramsegh.Eeta = Eeta; paramsegh.F0 = F0; paramsegh.sigma = sigma;
paramsegh.sigmaeta = sigmaeta; paramsegh.sigmaF = sigmaF; paramsegh.rho = rho;
paramsegh.RP = [R P];
clear B_QQ B_QE B_EE B_EF B_Qeta phi* gamma0 Q0 E0 Eeta sigma F0 rho...
    sigmaF sigma sigmaeta R P

% Create single-number tex files for paper
% Q0
fid = fopen(strcat(sndir, 'elec_Q0.tex'),'w');
fprintf(fid,'%8.0f', paramseg.Q0);
fclose(fid);
% E0
fid = fopen(strcat(sndir, 'elec_E0.tex'),'w');
fprintf(fid,'%8.0f', paramseg.E0);
fclose(fid);
% sigmaF
fid = fopen(strcat(sndir, 'elec_sigmaF.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', paramseg.sigmaF);
fclose(fid);
% F0
fid = fopen(strcat(sndir, 'elec_F0.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', paramseg.F0);
fclose(fid);
% P0
fid = fopen(strcat(sndir, 'elec_P0.tex'),'w');
fprintf(fid,'$\\$ %8.3f $', P0/1E9);
fclose(fid); clear P0
% rho
fid = fopen(strcat(sndir, 'elec_rho.tex'),'w');
fprintf(fid,'%8.2f', paramseg.RP(1));
fclose(fid);
% Elec demand elasticity
fid = fopen(strcat(sndir, 'elec_elastQ0.tex'),'w');
fprintf(fid,'%8.2f', elastQ0);
fclose(fid); clear elastQ0
% High elec demand elasticity
fid = fopen(strcat(sndir, 'elec_elastQh.tex'),'w');
fprintf(fid,'%8.2f', elastQ);
fclose(fid); clear elastQ
% Cullen-Mansur pct decrease in E
fid = fopen(strcat(sndir, 'elec_CullenMansur.tex'),'w');
fprintf(fid,'%8.1f\\%%', PctFallE_40perton);
fclose(fid); clear PctFallE_40perton
% dEdQ
fid = fopen(strcat(sndir, 'elec_dEdQ.tex'),'w');
fprintf(fid,'%8.2f', dEdQ);
fclose(fid); clear dEdQ
% Eng
fid = fopen(strcat(sndir, 'elec_Eng.tex'),'w');
fprintf(fid,'%8.2f', Eng);
fclose(fid); clear Eng
% Ecoal
fid = fopen(strcat(sndir, 'elec_Ecoal.tex'),'w');
fprintf(fid,'%8.2f', Ecoal);
fclose(fid); clear Ecoal
% PCoal0
fid = fopen(strcat(sndir, 'elec_PCoal0.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', PCoal0);
fclose(fid); clear PCoal0
% mmBtu_per_short
fid = fopen(strcat(sndir, 'elec_mmBtu_per_short.tex'),'w');
fprintf(fid,'%8.2f', mmBtu_per_short);
fclose(fid); clear mmBtu_per_short




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAPER-FORMATTED TABLES WITH INPUT VALUES FOR ALL THREE APPLICATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Footprint standards
filenameo = strcat(tdir,'tbl_footprintinputs.tex');
fileid = fopen(char(filenameo),'w');
fprintf(fileid,'$B_{QQ}$ & -\\$%3.0f / ft$^4$ & Ito and Sallee (2018), table 3, column (2)  \\\\ \n',-paramsfp.BQQ); 
fprintf(fileid,...
    '$B_{EE}$ & -\\$%4.0f / (gal/100mi)$^2$ & National Research Council (2015), tables 8.4a and 8.4b  \\\\ \n',-paramsfp.BEE); 
fprintf(fileid,...
    '$B_{QE}$ & \\$%4.1f / (gal$\\cdot$ft$^2$/100mi) & $B_{EE}$ times slope of U.S. footprint-based standard  \\\\ \n',paramsfp.BQE); 
fprintf(fileid,...
    '$B_{EF}$ & -%6.0f miles & Discounted vehicle lifetime from Busse {\\it et al.} (2013)  \\\\ \n',-paramsfp.BEF*100); 
fprintf(fileid,'$B_{Q\\eta}$ & 1 & normalization  \\\\ \n'); 
fprintf(fileid,'$\\sigma_F$ & \\$%3.2f / gal & historic volatility of gasoline prices  \\\\ \n',paramsfp.sigmaF); 
fprintf(fileid,'$\\sigma_\\eta$ & \\$%3.2f / ft$^2$ & historic volatility of vehicle footprints  \\\\ \n',paramsfp.sigmaeta); 
fprintf(fileid,...
    '$\\phi$ & \\$%3.0f  / (gal/100mi) & \\begin{minipage}[t]{0.5\\columnwidth} \n Kellogg (2018) externality of \\$0.42/gal, converted to vehicle lifetime using Busse {\\it et al.} (2013) \n \\end{minipage}  \\\\ \n',...
    paramsfp.phi); 
fclose('all');

% Endogenous miles traveled
filenameo = strcat(tdir,'tbl_VMTinputs.tex');
fileid = fopen(char(filenameo),'w');
fprintf(fileid,'$B_{QQ}$ & -\\$%3.2f / (1000mi)$^2$ & Gillingham (2020)  \\\\ \n',-paramsfe.BQQ*100); 
fprintf(fileid,...
    '$B_{EE}$ & -\\$%1.4f  / gallon$^2$ & $B_{EE}$ from table \\ref{tab:calfoot}, scaled by vehicle lifetime miles  \\\\ \n',-paramsfe.BEE); 
fprintf(fileid,...
    '$B_{QE}$ & \\$%1.3f / (gal$\\cdot$1000mi) & $B_{EE}$ from table \\ref{tab:calfoot}; Leard and McConnell (2017)  \\\\ \n',paramsfe.BQE*10); 
fprintf(fileid,...
    '$B_{EF}$ &  -1 & by construction  \\\\ \n'); 
fprintf(fileid,'$B_{Q\\eta}$ & 1 & normalization  \\\\ \n'); 
fprintf(fileid,'$\\sigma_F$ & \\$%3.2f / gal & historic volatility of gasoline prices  \\\\ \n',paramsfe.sigmaF); 
fprintf(fileid,'$\\sigma_\\eta$ & \\$%3.2f / 1000mi & historic volatility of miles traveled  \\\\ \n',paramsfe.sigmaeta*10); 
fprintf(fileid,...
    '$\\phi$ & \\$%3.2f  / gallon & Kellogg (2018) \\\\ \n',paramsfe.phi); 
fclose('all');

% Electricity
filenameo = strcat(tdir,'tbl_elecinputs.tex');
fileid = fopen(char(filenameo),'w');
fprintf(fileid,'$B_{QQ}$ & -\\$%6.0f / TWh$^2$ & Ito (2014), table 3  \\\\ \n',-paramseg.BQQ); 
fprintf(fileid,...
    '$B_{EE}$ & -\\$%6.0f / Mmton$^2$ & Cullen and Mansur (2017), table 2  \\\\ \n',-paramseg.BEE); 
fprintf(fileid,...
    '$B_{QE}$ & \\$%6.0f / (TWh$\\cdot$Mmton) & constant returns to scale assumption  \\\\ \n',paramseg.BQE); 
fprintf(fileid,...
    '$B_{EF}$ & -%2.1f million mmBtu per Mmton & Cullen and Mansur (2017) and EIA  \\\\ \n',-paramseg.BEF/1E6); 
fprintf(fileid,'$B_{Q\\eta}$ & 1 & normalization  \\\\ \n'); 
fprintf(fileid,'$\\sigma_F$ & \\$%3.2f / mmBtu & historic volatility of natural gas prices  \\\\ \n',paramseg.sigmaF); 
fprintf(fileid,'$\\sigma_\\eta$ & \\$%2.1f million / TWh & historic volatility of electricity consumption  \\\\ \n',paramseg.sigmaeta/1E6); 
fprintf(fileid,...
    '$\\phi$ & \\$%2.0f million / Mmton & Interagency Working Group (2013)  \\\\ \n',paramseg.phi/1E6); 
fclose('all');
clear filenameo fileid ans





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2
% Instantiate and run models with both eta and F uncertainty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FOOTPRINT MODEL
% Instantiate model
obj = indexmodelF(paramsfp);

% Run model
[EWfpo, EWfp0, EWfpi, EWfps, Qfp, Efp, EQfp, EEfp, MUfp, GAMMAfp] = resultsall(obj);
Q = Qfp; E = Efp; EQ = EQfp; EE = EEfp;

% Create arrow for eta direction
x0 = 56.2; y0 = 4.75;                % base of arrow
d = obj.BQQ * obj.BEE - obj.BQE^2;  % agents' SOC
xeta = -obj.BEE * obj.BQeta * obj.sigmaeta / d;
yeta = obj.BQE * obj.BQeta * obj.sigmaeta / d;

% Create arrow for F direction
xF = obj.BQE * obj.BEF * obj.sigmaF / d;
yF = -obj.BQQ * obj.BEF * obj.sigmaF / d;
clear d

% Plot results
plotgreen = 1/255*[0 153 51];
clf
plot(EQ,EE,'xk','MarkerSize',10); hold on       % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',2.5); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',2.5); % social opt
plot(Q.Q0,E.E0,':k','LineWidth',5);                 % flat std
plot(Q.Qi,E.Ei,'-.b','LineWidth',5);                % intensity std
plot(Q.Qs,E.Es,'-r','LineWidth',5);                 % optimal indexed std
xl = 56; xu = 57.75; yl = 3.4; yu = 5;               % axes
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20);   % eta label
text(x0+1.4*xF,y0+1.2*yF,'F','FontSize',20);            % F label
grid; xlabel('Footprint (square ft)');
ylabel('Gallons per 100 miles');
legend('Expected unconstrained choice','Unconstrained choice 95% c.i.',...
    'Social optimum 95% c.i.',...
    'Optimal flat standard','Optimal intensity standard',...
    'Optimal indexed standard','Location','southeast');
hold off
h = gcf;
set(gca,'FontSize',20);
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
outfile = strcat(fdir, 'FootprintPlot_doubleuncert.pdf');
print(outfile,'-dpdf','-fillpage');
clear h

% Plot results -- no optimal std
clf
plot(EQ,EE,'xk','MarkerSize',10); hold on       % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',2.5); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',2.5); % social opt
plot(Q.Q0,E.E0,':k','LineWidth',5);                 % flat std
plot(Q.Qi,E.Ei,'-.b','LineWidth',5);                % intensity std
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20);   % eta label
text(x0+1.4*xF,y0+1.2*yF,'F','FontSize',20);            % F label
grid; xlabel('Footprint (square ft)');
ylabel('Gallons per 100 miles');
legend('Expected unconstrained choice','Unconstrained choice 95% c.i.',...
    'Social optimum 95% c.i.',...
    'Optimal flat standard','Optimal intensity standard',...
    'Location','southeast');
hold off
h = gcf;
set(gca,'FontSize',20);
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
outfile = strcat(fdir, 'FootprintPlot_doubleuncert_noopt.pdf');
print(outfile,'-dpdf','-fillpage');
clear h

% Plot results -- no regs
clf
plot(EQ,EE,'xk','MarkerSize',10); hold on       % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',2.5); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',2.5); % social opt
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20);   % eta label
text(x0+1.4*xF,y0+1.2*yF,'F','FontSize',20);            % F label
grid; xlabel('Footprint (square ft)');
ylabel('Gallons per 100 miles');
legend('Expected unconstrained choice','Unconstrained choice 95% c.i.',...
    'Social optimum 95% c.i.',...
    'Location','southeast');
hold off
h = gcf;
set(gca,'FontSize',20);
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
outfile = strcat(fdir, 'FootprintPlot_doubleuncert_noregs.pdf');
print(outfile,'-dpdf','-fillpage');
clear h

% Plot results formatted for paper
clf
formatSpec = '$%2.0f';
str0 = [' Expected unconstrained choice: ', sprintf('$%1.0f',0)];
str1 = [' Unconstrained choices:               ', sprintf('$%1.0f',0)];
str2 = [' Socially optimal choices:            ',sprintf(formatSpec,EWfpo)];
str3 = [' Optimal flat std:                          ',sprintf(formatSpec,EWfp0)];
str4 = [' Optimal intensity std:                  ',sprintf(formatSpec,EWfpi)];
str5 = [' Optimal footprint-based std:       ',sprintf(formatSpec,EWfps)];
plotgreen = 1/255*[0 153 51];
plot(EQ,EE,'xk','MarkerSize',10,'DisplayName',str0); hold on       % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',1.5,'DisplayName',str1); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',1.5,'DisplayName',str2); % social opt
plot(Q.Q0,E.E0,':k','LineWidth',5,'DisplayName',str3);                 % flat std
plot(Q.Qi,E.Ei,'-.b','LineWidth',5,'DisplayName',str4);                % intensity std
plot(Q.Qs,E.Es,'-r','LineWidth',5,'DisplayName',str5);                 % optimal indexed std
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20,'FontName','Times New Roman');   % eta label
text(x0+1.4*xF,y0+1.2*yF,'F','FontSize',20,'FontName','Times New Roman');            % F label
grid; xlabel('Footprint (square ft)');
ylabel('Gallons per 100 miles');
legend('Location','southeast');
lgd = legend;
lgd.FontSize = 20; lgd.Title.String = 'Expected welfare ($/vehicle)';
hold off
h = gcf;
set(gca,'FontSize',20,'FontName','Times New Roman');
set(h,'PaperUnits','inches','PaperType','usletter')
set(h,'PaperOrientation','landscape','PaperPosition', [-0.3 -0.2 11.7 9.1]);
outfile = strcat(fdir, 'FootprintPlot_doubleuncert_paper.pdf');
print(h,outfile,'-dpdf');
clear h
clear str* formatSpec xl xu yl yu x0 y0 xeta yeta xF yF pos startx starty...
    endxe endye endxF endyF Q E EQ EE obj lgd

% Output single number tex
% EQfp
fid = fopen(strcat(sndir, 'foot_EQfp.tex'),'w');
fprintf(fid,'%8.1f', EQfp);
fclose(fid);
% EEfp
fid = fopen(strcat(sndir, 'foot_EEfp.tex'),'w');
fprintf(fid,'%8.2f', EEfp);
fclose(fid);
% EWfpo
fid = fopen(strcat(sndir, 'foot_EWfpo.tex'),'w');
fprintf(fid,'$\\$ %8.0f $', EWfpo);
fclose(fid);
% EWfp0
fid = fopen(strcat(sndir, 'foot_EWfp0.tex'),'w');
fprintf(fid,'$\\$ %8.0f $', EWfp0);
fclose(fid);
% EWfp0/EWfpo
fid = fopen(strcat(sndir, 'foot_ratio_flat_opt.tex'),'w');
fprintf(fid,'%8.0f\\%%', EWfp0 / EWfpo * 100);
fclose(fid);
% opt gamma
fid = fopen(strcat(sndir, 'foot_optgamma.tex'),'w');
fprintf(fid,'%8.3f', GAMMAfp(3));
fclose(fid);
% opt gamma to act gamma
fid = fopen(strcat(sndir, 'foot_ratio_optgamma_actgamma.tex'),'w');
fprintf(fid,'%8.0f\\%%', GAMMAfp(3) / paramsfp.gamma0 * 100);
fclose(fid);
% intensity gamma to act gamma
fid = fopen(strcat(sndir, 'foot_pctgain_ingamma_actgamma.tex'),'w');
fprintf(fid,'%8.0f\\%%', (GAMMAfp(2) / paramsfp.gamma0 - 1) * 100);
fclose(fid);



% FOOTPRINT MODEL WITH LARGE SIGMAF FROM MONTHLY GAS PRICE DATA
obj = indexmodelF(paramsfphf);
% Run model
[EWfphfo, EWfphf0, EWfphfi, EWfphfs, ~, ~, ~, ~, MUfphf, GAMMAfphf] = resultsall(obj);
% Output single number tex
% EWfphfs
fid = fopen(strcat(sndir, 'foot_EWfphfs.tex'),'w');
fprintf(fid,'$\\$ %8.0f $', EWfphfs);
fclose(fid);

% FOOTPRINT MODEL WITH RANDOM WALK ASSUMPTION ON Q
obj = indexmodelF(paramsfprw);
% Run model
[EWfprwo, EWfprw0, EWfprwi, EWfprws, ~, ~, ~, ~, MUfprw, GAMMAfprw] = resultsall(obj);
% Output single number tex
% opt gamma
fid = fopen(strcat(sndir, 'foot_optgamma_rw.tex'),'w');
fprintf(fid,'%8.3f', GAMMAfprw(3));
fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUEL ECONOMY MODEL
% Instantiate model
obj = indexmodelF(paramsfe);

% Run model
[EWfeo, EWfe0, EWfei, EWfes, Qfe, Efe, EQfe, EEfe, MUfe, GAMMAfe] = resultsall(obj);
Q = Qfe; E = Efe; EQ = EQfe; EE = EEfe;

% convert 100 miles to 1000 miles for plotting
Q.Qu = Q.Qu/10; Q.Qo = Q.Qo/10; Q.Q0=Q.Q0/10; Q.Qi=Q.Qi/10; Q.Qs=Q.Qs/10; EQ = EQ/10;

% Create arrow for eta direction
x0 = 104; y0 = 5400;                % base of arrow
d = obj.BQQ * obj.BEE - obj.BQE^2;  % agents' SOC
xeta = -obj.BEE * obj.BQeta * obj.sigmaeta / d / 10;
yeta = obj.BQE * obj.BQeta * obj.sigmaeta / d;

% Create arrow for F direction
xF = obj.BQE * obj.BEF * obj.sigmaF / d / 10;
yF = -obj.BQQ * obj.BEF * obj.sigmaF / d;
clear d

% Plot results
plotgreen = 1/255*[0 153 51];
clf
plot(EQ,EE,'xk','MarkerSize',10); hold on           % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',2.5); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',2.5); % social opt
plot(Q.Q0,E.E0,':k','LineWidth',5);                 % flat std
plot(Q.Qi,E.Ei,'-.b','LineWidth',5);                % intensity std
plot(Q.Qs,E.Es,'-r','LineWidth',5);                 % optimal indexed std
xl = 102; xu = 128; yl = 3500; yu = 5700;           % axes
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20);   % eta label
text(x0+1.4*xF,y0+1.2*yF,'F','FontSize',20);            % F label
grid; xlabel('Lifetime discounted thousands of miles driven per vehicle');
ylabel('Lifetime discounted gasoline use per vehicle (gal)');
legend('Expected unconstrained choice','Unconstrained choice 95% c.i.',...
    'Social optimum 95% c.i.',...
    'Optimal flat standard','Optimal intensity standard',...
    'Optimal indexed standard','Location','southeast');
hold off
h = gcf;
set(gca,'FontSize',20);
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
outfile = strcat(fdir, 'FuelEconomyPlot_doubleuncert.pdf');
print(outfile,'-dpdf','-fillpage');
clear h

% Plot results -- no regs
plotgreen = 1/255*[0 153 51];
clf
plot(EQ,EE,'xk','MarkerSize',10); hold on           % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',2.5); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',2.5); % social opt
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20);   % eta label
text(x0+1.4*xF,y0+1.2*yF,'F','FontSize',20);            % F label
grid; xlabel('Lifetime discounted thousands of miles driven per vehicle');
ylabel('Lifetime discounted gasoline use per vehicle (gal)');
legend('Expected unconstrained choice','Unconstrained choice 95% c.i.',...
    'Social optimum 95% c.i.',...
    'Location','southeast');
hold off
h = gcf;
set(gca,'FontSize',20);
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
outfile = strcat(fdir, 'FuelEconomyPlot_doubleuncert_noregs.pdf');
print(outfile,'-dpdf','-fillpage');
clear h

% Plot results formatted for paper
clf
formatSpec = '$%2.0f';
str0 = [' Expected unconstrained choice: ', sprintf('$%1.0f',0)];
str1 = [' Unconstrained choices:               ', sprintf('$%1.0f',0)];
str2 = [' Socially optimal choices:            ',sprintf(formatSpec,EWfeo)];
str3 = [' Optimal flat std:                          ',sprintf(formatSpec,EWfe0)];
str4 = [' Optimal intensity std:                  ',sprintf(formatSpec,EWfei)];
str5 = [' Optimal mileage-based std:         ',sprintf(formatSpec,EWfes)];
plotgreen = 1/255*[0 153 51];
plot(EQ,EE,'xk','MarkerSize',10,'DisplayName',str0); hold on       % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',1.5,'DisplayName',str1); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',1.5,'DisplayName',str2); % social opt
plot(Q.Q0,E.E0,':k','LineWidth',5,'DisplayName',str3);                 % flat std
plot(Q.Qi,E.Ei,'-.b','LineWidth',5,'DisplayName',str4);                % intensity std
plot(Q.Qs,E.Es,'-r','LineWidth',5,'DisplayName',str5);                 % optimal indexed std
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20,'FontName','Times New Roman');   % eta label
text(x0+1.4*xF,y0+1.2*yF,'F','FontSize',20,'FontName','Times New Roman');            % F label
grid; xlabel('Lifetime discounted thousands of miles driven per vehicle');
ylabel('Lifetime discounted gasoline use per vehicle (gal)');
legend('Location','southeast');
lgd = legend;
lgd.FontSize = 20; lgd.Title.String = 'Expected welfare ($/vehicle)';
hold off
h = gcf;
set(gca,'FontSize',20,'FontName','Times New Roman');
set(h,'PaperUnits','inches','PaperType','usletter')
set(h,'PaperOrientation','landscape','PaperPosition', [-0.3 -0.2 11.7 9.1]);
outfile = strcat(fdir, 'FuelEconomyPlot_doubleuncert_paper.pdf');
print(h,outfile,'-dpdf');
clear h
clear str* formatSpec xl xu yl yu x0 y0 xeta yeta xF yF pos startx starty...
    endxe endye endxF endyF Q E EQ EE obj lgd

% Output single number tex
% opt gamma
fid = fopen(strcat(sndir, 'vmt_optgamma.tex'),'w');
fprintf(fid,'%8.2f', GAMMAfe(3));
fclose(fid);
% opt intensity gamma
fid = fopen(strcat(sndir, 'vmt_ingamma.tex'),'w');
fprintf(fid,'%8.2f', GAMMAfe(2));
fclose(fid);
% EWfes
fid = fopen(strcat(sndir, 'vmt_EWfes.tex'),'w');
fprintf(fid,'$\\$ %8.1f $', EWfes);
fclose(fid);
% EWfe0
fid = fopen(strcat(sndir, 'vmt_EWfe0.tex'),'w');
fprintf(fid,'$\\$ %8.1f $', EWfe0);
fclose(fid);
% EWfei
fid = fopen(strcat(sndir, 'vmt_EWfei.tex'),'w');
fprintf(fid,'$\\$ %8.1f $', EWfei);
fclose(fid);


% FUEL ECONOMY MODEL WITH HIGH DEMAND ELASTICITY
% Instantiate model
obj = indexmodelF(paramsfeh);

% Run model
[EWfeho, EWfeh0, EWfehi, EWfehs, ~, ~, ~, ~, MUfeh, GAMMAfeh] = resultsall(obj);

% Output single number tex
% opt gamma
fid = fopen(strcat(sndir, 'vmt_optgammah.tex'),'w');
fprintf(fid,'%8.2f', GAMMAfeh(3));
fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% US ELECTRICITY MODEL
% Instantiate model
obj = indexmodelF(paramseg);

% Run model
[EWego, EWeg0, EWegi, EWegs, Qeg, Eeg, EQeg, EEeg, MUeg, GAMMAeg] = resultsall(obj);
Q = Qeg; E = Eeg; EQ = EQeg; EE = EEeg;

% Create arrow for eta direction
x0 = 3600; y0 = 2250;                % base of arrow
d = obj.BQQ * obj.BEE - obj.BQE^2;  % agents' SOC
xeta = -obj.BEE * obj.BQeta * obj.sigmaeta / d;
yeta = obj.BQE * obj.BQeta * obj.sigmaeta / d;

% Create arrow for F direction
xF = obj.BQE * obj.BEF * obj.sigmaF / d;
yF = -obj.BQQ * obj.BEF * obj.sigmaF / d;
clear d


% Plot results
plotgreen = 1/255*[0 153 51];
clf
plot(EQ,EE,'xk','MarkerSize',10); hold on           % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',2.5); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',2.5); % social opt
plot(Q.Q0,E.E0,':k','LineWidth',5);                 % flat std
plot(Q.Qi,E.Ei,'-.b','LineWidth',5);                % intensity std
plot(Q.Qs,E.Es,'-r','LineWidth',5);                 % optimal indexed std
xl = 3400; xu = 4700; yl = 1400; yu = 2500;         % axes
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20);   % eta label
text(x0+1.2*xF,y0+1.2*yF,'F','FontSize',20);            % F label
grid; xlabel('Annual electricity consumption (TWh per year)');
ylabel('Annual CO2 emissions (millions of mtons per year)');
legend('Expected unconstrained choice','Unconstrained choice 95% c.i.',...
    'Social optimum 95% c.i.',...
    'Optimal flat standard','Optimal intensity standard',...
    'Optimal indexed standard','Location','southeast');
hold off
h = gcf;
set(gca,'FontSize',20);
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
outfile = strcat(fdir, 'USElectricityPlot_doubleuncert.pdf');
print(outfile,'-dpdf','-fillpage');
clear h


% Plot results -- no regs
plotgreen = 1/255*[0 153 51];
clf
plot(EQ,EE,'xk','MarkerSize',10); hold on           % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',2.5); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',2.5); % social opt
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20);   % eta label
text(x0+1.2*xF,y0+1.2*yF,'F','FontSize',20);            % F label
grid; xlabel('Annual electricity consumption (TWh per year)');
ylabel('Annual CO2 emissions (millions of mtons per year)');
legend('Expected unconstrained choice','Unconstrained choice 95% c.i.',...
    'Social optimum 95% c.i.',...
    'Location','southeast');
hold off
h = gcf;
set(gca,'FontSize',20);
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
outfile = strcat(fdir, 'USElectricityPlot_doubleuncert_noregs.pdf');
print(outfile,'-dpdf','-fillpage');
clear h

% Plot results formatted for paper
clf
formatSpec = '$%2.2f';
str0 = [' Expected unconstrained choice: ', sprintf('$%1.0f',0)];
str1 = [' Unconstrained choices:               ', sprintf('$%1.0f',0)];
str2 = [' Socially optimal choices:            ',sprintf(formatSpec,EWego/1E9)];
str3 = [' Optimal flat std:                          ',sprintf(formatSpec,EWeg0/1E9)];
str4 = [' Optimal intensity std:                  ',sprintf(formatSpec,EWegi/1E9)];
str5 = [' Optimal output-based std:           ',sprintf(formatSpec,EWegs/1E9)];
plotgreen = 1/255*[0 153 51];
plot(EQ,EE,'xk','MarkerSize',10,'DisplayName',str0); hold on       % expected Q,E
plot(Q.Qu,E.Eu,'k','LineWidth',1.5,'DisplayName',str1); hold on      % unconstrained
plot(Q.Qo,E.Eo,'Color',plotgreen,'LineStyle','--','LineWidth',1.5,'DisplayName',str2); % social opt
plot(Q.Q0,E.E0,':k','LineWidth',5,'DisplayName',str3);                 % flat std
plot(Q.Qi,E.Ei,'-.b','LineWidth',5,'DisplayName',str4);                % intensity std
plot(Q.Qs,E.Es,'-r','LineWidth',5,'DisplayName',str5);                 % optimal indexed std
axis([xl xu yl yu]);
pos = get(gca, 'Position');                         % vector of window [left bottom width height]
startx = pos(1) + (x0-xl)/(xu-xl)*pos(3); starty = pos(2) + (y0-yl)/(yu-yl)*pos(4); % arrow start
endxe = pos(1) + (x0+xeta-xl)/(xu-xl)*pos(3); endye = pos(2) + (y0+yeta-yl)/(yu-yl)*pos(4); % eta arrow end
endxF = pos(1) + (x0+xF-xl)/(xu-xl)*pos(3); endyF = pos(2) + (y0+yF-yl)/(yu-yl)*pos(4);     % F arrow end
annotation('arrow', [startx endxe], [starty endye]);    % eta arrow
annotation('arrow', [startx endxF], [starty endyF]);    % F arrow
text(x0+1.05*xeta,y0+1.05*yeta,'\eta','FontSize',20,'FontName','Times New Roman');   % eta label
text(x0+1.2*xF,y0+1.2*yF,'F','FontSize',20,'FontName','Times New Roman');            % F label
grid; xlabel('Annual electricity consumption (TWh per year)');
ylabel('Annual CO2 emissions (millions of mtons per year)');
legend('Location','southeast');
lgd = legend;
lgd.FontSize = 20; lgd.Title.String = 'Expected welfare ($billion/yr)';
hold off
h = gcf;
set(gca,'FontSize',20,'FontName','Times New Roman');
set(h,'PaperUnits','inches','PaperType','usletter')
set(h,'PaperOrientation','landscape','PaperPosition', [-0.3 -0.2 11.7 9.1]);
outfile = strcat(fdir, 'USElectricityPlot_doubleuncert_paper.pdf');
print(h,outfile,'-dpdf');
clear h
clear str* formatSpec xl xu yl yu x0 y0 xeta yeta xF yF pos startx starty...
    endxe endye endxF endyF Q E EQ EE obj lgd

% Output single number tex
% EWegi
fid = fopen(strcat(sndir, 'elec_EWegi.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWegi/1e9);
fclose(fid);
% EWeg0
fid = fopen(strcat(sndir, 'elec_EWeg0.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWeg0/1e9);
fclose(fid);
% EWego
fid = fopen(strcat(sndir, 'elec_EWego.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWego/1e9);
fclose(fid);
% opt gamma
fid = fopen(strcat(sndir, 'elec_optgamma.tex'),'w');
fprintf(fid,'%8.2f', GAMMAeg(3));
fclose(fid);
% Ratio of opt gamma to intensity gamma
fid = fopen(strcat(sndir, 'elec_ratio_optgamma_intgamma.tex'),'w');
fprintf(fid,'%8.0f\\%%', GAMMAeg(3) / GAMMAeg(2) * 100);
fclose(fid);
% Ratio of sigma_eta to Eeta
fid = fopen(strcat(sndir, 'elec_ratio_sigmaeta_Eeta.tex'),'w');
fprintf(fid,'%8.1f\\%%', paramseg.sigmaeta / paramseg.Eeta * 100);
fclose(fid);



% US ELECTRICITY MODEL WITH HALF SIGMA ETA
% Instantiate model
obj = indexmodelF(paramseghalf);

% Run model
[EWeghalfo, EWeghalf0, EWeghalfi, EWeghalfs, ~, ~, ~, ~, MUeghalf, GAMMAeghalf] = resultsall(obj);



% US ELECTRICITY MODEL WITH DOUBLE SIGMA ETA
% Instantiate model
obj = indexmodelF(paramsegdoub);

% Run model
[EWegdoubo, EWegdoub0, EWegdoubi, EWegdoubs, ~, ~, ~, ~, MUegdoub, GAMMAegdoub] = resultsall(obj);



% US ELECTRICITY MODEL WITH HIGH DEMAND ELASTICITY
% Instantiate model
obj = indexmodelF(paramsegh);

% Run model
[EWegho, EWegh0, EWeghi, EWeghs, ~, ~, ~, ~, MUegh, GAMMAegh] = resultsall(obj);

% Output single number tex
% half eta, opt gamma
fid = fopen(strcat(sndir, 'elec_eghalf_optgamma.tex'),'w');
fprintf(fid,'%8.2f', GAMMAeghalf(3));
fclose(fid);
% double eta, opt gamma
fid = fopen(strcat(sndir, 'elec_egdoub_optgamma.tex'),'w');
fprintf(fid,'%8.2f', GAMMAegdoub(3));
fclose(fid);
% half eta, EWeg0
fid = fopen(strcat(sndir, 'elec_eghalf_EWeg0.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWeghalf0/1e9);
fclose(fid);
% half eta, EWegi
fid = fopen(strcat(sndir, 'elec_eghalf_EWegi.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWeghalfi/1e9);
fclose(fid);
% half eta, EWegs
fid = fopen(strcat(sndir, 'elec_eghalf_EWegs.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWeghalfs/1e9);
fclose(fid);
% doub eta, EWeg0
fid = fopen(strcat(sndir, 'elec_egdoub_EWeg0.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWegdoub0/1e9);
fclose(fid);
% doub eta, EWegi
fid = fopen(strcat(sndir, 'elec_egdoub_EWegi.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWegdoubi/1e9);
fclose(fid);
% doub eta, EWegs
fid = fopen(strcat(sndir, 'elec_egdoub_EWegs.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWegdoubs/1e9);
fclose(fid);
% high elast, EWeg0
fid = fopen(strcat(sndir, 'elec_egh_EWeg0.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWegh0/1e9);
fclose(fid);
% high elast, EWegi
fid = fopen(strcat(sndir, 'elec_egh_EWegi.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWeghi/1e9);
fclose(fid);
% high elast, EWego
fid = fopen(strcat(sndir, 'elec_egh_EWego.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWegho/1e9);
fclose(fid);
% high elast, opt gamma
fid = fopen(strcat(sndir, 'elec_egh_optgamma.tex'),'w');
fprintf(fid,'%8.2f', GAMMAegh(3));
fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODELS WITH NON-ZERO CORRELATION BETWEEN ETA AND F
% Footprint
obj = indexmodelF(paramsfpr);
[EWfpro, EWfpr0, EWfpri, EWfprs, ~, ~, ~, ~, MUfpr, GAMMAfpr] = resultsall(obj);
% VMT
obj = indexmodelF(paramsfer);
[EWfero, EWfer0, EWferi, EWfers, ~, ~, ~, ~, MUfer, GAMMAfer] = resultsall(obj);
% Electricity
obj = indexmodelF(paramsegr);
[EWegro, EWegr0, EWegri, EWegrs, ~, ~, ~, ~, MUegr, GAMMAegr] = resultsall(obj);
clear obj
% Electricity, high rho
obj = indexmodelF(paramsegrh);
[EWegrho, EWegrh0, EWegrhi, EWegrhs, ~, ~, ~, ~, MUegrh, GAMMAegrh] = resultsall(obj);
clear obj

% Output single number tex
% elec, high rho, EWegrh0
fid = fopen(strcat(sndir, 'elec_egrh_EWeg0.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWegrh0/1e9);
fclose(fid);
% elec, high rho, EWegrhi
fid = fopen(strcat(sndir, 'elec_egrh_EWegi.tex'),'w');
fprintf(fid,'$\\$ %8.2f $', EWegrhi/1e9);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3
% Instantiate and run models with eta uncertainty only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FOOTPRINT MODEL
obj = indexmodel(paramsfp);
[eta_EWfpo, eta_EWfp0, eta_EWfpi, eta_EWfps, ~, ~, eta_MUfp, eta_GAMMAfp]...
    = resultsall(obj);


% FUEL ECONOMY MODEL
obj = indexmodel(paramsfe);
[eta_EWfeo, eta_EWfe0, eta_EWfei, eta_EWfes, ~, ~, eta_MUfe, eta_GAMMAfe]...
    = resultsall(obj);


% US ELECTRICITY MODEL
obj = indexmodel(paramseg);
[eta_EWego, eta_EWeg0, eta_EWegi, eta_EWegs, ~, ~, eta_MUeg, eta_GAMMAeg]...
    = resultsall(obj);

% Output single number tex
% fp opt gamma
fid = fopen(strcat(sndir, 'foot_eta_optgamma.tex'),'w');
fprintf(fid,'%8.3f', eta_GAMMAfp(3));
fclose(fid);
% fe opt gamma
fid = fopen(strcat(sndir, 'vmt_eta_optgamma.tex'),'w');
fprintf(fid,'%8.2f', eta_GAMMAfe(3));
fclose(fid);
% eg opt gamma
fid = fopen(strcat(sndir, 'elec_eta_optgamma.tex'),'w');
fprintf(fid,'%8.2f', eta_GAMMAeg(3));
fclose(fid);
% EWfps
fid = fopen(strcat(sndir, 'foot_eta_EWfps.tex'),'w');
fprintf(fid,'$\\$ %8.0f $', eta_EWfps);
fclose(fid);
% EWfes
fid = fopen(strcat(sndir, 'vmt_eta_EWfes.tex'),'w');
fprintf(fid,'$\\$ %8.0f $', eta_EWfes);
fclose(fid);
% EWegs
fid = fopen(strcat(sndir, 'elec_eta_EWegs.tex'),'w');
fprintf(fid,'$\\$ %8.1f $', eta_EWegs/1e9);
fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4
% Export csv with results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Table of gammas
filenameo = strcat(tdir,'outputs_gamma.csv');
fileid = fopen(char(filenameo),'w');
fprintf(fileid,'Application, Flat, Intensity, Optimal \n');
fprintf(fileid,'MODELS WITH BOTH ETA AND F VOLATILITY, , , \n');
fprintf(fileid,'Footprint, %1.4f, %1.4f, %1.4f \n', GAMMAfp');
fprintf(fileid,'Miles driven, %1.4f, %1.4f, %1.4f \n', GAMMAfe');
fprintf(fileid,'Elec, %1.4f, %1.4f, %1.4f \n', GAMMAeg');
fprintf(fileid,'Footprint with high sigmaF, %1.4f, %1.4f, %1.4f \n', GAMMAfphf');
fprintf(fileid,'Elec with half sigmaeta, %1.4f, %1.4f, %1.4f \n', GAMMAeghalf');
fprintf(fileid,'Elec with double sigmaeta, %1.4f, %1.4f, %1.4f \n', GAMMAegdoub');
fprintf(fileid,'Elec with high demand elast, %1.4f, %1.4f, %1.4f \n', GAMMAegh');
fprintf(fileid,'MODELS WITH NONZERO CORR BETWEEN ETA AND F, , , \n');
fprintf(fileid,'Footprint, %1.4f, %1.4f, %1.4f \n', GAMMAfpr');
fprintf(fileid,'Miles driven, %1.4f, %1.4f, %1.4f \n', GAMMAfer');
fprintf(fileid,'Elec, %1.4f, %1.4f, %1.4f \n', GAMMAegr');
fprintf(fileid,'MODELS WITH ONLY ETA VOLATILITY, , , \n');
fprintf(fileid,'Footprint, %1.4f, %1.4f, %1.4f \n', eta_GAMMAfp');
fprintf(fileid,'Miles driven, %1.4f, %1.4f, %1.4f \n', eta_GAMMAfe');
fprintf(fileid,'Elec, %1.4f, %1.4f, %1.4f \n', eta_GAMMAeg');
fclose('all');

% Table of expected welfare
filenameo = strcat(tdir,'outputs_welfare.csv');
fileid = fopen(char(filenameo),'w');
fprintf(fileid,'Application, Pigou, Flat, Intensity, Optimal \n');
fprintf(fileid,'MODELS WITH BOTH ETA AND F VOLATILITY, , , \n');
fprintf(fileid,'Footprint, %2.4f, %2.4f, %2.4f, %2.4f \n',...
    [EWfpo EWfp0 EWfpi EWfps]');
fprintf(fileid,'Miles driven, %2.4f, %2.4f, %2.4f, %2.4f \n',...
    [EWfeo EWfe0 EWfei EWfes]');
fprintf(fileid,'Elec, %1.4e, %1.4e, %1.4e, %1.4e \n',...
    [EWego EWeg0 EWegi EWegs]');
fprintf(fileid,'Footprint with high sigmaF, %2.4f, %2.4f, %2.4f, %2.4f \n',...
    [EWfphfo EWfphf0 EWfphfi EWfphfs]');
fprintf(fileid,'Elec with half sigmaeta, %1.4e, %1.4e, %1.4e, %1.4e \n',...
    [EWeghalfo EWeghalf0 EWeghalfi EWeghalfs]');
fprintf(fileid,'Elec with double sigmaeta, %1.4e, %1.4e, %1.4e, %1.4e \n',...
    [EWegdoubo EWegdoub0 EWegdoubi EWegdoubs]');
fprintf(fileid,'Elec with high demand elast, %1.4e, %1.4e, %1.4e, %1.4e \n',...
    [EWegho EWegh0 EWeghi EWeghs]');
fprintf(fileid,'MODELS WITH NONZERO CORR BETWEEN ETA AND F, , , \n');
fprintf(fileid,'Footprint, %2.4f, %2.4f, %2.4f, %2.4f \n',...
    [EWfpro EWfpr0 EWfpri EWfprs]');
fprintf(fileid,'Miles driven, %2.4f, %2.4f, %2.4f, %2.4f \n',...
    [EWfero EWfer0 EWferi EWfers]');
fprintf(fileid,'Elec, %1.4e, %1.4e, %1.4e, %1.4e \n',...
    [EWegro EWegr0 EWegri EWegrs]');
fprintf(fileid,'MODELS WITH ONLY ETA VOLATILITY, , , \n');
fprintf(fileid,'Footprint, %2.4f, %2.4f, %2.4f, %2.4f \n',...
    [eta_EWfpo eta_EWfp0 eta_EWfpi eta_EWfps]');
fprintf(fileid,'Miles driven, %2.4f, %2.4f, %2.4f, %2.4f \n',...
    [eta_EWfeo eta_EWfe0 eta_EWfei eta_EWfes]');
fprintf(fileid,'Elec, %1.4e, %1.4e, %1.4e, %1.4e \n',...
    [eta_EWego eta_EWeg0 eta_EWegi eta_EWegs]');
fclose('all');
clear filenameo fileid ans

        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save workspace
outfile = strcat(workspacedir,'matlabworkspace');   
save(outfile);

