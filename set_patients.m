%AUTHOR: BLANCHE MONGEON @ UNIVERSITÉ DE MONTRÉAL
%ARTICLE: Virtual clinical trial reveals significant clinical potential of targeting tumour-associated macrophages and microglia to treat glioblastoma  
%DATE: DECEMBER 3RD, 2024
%File to set general parameter values for the model - see set_treatment.m
%for treatment administration details

%clear p

%NOTE: if ~isfield(p,paramName) assigns a value to paramName if it does not
%already have a value in the parameter structure p.

%% General parameters
p.K = 158040;%carrying capacity of tumour cells (10^6 cells)
p.Vol = p.K.*10^(6)./(0.8.*10.^(9));%Tumour microenvironment volume (mL)
p.volDeath=(4/3).*(30^3).*pi;%Lethal volume (10^6 cells)
p.volDiagnosis = 17700;%Average volume at diagnosis (10^6 cells)

%% Parameters for tumour equation
if ~isfield(p,'beta')
    p.beta = 0.0134;%Tumour intrinsic growth rate (1/day)
end
if ~isfield(p,'alpha')
    p.alpha = 0.0261;%Maximal killing rate of the GBM cells by the CD8+ T cells (1/day)
end

%% Parameters for TAM equations
if ~isfield(p,'lambdaM')
    p.lambdaM = 0.0248;%Daily source of resident TAMs (10^6 cells)
end
p.deltaM = 0.2;%natural death rate of M1 and M2 TAMs (1/day)
p.deltaMR = log(2)./35;%Natural death rate of resident TAMs (1/day)
if ~isfield(p,'aDM')
    p.aDM = 1.1*10^(3)*10^(-3)*(1/p.Vol);%activation rate of resident TAMs by dead (and damaged) tumour cells (1/(10^6 cells * day))
end
p.aMM = 1.1*10^(3)*10^(-3)*(1/p.Vol)*(1/10);%activation rate of resident TAMs by M1 and M2 TAMs (1/(10^6 cells * day))
if ~isfield(p,'muM2')
    p.muM2 = 0.05; %M2-to-M1 repolarisation rate (1/day)
end
p.muM1 = 0.075;%M1-M2 polarisation rate (1/day)
if ~isfield(p,'q')
    p.q = 151/201; %Fraction of TAMs that acquire M2 phenotype upon activate (unit-less)
end

%set to 0 if you want to simulate without inhibition of the T cells, and to 1 to simulate with inhibition
if ~isfield(p,'inh')
    p.inh = 1; 
end

if ~isfield(p,'lambdaT')
    p.lambdaT=0; %CD8+ T cells source - used to increase CD8+ T cell recruitment during treatment (10^6 cells/day)
end

%% Parameters for CD8+ T cells and PD-1/PD-L1
p.kYQ = 1.365.*10.^(-18)*10^(24);%Inhibition of the CD8+ T cell activation by the PD-1/PD-L1 signalling axis (pg^2/ml^2)
p.rhoG = 10*9282*5.8*10^(-20)*(1/p.Vol)*10^(6)*10^(12);%PD-L1 concentration per 10^6 G cells (pg/ml/10^6 cells)
p.rhoM = 10*9282*5.8*10^(-20)*(1/p.Vol)*10^(6)*10^(12);%PD-L1 concentration per 10^6 TAMs (pg/ml/10^6 cells)
p.rhoT = 3096*8.3*10^(-20)*(1/p.Vol)*10^(6)*10^(12);%PD-1 concentration per 10^6 T cells (pg/ml/10^6 cells)
p.deltaT = 0.4;%natural death rate of T cells (1/day)
p.n2=0.1;%Half-saturating constant of the activation of the CD8+ T cells by the innate immune system (unit-less)
p.n1 = 0.0464;%Maximal activation rate of the CD8+ T cells by the innate immune system (10^6 cells/day)
p.hT = 40*p.K./(10^(11));%Half-saturating constant of the activation of the CD8+ T cells by the GBM cells (10^6 cells)
if ~isfield(p,'aAT')
    p.aAT = 0.0375;%Maximal activation rate of the CD8+ T cells by the GBM cells (1/day)
end

%% Parameters for dead and damaged cells
p.dDd = 8;%Dead cells death rate (1/day)
p.deltaDM = 8.03*10^(-3)*(1/p.Vol);%rate TAMs phagocytose dead cells (1/(10^6 cells * day))

%% Parameters for ICI treatment equations
p.deltaI = log(2)./27;%natural clearance of anti-PD1-1 (1/day)
p.muI = p.deltaI/(9*747);%Binding rate of anti-PD-1 to PD-1 receptors (mL/pg/day)
p.muP = p.muI*10^(9);%anti-PD-1 blocking rate of PD-1 receptors (mL/pg/day)

%% Parameters TAM-targeting strategies
p.lambdaS1 = 1;
p.lambdaS2 = 1;
p.lambdaS3 = 0;
if ~isfield(p,'dPMac')
    p.dPMac = 2;%Maximal phagocytosis rate of tumour cells by M1 TAMs (1/(10^6 cells * day))
end
p.lambdaS4v1 = 1;
p.lambdaS4v2 = 1;
p.newLambdaMacValues = [p.lambdaS1 p.lambdaS2 p.lambdaS3 p.lambdaS4v1 p.lambdaS4v2];

%Set to 'true' if you want to model on/off strategies (see set_treatment.m)
if ~isfield(p,'breakMac')
    p.breakMac = false;
end

%% Initial conditions
p.G0 = 5;%initial tumor volume (10^6 cells)
p.totalM = p.G0*391/1562;%total initial macrophages count (10^6 cells) 

if isfield(p,'M1ko') %M1 knockout-out scenario -- all activated TAMs are M2
    p.M20 = max(1*10^(-6),0.0000032*p.totalM);%initial M2 TAMs (10^6 cells)
    p.M10 = 0;%initial M1 TAMs (10^6 cells)
    p.muM2 = 0; %M2 TAMs dont't change phenotype
    p.q = 1; %all resident TAM acquire M2 phenotype upon activation
    p.M2ko = false; %can't simulate both knock-out at the same time
elseif isfield(p,'M2ko') %M2 knockout-out scenario -- all activated TAMs are M1
    p.M1ko = false; %can't simulate both knock-out at the same time
    p.M20 = 0;
    p.M10 = max(1*10^(-6),0.0000032*p.totalM);%
    p.q = 0; %no resident TAM acquire M2 phenotype upon activation
    p.muM1 = 0; %M1 TAMs dont't change phenotype
else
    if(p.q==0) %-> M2 knockout
        p.M20 = 0;
        p.M10 = max(1*10^(-6),0.000001*p.totalM);%
        p.M1ko = false;
        p.M2ko = true;
        p.muM1 = 0;
    elseif(p.q<1)
        p.M20 = max(1*10^(-6),0.000001*(p.q/(1-p.q))*p.totalM);
        p.M10 = max(1*10^(-6),0.000001*p.totalM);
        p.M1ko = false;
        p.M2ko = false;
    else %-> M1 knockout
        p.M10 = 0;
        p.M20 = max(1*10^(-6),0.0000032*p.totalM);
        p.muM2 = 0;
        p.M1ko = true;
        p.M2ko = false;
    end
end
p.M0 = p.totalM-p.M20-p.M10;%initial count of resident TAMs (10^6 cells)
p.T0 = p.G0*10/1562;%initial T cells count (10^6 cells)
p.P0 = p.rhoT*p.T0;%initial PD-1 concentration (pg/ml)
p.D0 = p.rhoG*p.G0+p.rhoM*(p.totalM);%initial PD-L1 concentration (pg/ml)
p.I0 = 0;%initial ICI concentration
p.A10 = 0;%initial concentration of TMZ in GIT
p.A20 = 0;%initial concentration of TMZ in blood
p.A30 = 0;%initial concentration of TMZ
p.Dd0 = 0;%initial dead and damaged tumor cells count
p.Tb0 = 0;%Necessary to keep track of fraction of tumour cells killed by treatment (excluding resection) or by the immune system at every time point
p.Z0=0;%anti-CD47 efficacy

p.initialcondition = [p.M10,p.M20,p.G0,p.T0,p.P0,p.I0,p.D0,p.A10,p.A20,p.A30,p.M0,p.Dd0,p.Tb0,p.Z0,0];

%% Parameters calculated at homeostasis and/or that need the initial conditions
p.k1 = 40.*p.T0;%Half-saturating constant of the cytotoxic activity of the CD8+ T cells (10^6 cells)

%% Parameters for TMZ treatment equations
p.CL = 10*24;%TMZ clearance rate (L/day)
p.Vd = 30.3;%%TMZ Volume of distribution (L)
p.Vp = 140;%Volume of distribution in CSF (brain) (mL);
p.ka = 5.8*24;%TMZ absorption rate (1/day)
p.k23 = 7.2*10^(-4)*24;%TMZ transfer rate from blood to CSF (1/day)
p.k32 = 0.76*24;%TMZ transfer rate from CSF to blood (1/day)
p.ImaxTMZ = 0.30001;%TMZ Maximal effect (1/day)
p.hTMZ = 0.5682;%TMZ Hill coefficient (unit-less)
p.IC50TMZ = 0.0011;%TMZ IC50 (mg/mL)

%% Set time info and variables
p.timepoints = linspace(0,10*365,12000); %5 years
p.tspan = [0 p.timepoints(end)];%set timespan for simulations

%% Parameters to model no treatment administration
%Unless already specified, we simulate lambdas with constant values
if(~isfield(p,'TAMT3decayingEfficacy'))
    p.TAMT3decayingEfficacy=false;
end
%Unless already specified, we decrease death rate of all TAMs
if(~isfield(p,'decreaseDeathRateOnlyM2'))
    p.decreaseDeathRateOnlyM2=false;
end

%These will be overwritten if you run set_treatment;
p.AdministrationTimesICI = [];
p.NumberAdmins = 1;
p.treatmentBoolean = 0;
p.treatmentICI = 0;
p.treatmentRT = 0;
p.treatmentTMZ = 0;
p.treatmentMac = 0;

%% Keep track of original TAM dynamics values
%Useful when doing M1ko or M2 ko during treatment only
p.muM1_baseline = p.muM1;
p.muM2_baseline = p.muM2;
p.q_baseline = p.q; 

