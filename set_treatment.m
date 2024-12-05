%AUTHOR: BLANCHE MONGEON @ UNIVERSITÉ DE MONTRÉAL
%ARTICLE: Virtual clinical trial reveals significant clinical potential of targeting tumour-associated macrophages and microglia to treat glioblastoma  
%DATE: DECEMBER 3RD, 2024
%File to set treatment administration details - see set_patients.m
%for general parameter values

%% General information
%p.treatmentBoolean: 1 if we resect, 0 if not
%p.treatmentRT: 1 if we simulate radiotherapy, 0 if not
%p.treatmentTMZ: 1 if we simulate TMZ chemotherapy, 0 if not
%p.treatmentICI: 1 if we simulate ICI with nivolumab, 0 if not
%p.treatmentMac: 1 if we simulate a TAM-targeting strategy, 0 if not
%p.NumberAdmins: Number of administration - all treatments considered

%% Set treatment administration
% We need to find the time at diagnosis - let the tumour grow without
% treatment and retrieve the first time it reaches the volume at diagnosis
p.treatmentBoolean=0;
p.treatmentICI=0;
p.treatmentRT=0;
p.treatmentTMZ=0;
p.treatmentMac=0;
p.NumberAdmins=1;
p.AdministrationTimesICI = [];

[p.solNo,~] = solver_gbm(p);

%XX
%p.volDiagnosis = 17700;%Average volume at diagnosis (10^6 cells)
%p.volDeath = (4/3)*pi*(30^3);%Lethal volume (10^6 cells)

p.all_popNo = real(deval(p.solNo,p.timepoints));
p.ind_time_diagnosis = find(p.all_popNo(3,:)>=p.volDiagnosis,1,'first');
if isempty(p.ind_time_diagnosis)
    %the tumor never reaches the diagnostic volume: there's no given
    %treatment and so we don't change anything in the parameters
    disp("no treatment will be given")
else
    %the tumor reaches the diagnostic volume: we operate
    p.time_diagnostic = round(p.timepoints(p.ind_time_diagnosis));
    %retrieve info for the initial conditions -- values at diagnosis
    p.initialInfo = p.all_popNo(:,p.ind_time_diagnosis);%
    p.initialcondition = p.initialInfo;%

    %resection: 3 weeks=21 days after diagnosis
    p.resectionTime = 21;%
    p.resectionExtent = 0.8;

    %TMZ Administration
    %Cycle 0 to the TMZ administration concumitant with radiotherapy, when TMZ
    %is administered daily at a dose of 75mg/m^2. It is followed by a four week
    %break.
    %TMZ is then administered on days 1-5 of 6 28-days cycles, at a dose of
    %150mg/m^2 for the first cycle, then 200mg/m^2 for the following cycles.
    
    %First chemotherapy: 21 days after resection
    p.TimeFirstAdminTMZ = p.resectionTime+21;
    p.cycle0Length = 6;%in weeks
    p.cycle0Length = p.cycle0Length*7;%in days;
    %daily administration during cycle 0
    p.AdministrationTimesTMZ = linspace(p.TimeFirstAdminTMZ,p.TimeFirstAdminTMZ+p.cycle0Length-1,6*7);
    p.TMZBreakLength= 4;%break between cycle 0 and cycle 1 (weeks)
    p.lastAdminTMZ = p.AdministrationTimesTMZ(end);
    p.numberCyclesTMZ = 6;
    p.numberAdminPerTMZCycle = 5;
    p.i=1;
    while(p.i<=p.numberCyclesTMZ)
        p.j=1;
        while(p.j<=p.numberAdminPerTMZCycle)
            p.AdministrationTimesTMZ(end+1) = p.lastAdminTMZ+7*p.TMZBreakLength+(p.i-1)*28+p.j;
            p.j = p.j + 1;
        end
        p.i = p.i+1;
    end
    
    p.BSA = 1.7; %average body surface area (m^2)
    p.DoseTMZ = 75*p.BSA; %total dose in mg
    %Following two parameters are necessary to increase TMZ dose during
    %treatment
    p.currentCycle = 0; 
    p.currentAdminInTMZCycle = 0;
    if any(p.TimeFirstAdminTMZ == 0)
        p.A10 = p.DoseTMZ;
    end
    
    % Dosing information -- ICI
    p.TimeFirstAdminICI = p.resectionTime+21; %time of first administration - same as TMZ
    p.DoseICI = 0.02*3+0.0092;
    p.AdministrationTimesICI = [];
    %First we have 8 doses administered every two weeks
    p.numberFirstDosesICI = 8;
    p.timeBetweenICI = 2;%in weeks
    p.i = 1;
    while(p.i <= p.numberFirstDosesICI)
        p.AdministrationTimesICI(end+1) = p.TimeFirstAdminICI+(p.i-1)*7*p.timeBetweenICI;
        p.i = p.i+1;
    end
    %Doses are then increased and administered every 4 weeks
    p.timeBetweenICI = 4;%in weeks
    while(p.AdministrationTimesICI(end)+7*p.timeBetweenICI+1 < p.AdministrationTimesTMZ(end))
        p.AdministrationTimesICI(end+1) = p.AdministrationTimesICI(end)+7*p.timeBetweenICI+1;
    end
    p.currentAdminICI = 1;
    
    % Radiotherapy regimen
    p.DoseRT = 2; %in Gy
    p.alphaRT = 0.19;%
    p.betaRT = 0.032;%
    p.RTDuration = 6;%in weeks
    p.TimeFirstAdminRT = p.resectionTime+21;
    p.AdministrationTimesRT = [];
    p.i=1;
    while(p.i<=p.RTDuration)
        p.j=1;
        while(p.j<=5)
            p.AdministrationTimesRT(end+1) = p.TimeFirstAdminRT+(p.i-1)*7+(p.j-1);
            p.j = p.j+1;
        end
        p.i = p.i+1;
    end
    p.givenRTDose = 0;
    
    %Targeting macrophages
    p.TAMT3efficacy = 1; %Maximal efficacy of anti-CD47
    %WITH LAMBDAS - CONSTANT EFFICACY VALUES 
    %let i, j and k be 3 times in p.AdministrationTimesMacrophage with i
    %and k even, an j odd. Then [i,j] is a 'on-cycle' during which lambda
    %is set to its treatment value and [j,k] is a 'off-cycle' during which
    %lambda is set to its no treatment value. 
    %WITH DECAYING EFFICACY (FOR CLINICAL POTENTIAL)
    %p.AdministrationTimesMacrophages contains the time at which we
    %simulate administration of anti-CD47 and reset the efficacy to 100%

    %start during chemo and radio
    %p.AdministrationTimesMacrophages = [p.AdministrationTimesTMZ(1)];
    %start during adjuvant chemo
    p.AdministrationTimesMacrophages = [p.AdministrationTimesTMZ(43)];

    %If p.breakMac=false, then we only reset lambda to its no treatment
    %value after TMZ has ended (not applicable with decaying efficacy)
    if(~p.breakMac)
        p.nPeriodMac = p.AdministrationTimesTMZ(end)-p.AdministrationTimesTMZ(43);
    elseif ~isfield(p,'nPeriodMac')
        %nPeriodMac: length of on- and off-cycles (constant values) or
        %length between two administrations (decaying efficacy)
        p.nPeriodMac = 7;
    end
    if ~isfield(p,'addMonthsMacTreatment')
        %addMonthsMacTreatment: number of months after SOC for which we
        %simulate the TAM-targeting treatment going on
        p.addMonthsMacTreatment = 0;
    end

    if(p.breakMac)
        while(p.AdministrationTimesMacrophages(end)<(p.AdministrationTimesTMZ(end)+p.addMonthsMacTreatment*30.417))
            p.AdministrationTimesMacrophages(end+1) = p.AdministrationTimesMacrophages(end)+p.nPeriodMac;
        end
    else
        p.AdministrationTimesMacrophages(end+1) = p.AdministrationTimesMacrophages(end)+p.nPeriodMac+p.addMonthsMacTreatment*30.417;
    end

    %XX
    %make sure we end on an off week which is equivalent to having an even
    %number of weeks in the administration times vector
    p.nWeekMacs = length(p.AdministrationTimesMacrophages);
    if(mod(p.nWeekMacs,2)==1)
        %add off cycle at the end
        %p.AdministrationTimesMacrophages(end+1) = p.AdministrationTimesMacrophages(end)+p.nPeriodMac;
        %remove last on cycle
        %p.AdministrationTimesMacrophages(end) = []; 
    end
    p.CurrentWeekMacrophages = 1;

    % Combining the treatments
    p.AdministrationTimesResection = [p.resectionTime];
    p.AdministrationTimes = unique(cat(2,0,p.AdministrationTimesTMZ,p.AdministrationTimesICI,...
        p.AdministrationTimesRT,p.AdministrationTimesResection,p.AdministrationTimesMacrophages),'sorted');
    p.NumberAdmins = length(p.AdministrationTimes);
    p.treatmentBoolean =1;
    p.treatmentICI=0;
    p.treatmentRT=1;
    p.treatmentTMZ=1;
    p.treatmentMac = 0;
  
end


