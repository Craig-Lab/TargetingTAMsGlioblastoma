%AUTHOR: BLANCHE MONGEON @ UNIVERSITÉ DE MONTRÉAL
%ARTICLE: Virtual clinical trial reveals significant clinical potential of targeting tumour-associated macrophages and microglia to treat glioblastoma  
%DATE: DECEMBER 3RD, 2024
%File to reproduce figures in the article

%Run sections 1 to 4, then all the subsequent sections can be run
%independently 

%% SECTION 1
close all
clearvars

%% SECTION 2: Get folders

folder_code = cd('VCT');
folder_VCT = cd(folder_code);

cd('eFAST_GBM/')
folder_efast = cd(folder_code);

cd('MSurv/')
folder_MatSurv = cd(folder_code);

%% SECTION 3: set colors

color_paleblue = '#a6cee3';%
color_blue = '#1f78b4';%
color_palegreen = '#b2df8a';%
color_green = '#33a02c';%
color_palered = '#fb9a99';
color_red = '#e31a1c';%
color_paleorange = '#fdbf6f';
color_orange = '#ff7f00';
color_palepurple = '#cab2d6';%
color_purple = '#6a3d9a';%
color_yellow = '#ffff99';%'
color_brown = '#b15928';
color_black = '#000000';
color_white = '#ffffff';
color_immunotherapy = '#C98473';

%% SECTION 4: set constants used for all panels

time = linspace(0,2*365,10000);
timeVCT = linspace(0,2*365,20000);
options.alpha      = 0.5;
options.line_width = 3;
options.error      = 'std';
options.x_axis = time;

col1 = hex2rgb('#b2182b');%red
col2 = hex2rgb('#2166ac');%blue
colours = interp1( [0 1], [col1;col2], linspace(0,1,4) );

%smallest response/benefit group
optionscat1 = options;
optionscat1.color_line = colours(1,:);%red

%second smallest response/benefit group
optionscat2 = options;
optionscat2.color_line = colours(2,:);%

%second best response/benefit group
optionscat3 = options;
optionscat3.color_line = colours(3,:);%

%best response/benefit group
optionscat4 = options;
optionscat4.color_line = colours(4,:);%blue

cd(folder_code)
clear p;
set_patients;
set_treatment;
lastTMZadmin = p.AdministrationTimesTMZ(end);
startTxTime = p.AdministrationTimes(3);
endTxTime = p.AdministrationTimes(end);
cd(folder_code)
baselineBeta = p.beta;%
baselineaAT = p.aAT;%
baselineq = p.q;%
baselineratio = p.q/(1-p.q);%
baselineCL = p.CL;
baselinek23 = p.k23;
baselineka = p.ka;
baselineParam = [baselineCL baselinek23 baselineka baselineBeta baselineaAT 100.*baselineq 0];
volDeath = p.volDeath;

axis_fontsize = 10;%32;
axis_linewidth = 1.5;%3;
plot_linewidth = 2;
legend_fontsize = 8;

cd(folder_VCT)
load('VCT_SOC.mat')
tm = squeeze(VCT_SOC(:,3,:));
%get tumour cell count after SOC
tmpostSOC = tm(:,find(timeVCT>=(p.AdministrationTimesTMZ(end)+7),1,'first'))';

%get tumour growth during treatment period but without treatment
load('VCT_none.mat')
growth = squeeze(VCT_none(:,3,find(timeVCT>=(p.AdministrationTimesTMZ(end)+7),1,'first')))-squeeze(VCT_none(:,3,1));
cd(folder_code)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PANEL 2 and Supplementary Figure 2
        
    clear p;
    set_patients;
    set_treatment;
    [solSOC,~] = solver_gbm(p);
    allpopSOC = deval(solSOC,time);
    tumourSOC = allpopSOC(3,:);
    
    clear p;
    set_patients;
    set_treatment;
    p.treatmentICI = 1;
    [solICI,~] = solver_gbm(p);
    allpopICI = deval(solICI,time);
    tumourICI = allpopICI(3,:);

    clear p;
    set_patients;
    lambdaT = 10.*p.initialcondition(4);
    set_treatment;
    p.lambdaT = lambdaT;
    p.treatmentICI = 1;
    [solICIincreasedCD8,pincreasedCD8] = solver_gbm(p);
    allpopICIincreasedCD8 = deval(solICIincreasedCD8,time);
    tumourICIincreasedCD8 = allpopICIincreasedCD8(3,:);
    
    %Figure 2B
    h_fig = figure;
    hold on
    plot(time./30.417,tumourSOC,'LineWidth',plot_linewidth,'Color',color_blue)
    plot(time./30.417,tumourICI,'LineWidth',plot_linewidth,'Color',color_green)
    plot(time./30.417,tumourICIincreasedCD8,'LineWidth',plot_linewidth,'Color',color_purple)
    yline(p.volDeath,'LineWidth',plot_linewidth,'LineStyle','-','Color',color_red)
    hold off
    xlabel('Time (months)')
    xlim([0 24])
    ylabel('Tumour cells (10^6 cells)')
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    legend({'SOC','SOC+ICI','SOC+ICI+ \uparrow CD8+ T cells'},'Location','northoutside')

    %close up
    h_fig = figure;
    hold on
    plot(time./30.417,tumourSOC,'LineWidth',plot_linewidth,'Color',color_blue)
    plot(time./30.417,tumourICI,'LineWidth',plot_linewidth,'Color',color_green)
    yline(p.volDeath,'LineWidth',plot_linewidth,'LineStyle','-','Color',color_red)
    hold off
    yticks([])
    xlabel('')
    xlim([14 14.25])
    xticks([14 14.25])
    ylabel('')
    set(gca,'FontSize',10,'LineWidth',axis_linewidth)

    %find times of death (in months after diagnosis)
    deathSOC = time(find(tumourSOC>=p.volDeath,1,'first'))./30.417
    deathICI = time(find(tumourICI>=p.volDeath,1,'first'))./30.417
    deathdiff = deathICI-deathSOC
    
    % FIGURE 2C

    clear p;
    set_patients;
    set_treatment;
    p.treatmentMac = 0;
    p.treatmentICI = 0;
    [solSOC,~] = solver_gbm(p);
    PUSOC = deval(solSOC,time,5);
    PBSOC = deval(solSOC,time,15);
    ISOC = deval(solSOC,time,6);
    DSOC = deval(solSOC,time,7);
    TSOC = deval(solSOC,time,4);
    inhSOC = 1./(1+(PUSOC.*DSOC)./p.kYQ);
    
    clear p;
    set_patients;
    set_treatment;
    p.treatmentMac = 0;
    p.treatmentICI = 1;
    [solSOCICI,~] = solver_gbm(p);
    PUSOCICI = deval(solSOCICI,time,5);
    PBSOCICI = deval(solSOCICI,time,15);
    ISOCICI = deval(solSOCICI,time,6);
    DSOCICI = deval(solSOCICI,time,7);
    TSOCICI = deval(solSOCICI,time,4);
    inhSOCICI = 1./(1+(PUSOCICI.*DSOCICI)./p.kYQ);
    
    startICI = p.AdministrationTimesICI(1);
    endICI = p.AdministrationTimesICI(end);
    
    h_fig = figure;
    hold on
    plot(time./30.417,inhSOC.*100,'LineWidth',plot_linewidth,'Color',color_blue)
    plot(time./30.417,inhSOCICI.*100,'LineWidth',plot_linewidth,'Color',color_green)
    ylabel('T cell activation (%)')
    x_points = [(p.time_diagnostic+startICI)./30.417 (p.time_diagnostic+startICI)./30.417 (p.time_diagnostic+endICI)./30.417 (p.time_diagnostic+endICI)./30.417];
    yl = ylim;
    y_points = [yl(1) yl(2) yl(2) yl(1)]; 
    fill(x_points, y_points,hex2rgb(color_green),'FaceAlpha',0.4,'LineStyle','none');
    xlim([0 24])
    xlabel('Time (months)')
    legend({'SOC','SOC+ICI','ICI treatment period'},'Location','best')
    h = get(gca,'Children');
    newPlotOrder = h(2:end);
    newPlotOrder = [newPlotOrder ; h(1)];
    set(gca,'Children',newPlotOrder)
    ylim(yl)
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)

    %Effect of TME bias towards M2 phenotype on T cell activation
    highratio = 1;%all resident TAMs get M2 phenotype upon activation
    midratio = 3.2./(1+3.2);%average bias
    lowratio = 0;%no resident TAMs get M2 phenotype upon activation
    ratio = [highratio midratio lowratio];
    highRatioColor = color_purple;
    midRatioColor = color_green;
    lowRatioColor = color_blue;
    colors = [hex2rgb(highRatioColor);hex2rgb(midRatioColor);hex2rgb(lowRatioColor)]; 
    
    %initialise metric arrays 
    tempTime = linspace(0,3*365,20000);
    resultsSOCICI = zeros(length(ratio),14,length(tempTime));
    resultsSOC = zeros(length(ratio),14,length(tempTime));
    actCD8byMac = zeros(length(ratio),length(tempTime));
    inhibition = zeros(length(ratio),length(tempTime));
    actCD8byMac_SOC = zeros(length(ratio),length(tempTime));
    inhibition_SOC = zeros(length(ratio),length(tempTime));
    Dkill_SOCICI = zeros(length(ratio),length(tempTime));
    seICI = zeros(length(ratio),1);
    survivalSOCICI = zeros(length(ratio),1);

    %set to 'false' for complete KO starting at cancer initiation
    koDuringTreatmentOnly = true; 
    for nQ=1:length(ratio)
        clear p

        if(koDuringTreatmentOnly)
            set_patients;
            set_treatment;
        end
        if(nQ==1)
            if(koDuringTreatmentOnly)
                p.M1koTX = true; 
                p.q = 1;
                p.muM2 = 0;
                p.initialcondition(7) = p.initialcondition(7)-p.rhoM.*p.initialcondition(1);
                p.initialcondition(1) = 0;
            end
            p.M1ko = true;
        elseif(nQ==3)
            if(koDuringTreatmentOnly)
                p.M2koTX = true; 
                p.q = 0;
                p.muM1 = 0;
                p.initialcondition(7) = p.initialcondition(7)-p.rhoM.*p.initialcondition(2);
                p.initialcondition(2) = 0;
            end
            p.M2ko = true;
        end
        if(~koDuringTreatmentOnly)
            set_patients;
            set_treatment;
        end

        p.treatmentICI = 1;
        p.treatmentMac = 0;
        [sol,p] = solver_gbm(p);
        allpop = deval(sol,tempTime);
        burdenSOCICI = allpop(3,find(tempTime>=256+7,1,'first'));
        allpop(14,:) = [];
        resultsSOCICI(nQ,:,:) = allpop;
        ratioM1 = allpop(1,:)./(allpop(1,:)+allpop(2,:));
        ratioM2 = allpop(2,:)./(allpop(1,:)+allpop(2,:));
        inh = (1./(1+((allpop(5,:).*(allpop(7,:)))./p.kYQ)));
        actCD8byMac_temp = inh.*p.n1.*(ratioM1./(p.n2+ratioM2));
        actCD8byMac(nQ,:) = actCD8byMac_temp;
        inhibition(nQ,:) = inh;
        Dkill_SOCICI(nQ,:) = p.alpha.*allpop(4,:).*allpop(3,:)./(allpop(4,:)+p.k1);
        indDeath = find(allpop(3,:)>=p.volDeath,1,'first');
        if(~isempty(indDeath))
            survivalSOCICI(nQ,1) = tempTime(indDeath)./30.417;
        else
            survivalSOCICI(nQ,1) = -1;
        end

        clear p

        if(koDuringTreatmentOnly)
            set_patients;
            set_treatment;
        end
        if(nQ==1)
            if(koDuringTreatmentOnly)
                p.M1koTX = true; 
                p.q = 1;
                p.muM2 = 0;
                p.initialcondition(7) = p.initialcondition(7)-p.rhoM.*p.initialcondition(1);
                p.initialcondition(1) = 0;
            end
            p.M1ko = true;
        elseif(nQ==3)
            if(koDuringTreatmentOnly)
                p.M2koTX = true; 
                p.q = 0;
                p.muM1 = 0;
                p.initialcondition(7) = p.initialcondition(7)-p.rhoM.*p.initialcondition(2);
                p.initialcondition(2) = 0;
            end
            p.M2ko = true;
        end
        if(~koDuringTreatmentOnly)
            set_patients;
            set_treatment;
        end

        p.treatmentICI = 0;
        p.treatmentMac = 0;
        [sol,~] = solver_gbm(p);
        allpop = deval(sol,tempTime);
        allpop(14,:) = [];
        resultsSOC(nQ,:,:) = allpop;
        burdenSOC = allpop(3,find(tempTime>=256+7,1,'first'));
        seICI(nQ,1) = -100.*(burdenSOCICI-burdenSOC)./burdenSOC;
        ratioM1 = allpop(1,:)./(allpop(1,:)+allpop(2,:));
        ratioM2 = allpop(2,:)./(allpop(1,:)+allpop(2,:));
        inh = (1./(1+((allpop(5,:).*(allpop(7,:)))./p.kYQ)));
        actCD8byMac_temp = inh.*p.n1.*(ratioM1./(p.n2+ratioM2));
        actCD8byMac_SOC(nQ,:) = actCD8byMac_temp;
        inhibition_SOC(nQ,:) = inh;
    end

    %SUPPLEMENTARY FIGURE 3
    inset=true;
    h_fig = figure;
    hold on
    if(~inset)
        for nQ=1:length(ratio)
            %PLOT TUMOUR CELLS
            %plot(tempTime./30.417,squeeze(resultsSOCICI(nQ,3,:)),'LineWidth',plot_linewidth,'Color',colors(nQ,:),'LineStyle','-')
            %PLOT ACTIVATION OF CD8+ T CELLS BY TAMs
            %plot(tempTime./30.417,actCD8byMac(nQ,:),'LineWidth',3,'Color',colors(nQ,:),'LineStyle','-')
            %PLOT CYTOTOXIC ACTIVITY OF CD8+ T CELLS
            plot(tempTime./30.417,Dkill_SOCICI(nQ,:),'LineWidth',plot_linewidth,'Color',colors(nQ,:),'LineStyle','-')
        end
    else
        for nQ=1:1
            %PLOT TUMOUR CELLS
            %plot(tempTime./30.417,squeeze(resultsSOCICI(nQ,3,:)),'LineWidth',plot_linewidth,'Color',colors(nQ,:),'LineStyle','-')
            %PLOT ACTIVATION OF CD8+ T CELLS BY TAMs
            %plot(tempTime./30.417,actCD8byMac(nQ,:),'LineWidth',3,'Color',colors(nQ,:),'LineStyle','-')
            %PLOT CYTOTOXIC ACTIVITY OF CD8+ T CELLS
            plot(tempTime./30.417,Dkill_SOCICI(nQ,:),'LineWidth',plot_linewidth,'Color',colors(nQ,:),'LineStyle','-')
        end
    end

    %PLOT TUMOUR CELLS
    %ylabel({'Tumour cells', '(10^6 cells)'})
    %PLOT ACTIVATION OF CD8+ T CELLS BY TAMs
    %ylabel({'CD8+ cells', '(10^6 cells)'})
    %PLOT CYTOTOXIC ACTIVITY OF CD8+ T CELLS
    ylabel({'Tumour cells', '(10^6 cells)'})

    %uncomment following four lines if you plot tumour cells -- comment out
    %if not
    %x2 = [tempTime./30.417 fliplr(tempTime./30.417)];
    %inBetween = [squeeze(resultsSOCICI(3,3,:))', fliplr(squeeze(resultsSOCICI(1,3,:))')];
    %yline(p.volDeath,'LineWidth',plot_linewidth,'LineStyle','-','Color',color_red)
    %fill(x2, inBetween, hex2rgb('#808080'),'FaceAlpha',0.15,'EdgeColor','none');

    xlabel('Time (months)')
    %legend({'M1 TAM knockout','GBM Average M2:M1 ratio','M2 TAM knockout','Lethality threshold'},'Location','northoutside','NumColumns',4,'FontSize',legend_fontsize)
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    xlim([0 36])
    h = get(gca,'Children');
    newPlotOrder = [];
    for nPlot=1:length(h)
        newPlotOrder(end+1) = h(end-nPlot+1);
    end
    set(gca,'Children',newPlotOrder)

    % FIGURE 2D
    hfig = figure;
    barh(log2(cat(1,seICI(1,1),seICI(3,1))./(seICI(2))),'FaceColor',color_blue)
    xline(0,'LineWidth',plot_linewidth)
    box off
    xlabel({'log_2 fold-change in','ICI efficacy'})
    set(gca,'ytick',[])
    set(gca,'ycolor',[1 1 1])
    xlim([-7 7])
    set(hfig,'Units','Centimeters')
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
 
    % find time of death
    deathOnlyM1 = time(find(squeeze(resultsSOC(3,3,:))>=p.volDeath,1,'first'))./30.417
    deathOnlyM2 = time(find(squeeze(resultsSOC(1,3,:))>=p.volDeath,1,'first'))./30.417
    deathDiffM2M1 = deathOnlyM1-deathOnlyM2
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOC+ICI VCT

    cd(folder_VCT)

    load('VCTPatients.mat')%
    N = size(VCTPatients,1); %nbr of patients
    VCTPatients = VCTPatients(1:N,:);
    VCTPatients(:,6) = VCTPatients(:,6).*100;

    load('VCT_SOC_ICI.mat')%
    load('VCT_SOC.mat')%
    load('deathsSOCICIVCT.mat')%
    
    benefits = deathsVCT(3,:)-deathsVCT(2,:);
    benefits = benefits./30.417; %in months

    cd(folder_code)

    Cl = VCTPatients(:,1);
    k23 = VCTPatients(:,2);
    ka = VCTPatients(:,3);
    beta = VCTPatients(:,4);
    aAT = VCTPatients(:,5);
    q = VCTPatients(:,6);
    resection = VCTPatients(:,7);
    M2M1ratio = q./(1-q);
    
    Variables = {'M1','M2','Cancer','T Cells','Unbound PD-1','ICI','PD-L1','GIT TMZ',...
        'Plasma TMZ','CSF TMZ','Resident Macrophages','Dead cells','Fraction of Tumour Cells Killed','Bounded PD-1'};
    Variables_short = {'M1','M2','G','T','PU','I','D','A1','A2','A3','M','Dd','Tb','PB'};
    
    %Divide virtual patients in 4 groups depending on their benefit gained
    %with SOC+ICI combination treatment
    Q1B = quantile(benefits,0.25);
    Q2B = quantile(benefits,0.5);
    Q3B = quantile(benefits,0.75);
    
    indCat1B = find(benefits(1,:)<=Q1B); 
    indCat2B = find(benefits(1,:)<=Q2B & benefits(1,:)>Q1B); 
    indCat3B = find(benefits(1,:)<=Q3B & benefits(1,:)>Q2B); 
    indCat4B = find(benefits(1,:)>Q3B); 

    VCT_cat1B = VCT_SOC_ICI(indCat1B,:,:);
    VCT_cat2B = VCT_SOC_ICI(indCat2B,:,:);
    VCT_cat3B = VCT_SOC_ICI(indCat3B,:,:);
    VCT_cat4B = VCT_SOC_ICI(indCat4B,:,:);
    
    patients1B = VCTPatients(indCat1B,:);
    patients2B = VCTPatients(indCat2B,:);
    patients3B = VCTPatients(indCat3B,:);
    patients4B = VCTPatients(indCat4B,:);
    
    parameters = ["CL","k_{23}","k_a","\beta","\alpha_{AT}","M2:M1 ratio","Resection extent"];
    labelsParam = {'CL (L/day)','k_{23} (1/day)','k_a (1/day)',...
        '\beta (1/day)','\alpha_{AT} (1/day)','q (unit-less)',''};
    titlesParam = {'TMZ Clearance','TMZ Plasma-CSF transfer rate','TMZ Absorption rate',...
        'Intrinsic tumour growth rate','Antigenicity of tumour cells','TME bias towards M2 phenotype',''};
    
    clear p;
    set_patients;

    nFig = 1;
    %p-values for the two-sided Wilcoxon rank sum test that compares worst
    %and best responders
    pvalsSOCICIVCT = zeros(1,(length(parameters)-1));
    
    for nParam=1:(length(parameters)-1)
        [pval,~] = ranksum(patients1B(:,nParam),patients4B(:,nParam));
        pvalsSOCICIVCT(1,nParam) = pval;
        h_fig = figure;
        hold on

        hgCat1 = histogram(log(patients1B(:,nParam)));
        hgCat1.FaceColor = optionscat1.color_line;
        hgCat1.LineWidth = 1;
        hgCat1.FaceAlpha = 0.9;
        hgCat1.EdgeAlpha = 0.5;

        hgCat4 = histogram(log(patients4B(:,nParam)));
        hgCat4.FaceColor = optionscat4.color_line;
        hgCat4.LineWidth = 1;
        hgCat4.FaceAlpha = 0.6;
        hgCat4.EdgeAlpha = 0.5;

        yl = ylim;
        if(pval<0.01)
            plot_pvalueline(hgCat1,hgCat4,yl,plot_linewidth)
        end

        xlabel(labelsParam{nParam})
        ylabel('Number of patients')
        title(titlesParam{nParam})
        leg = legend({'Worst responders','Best responders'});
        set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
        nFig = nFig+1;

    end

    %Kaplan-Meier curves for four surviving groups
    fig = figure;
    hold on
    KM_curve(deathsVCT(3,indCat1B)./30.417,'LineColor',optionscat1.color_line,'LineWidth',optionscat1.line_width,'ShowMedian',false)
    KM_curve(deathsVCT(3,indCat2B)./30.417,'LineColor',optionscat2.color_line,'LineWidth',optionscat2.line_width,'ShowMedian',false)
    KM_curve(deathsVCT(3,indCat3B)./30.417,'LineColor',optionscat3.color_line,'LineWidth',optionscat3.line_width,'ShowMedian',false)
    KM_curve(deathsVCT(3,indCat4B)./30.417,'LineColor',optionscat4.color_line,'LineWidth',optionscat4.line_width,'ShowMedian',false)
    hold off
    xlabel('Time post diagnosis (months)')
    ylabel('Percent Survival (%)')
    xl = xlim;
    xlim([12 24])%xl(2)])
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)

    %C CD8+ T cells activation by macrophages - Best versus worst
    %responders
    clear p;
    set_patients;
    timePoints = linspace(0,2*365,20000);
    optionscat1.x_axis = timePoints;
    optionscat4.x_axis = timePoints;
    
    actMac = zeros(N,length(timePoints));
    actTumour = zeros(N,length(timePoints));
    
    for nPx=1:N
        M1 = squeeze(VCT_SOC_ICI(nPx,1,:));
        M2 = squeeze(VCT_SOC_ICI(nPx,2,:));
        ratioM1 = M1./(M1+M2);
        ratioM2 = M2./(M1+M2);
        PU = squeeze(VCT_SOC_ICI(nPx,5,:));
        D = squeeze(VCT_SOC_ICI(nPx,7,:));
        G = squeeze(VCT_SOC_ICI(nPx,3,:));
        T = squeeze(VCT_SOC_ICI(nPx,4,:));
        inh = 1./(1+(PU.*D)./p.kYQ);
        actMacTEMP = inh.*p.n1.*(ratioM1./(p.n2+ratioM2));
        actMac(nPx,:) = actMacTEMP;
        actTumourTEMP = T.*inh.*aAT(nPx).*G./(G+p.hT);
        actTumour(nPx,:) = actTumourTEMP;
    end
    
    actTumourCat1B = actTumour(indCat1B,:);
    actTumourCat4B = actTumour(indCat4B,:);
    
    actMacCat1B = actMac(indCat1B,:);
    actMacCat4B = actMac(indCat4B,:);
    
    fig = figure;
    hold on
    plot_areaerrorbar(actMacCat1B,optionscat1)
    plot_areaerrorbar(actMacCat4B,optionscat4)
    xlabel('Time post diagnosis (months)')
    ylabel({'CD8+ T cells (10^6 cells)'})
    title({'Activation of CD8+ T cells','by macrophages'})
    xlim([0 24*30.417])
    xticks([0 6 12 18 24].*30.417)
    xticklabels({'0','6','12','18','24'})
    legend({'','Worst responders','','Best responders'})
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUUPLEMENTARY FIGURE 6 -- Heatmaps 

    %First row: SOC
    %Second row: SOC+ICI
    %Left: No TAM-targeting treatment
    %right: Maximal efficacy

    %Strategy 1 - Increase death rate of TAMs
    clear p;
    set_patients;
    deltaM = p.deltaM;
    lambdaS1Values = linspace(1,1./deltaM,20);
    load('strategy1Results.mat')
    strategy1Results(2,:) = [];
    strategy1Results_normd = plot_HeatmapStrategy(strategy1Results,1,lambdaS1Values,...
        '1',axis_linewidth,axis_fontsize,'Strategy1');

    %Strategy 2 - Inhibit recruitment of TAMs
    lambdaS2Values = linspace(0,1,20);
    load('strategy2Results.mat')
    strategy2Results(2,:) = [];
    strategy2Results_normd = plot_HeatmapStrategy(strategy2Results,length(lambdaS2Values),lambdaS2Values,'2',...
        axis_linewidth,axis_fontsize,'Strategy2');

    %colormap Strategies 1 and 2
    minLimColorbar = -0.3;%
    maxLimColorbar = 1.55;%
    breakpoint1 = -10^(-6);
    breakpoint2 = 10^(-6);
    largerAbsoluteValue = 1.55;%
    cstart = hex2rgb('#08519c');
    cavg1 = hex2rgb('#f7fbff');
    cavg2 = hex2rgb('#fff5f0');
    cend = hex2rgb('#a50f15');
    % Set up an interpolation, from our non-uniform colour array including the
    % breakpoint to a nice evenly spaced colour map which changes the same
    colours = [cstart; cavg1; cavg2; cend];
    breakpoints = [-largerAbsoluteValue; breakpoint1; breakpoint2; largerAbsoluteValue];
    colours = interp1( breakpoints, colours, linspace(minLimColorbar,maxLimColorbar,100) );

    hf = figure; 
    colormap(colours)
    hCB = colorbar('southoutside');
    set(gca,'Visible',false)
    ax = gca;
    ax.CLim = [minLimColorbar maxLimColorbar];
    set(gcf, 'Position',  [0, 0, 996, 150])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 3 - Heatmaps for Strategies 3, 4 and 5
    %First row: SOC
    %Second row: SOC+ICI
    %Left: No TAM-targeting treatment
    %right: Maximal efficacy

    %Strategy 3 - Restore M1 phagocytic activity
    clear p;
    set_patients;
    deltaM = p.deltaM;
    lambdaS3Values = linspace(0,1,20);
    load('strategy3Results.mat')
    strategy3Results(2,:) = [];
    strategy3Results_normd = plot_HeatmapStrategy(strategy3Results,1,lambdaS3Values,'3',...
        axis_linewidth,axis_fontsize,'Strategy3');

    %Strategy 4 - Increase M2-to-M1 repolarisation rate
    clear p;
    set_patients;
    muM2 = p.muM2;
    lambdaS4v1Values = linspace(1,1/muM2,20);
    load('strategy4Results.mat')
    strategy4v1Results(2,:) = []; 
    strategy4v1Results_normd = plot_HeatmapStrategy(strategy4v1Results,1,lambdaS4v1Values,...
        '4',axis_linewidth,axis_fontsize,'Strategy4v1');  
    
    %Strategy 5 - Decrease fraction of resident TAMs that get M2 phenotype
    %upon activation
    clear p;
    lambdaS4v2Values = linspace(0,1,20);
    load('strategy5Results.mat')
    strategy4v2Results(2,:) = [];
    strategy4v2Results_normd = plot_HeatmapStrategy(strategy4v2Results,length(lambdaS4v2Values),...
        lambdaS4v2Values,'5',axis_linewidth,axis_fontsize,'Strategy4v2');

    %colormap Strategies 3, 4  and 5 
    minLimColorbar = -100;
    maxLimColorbar = 0;
    breakpoint = -50;
    % Set the breakpoint value
    cend = hex2rgb('#f7fbff');
    cstart = hex2rgb('#08519c');
    cavg = mean( [cstart; cend] );
    % Set up an interpolation, from our non-uniform colour array including the
    % breakpoint to a nice evenly spaced colour map which changes the same
    colours = [cstart; cavg; cend];
    breakpoints = [minLimColorbar; breakpoint; maxLimColorbar];
    colours = interp1( breakpoints, colours, linspace(minLimColorbar,maxLimColorbar,100) );

    hf = figure; 
    colormap(colours)
    hCB = colorbar('southoutside');
    set(gca,'Visible',false)
    ax = gca;
    ax.CLim = [minLimColorbar maxLimColorbar];
    set(gcf, 'Position',  [0, 0, 996, 150])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PANEL 7 -- TAM Targeting VCT : Enhancing phagocytosis STRATEGY 3 

    cd(folder_VCT)

    load("VCTPatients.mat")%
    N = size(VCTPatients,1); %nbr of patients

    load('VCT_SOC.mat')%
    tm = squeeze(VCT_SOC(:,3,:));
    timeVCT = linspace(0,2*365,20000);
    tmpostSOC = tm(:,find(timeVCT>=(lastTMZadmin+7),1,'first'))';
    
    load('VCT_TAM_Strategy3.mat')%
    %Percent decrease (between 0/no change and 100/total eradication) in
    %cancer cells after treatment
    maxChangeSOCTAM = 100.*(tmpostSOC(1,:)-VCT_TAM_Strategy3(1,:))./tmpostSOC(1,:);%
    maxChangeSOCTAMICI = 100.*(tmpostSOC(1,:)-VCT_TAM_Strategy3(2,:))./tmpostSOC(1,:);%
    
    %Normalize percent change with tumour growth over same period when
    %untreated (see Section 4 for 'growth' definition)
    dataSOCTAM = maxChangeSOCTAM./growth';
    dataSOCTAMICI = maxChangeSOCTAMICI./growth';
    
    cd(folder_code)
    h_fig = figure;
    hold on
    hgSOCTAM = histogram(-maxChangeSOCTAMICI);%percent change -> decrease is negative
    hgSOCTAM.FaceColor = color_green;
    hgSOCTAM.LineWidth = 1;
    hgSOCTAM.FaceAlpha = 0.75;
    hgSOCTAMICI = histogram(-maxChangeSOCTAM);%
    hgSOCTAMICI.FaceColor = color_blue;
    hgSOCTAMICI.LineWidth = 1;
    hgSOCTAMICI.FaceAlpha = 0.75;
    xlabel({'Percent change with SOC in glioblastoma', 'cell count after treatment'})
    ylabel('Number of patients')
    %legend({'SOC+ICI+Restoring phagocytosis','SOC+Restoring phagocytosis'},'Location','southoutside','FontSize',legend_fontsize)
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)

    %dataSOCTAM is lower bounded by 0 -- the larger its value, the larger
    %the decrease (see above)
    Q3 = quantile(dataSOCTAM,0.75);
    Q2 = quantile(dataSOCTAM,0.5);
    Q1 = quantile(dataSOCTAM,0.25);
    
    indCat1 = find(dataSOCTAM(:)<=Q1); 
    indCat2 = find(dataSOCTAM(:)<=Q2 & dataSOCTAM(:)>Q1); 
    indCat3 = find(dataSOCTAM(:)<=Q3 & dataSOCTAM(:)>Q2); 
    indCat4 = find(dataSOCTAM(:)>Q3);

    patientsInfo = VCTPatients;
    patientsInfo = log(patientsInfo);%natural logarithm
    baselineParam = log(baselineParam);

    patients1 = patientsInfo(indCat1,:);
    patients2 = patientsInfo(indCat2,:);
    patients3 = patientsInfo(indCat3,:);
    patients4 = patientsInfo(indCat4,:);
    patientCategory = zeros(size(patientsInfo,1),1);
    patientCategory(indCat1,1) = 1;
    patientCategory(indCat2,1) = 2;
    patientCategory(indCat3,1) = 3;
    patientCategory(indCat4,1) = 4;
    
    parameters = ["CL","k_{23}","k_a","\beta","\alpha_{AT}","M2:M1 ratio","Resection extent"];
    %get min and  max for each parameters to ensure that figures have the same
    %axis limits
    pmin= min(patientsInfo);
    pmax= max(patientsInfo);
    
    labelsParam = {'ln(CL)','ln(k_{23})',...
        'ln(k_a)','ln(\beta )',...
        'ln(\alpha )','ln(q )',''};
    titlesParam = {'TMZ Clearance','TMZ Plasma to CSF transfer rate','TMZ Absorption rate',...
        'Intrinsic tumour growth rate','Antigenicity of tumour cells','TME bias towards M2 phenotype',''};

    figTitles = {strcat('VCTStrategy','3','_CLhist_compared'),...
        strcat('VCTStrategy','3','_k23hist_compared'),...
        strcat('VCTStrategy','3','_kahist_compared'),...
        strcat('VCTStrategy','3','_betahist_compared'),...
        strcat('VCTStrategy','3','_aAThist_compared'),...
        strcat('VCTStrategy','3','_qhist_compared'),...
        strcat('VCTStrategy','3','_change')};

    nFig = 1;
    pvals_TAMT3_Wilcoxon = zeros(1,(length(parameters)-1));
    %FIGURES 4C-E and similar figures for TMZ PK parameters 
    for nParam=1:(length(parameters)-1)
        [pvalWilcoxon,~] = ranksum(patients1(:,nParam),patients4(:,nParam));
        pvals_TAMT3_Wilcoxon(nParam) = pvalWilcoxon;
        h_fig = figure;
        hold on
        hgCat1 = histogram(patients1(:,nParam));
        hgCat1.FaceColor = optionscat1.color_line;
        hgCat1.LineWidth = 1;
        hgCat1.FaceAlpha = 0.9;
        hgCat1.EdgeAlpha = 0.5;

        hgCat4 = histogram(patients4(:,nParam));
        hgCat4.FaceColor = optionscat4.color_line;
        hgCat4.LineWidth = 1;
        hgCat4.FaceAlpha = 0.6;
        hgCat4.EdgeAlpha = 0.5;

        yl = ylim;
        %we plot a line between two median if distributions are
        %significantly different with two-sided Wilcoxon test
        if(pvalWilcoxon<0.01)
            plot_pvalueline(hgCat1,hgCat4,yl,plot_linewidth)
        end
        xlabel(labelsParam{nParam})
        ylabel('Number of patients')
        %leg = legend({'Worst responders','Best responders'});
        set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
        nFig = nFig+1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Strategy 3 VCT: KM curves

    cd(folder_VCT)
    load("deathsStrategy3.mat")%
    load("deathsSOCICIVCT.mat")%

    tSOC = deathsVCT(2,:)./30.417;
    eSOC = ones(1,length(tSOC));
    eSOC(tSOC<0) = 0; %alive - censored
    tSOC(tSOC<0) = 10^6; 
    tSOC = tSOC-endTxTime./30.417; %0 = end of treatment instead of diagnosis
    grSOC = cell(1,length(tSOC));
    grSOC(1,:) = {'SOC'};

    tSOCICI = deathsVCT(3,:)./30.417;
    eSOCICI = ones(1,length(tSOCICI));
    eSOCICI(tSOCICI<0) = 0; %alive - censored
    tSOCICI(tSOCICI<0) = 10^6; 
    tSOCICI = tSOCICI-endTxTime./30.417; %0 = end of treatment instead of diagnosis
    grSOCICI = cell(1,length(tSOCICI));
    grSOCICI(1,:) = {'SOC+ICI'};

    tSOCTAMT3 = deathsStrategy3(1,:)./30.417;
    eSOCTAMT3 = ones(1,length(tSOCTAMT3));
    eSOCTAMT3(tSOCTAMT3<0) = 0; %alive - censored
    tSOCTAMT3(tSOCTAMT3<0) = 10^6; 
    tSOCTAMT3 = tSOCTAMT3-endTxTime./30.417; %0 = end of treatment instead of diagnosis
    grSOCTAMT3 = cell(1,length(tSOCTAMT3));
    grSOCTAMT3(1,:) = {'SOC+TAMT3'};

    tSOCTAMT3ICI = deathsStrategy3(2,:)./30.417;
    eSOCTAMT3ICI = ones(1,length(tSOCTAMT3ICI));
    eSOCTAMT3ICI(tSOCTAMT3ICI<0) = 0; %alive - censored
    tSOCTAMT3ICI(tSOCTAMT3ICI<0) = 10^6; 
    tSOCTAMT3ICI = tSOCTAMT3ICI-endTxTime./30.417; %0 = end of treatment instead of diagnosis
    grSOCTAMT3ICI = cell(1,length(tSOCTAMT3ICI));
    grSOCTAMT3ICI(1,:) = {'SOC+TAMT3+ICI'};

    adjustLW = -0.75;
    %Figure 4B
    fig = figure;
    hold on
    KM_curve(tSOC,'LineColor',hex2rgb(color_black),'LineStyle','-','ShowMedian',false,'LineWidth',plot_linewidth+adjustLW)
    KM_curve(tSOCICI,'LineColor',hex2rgb(color_black),'LineStyle',':','ShowMedian',false,'LineWidth',plot_linewidth+adjustLW)
    KM_curve(tSOCTAMT3,'LineColor',hex2rgb(color_blue),'LineStyle','-','ShowMedian',false,'LineWidth',plot_linewidth+adjustLW)
    KM_curve(tSOCTAMT3ICI,'LineColor',hex2rgb(color_blue),'LineStyle',':','ShowMedian',false,'LineWidth',plot_linewidth+adjustLW)
    xlim([0 (2*365-256)/30.417])
    xlabel({'Time after treatment (months)',''})
    ylabel('SP (%)')
    %legend({'SOC','SOC+ICI','SOC+Restoring phagocytosis','SOC+ICI+Restoring phagocytosis'},'Location','southoutside','FontSize',legend_fontsize)
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Decreasing M2:M1 ratio of 25% patients that respond the least

    cd(folder_code)
    clear p
    set_patients;
    set_treatment;
    
    %FIGURE 4 H
    %plot final cancer cell count
    %fcc: final cancer cell count (Modified is SOC+S3+S4, nonModified is SOC+S3)
    cd(folder_VCT)
    load('fccBestResponders_nonModified.mat')
    load('fccWorstResponders_ModifiedS4.mat')
    load('fccWorstResponders_nonModified.mat')
    fccWorstResponders_Modified = fccWorstResponders_Modified.*growth(indCat1);
    fccWorstResponders_nonModified = fccWorstResponders_nonModified.*growth(indCat1);
    fccBestResponders_nonModified = fccBestResponders_nonModified.*growth(indCat4);
    cd(folder_code)
    violinplot(cat(2,-fccWorstResponders_nonModified,-fccWorstResponders_Modified,-fccBestResponders_nonModified),...
        {'','',''},'ViolinColor',cat(1,optionscat1.color_line,optionscat1.color_line,optionscat4.color_line),'MarkerSize',6,'MedianMarkerSize',10)
    xticks([1 2 3])
    xticklabels({'','',''})
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    
    %FIGURE 4G
    % plotM2:M1 ratio
    cd(folder_VCT)
    load('M2M1ratioWorstResponders_nonModified')
    load('M2M1ratioWorstResponders_Modified')
    load('M2M1ratioBestResponders_nonModified')
    cd(folder_code)
    figure
    hold on
    violinplot(cat(2,M2M1ratioWorstResponders_nonModified(:,2),M2M1ratioWorstResponders_Modified(:,2),...
        M2M1ratioBestResponders_nonModified(:,2)),{'','',''},'ViolinColor',cat(1,optionscat1.color_line,optionscat1.color_line,optionscat4.color_line),...
        'MarkerSize',6,'MedianMarkerSize',10)
    xticks([1 2 3])
    xticklabels({'','',''})
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 4F: M1-biased vs M2-biased population

    cd(folder_VCT)
    load("deathsStrategy3_M1biased.mat")
    deathsM1 = deathsStrategy3;
    load("deathsStrategy3_M2biased.mat")
    deathsM2 = deathsStrategy3;
    load("deathsStrategy3.mat")
    figure
    hold on
    KM_curve(deathsM1(1,:)./30.417-endTxTime./30.417,'LineColor',color_green,'LineWidth',plot_linewidth,'ShowMedian',false)
    KM_curve(deathsM1(2,:)./30.417-endTxTime./30.417,'LineColor',color_green,'LineWidth',plot_linewidth,'ShowMedian',false,'LineStyle','--')
    KM_curve(deathsM2(1,:)./30.417-endTxTime./30.417,'LineColor',color_purple,'LineWidth',plot_linewidth,'ShowMedian',false)
    KM_curve(deathsM2(2,:)./30.417-endTxTime./30.417,'LineColor',color_purple,'LineWidth',plot_linewidth,'ShowMedian',false,'LineStyle','--')
    xlabel('Time after treatment (months)')
    ylabel('SP (%)')
    legend({'M1-biased','','M2-biased',''})
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    cd(folder_code)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 5: Potential in the clinic
    
    cd(folder_code)
    
    clear p;
    p.TAMT3decayingEfficacy = true;
    p.breakMac = 1;
    set_patients;
    p.nPeriodMac = 14;
    set_treatment;
    p.treatmentMac = 0;
    [solSOC,pSOC] = solver_gbm(p);
    allpopSOC = deval(solSOC,time);
    tumourSOC = allpopSOC(3,:);
    
    %Administration every week
    clear p;
    p.TAMT3decayingEfficacy = true;
    p.breakMac = 1;
    set_patients;
    p.nPeriodMac = 7;
    set_treatment;
    p.treatmentMac = 1;
    [solTAM_1wk,p_1wk] = solver_gbm(p);
    allpopTAM_1wk = deval(solTAM_1wk,time);
    tumourTAM_1wk = allpopTAM_1wk(3,:);
    
    p.addMonthsMacTreatment = 12;
    set_treatment;
    p.treatmentMac = 1;
    [solTAM_1wk_longer,p_1wk_longer] = solver_gbm(p);
    allpopTAM_1wk_longer = deval(solTAM_1wk_longer,time);
    tumourTAM_1wk_longer = allpopTAM_1wk_longer(3,:);
    
    %Administration every two weeks
    clear p;
    p.TAMT3decayingEfficacy = true;
    p.breakMac = 1;
    set_patients;
    p.nPeriodMac = 14;
    set_treatment;
    p.treatmentMac = 1;
    [solTAM_2wk,p_2wk] = solver_gbm(p);
    allpopTAM_2wk = deval(solTAM_2wk,time);
    tumourTAM_2wk = allpopTAM_2wk(3,:);
    
    p.addMonthsMacTreatment = 12;
    set_treatment;
    p.treatmentMac = 1;
    [solTAM_2wk_longer,p_2wk_longer] = solver_gbm(p);
    allpopTAM_2wk_longer = deval(solTAM_2wk_longer,time);
    tumourTAM_2wk_longer = allpopTAM_2wk_longer(3,:);
    
    detectionVolumeThreshold = (4/3).*(15^3).*pi;%
    
    %Figure 5B
    h_fig = figure;
    hold on
    plot(time./30.417,tumourSOC,'Color',color_blue,'LineWidth',plot_linewidth)
    plot(time./30.417,tumourTAM_1wk,'Color',color_green,'LineWidth',plot_linewidth)
    plot(time./30.417,tumourTAM_2wk,'Color',color_orange,'LineWidth',plot_linewidth)
    plot(time./30.417,tumourTAM_1wk_longer,'Color',color_green,'LineStyle',':','LineWidth',plot_linewidth)
    plot(time./30.417,tumourTAM_2wk_longer,'Color',color_orange,'LineStyle',':','LineWidth',plot_linewidth)
    yline(p.volDeath,'LineWidth',plot_linewidth,'Color',color_red)
    xline((p.AdministrationTimesTMZ(end)+7)./30.417,'LineWidth',2,'Color',color_black)
    xa = linspace(0, 24, 25);                             
    ya = detectionVolumeThreshold.*ones(1,25);                                  
    patch([xa fliplr(xa)], [zeros(size(ya)) fliplr(ya)],[0.7 0.7 0.7],'FaceAlpha',0.3,'EdgeColor','none')
    hold off
    xlabel('Time after diagnosis (months)')
    xlim([0 24])
    ylabel('Cancer cells (10^6 cells)')
    legend({'SOC','SOC + anti-CD47 every week','SOC + anti-CD47 every 2 weeks',...
        'SOC + anti-CD47 every week until 1 year post-SOC','SOC + anti-CD47 every 2 weeks until 1 year post-SOC'}, ...
        'Location','southoutside','NumColumns',1,'FontSize',legend_fontsize)
    set(gca,'Fontsize',axis_fontsize,'LineWidth',axis_linewidth)

    %close up
    h_fig = figure;
    hold on
    plot(time./30.417,tumourSOC,'Color',color_blue,'LineWidth',plot_linewidth)
    plot(time./30.417,tumourTAM_1wk,'Color',color_green,'LineWidth',plot_linewidth)
    plot(time./30.417,tumourTAM_2wk,'Color',color_orange,'LineWidth',plot_linewidth)
    plot(time./30.417,tumourTAM_1wk_longer,'Color',color_green,'LineStyle',':','LineWidth',plot_linewidth)
    plot(time./30.417,tumourTAM_2wk_longer,'Color',color_orange,'LineStyle',':','LineWidth',plot_linewidth)
    xline((p.AdministrationTimesTMZ(end)+7)./30.417,'LineWidth',1,'Color',color_black)
    xa = linspace(0, 24, 25);                             
    ya = detectionVolumeThreshold.*ones(1,25);                                  
    patch([xa fliplr(xa)], [zeros(size(ya)) fliplr(ya)],[0.7 0.7 0.7],'FaceAlpha',0.3,'EdgeColor','none')
    hold off
    xlabel('')
    xlim([8 9])
    ylabel('')
    yticks([])
    set(gca,'Fontsize',10,'LineWidth',axis_linewidth)

    %Figure 5A
    h_fig = figure;
    hold on
    plot(time./30.417,100.*allpopTAM_1wk_longer(14,:),'Color',color_green,'LineWidth',plot_linewidth)
    plot(time./30.417,100.*allpopTAM_2wk_longer(14,:),'Color',color_orange,'LineWidth',plot_linewidth)
    hold off
    xlabel({'Time after first administration','(months)'})
    xlim([p.AdministrationTimesMacrophages(1) (p.AdministrationTimesMacrophages(1)+60)]./30.417)
    ylabel({'Treatment', 'efficacy (%)'})
    ylim([40 100])
    set(gca,'Fontsize',axis_fontsize,'LineWidth',axis_linewidth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TMZ PD model - Fitting

    load('lodicoDataTransformed.mat')
    C = exp(lodicoDataTransformed(:,1));%lodicoDataTransformed(:,1) was the natural logarithm of the concentrations
    effect = lodicoDataTransformed(:,2);
    
    Emax = 90.0031; %Maximal effect over 72 hours
    IC50 = 0.0011;% drug concentration at which the effect on cell viability if half of it's possible effect  
    h    = 0.5682; %hill coefficient
    
    x = linspace(C(1),C(end),1000);
    y = (Emax.*x.^h)./(IC50.^h+x.^h);
    
    h_fig = figure;
    hold on
    plot(x,y./3,'LineStyle','-','LineWidth',plot_linewidth,'Color',color_blue)%
    scatter(C,effect./3,30,hex2rgb(color_red),"filled")
    ylabel({'TMZ daily','cytotoxicity (%)'})
    xlabel('TMZ Concentration (mg/ml)')
    set(gca,'XScale','log') %log-scale base 10
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    legend({'Fit','Data'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Figures 7, 10 and 11:  Visualise Virtual patients distributions

    cd(folder_VCT)
    %You can change the data set you want to visualise
    load('VCTPatients_M2biased.mat')
    cd(folder_code)
    
    CL = VCTPatients(:,1);
    k23 = VCTPatients(:,2);
    ka = VCTPatients(:,3);
    beta = VCTPatients(:,4);
    aAT = VCTPatients(:,5);
    q = VCTPatients(:,6);
    
    %General population
    %muTumour = [0.013 0.0375 0.7512];
    %stdTumour = [0.0004 0.01127 0.0783];

    %M1-biased
    %muTumour = [0.013 0.0375 0.25];
    %stdTumour = [0.0004 0.01127 0.2310];

    %M2-biased
    muTumour = [0.013 0.0375 0.75];
    stdTumour = [0.0004 0.01127 0.0959];

    muPK = [10 7.2*10^(-4) 5.8].*24;
    stdPK = [0.0470 0.8961 0.1649];
    
    %TMZ Clearance
    h_fig = figure;
    hold on
    histogram(log(CL),'FaceColor',color_green)
    xline(log(muPK(1)),'LineWidth',plot_linewidth)
    hold off
    xlabel('ln(k_{CL})')
    ylabel('Number of patients')
    %title('TMZ Clearance')
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    
    %TMZ Plasma-to-CSF transfer rate
    h_fig = figure;
    hold on
    histogram(log(k23),'FaceColor',color_green)
    xline(log(muPK(2)),'LineWidth',plot_linewidth)
    xlabel('ln(k_23 )')
    ylabel('Number of patients')
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    
    %TMZ Absorption rate
    h_fig = figure;
    hold on
    histogram(log(ka),'FaceColor',color_green)
    xline(log(muPK(3)),'LineWidth',plot_linewidth)
    xlabel('ln(k_a )')
    ylabel('Number of patients')
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    
    %beta
    h_fig = figure;
    hold on
    histogram(log(beta),'FaceColor',color_red)
    xline(log(muTumour(1)),'LineWidth',plot_linewidth)
    hold off
    xlabel('ln(\beta )')
    ylabel('Number of patients')
    %title('Tumour Intrinsinc Growth Rate')
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    
    %aAT
    h_fig = figure;
    hold on
    histogram(log(aAT),'FaceColor',color_red)
    xline(log(muTumour(2)),'LineWidth',plot_linewidth)
    hold off
    xlabel('ln(\alpha )')
    ylabel('Number of patients')
    %title('Antigenicity of Tumour Cells')
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    
    %q
    h_fig = figure;
    hold on
    histogram(log(q),'FaceColor',color_red)
    xline(log(muTumour(3)),'LineWidth',plot_linewidth)
    hold off
    xlabel('ln(q)')
    ylabel('Number of patients')
    %title('Fraction of macrophages with M2 phenotype')
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Figure 9: Kaplan-meier curves depending on median of parameters

    cd(folder_VCT)
    load('VCTPatients.mat')
    load('VCT_SOC_ICI.mat')
    load('VCT_SOC.mat')
    N = size(VCTPatients,1);
    load('deathsSOCICIVCT.mat')
    deathsSOC = deathsVCT(2,:);
    deathsICI = deathsVCT(3,:);
    tSOC = deathsSOC./30.417;
    tSOC = tSOC-(256/30.417);
    eSOC = ones(1,length(tSOC));
    eSOC(tSOC<0) = 0; %alive - censored
    tSOC(tSOC<0) = 10^6;

    tICI = deathsICI./30.417;
    tICI = tICI-(256/30.417);
    eICI = ones(1,length(tICI));
    eICI(tICI<0) = 0; %alive - censored
    tICI(tICI<0) = 10^6;

    gr = cell(1,length(tSOC)+length(tICI));
    gr(1,1:length(tSOC)) = {'SOC'};
    gr(1,length(tSOC)+1:end) = {'ICI'};

    cd(folder_code)
    %Supplementary figure B
    figure
    hold on
    KM_curve(tSOC,'ShowMedian',false,'Linewidth',plot_linewidth,'Linecolor',color_blue)
    KM_curve(tICI,'ShowMedian',false,'Linewidth',plot_linewidth,'Linecolor',color_green)
    xlabel('Time after treatment (months)')
    ylabel({'Survival probability', '(%)'})
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    
    param = ["k_{CL}","k_{23}","k_a","\beta","\alpha_{AT}","q"];
    %Supplementary figure 8 C-H
    for i=1:6
        figure
        hold on
        KM_curve(tSOC(VCTPatients(:,i)<median(VCTPatients(:,i))),'ShowMedian',false,'Linewidth',plot_linewidth,'Linecolor',color_blue)
        KM_curve(tSOC(VCTPatients(:,i)>=median(VCTPatients(:,i))),'ShowMedian',false,'Linewidth',plot_linewidth,'Linecolor',color_green)
        %xlabel({'Time after','treatment (months)'})
        ylabel('SP (%)')
        set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
        lh = legend({strcat(param(i),'<',num2str(median(VCTPatients(:,i)))),...
            strcat(param(i),'>=',num2str(median(VCTPatients(:,i))))},'FontSize',10,'Location','southwest');
        lh.ItemTokenSize = [5, 5];
    end

    %Supplementary figure A
    [minSurvival,indMinSurvival] = min(tICI);
    [maxSurvival,indMaxSurvival] = max(tICI);
    clear p
    set_patients;
    figure
    hold on
    plot(timeVCT./30.1417,squeeze(VCT_SOC_ICI(indMinSurvival,3,:)),'LineWidth',plot_linewidth,'Color',color_blue)
    plot(timeVCT./30.1417,squeeze(VCT_SOC_ICI(indMaxSurvival,3,:)),'LineWidth',plot_linewidth,'Color',color_green)
    yline(p.volDeath,'LineWidth',plot_linewidth,'Color',color_red)
    x2 = [timeVCT./30.417 fliplr(timeVCT./30.417)];
    inBetween = [squeeze(VCT_SOC_ICI(indMinSurvival,3,:))', fliplr(squeeze(VCT_SOC_ICI(indMaxSurvival,3,:))')];
    yline(p.volDeath,'LineWidth',plot_linewidth,'LineStyle','-','Color',color_red)
    fill(x2, inBetween, hex2rgb('#808080'),'FaceAlpha',0.15,'EdgeColor','none');
    ylabel({'Tumour cells','(10^6 cells)'})
    xlabel('Time after diagnosis (months)')
    xlim([0 24])
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    legend({'Minimum survival','Maximum survival','Death threshold'},'FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Figure 5: eFAST: Plot Sti for tumour burden - beta
    cd(folder_efast)
    load('Model_efastNs1000Nr5_tumourburden.mat')
    cd(folder_code)
    
    test_pvalue = 0.01;% 
    
    tempSti = s_Y.Sti;
    %Reordering the variables
    tempSti(5,:) = s_Y.Sti(7,:); %n1
    tempSti(7,:) = s_Y.Sti(9,:); %alpha
    tempSti(8,:) = s_Y.Sti(5,:); %beta
    tempSti(9,:) = s_Y.Sti(8,:); %mu_P
    
    temp_pval = s_Y.p_Sti;
    temp_pval(5,:) = s_Y.p_Sti(7,:);
    temp_pval(7,:) = s_Y.p_Sti(9,:);
    temp_pval(8,:) = s_Y.p_Sti(5,:);
    temp_pval(9,:) = s_Y.p_Sti(8,:);

    tempSti(end-1,:) = [];
    temp_pval(end,:) = [];
    
    temp_efast_var={'\lambda_M','q','\mu_{M1}','\mu_{M2}','\eta_1','\alpha_{AT}','\alpha_2','\beta','\mu_P','dummy'};%,
    temp_efast_var(end-1) = '';

    nFig = 1;
    significant = zeros(0,2);
    p_values = squeeze(temp_pval(:,1,:,nFig));
    for kk=1:size(p_values,1)
        for ll=1:size(p_values,2)
            p_to_test = p_values(kk,ll);
            if (p_to_test <= test_pvalue) && (tempSti(kk,1)>tempSti(end,1))
                significant(end+1,1) = ll;
                significant(end,2) = kk;
            end
        end
    end
    fig = figure;
    dataEfast = tempSti(end-1,:,nFig);
    hold on
    b = bar(dataEfast,'FaceColor','flat');
    yl = ylim;
    b(1).CData = hex2rgb(color_palegreen);
    xticks([1])
    ylabel('S_{Ti}')
    set(gca,'XTickLabel',temp_efast_var,'LineWidth',axis_linewidth,'FontSize',axis_fontsize)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig)

%% Supplementary Figure 5: eFAST: Plot Sti for tumour burden - without beta
    cd(folder_efast)
    load('Model_efastNs1000Nr5_tumourburden_withoutBeta.mat')
    cd(folder_code)
    
    test_pvalue = 0.01;%
    color_pvalue = [251 128 114]./255;
    
    tempSti = s_Y.Sti;
    tempSti(5,:) = s_Y.Sti(7,:);
    tempSti(7,:) = s_Y.Sti(9,:);
    tempSti(8,:) = s_Y.Sti(5,:);
    tempSti(9,:) = s_Y.Sti(8,:);
    
    temp_pval = s_Y.p_Sti;
    temp_pval(5,:) = s_Y.p_Sti(7,:);
    temp_pval(7,:) = s_Y.p_Sti(9,:);
    temp_pval(8,:) = s_Y.p_Sti(5,:);
    temp_pval(9,:) = s_Y.p_Sti(8,:);
    
    tempSti(end-1,:) = [];
    temp_pval(end,:) = [];
    
    temp_efast_var={'\lambda_M','q','\mu_{M1}','\mu_{M2}','\eta_1','\alpha_{AT}','\alpha_2','\beta','\mu_P','dummy'};%,
    temp_efast_var(end-1) = '';


    nFig = 1;
    significant = zeros(0,2);
    p_values = squeeze(temp_pval(:,1,:,nFig));
    for kk=1:size(p_values,1)
        for ll=1:size(p_values,2)
            p_to_test = p_values(kk,ll);
            if (p_to_test <= test_pvalue) && (tempSti(kk,1)>tempSti(end,1))
                significant(end+1,1) = ll;
                significant(end,2) = kk;
            end
        end
    end
    
    fig = figure;
    dataEfast = tempSti(:,:,nFig);
    dataEfast(end-2,:) = []; %remove beta
    hold on
    b = bar(dataEfast,'FaceColor','flat');
    yl = ylim;
    if size(significant,1)>0
        for gg=1:size(significant,1)
            tempTimePoint = significant(gg,1);
            tempParam = significant(gg,2);
            xpt = significant(gg,2);
            ypt = b.YData(xpt);
            ypt = ypt+0.01*(yl(2)-yl(1));
        end
    end
    for nParam=1:size(dataEfast,1)
        if nParam<6
            b(1).CData(nParam,:) = hex2rgb(color_paleblue);
        elseif nParam<8
            b(1).CData(nParam,:) = hex2rgb(color_palered);
        elseif nParam<size(dataEfast,1)
            b(1).CData(nParam,:) = hex2rgb(color_paleorange);
        else
            b(1).CData(nParam,:) = hex2rgb(color_black);
        end
    end
    xticks([1:size(tempSti,1)])
    ylabel('S_{Ti}')
    set(gca,'XTickLabel',temp_efast_var,'LineWidth',axis_linewidth,'FontSize',axis_fontsize)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Supplementary Figure 5: eFAST: Plot Sti for SE(SOC+ICI) - beta
    cd(folder_efast)
    load('Model_efastNs1000Nr5_nivolumab_withBeta.mat')
    cd(folder_code)
    
    test_pvalue = 0.01;%
    color_pvalue = [251 128 114]./255;
    
    tempSti = s_Y.Sti;
    tempSti(5,:) = s_Y.Sti(7,:);
    tempSti(7,:) = s_Y.Sti(9,:);
    tempSti(8,:) = s_Y.Sti(5,:);
    tempSti(9,:) = s_Y.Sti(8,:);
    
    temp_pval = s_Y.p_Sti;
    temp_pval(5,:) = s_Y.p_Sti(7,:);
    temp_pval(7,:) = s_Y.p_Sti(9,:);
    temp_pval(8,:) = s_Y.p_Sti(5,:);
    temp_pval(9,:) = s_Y.p_Sti(8,:);
    
    temp_efast_var={'\lambda_M','q','\mu_{M1}','\mu_{M2}','\eta_1','\alpha_{AT}','\alpha_2','\beta','\mu_P','dummy'};%,
    nFig = 1;
    significant = zeros(0,2);
    p_values = squeeze(temp_pval(:,1,:,nFig));
    for kk=1:size(p_values,1)
        for ll=1:size(p_values,2)
            p_to_test = p_values(kk,ll);
            if (p_to_test <= test_pvalue) && (tempSti(kk,1)>tempSti(end,1))
                significant(end+1,1) = ll;
                significant(end,2) = kk;
            end
        end
    end
    
    fig = figure;
    dataEfast = tempSti(end-2,:,nFig);
    hold on
    b = bar(dataEfast,'FaceColor','flat');
    yl = ylim;
    b(1).CData = hex2rgb(color_palegreen);
    xticks([1])
    xtickangle(45)
    ylabel('S_{Ti}')
    set(gca,'XTickLabel',{'dummy'},'LineWidth',axis_linewidth,'FontSize',axis_fontsize)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig)

%% Supplementary Figure 5: eFAST: Plot Sti for SE(SOC+ICI) - without beta
    cd(folder_efast)
    load('Model_efastNs1000Nr5_nivolumab_withoutBeta.mat')
    cd(folder_code)
    
    test_pvalue = 0.01;%
    color_pvalue = [251 128 114]./255;
    
    tempSti = s_Y.Sti;
    tempSti(5,:) = s_Y.Sti(7,:);
    tempSti(7,:) = s_Y.Sti(9,:);
    tempSti(8,:) = s_Y.Sti(5,:);
    tempSti(9,:) = s_Y.Sti(8,:);
    
    temp_pval = s_Y.p_Sti;
    temp_pval(5,:) = s_Y.p_Sti(7,:);
    temp_pval(7,:) = s_Y.p_Sti(9,:);
    temp_pval(8,:) = s_Y.p_Sti(5,:);
    temp_pval(9,:) = s_Y.p_Sti(8,:);
    
    temp_efast_var={'\lambda_M','q','\mu_{M1}','\mu_{M2}','\eta_1','\alpha_{AT}','\alpha_2','\beta','\mu_P','dummy'};%,
    nFig = 1;
    significant = zeros(0,2);
    p_values = squeeze(temp_pval(:,1,:,nFig));
    for kk=1:size(p_values,1)
        for ll=1:size(p_values,2)
            p_to_test = p_values(kk,ll);
            if (p_to_test <= test_pvalue) && (tempSti(kk,1)>tempSti(end,1))
                significant(end+1,1) = ll;
                significant(end,2) = kk;
            end
        end
    end
    
    fig = figure;
    dataEfast = tempSti(:,:,nFig);
    dataEfast(end-2,:) = []; %remove beta
    hold on
    b = bar(dataEfast,'FaceColor','flat');
    yl = ylim;
    if size(significant,1)>0
        for gg=1:size(significant,1)
            tempTimePoint = significant(gg,1);
            tempParam = significant(gg,2);
            xpt = significant(gg,2);
            ypt = b.YData(xpt);
            ypt = ypt+0.01*(yl(2)-yl(1));
        end
    end
    for nParam=1:size(dataEfast,1)
        if nParam<6
            b(1).CData(nParam,:) = hex2rgb(color_paleblue);
        elseif nParam<8
            b(1).CData(nParam,:) = hex2rgb(color_palered);
        elseif nParam<size(dataEfast,1)
            b(1).CData(nParam,:) = hex2rgb(color_paleorange);
        else
            b(1).CData(nParam,:) = hex2rgb(color_black);
        end
    end
    xticks([1:size(tempSti,1)])
    ylabel('S_{Ti}')
    set(gca,'XTickLabel',temp_efast_var,'LineWidth',axis_linewidth,'FontSize',axis_fontsize)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Figure 6B - Only decreasing death rates of M2 macrophages

    clear p;
    set_patients;
    set_treatment;
    l = 10;
    lambdaS1 = linspace(1,1/max(p.deltaM,p.deltaMR),l);
    
    %get color scale
    cstart = [253 208 162]./255;
    cend = [166 54 3]./255;
    cavg = mean( [cstart; cend] );
    colours = [cstart; cavg; cend];
    breakpoints = [0; 0.5; 1];
    colours = interp1( breakpoints, colours,linspace(0,1,10) );
    S1cells = zeros(14,l,length(time));

    for nLambda=1:l
        nLambda
        %Strategy 1
        clear p;
        set_patients;
        set_treatment;
        p.treatmentICI = 1;
        p.treatmentMac = 1;
        p.newLambdaMacValues = [lambdaS1(nLambda) p.lambdaS2 p.lambdaS3 p.lambdaS4v1 p.lambdaS4v2];
        [sol,~] = solver_gbm(p);
        pop = deval(sol,time);
        pop(14,:) = [];
        S1cells(:,nLambda,:) = pop;
    end
    
    % SOC
    SOCcells = zeros(14,l,length(time));
    for nLambda=1:l
        nLambda
        %Strategy 1
        clear p;
        set_patients;
        set_treatment;
        p.treatmentICI = 1;
        [sol,~] = solver_gbm(p);
        pop = deval(sol,time);
        pop(14,:) = [];
        SOCcells(:,nLambda,:) = pop;
    end

    %DECREASE ONLY DEATH RATE OF M2 MACROPHAGES 
    S1cells_v2 = zeros(14,l,length(time));
    for nLambda=1:l
        nLambda
        clear p;
        p.decreaseDeathRateOnlyM2 = true;
        set_patients;
        set_treatment;
        p.treatmentICI = 1;
        p.treatmentMac = 1;
        p.newLambdaMacValues = [lambdaS1(nLambda) p.lambdaS2 p.lambdaS3 p.lambdaS4v1 p.lambdaS4v2];
        [sol,~] = solver_gbm(p);
        pop = deval(sol,time);
        pop(14,:) = [];
        S1cells_v2(:,nLambda,:) = pop;
    end

    optionsS1 = optionscat1;
    optionsS1.linewidth = plot_linewidth;
    optionsS1.color_line = hex2rgb(color_red);
    optionsS1.x_axis = optionsS1.x_axis./30.417;
    optionsS1v2 = optionsS1;
    optionsS1v2.color_line = hex2rgb(color_blue);
    optionsSOC = optionsS1;
    optionsSOC.color_line = hex2rgb(color_green);
    
    %S1 Tumour burden
    hfig= figure;
    hold on
    plot(time./30.417,squeeze(S1cells(3,end,:)),'LineWidth',plot_linewidth,'Color',color_red)
    plot(time./30.417,squeeze(S1cells_v2(3,end,:)),'LineWidth',plot_linewidth,'Color',color_blue)
    plot(time./30.417,squeeze(SOCcells(3,end,:)),'LineWidth',plot_linewidth,'Color',color_green)
    xlim([0 10])
    xl = xlim;
    yl = ylim;
    h=fill([startTxTime,endTxTime,endTxTime,startTxTime]./30.417,[yl(1),yl(1),yl(2),yl(2)],[128 128 128]/255,'EdgeColor','none');
    h.FaceAlpha=0.3;
    xlabel('Time after diagnosis (months)')
    ylabel('Tumour cells (10^6 cells)')
    legend({'Target all TAMs','Target only M2 TAMs','SOC','Treatment period'},'Location','best','FontSize',legend_fontsize)
    set(gca,'LineWidth',axis_linewidth,'FontSize',axis_fontsize)

    indEndTx = find(time>=endTxTime,1,'first');
    changeS1v2 = 100.*(S1cells_v2(3,end,indEndTx)-SOCcells(3,end,indEndTx))./SOCcells(3,end,indEndTx);

    %S1 Tumour burden - close up
    fig = figure;
    hold on
    plot(time./30.417,squeeze(S1cells(3,end,:)),'LineWidth',1,'Color',color_red)
    plot(time./30.417,squeeze(S1cells_v2(3,end,:)),'LineWidth',1,'Color',color_blue)
    plot(time./30.417,squeeze(SOCcells(3,end,:)),'LineWidth',1,'Color',color_green)
    xlim([8.25 8.5])
    xl = xlim;
    yl = ylim;
    ylim([yl(1) yl(2)])
    h=fill([startTxTime,endTxTime,endTxTime,startTxTime]./30.417,[yl(1),yl(1),yl(2),yl(2)],[128 128 128]/255,'EdgeColor','none');
    h.FaceAlpha=0.3;
    yticks([])
    set(gca,'LineWidth',1,'FontSize',12) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Figure 2: fitting alpha and beta

    %no treatment
    clear p;
    set_patients;
    set_treatment;
    p.treatmentBoolean = 0;
    p.treatmentRT = 0;
    p.treatmentTMZ = 0;
    p.treatmentICI = 0;
    p.treatmentMac = 0;
    [solNo,~] = solver_gbm(p);
    tumourNo = deval(solNo,time,3);
    
    %RT
    clear p;
    set_patients;
    set_treatment;
    p.treatmentBoolean = 1;
    p.treatmentRT = 1;
    p.treatmentTMZ = 0;
    p.treatmentICI = 0;
    p.treatmentMac = 0;
    [solRT,~] = solver_gbm(p);
    tumourRT = deval(solRT,time,3);
    
    %SOC
    clear p;
    set_patients;
    set_treatment;
    p.treatmentBoolean = 1;
    p.treatmentRT = 1;
    p.treatmentTMZ = 1;
    p.treatmentICI = 0;
    p.treatmentMac = 0;
    [solSOC,~] = solver_gbm(p);
    tumourSOC = deval(solSOC,time,3);
    
    h_fig = figure;
    hold on
    plot(time./30.417,tumourNo,'LineWidth',plot_linewidth,'Color',color_blue)
    plot(time./30.417,tumourRT,'LineWidth',plot_linewidth,'Color',color_green)
    plot(time./30.417,tumourSOC,'LineWidth',plot_linewidth,'Color',color_orange)
    yline(p.volDeath,'LineStyle','-','LineWidth',plot_linewidth,'Color',color_red)
    yl = ylim;
    CIno=fill([3,5,5,3],[yl(1),yl(1),yl(2),yl(2)],hex2rgb(color_blue));
    CIno.FaceAlpha=0.3;
    CIno.EdgeColor = 'none';
    CIrt=fill([11.2,13,13,11.2],[yl(1),yl(1),yl(2),yl(2)],hex2rgb(color_green));
    CIrt.FaceAlpha=0.3;
    CIrt.EdgeColor = 'none';
    CIsoc=fill([13.2,16.8,16.8,13.2],[yl(1),yl(1),yl(2),yl(2)],hex2rgb(color_orange));
    CIsoc.FaceAlpha=0.3;
    CIsoc.EdgeColor = 'none';
    ylabel('Glioblastoma cells (10^6 cells)')
    xlabel('Time after diagnosis (months)')
    legend({'No treatment','Surgery+RT','SOC'},'Location','best','FontSize',legend_fontsize)
    set(gca,'LineWidth',axis_linewidth,'FontSize',axis_fontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Figure 8: Validity of virtual populations

    optionsVP = options;
    optionsVP.x_axis = timeVCT./30.417;
    optionsVP.color_line = hex2rgb(color_blue);
    optionsVP.color_area = hex2rgb(color_paleblue);

    cd(folder_VCT)
    load('M2M1ratiodata.mat')
    load('karimiCellFrequency.mat')
    cd(folder_code)

    tNoTx = zeros(1,size(VCT_none,1));
    for nP=1:length(tNoTx)
        tNoTx(1,nP) = timeVCT(find(squeeze(VCT_none(nP,3,:))>=volDeath,1,'first'));
    end
    [minSurvival,indMinSurvival] = min(tNoTx);
    [maxSurvival,indMaxSurvival] = max(tNoTx);

    %Panel A: Tumour evolution without treatment
    fig = figure;
    hold on
    plot(timeVCT./30.1417,squeeze(VCT_none(indMinSurvival,3,:)),'LineWidth',plot_linewidth,'Color',color_blue)
    plot(timeVCT./30.1417,squeeze(VCT_none(indMaxSurvival,3,:)),'LineWidth',plot_linewidth,'Color',color_green)
    yline(volDeath,'LineWidth',plot_linewidth,'Color',color_red)
    x2 = [timeVCT./30.417 fliplr(timeVCT./30.417)];
    inBetween = [squeeze(VCT_none(indMinSurvival,3,:))', fliplr(squeeze(VCT_none(indMaxSurvival,3,:))')];
    yline(p.volDeath,'LineWidth',plot_linewidth,'LineStyle','-','Color',color_red)
    fill(x2, inBetween, hex2rgb('#808080'),'FaceAlpha',0.15,'EdgeColor','none');
    ylabel({'Tumour cells','(10^6 cells)'})
    xlabel({'Time after diagnosis','(months)'})
    xlim([0 10])
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    legend({'Minimum survival','Maximum survival','Death threshold'},'FontSize',legend_fontsize)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig)
    
    M1cells = squeeze(VCT_none(:,1,1));
    M2cells = squeeze(VCT_none(:,2,1));
    M0cells = squeeze(VCT_none(:,11,1));
    CD8cells = squeeze(VCT_none(:,4,1));
    TMEcells = M1cells+M2cells+M0cells+CD8cells;
    CD8percent = 100.*CD8cells./TMEcells;
    M2M1ratio = M2cells./M1cells;
    
    %Panel B: CD8+ T cell infiltration at diagnosis
    fig1 = figure;
    hold on
    violinplot(CD8percent,{''},'ViolinColor',hex2rgb(color_blue),'MarkerSize',6,'MedianMarkerSize',10,'ShowBox',false)
    plot([0.8 1.2],100.*[mean(karimiCellFrequency(:,1)) mean(karimiCellFrequency(:,1))],'LineWidth',3,'LineStyle',':','Color','black')
    plot([0.8 1.2],[mean(CD8percent) mean(CD8percent)],'LineWidth',3,'LineStyle',':','Color',color_red)
    xticks([1])
    xticklabels({''})
    ylabel({'CD8+ T cells','(% of TME cells)'})
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    set(fig1,'Units','Inches');
    pos = get(fig1,'Position');
    set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig1)

    %Panel C: M2M1 ratio at diagnosis
    fig3 = figure;
    hold on
    violinplot(M2M1ratio,{''},'ViolinColor',hex2rgb(color_blue),'MarkerSize',6,'MedianMarkerSize',10,'ShowBox',false)
    plot([0.8 1.2],[3.02 3.02],'LineWidth',3,'LineStyle',':','Color','black')
    plot([0.8 1.2],[mean(M2M1ratio) mean(M2M1ratio)],'LineWidth',3,'LineStyle',':','Color',color_red)
    xticks([1])
    xticklabels({''})
    ylabel('M2:M1 ratio')
    ylim([0 5])
    set(gca,'FontSize',axis_fontsize,'LineWidth',axis_linewidth)
    set(fig3,'Units','Inches');
    pos = get(fig3,'Position');
    set(fig3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    printFig(fig3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS

function results_normd = plot_HeatmapStrategy(results,indNoTAMtx,lambdaVal,strategyNumber,lw,fs,figTitle)

    %normalize
    results_normd = ((results./results(1,indNoTAMtx))-1).*100;

    if(indNoTAMtx>1)
        results_normd = flip(results_normd,2); %no treatment column to the left
    end

    %Create color scale
    % Get the min and max values of the mesh, need this for scaling
    if(str2num(strategyNumber)>2)
        minLimColorbar = -100;%
        maxLimColorbar = 0;%
    else
        minLimColorbar = -0.3;%
        maxLimColorbar = 1.55;%
    end

    if(maxLimColorbar>0)
        breakpoint1 = -10^(-6);
        breakpoint2 = 10^(-6);
        largerAbsoluteValue = 1.55;%
        cstart = hex2rgb('#08519c');
        cavg1 = hex2rgb('#f7fbff');
        cavg2 = hex2rgb('#fff5f0');
        cend = hex2rgb('#a50f15');
        % Set up an interpolation, from our non-uniform colour array including the
        % breakpoint to a nice evenly spaced colour map which changes the same
        colours = [cstart; cavg1; cavg2; cend];
        %breakpoints = [minLimColorbar; breakpoint1; breakpoint2; maxLimColorbar];
        breakpoints = [-largerAbsoluteValue; breakpoint1; breakpoint2; largerAbsoluteValue];
        colours = interp1( breakpoints, colours, linspace(minLimColorbar,maxLimColorbar,100) );
    else 
        breakpoint = (minLimColorbar+maxLimColorbar)./2;
        % Set the breakpoint value
        cend = hex2rgb('#f7fbff');
        cstart = hex2rgb('#08519c');
        cavg = mean( [cstart; cend] );
        % Set up an interpolation, from our non-uniform colour array including the
        % breakpoint to a nice evenly spaced colour map which changes the same
        colours = [cstart; cavg; cend];
        breakpoints = [minLimColorbar; breakpoint; maxLimColorbar];
        colours = interp1( breakpoints, colours, linspace(minLimColorbar,maxLimColorbar,100) );
    end
    
    h_fig = figure;
    imagesc(results_normd)
    colormap( colours );
    colorbar('off')%'southoutside')
    xt = [1,length(lambdaVal)]; 
    xtlbl = {};
    strategyString = strcat('S',strategyNumber);

    yt = [1:size(results_normd,1)]; 
    ytlbl = {};
    set(gca,'LineWidth',lw,'FontSize',fs, 'XTick',xt, 'YTick',yt,'XTickLabel',xtlbl,'YTickLabel',ytlbl)
    ax = gca;
    ax.CLim = [minLimColorbar maxLimColorbar];
    if(str2num(strategyNumber)>2)
        set(gcf, 'Position',  [0,   0, 996, 179]) %S3, S4 and S5   255.5
    else
        set(gcf, 'Position',  [0,   0, 943.5, 232.5]) %S1 and S2
    end
    %colorbar
    set(h_fig,'Units','Inches');
    pos = get(h_fig,'Position');
    set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

end

function plot_pvalueline(hg1,hg2,yaxislim,lw)
    median1 = median(hg1.Data);
    median2 = median(hg2.Data);

    %find the bin the median1 falls in
    searchBin = true;
    idxBin = 1;
    while(searchBin)
        %we are in the last bin
        if(idxBin==length(hg1.BinEdges))
            searchBin = false;
        elseif((hg1.BinEdges(idxBin)<=median1)&&(median1<hg1.BinEdges(idxBin+1)))
            searchBin = false;
        else
            idxBin = idxBin+1;
        end
    end
    %get middle coordinate of the bin the median falls in
    median1_x = hg1.BinEdges(idxBin)+hg1.BinWidth/2;
    %get y coordinate of the bin the median falls in
    median1_y = hg1.Values(idxBin);

    %find the bin the median2 falls in
    searchBin = true;
    idxBin = 1;
    while(searchBin)
        %we are in the last bin
        if(idxBin==length(hg2.BinEdges))
            searchBin = false;
        elseif((hg2.BinEdges(idxBin)<=median2)&&(median2<hg2.BinEdges(idxBin+1)))
            searchBin = false;
        else
            idxBin = idxBin+1;
        end
    end
    %get middle coordinate of the bin the median falls in
    median2_x = hg2.BinEdges(idxBin)+hg2.BinWidth/2;
    %get y coordinate of the bin the median falls in
    median2_y = hg2.Values(idxBin);

    height_pvalue = (yaxislim(2)-yaxislim(1))./20;
    plot([median2,median2],[median2_y,(max(max(hg1.Values),max(hg2.Values))+height_pvalue)],...
        'LineWidth',lw,'Color',hg2.FaceColor)
    plot([median1,median1],[median1_y,(max(max(hg1.Values),max(hg2.Values))+height_pvalue)],...
        'LineWidth',lw,'Color',hg1.FaceColor)
    plot([min(median1,median2),max(median1,median2)],[(max(max(hg1.Values),max(hg2.Values))+height_pvalue),...
        (max(max(hg1.Values),max(hg2.Values))+height_pvalue)],...
        'LineWidth',lw,'Color','black')
end

