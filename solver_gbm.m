%AUTHOR: BLANCHE MONGEON @ UNIVERSITÉ DE MONTRÉAL
%ARTICLE: Virtual clinical trial reveals significant clinical potential of targeting tumour-associated macrophages and microglia to treat glioblastoma  
%DATE: DECEMBER 3RD, 2024
%Solver file 
%We use fullmodelICIadmin function when we administer ICI as it considers
%the time of administration which enables continuous administration of
%ICI during 1 hour
  
function [sol p] = solver_gbm(p)

if p.NumberAdmins==1   
    %Call the solver to calculate solution over initial and only administration
    tstart = tic;
    if any(p.AdministrationTimesICI(:)== p.tspan(1)) && p.treatmentICI == 1
        sol = ode23t(@(t,y)fullmodelICIadmin(t,y,p,p.tspan(1)),p.tspan,p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
    else
        sol = ode23t(@(t,y)fullmodel(t,y,p),p.tspan,p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
    end
else
    %Call the solver to calculate solution over initial administration
    tstart = tic;
    if any(p.AdministrationTimesICI(:)== p.AdministrationTimes(1)) && p.treatmentICI == 1
        sol = ode23t(@(t,y)fullmodelICIadmin(t,y,p,p.AdministrationTimes(1)),[p.AdministrationTimes(1) p.AdministrationTimes(2)],p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
    else
        sol = ode23t(@(t,y)fullmodel(t,y,p),[p.AdministrationTimes(1) p.AdministrationTimes(2)],p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
    end
    %There are more than 2 administrations
    if numel(p.AdministrationTimes) > 2
        %Extend the solution for each subsequent administration
        for nn=3:1:numel(p.AdministrationTimes)
            p.CurrentAdmin=p.AdministrationTimes(nn-1);
            %Set ICs to be final solution in last interval
            p.initialcondition = sol.y(:,end);
            %Resection
            if p.CurrentAdmin == p.resectionTime && p.treatmentBoolean == 1
                tempKilledTumourCells = p.initialcondition(13);
                p.initialcondition = p.initialcondition.*(1-p.resectionExtent);
                %we exclude resection wehn looking at fraction of tumour
                %cells killed by treatment or by the immune systems at
                %every time point
                p.initialcondition(13) = tempKilledTumourCells;
            end
            %Immunotherapy
            if any(p.AdministrationTimesICI(:)== p.CurrentAdmin) && p.treatmentICI == 1
                p.initialcondition(6) = p.initialcondition(6)+p.DoseICI;
                if p.currentAdminICI > p.numberFirstDosesICI
                    %Doses are increased
                    p.DoseICI = 0.02*6+0.0092;
                end
                p.currentAdminICI = p.currentAdminICI+1;
            end

            %TAM targeting 
            if any(p.AdministrationTimesMacrophages(:)== p.CurrentAdmin) && p.treatmentMac == 1
                %decaying efficacy
                if(p.TAMT3decayingEfficacy)
                    p.initialcondition(14) = p.TAMT3efficacy;
                else
                    if(mod(p.CurrentWeekMacrophages,2)==1)
                        %week on: we change the value of the parameter of
                        %macrophages targeting
                        p.lambdaS1 = p.newLambdaMacValues(1);
                        p.lambdaS2 = p.newLambdaMacValues(2);
                        p.lambdaS3 = p.newLambdaMacValues(3);
                        p.lambdaS4v1 = p.newLambdaMacValues(4);
                        p.lambdaS4v2 = p.newLambdaMacValues(5);
                    else
                        %week off: no treatment for the next week we set back
                        %the values to that of no treatment
                        p.lambdaS1 = 1;
                        p.lambdaS2 = 1;
                        p.lambdaS3 = 0;
                        p.lambdaS4v1 = 1;
                        p.lambdaS4v2 = 1;
                    end
                end
                p.CurrentWeekMacrophages = p.CurrentWeekMacrophages+1;
            end

            %Chemotherapy
            if any(p.AdministrationTimesTMZ(:)== p.CurrentAdmin) && p.treatmentTMZ == 1
                p.initialcondition(8) = p.initialcondition(8)+p.DoseTMZ;
                p.currentAdminInTMZCycle = p.currentAdminInTMZCycle+1;
                if p.currentCycle == 0 && p.currentAdminInTMZCycle == p.cycle0Length
                    %we're then after in cycle 1
                    p.DoseTMZ = 150*p.BSA;
                    p.currentCycle = 1;
                    p.currentAdminInTMZCycle = 0;
                elseif p.currentCycle == 1 && p.currentAdminInTMZCycle == p.numberAdminPerTMZCycle
                    p.DoseTMZ = 200*p.BSA;
                    p.currentCycle = 2;
                    p.currentAdminInTMZCycle = 0;
                end
                
            end
            %Radiotherapy : direct tumor volume reduction
            if any(p.AdministrationTimesRT(:) == p.CurrentAdmin) && p.treatmentRT == 1
                volumeBefore = p.initialcondition(3);
                cellSurvival = exp(-p.DoseRT.*(p.alphaRT+p.betaRT*p.DoseRT));
                volumeAfter = volumeBefore*cellSurvival;
                p.initialcondition(3) = volumeAfter;
                deadVol = volumeBefore-volumeAfter;
                p.initialcondition(7) = p.initialcondition(7) - deadVol.*p.rhoG;%remove PD-L1 concentration due to cells that are now dead
                p.initialcondition(12) = p.initialcondition(12)+deadVol;%dead and damaged cell count increases
                p.initialcondition(13) = 1-cellSurvival;%Fraction of tumour cells killed by radiotherapy
            end
            if any(p.AdministrationTimesICI(:)== p.CurrentAdmin) && p.treatmentICI == 1
                sol = odextend(sol,@(t,y)fullmodelICIadmin(t,y,p,p.CurrentAdmin),p.AdministrationTimes(nn),p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
            else
                sol = odextend(sol,@(t,y)fullmodel(t,y,p),p.AdministrationTimes(nn),p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
            end
        end
        %Calculate solution from end of last admin to final time
        p.CurrentAdmin=p.AdministrationTimes(end);
        p.initialcondition=sol.y(:,end);
        %ICI
        if any(p.AdministrationTimesICI(:)== p.CurrentAdmin) && p.treatmentICI == 1
                p.initialcondition(6) = p.initialcondition(6)+p.DoseICI*p.treatmentICI;
        end
        %CT
        if any(p.AdministrationTimesTMZ(:)== p.CurrentAdmin)
                p.initialcondition(8) = p.initialcondition(8)+p.DoseTMZ.*p.treatmentTMZ;
        end
        %RT
        if any(p.AdministrationTimesRT(:) == p.CurrentAdmin) && p.treatmentRT == 1
                volumeBefore = p.initialcondition(3);
                cellSurvival = exp(-p.DoseRT.*(p.alphaRT+p.betaRT*p.DoseRT));
                volumeAfter = volumeBefore*cellSurvival;
                p.initialcondition(3) = volumeAfter;
                deadVol = volumeBefore-volumeAfter;
                p.initialcondition(7) = p.initialcondition(7) - deadVol.*p.rhoG;%remove PD-L1 concentration due to cells that are now dead
                p.initialcondition(12) = p.initialcondition(12)+deadVol;%dead and damaged cell count increases
                p.initialcondition(13) = 1-cellSurvival;%
        end
        %TAM targeting 
        if any(p.AdministrationTimesMacrophages(:)== p.CurrentAdmin) && p.treatmentMac == 1
            %decaying efficacy
            if(p.TAMT3decayingEfficacy)
                p.initialcondition(14) = p.TAMT3efficacy;
            else
                if(mod(p.CurrentWeekMacrophages,2)==1)
                    %week on: we change the value of the parameter of
                    %macrophages targeting
                    p.lambdaS1 = p.newLambdaMacValues(1);
                    p.lambdaS2 = p.newLambdaMacValues(2);
                    p.lambdaS3 = p.newLambdaMacValues(3);
                    p.lambdaS4v1 = p.newLambdaMacValues(4);
                    p.lambdaS4v2 = p.newLambdaMacValues(5);
                else
                    %week off: no treatment for the next week we set back
                    %the values to that of no treatment
                    p.lambdaS1 = 1;
                    p.lambdaS2 = 1;
                    p.lambdaS3 = 0;
                    p.lambdaS4v1 = 1;
                    p.lambdaS4v2 = 1;
                end
            end
            p.CurrentWeekMacrophages = p.CurrentWeekMacrophages+1;
            
        end
        %Stop M1 or M2 knockout if only simulated during treatment
        if isfield(p,'M1koTX')
            p.M1ko = false;
            p.q = p.q_baseline;
            p.muM1 = p.muM1_baseline;
            p.muM2 = p.muM2_baseline;
        elseif isfield(p,'M2koTX')
            p.M2ko = false;
            p.q = p.q_baseline;
            p.muM1 = p.muM1_baseline;
            p.muM2 = p.muM2_baseline;
        end
        if any(p.AdministrationTimesICI(:)== p.CurrentAdmin) && p.treatmentICI == 1
            p.lambdaT = 0;%Stop extra-recruitment of CD8+T cells if lambdaT was set to not 0
            sol = odextend(sol,@(t,y)fullmodelICIadmin(t,y,p,p.CurrentAdmin),p.tspan(2),p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
        else
            p.lambdaT = 0;
            sol = odextend(sol,@(t,y)fullmodel(t,y,p),p.tspan(2),p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
        end
    else
        %Set ICs to be final solution in last interval
        p.initialcondition=sol.y(:,end);
        p.CurrentAdmin=p.AdministrationTimes(end);
        %ICI
        if any(p.AdministrationTimesICI(:)== p.CurrentAdmin) && p.treatmentICI == 1
                p.initialcondition(6) = p.initialcondition(6)+p.DoseICI;
        end
        %CT
        if any(p.AdministrationTimesTMZ(:)== p.CurrentAdmin)
                p.initialcondition(8) = p.initialcondition(8)+p.DoseTMZ.*p.treatmentTMZ;
        end
        %RT
        if any(p.AdministrationTimesRT(:) == p.CurrentAdmin) && p.treatmentRT == 1
                volumeBefore = p.initialcondition(3);
                cellSurvival = exp(-p.DoseRT.*(p.alphaRT+p.betaRT*p.DoseRT));
                volumeAfter = volumeBefore*cellSurvival;
                p.initialcondition(3) = volumeAfter;
                deadVol = volumeBefore-volumeAfter;
                p.initialcondition(7) = p.initialcondition(7) - deadVol.*p.rhoG;%remove PD-L1 concentration due to cells that are now dead
                p.initialcondition(12) = p.initialcondition(12)+deadVol;%dead and damaged cell count increases
                p.initialcondition(13) = 1-cellSurvival;%
        end
        %Calculate solution from end of last admin to final time
        if any(p.AdministrationTimesICI(:)== p.CurrentAdmin) && p.treatmentICI == 1
            p.lambdaT = 0; %end additionnal recruitment of CD8+ T cells
            sol = odextend(sol,@(t,y)fullmodelICIadmin(t,y,p,p.CurrentAdmin),p.tspan(2),p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
        else
            p.lambdaT = 0; %end additionnal recruitment of CD8+ T cells
            sol = odextend(sol,@(t,y)fullmodel(t,y,p),p.tspan(2),p.initialcondition,odeset('Events',@(t,y) myevent(t,y,tstart),'RelTol',1e-4,'AbsTol',1e-6));
        end
    end

end


function [dydt] = fullmodel(t,y,p)

    M1 = y(1); %M1 TAMs
    if(M1<10^-6)
        M1=0;
    end
    M2 = y(2); %M2 TAMs
    G = y(3); %tumour cells
    T = y(4); %CD8+ T cells
    if(T<10^-6)
        T=0;
    end
    PU = y(5); %unbound PD-1 concentration
    I = y(6); %ICI concentration
    D = y(7); %PD-L1 concentration
    A1 = real(y(8)); %drug amount in GIT
    A2 = real(y(9)); %drug amount in blood
    A3 = real(y(10)); %drug (TMZ) amount in brain
    M = y(11); %Resident TAMs
    Dd = y(12); %Dead and damaged tumour cells 
    Tb = y(13);%Fraction of tumour cells killed by treatment (excluding resection) or by immune cells
    Z = y(14);%Decaying efficacy of enhancement of phagocytosis
    PB = y(15); %ICI-Bound PD-1 concentration
    
    dA1 = -p.ka.*A1;
    dA2 = (p.ka.*A1)-(p.CL.*A2./p.Vd)-(p.k23.*A2)+(p.k32.*A3);
    dA3 = (p.k23.*A2)-(p.k32.*A3);
    
    dZ = (log(0.5)./13).*Z;%anti-CD47 treatment
    reprogramRate = p.lambdaS4v1;
    biasFraction = p.q;
    
    TMZ =real(A3./p.Vp);%drug concentration in mg/ml
    if TMZ <0
        TMZ=0;
    end
   
    %ODE: Resident TAMs
    if(p.decreaseDeathRateOnlyM2)
        dM = p.lambdaS2.*p.lambdaM - p.aDM.*M.*(Dd)-(p.aMM).*(M1+M2).*M-p.deltaMR.*M;%
    else
        dM = p.lambdaS2.*p.lambdaM - p.aDM.*M.*(Dd)-(p.aMM).*(M1+M2).*M-p.lambdaS1*p.deltaMR.*M;%
    end
    
    %ODE: M1 TAMs
    if(p.M1ko)
        dM1 = 0;
    else
        if(p.decreaseDeathRateOnlyM2)
            dM1 = -(p.deltaM*M1)-p.muM1*M1+(p.aDM*M*(Dd)+(p.aMM).*(M1+M2).*M)*(1-(p.lambdaS4v2.*p.q))+p.lambdaS4v1.*p.muM2*M2;%
        else
            dM1 = -p.lambdaS1.*(p.deltaM*M1)-p.muM1*M1+(p.aDM*M*(Dd)+(p.aMM).*(M1+M2).*M)*(1-(p.lambdaS4v2.*biasFraction))+reprogramRate.*p.muM2*M2;%
        end
    end
    
    %ODE: M2 TAMs
    if(p.M2ko)
        dM2 = 0;
    else
        dM2 = p.muM1*M1 - p.lambdaS1*p.deltaM*M2+(p.aDM*M*(Dd)+(p.aMM).*(M1+M2).*M)*p.lambdaS4v2.*biasFraction-reprogramRate.*p.muM2*M2;%
    end

    %ODE: Tumour cells
    if(~p.TAMT3decayingEfficacy)
        phagocytosisMac = p.lambdaS3.*p.dPMac.*G.*(M1);%;constant value
    else
        phagocytosisMac = Z.*p.dPMac.*G.*(M1);%decaying efficacy of anti-CD47 treatment
    end
    Dkill = p.alpha.*((T)./(T+p.k1)).*G;%Tumour cells killed by CD8+ T cells
    ETMZ = (p.ImaxTMZ*TMZ^p.hTMZ)/(TMZ^p.hTMZ+p.IC50TMZ^p.hTMZ)*p.treatmentTMZ;%Effect of TMZ
    dG = G.*p.beta.*log(p.K./G)-Dkill-phagocytosisMac-ETMZ*G;

    %ODE: PD-L1 Concentration
    dD = p.rhoG*dG + p.rhoM*(M1+M2+M);%

    %fraction of tumour cells
    %We retrieve the last value (Tb) so that we'll obtain the fraction at
    %every timepoint when evaluating the solution.
    dTb = -Tb+(Dkill+ETMZ*G+phagocytosisMac)./G;
    
    %ODE: CD8+ T cells
    Dinh = (1./(1+p.inh.*((PU.*D)./p.kYQ)));%Inhibition of CD8+ T cell activation due to PD-1/PD-L1 complexes
    ratioM1=M1/(M2+M1); %M1:total activated TAMs ratio
    ratioM2=M2/(M1+M2); %M2:total activated TAMs ratio
    %CD8+ T cells activated by TAMs
    if(p.M1ko)
        activatedTcells_mac = 0;
    else
        activatedTcells_mac =  p.n1.*(ratioM1/(p.n2+ratioM2));
    end
    %CD8+ T cells activated by tumour cells
    activatedTcells_tumour =  T.*p.aAT.*G./(G+p.hT);%
    activatedTcells = (activatedTcells_mac+activatedTcells_tumour).*Dinh;
    deadTcells = (p.deltaT*T);
    dT = p.lambdaT-deadTcells+activatedTcells;%
    
    %ODE: Unbound PD-1 receptors
    dPU = p.rhoT.*activatedTcells-p.muP.*PU.*I-PU/(PU+PB)*p.rhoT*deadTcells;

    %ODE: Bound PD-1 receptors
    dPB = p.muP.*PU.*I - PB/(PU+PB)*p.rhoT*deadTcells;
    
    %ODE: ICI
    dI =  - p.deltaI*I-p.muI*PU*I;
    
    %ODE: Dead and damaged cells
    dDd = -p.deltaDM.*Dd*(M1+M2+M)-p.dDd.*Dd+Dkill+ETMZ*G;%
    
    %        1   2   3  4  5  6  7  8   9   10 11  12  13  14  15
    dydt = [dM1;dM2;dG;dT;dPU;dI;dD;dA1;dA2;dA3;dM;dDd;dTb;dZ;dPB];

   
end

function [dydt] = fullmodelICIadmin(t,y,p,initialTime)

    M1 = y(1); %M1 TAMs
    if(M1<10^-6)
        M1=0;
    end
    M2 = y(2); %M2 TAMs
    G = y(3); %tumour cells
    T = y(4); %CD8+ T cells
    if(T<10^-6)
        T=0;
    end
    PU = y(5); %unbound PD-1 concentration
    I = y(6); %ICI concentration
    D = y(7); %PD-L1 concentration
    A1 = real(y(8)); %drug amount in GIT
    A2 = real(y(9)); %drug amount in blood
    A3 = real(y(10)); %drug (TMZ) amount in brain
    M = y(11); %Resident macrophages
    Dd = y(12); %Dead and damaged cells 
    Tb = y(13);%Fraction of tumour cells killed by treatment (excluding resection) or by immune cells
    Z = y(14);%Decaying efficacy of enhancement of phagocytosis
    PB = y(15); %ICI-Bound PD-1 concentration
    
    %ODE: TMZ amount in GIT (A1), Plasma (A2) and CSF (A3)
    dA1 = -p.ka.*A1;
    dA2 = (p.ka.*A1)-(p.CL.*A2./p.Vd)-(p.k23.*A2)+(p.k32.*A3);
    dA3 = (p.k23.*A2)-(p.k32.*A3);
    
    %ODE: anti-CD47 treatment
    dZ = (log(0.5)./13).*Z;
    reprogramRate = p.lambdaS4v1;
    biasFraction = p.q;
    
    TMZ =real(A3./p.Vp);%drug concentration in mg/ml
    if TMZ <0
        TMZ=0;
    end
   
    %ODE: Resident TAMs
    if(p.decreaseDeathRateOnlyM2)
        dM = p.lambdaS2.*p.lambdaM - p.aDM.*M.*(Dd)-(p.aMM).*(M1+M2).*M-p.deltaMR.*M;%
    else
        dM = p.lambdaS2.*p.lambdaM - p.aDM.*M.*(Dd)-(p.aMM).*(M1+M2).*M-p.lambdaS1*p.deltaMR.*M;%
    end
    
    %ODE: M1 TAMs
    if(p.M1ko)
        dM1 = 0;
    else
        if(p.decreaseDeathRateOnlyM2)
            dM1 = -(p.deltaM*M1)-p.muM1*M1+(p.aDM*M*(Dd)+(p.aMM).*(M1+M2).*M)*(1-(p.lambdaS4v2.*p.q))+p.lambdaS4v1.*p.muM2*M2;%
        else
            dM1 = -p.lambdaS1.*(p.deltaM*M1)-p.muM1*M1+(p.aDM*M*(Dd)+(p.aMM).*(M1+M2).*M)*(1-(p.lambdaS4v2.*biasFraction))+reprogramRate.*p.muM2*M2;%
        end
    end
    
    %ODE: M2 TAMs
    if(p.M2ko)
        dM2 = 0;
    else
        dM2 = p.muM1*M1 - p.lambdaS1*p.deltaM*M2+(p.aDM*M*(Dd)+(p.aMM).*(M1+M2).*M)*p.lambdaS4v2.*biasFraction-reprogramRate.*p.muM2*M2;%
    end

    %ODE: Tumour cells
    if(~p.TAMT3decayingEfficacy)
        phagocytosisMac = p.lambdaS3.*p.dPMac.*G.*(M1);%;constant value
    else
        phagocytosisMac = Z.*p.dPMac.*G.*(M1);%decaying efficacy of anti-CD47 treatment
    end
    Dkill = p.alpha.*((T)./(T+p.k1)).*G;%Tumour cells killed by CD8+ T cells
    ETMZ = (p.ImaxTMZ*TMZ^p.hTMZ)/(TMZ^p.hTMZ+p.IC50TMZ^p.hTMZ)*p.treatmentTMZ;%Effect of TMZ
    dG = G.*p.beta.*log(p.K./G)-Dkill-phagocytosisMac-ETMZ*G;

    %ODE: PD-L1 Concentration
    dD = p.rhoG*dG + p.rhoM*(M1+M2+M);%

    %fraction of tumour cells
    %We retrieve the last value (Tb) so that we'll obtain the fraction at
    %every timepoint when evaluating the solution.
    dTb = -Tb+(Dkill+ETMZ*G+phagocytosisMac)./G;
    
    %ODE: CD8+ T cells
    Dinh = (1./(1+p.inh.*((PU.*D)./p.kYQ)));%Inhibition of CD8+ T cell activation due to PD-1/PD-L1 complexes
    ratioM1=M1/(M2+M1); %M1:total activated TAMs ratio
    ratioM2=M2/(M1+M2); %M2:total activated TAMs ratio
    %CD8+ T cells activated by TAMs
    if(p.M1ko)
        activatedTcells_mac = 0;
    else
        activatedTcells_mac =  p.n1.*(ratioM1/(p.n2+ratioM2));
    end
    %CD8+ T cells activated by tumour cells
    activatedTcells_tumour =  T.*p.aAT.*G./(G+p.hT);%
    activatedTcells = (activatedTcells_mac+activatedTcells_tumour).*Dinh;
    deadTcells = (p.deltaT*T);
    dT = p.lambdaT-deadTcells+activatedTcells;%
    
    %ODE: Unbound PD-1 receptors
    dPU = p.rhoT.*activatedTcells-p.muP.*PU.*I-PU/(PU+PB)*p.rhoT*deadTcells;

    %ODE: Bound PD-1 receptors
    dPB = p.muP.*PU.*I - PB/(PU+PB)*p.rhoT*deadTcells;
    
    %ODE: ICI
    %time units: days, we want to administer ICI during one hour
    if(t>=initialTime && t<(initialTime+(1/24)))
        dI =  p.DoseICI - p.deltaI*I-p.muI*PU*I;
    else
        dI =  - p.deltaI*I-p.muI*PU*I;
    end
    
    %ODE: Dead and damaged cells
    dDd = -p.deltaDM.*Dd*(M1+M2+M)-p.dDd.*Dd+Dkill+ETMZ*G;%
    
    %        1   2   3  4  5  6  7  8   9   10 11  12  13  14  15
    dydt = [dM1;dM2;dG;dT;dPU;dI;dD;dA1;dA2;dA3;dM;dDd;dTb;dZ;dPB];

   
end

end