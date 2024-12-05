%deaths must be a 1xN array where N is the number of patients
%deaths(1,i) indicates the time of death of the ith patient
%if deaths(1,i) = -1, then the patient is not dead by the end of the
%simulations
%options: options.color_line and options.line_width
function KM_curve(deaths,varargin)

    options = KMParseInput(varargin{:});

    N = size(deaths,2);
    init = N;
    deathTimes = zeros(1,1);
    survival = ones(1,1);%everyone is alive at time 0
    temp = sort(deaths(deaths>0));
    for nnDeath=1:length(temp)
        deathTimes(1,end+1) = temp(nnDeath);
        init = init-1;
        survival(1,end+1) = init/N;
    end
    hold on
    survival = survival.*100;
    stairs(deathTimes,survival,'color',options.LineColor,'LineWidth',options.LineWidth,'LineStyle',options.LineStyle)
    if(options.ShowMedian)
        medianSurvival = deathTimes(find(survival<=50,1,'first'));
        if ~(isempty(medianSurvival))
            plot([medianSurvival medianSurvival],[0 50],'color',options.LineColor,'LineWidth',options.LineWidth-1)
        else
            plot(0,0)
        end
    end
    hold off

    function params = KMParseInput(varargin)
        %Parse input and set defualt values
        p = inputParser;
        
        p.addParameter('ShowMedian',true);  
        p.addParameter('LineStyle','-')
        p.addParameter('LineWidth',2)
        p.addParameter('LineColor',[0 0 0])
   
        parse(p,varargin{:});
        params = p.Results;

    end

end