function printFig(figHandle)
    if nargin<1
        figHandle = gcf;
    end
    set(figHandle,'Units','Inches');
    pos = get(figHandle,'Position');
    set(figHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
end