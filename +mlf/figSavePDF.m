function figSavePDF(name,factorHeight)

    L = 35;
    if nargin == 2
        H = L*factorHeight;
    else
        H = L*.75;
    end
    set(gcf,'PaperSize',[L H]*1.02,'PaperPosition',[0 0 L H])
    print('-dpdf',name)

end