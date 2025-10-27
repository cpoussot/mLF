function figSavePDF(name,factorHeight)

    L = 35;
    if nargin == 2
        H = L*factorHeight;
    else
        H = L*.75;
    end
    set(gcf,'PaperSize',[L H]*1,'PaperPosition',[0 0 L H])
    %set(gcf,'PaperSize',[L H]*1)
    % L = 35;
    % H = 20;
    %set(gcf,'PaperSize',[L H],'PaperPosition',[0 0 L H])
    print('-dpdf',name)

end