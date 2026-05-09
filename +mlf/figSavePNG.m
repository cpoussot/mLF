function figSavePNG(name,factorHeight)

L = 35;
if nargin == 2
    H = L*factorHeight;
else
    H = L*.75;
end
set(gcf,'PaperSize',[L H]*1,'PaperPosition',[0 0 L H])
set(gcf,'Color','none')
set(gca,'Color','none')
saveas(gcf,name,'png')