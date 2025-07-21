function figSavePNG(name)

set(gcf,'Color','none')
set(gca,'Color','none')
saveas(gcf,name,'png')