fig=figure(4)
set(fig,'Units','Inches');    
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'11_13_23_mccs2.pdf','-dpdf','-r1200')