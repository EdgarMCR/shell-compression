x = 0:0.01:20;
y1 = 200*exp(-0.05*x).*sin(x);
y2 = 0.8*exp(-0.5*x).*sin(10*x);
[AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
set(get(AX(1),'Ylabel'),'String','Slow Decay')
set(get(AX(2),'Ylabel'),'String','Fast Decay')
xlabel('Time (\musec)')
title('Multiple Decay Rates')
set(H1,'LineStyle','--')
set(H2,'LineStyle',':')

% axes(AX(1))
% hold on
% plot(randn(20,1)*100)
% 
% axes(AX(2))
% hold on
% plot(rand(20,1),'y--')
% plot(-rand(20,1),'k.')
% plot(rand(20,1)-0.5,'ro')