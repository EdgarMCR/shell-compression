function[] = runScript(Np, H)
% fprintf('%.4f \n', H);
global counter
counter=1;
leng=100;
Step=(H-0.01)/leng;
for ii=1:leng
    [X,p]=initialGuess(Np, H);
    problem(Np,  H, 0.00001, X, 'sol.txt', 1, 'sol.txt');    
    H=H-Step;
    counter=counter+1;
end

% [X,p]=initialGuess_closeToSpherical(Np, H);


% function[] = runScript(leng, H)
% len=leng;
% N=20;
% 
% r=zeros([len, 1]);
% nPoints=zeros([len, 1]);
% for ii=1:len
%     Na=N*2^(ii-1);
%     fprintf('Step %d with %d points \n', ii, Na);
%     nPoints(ii,1)=Na;
%     [r(ii,1)] = testResidue(Na);
% end
% 
% h=figure(); %set(gcf,'Visible', 'off'); 
% plot(nPoints, r,  'ro', 'MarkerFaceColor', 'r')
% titleString=sprintf('Norm of Residue with guessed Solution');
% title(titleString,'FontSize', 16)
% xlabel('Number of points','FontSize', 16) 
% ylabel('Norm of Residue','FontSize', 16)
% 
% filename=sprintf('residuetest\\H=%.2f_len=%d.png',H, len);
% print(h,filename, '-dpng')
% 
% h2=figure(); %set(gcf,'Visible', 'off'); 
% loglog(nPoints, r,  'ro', 'MarkerFaceColor', 'r')
% titleString=sprintf('Norm of Residue with guessed Solution');
% title(titleString,'FontSize', 16)
% xlabel('Number of points','FontSize', 16)
% ylabel('Norm of Residue','FontSize', 16)
% filename=sprintf('residuetest\\Loglog_H=%.2f_len=%d.png',H, len);
% print(h2,filename, '-dpng')
fprintf('\n');
end