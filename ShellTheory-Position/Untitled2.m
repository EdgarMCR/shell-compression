theta=-2*pi:pi/512:2*pi;
length(theta)

figure()
plot(theta, sin(theta), 'r');
hold on;
plot(theta, cos(theta), 'b');

% r=abs(cos(5*(theta-pi/2)));

% for i=1:length(theta)
%     x(i)=r(i)*cos(theta(i));
%     y(i)=r(i)*sin(theta(i));
% end
% 
% plot(x,y, 'o')