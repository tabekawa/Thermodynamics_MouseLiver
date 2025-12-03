function [b] = BarGraphWithErrors(x, y, By)

ci = conf_int(By,0.05);
% Calculate upper and lower error bounds
upperBound = ci(:,2) - y;
lowerBound = y - ci(:,1);

% Plot the bar graph
b = bar(x, y);
b.EdgeColor = 'none';
hold on;

% Plot the error bars
errorbar(x, y, lowerBound, upperBound,'k','LineStyle','none');
% errorbar(x, y, lowerBound, upperBound,'k','LineStyle','none', 'CapSize',10); % EC ratio of pathway

% For ddrG
set(gca, 'ylim',[-15,15], 'Ygrid', 'on', 'YMinorGrid', 'on');

% % For total EC
% set(gca, 'Ygrid', 'on', 'YMinorGrid', 'on', 'Yscale', 'log');

hold off;

end