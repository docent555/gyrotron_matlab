function [lHandleFr, lHandleFa, hF] = makeFig(z, t)

hF = figure;
axFr = subplot(2,1,1);

axFr.XLim = [0 t(end)];
lHandleFr = line(axFr, nan, nan); %# Generate a blank line and return the line handle
lHandleFr.Color = 'black';
lHandleFr.LineWidth = 1;
axFr.FontSize = 12;
axFr.XLabel.String = 't';
axFr.YLabel.String = '|F|_{max}';

axFa = subplot(2,1,2);

lHandleFa = line(axFa, nan, nan); %# Generate a blank line and return the line handle
lHandleFa.Color = 'black';
lHandleFa.LineWidth = 1;
axFa.FontSize = 12;
axFa.XLabel.String = 'z';
axFa.YLabel.String = '|F|';

lHandleFa.XData = z;
end