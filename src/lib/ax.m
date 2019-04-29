function ax(xlab, ylab, xlim, ylim)
set(gcf, 'Position', [1     1   473   261]);
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 8)
grid on
if nargin > 0
	xlabel(xlab);
end

if nargin > 1
	ylabel(ylab);
end

if nargin > 2
	set(gca, 'XLim', xlim);
end

if nargin > 3
	set(gca, 'YLim', ylim);
end

