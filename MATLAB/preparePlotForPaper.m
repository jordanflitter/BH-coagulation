function preparePlotForPaper(fontSize, lineWidth, AR, fig)

if nargin < 1, fontSize = 25; end
if nargin < 2, lineWidth = 2; end
if nargin < 3, AR = 3/4;     end
if nargin < 4, fig = gcf;     end

%AR = 9/16;%8/6;%1;%9/16;
W = 36;
set(fig, 'units', 'centimeters', 'position', [5 5 W W*AR]);

% Get widest axes
ax = findobj(fig,'type','axes');

p = [];

for i=1:length(ax)
    p(i,:) = get(ax(i),'position');
end

[~,I] = max(p(:,3));
ax = ax(I);
%%%%

set(ax, 'fontsize', fontSize);

lines = findobj(ax,'Type','line');

for i = 1:length(lines)
    set(lines(i), 'linewidth', lineWidth);
end

t = get(ax, 'title');
set(t, 'fontsize', fontSize);

xl = get(ax, 'xlabel');
set(xl, 'fontsize', fontSize);

yl = get(ax, 'ylabel');
set(yl, 'fontsize', fontSize);

leg = findobj(fig,'Type','axes','Tag','legend');
set(leg, 'fontsize', fontSize);

end

