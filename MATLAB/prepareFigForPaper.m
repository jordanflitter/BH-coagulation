function prepareFigForPaper(correctYLabel, fig)

if nargin < 1, correctYLabel = 0; end
if nargin < 2, fig = gcf;         end

% Get widest axes
ax = findobj(fig,'type','axes');

p = [];

for i=1:length(ax)
    p(i,:) = get(ax(i),'position');
end

[~,I] = max(p(:,3));
ax = ax(I);
%%%%

ax_u = get(ax, 'units');

set(ax, 'units',      'centimeters', ...
        'box',        'on');

p  = get(ax, 'Position');   % [left,    bottom,    width,    height]
ti = get(ax, 'TightInset'); % [dx_left, dy_bottom, dx_right, dy_top]

w = p(3)+ti(3)+ti(1);
h = p(4)+ti(4)+ti(2);

li = ti + repmat(0.2, size(ti));

if (correctYLabel)
    yl = get(ax, 'ylabel');
    set(yl, 'units', 'centimeters');
    ex = get(yl, 'extent');
    
    f = 1.5;
    w = w + ex(4)*f;
    li(1) = li(1) + ex(4)*f;
end

set(ax, 'LooseInset', li);
        
set(fig, 'PaperUnits',        'centimeters', ...
	     'PaperSize',         [w h],         ...
         'PaperPositionMode', 'manual',      ...
         'PaperPosition',     [0 0 w h]);

set(ax, 'units', ax_u);

end

