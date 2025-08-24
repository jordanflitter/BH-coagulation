function out = MakeStairPlot(C,par,names,legitems,sc,fsize,c)
% out = MakeStairPlot(C,par,names,sc,c)
% Makes stairstep Fisher matrix plot
%
% INPUTS:
% C          Covariance matrix
%
% par        Maximum likelihood parameter values from Fisher matrix
%
% labels     cell array of parameter names for axis labels
%
% sc         Scale factor for each parameter.  Parameter values and
%            uncertanties are all multiplied by sc. Defaults to 1 for all
%            parameters
%
% colorfile  Filename for array of line colors. Defaults to
%            'DefaultColors.dat'.  Should be formatted as a (number of
%            colors)x6 matrix where the first three columns give rgb values
%            for 2-sigma ellipses and the second three give rgb values for
%            the 1-sigma ellipses.
%
% legitems   Cell array of legend names.  Only displays for 1-sigma
%            ellipses.
%
% OUTPUT:
% out       vector containing handle for figure followed by handles for
%           each axis

if nargin<7
    colorfile = 'DefaultColors.dat';
    c = dlmread(colorfile);
end
c1 = c(:,1:3);
c2 = c(:,4:6);

FillAlpha = 1;

if nargin<5
    sc = ones(size(par));
end

if nargin<6
    fsize=12;
end

[N,~,N2] = size(C);

if N~=numel(par)
    error('Number of paramter values does not match covariance matrix')
end

[nc,~] = size(c);

if nc<N2
    c1 = repmat(c1,[ceil(N2/nc),1]);
    c2 = repmat(c2,[ceil(N2/nc),1]);
else

f = figure;

% Set margins on sides of figure
tm = .08;
bm = .12;
lm = .13;
rm = .05;

% total figure size
totalheight = 1-tm-bm;
totalwidth = 1-lm-rm;

% Number of rows and columns
Nrc = N-1;

% Size of each plot
ht = totalheight/Nrc;
wd = totalwidth/Nrc;

% Number of plots
Nplot = N*(N+1)/2;

% Set positions of each plot
% xpos = zeros(size(Nplot));
% ypos = xpos;
ind = 1;
s = zeros(size(Nplot));

% Get error on each parameter to set axis ranges
err = ones(1,N);
for ii=1:N
    if N2>1
        err(ii) = max(sqrt(C(ii,ii,:)));
    else
        err(ii) = sqrt(C(ii,ii));
    end
end

for ii=1:Nrc % count across rows
    for jj=1:(Nrc-ii+1) % count across columns
        
        % Set axis positions
        xpos = lm+wd*(ii-1)+.000001;
        ypos = bm+ht*(jj-1)+.000001;
        
        nj = N-jj+1;
        
        C2 = C([ii,nj],[ii,nj],1); % Marginalized covariance matrix
        
        % Create axes
        s(ind) = subplot('position',[xpos,ypos,wd,ht]);
        
        % Plot Fiducial parameters
%         fx = plot([par(ii),par(ii)],[par(nj)-10*err(nj),...
%             par(nj)+10*err(nj)],':');
%         set(fx,'LineWidth',1,'Color',[.2,.2,.2]);
%         hold on
%         fy = plot([par(ii)-10*err(ii),par(ii)+10*err(ii)],...
%             [par(nj),par(nj)],':');
%         set(fy,'LineWidth',1,'Color',[.2,.2,.2]);
        
        sc2 = [sc(ii),sc(nj)];
        
        % Plot error ellipses
        %a1 = error_ellipse(C2,[par(ii),par(nj)],'conf',.68,'style','b');
        a1 = error_ellipse2(C2,[par(ii),par(nj)],'conf',.95,...
                    'scale',sc2,'Fill',1,'FillColor',c1(1,:),...
                    'EdgeColor',c2(1,:),'FillAlpha',FillAlpha);
        hold on
        %a2 = error_ellipse(C2,[par(ii),par(nj)],'conf',.95,'style','b');
        a2 = error_ellipse2(C2,[par(ii),par(nj)],'conf',.68,...
                    'scale',sc2,'Fill',1,'FillColor',c2(1,:),...
                    'EdgeColor',c2(1,:),'FillAlpha',FillAlpha);
        legentries = a2;
        hold off
        set(a1,'LineWidth',1);
        set(a2,'LineWidth',1);
        
        % Multiple ellipses?
        if N2 > 1
            hold on
            for kk=2:N2
                C2 = C([ii,nj],[ii,nj],kk); % Marginalized covariance matrix
                a1 = error_ellipse2(C2,[par(ii),par(nj)],'conf',.95,...
                    'scale',sc2,'Fill',1,'FillColor',c1(kk,:),...
                    'EdgeColor',c2(kk,:),'FillAlpha',FillAlpha);
                %a2 = error_ellipse(C2,[par(ii),par(nj)],'conf',.95);
                a2 = error_ellipse2(C2,[par(ii),par(nj)],'conf',.68,...
                    'scale',sc2,'Fill',1,'FillColor',c2(kk,:),...
                    'EdgeColor',c2(kk,:),'FillAlpha',FillAlpha);
                
                legentries = [legentries,a2];
                
                set(a2,'LineWidth',1);
                set(a1,'LineWidth',1);
            end
            hold off
            
        end
        
        % Set axis properties
        set(gca,'FontSize',fsize);
        set(gca,'XLim',sc(ii)*[par(ii)-4*err(ii),par(ii)+4*err(ii)]);
        set(gca,'YLim',sc(nj)*[par(nj)-4*err(nj),par(nj)+4*err(nj)]);
        box on
        
        if ii ~= 1
            set(gca,'YTick',[]);
        else
            ylabel(names(nj));
        end
        
        if jj ~= 1
            set(gca,'XTick',[]);
        else
            xlabel(names(ii));
        end
        
        set(s(ind),'Layer','top');
        ind = ind+1;
    end
end

out = [f,s];

if nargin>=4
    leg = legend(fliplr(legentries),legitems);
    legPos = leg.Position;
    legwd = legPos(3);
    leght = legPos(4);
    legx = 1-rm-legwd;
    legy = 1-tm-leght;
    leg.Position = [legx,legy,legwd,leght];
    leg.Box = 'off';
    set(leg,'interpreter','latex')
    set(leg,'FontSize',fsize)
    out = [out,leg];
end


end
        
        
        
