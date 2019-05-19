
function fh = schedulePlot(schedule,C,xlabels,ylabels)
%SCHEDULEPLOT Visualize a binary schedule
%   Given a matrix SCHEDULE, create a grid visualization.

if nargin < 3 || isempty(C)
    C = [0.8 1 0.8];
end

nx = size(schedule,2)+1;
ny = size(schedule,1)+1;
tol = 1E-3;
schedule(schedule<=tol) = 0;

[ys,xs] = ind2sub(size(schedule),find(schedule));

Y = [ys'; ys'; ys'-1; ys'-1];
X = [xs'; xs'-1; xs'-1; xs'];

Y = ny - 1 - Y;

fh = figure('Renderer','painters','Position',[100 100 1000 750]);

C = schedule(schedule>tol)';
fill(X,Y,C); % Visualization of on/off
colormap winter

ax1 = gca;
set(ax1,'XLim',[0 nx-1],'YLim',[0 ny-1],...
        'XTick',0:nx,'YTick',0:ny,...
        'XGrid','on','YGrid','on',...
        'XTickLabel',[],'YTickLabel',[],...
        'TickLength',[0 0],...
        'GridLineStyle','-');

ax2 = axes;
xlim = get(ax1,'XLim');
ylim = get(ax1,'YLim');
xt = get(ax1,'XTick');
yt = get(ax1,'YTick');
set(ax2,'Color','none',...
        'XLim',xlim,'YLim',ylim,...
        'XTick',xt(1:end-1)+diff(xt)/2,'YTick',yt(1:end-1)+diff(yt)/2,...
        'TickLength',[0 0]);
    
if nargin > 2 && ~isempty(xlabels)
    set(ax2,'XTickLabel',xlabels);
else
    nxlabels = 50;
    allxlabels = cellstr(num2str((1:(nx-1))'));
    if length(allxlabels)>nxlabels
        spacing = ceil(length(allxlabels)/50);
        idxs = ones(length(allxlabels),1);
        keepers = 1:spacing:(nx-1);
        idxs(keepers) = 0;
        idxs = logical(idxs);
        allxlabels(idxs) = repmat({''},sum(idxs),1);
    end
    set(ax2,'XTickLabel',allxlabels);
end
ax2.XTickLabelRotation = 90;

if nargin > 3 && ~isempty(ylabels)
    set(ax2,'YTickLabel',ylabels);
else
    set(ax2,'YTickLabel',num2str((1:(ny-1))'));
end

xlabel('Time (hrs)');
ylabel('Generator');
title('Generator Schedule','FontSize',12,'FontWeight','bold');

c = colorbar;
set(c,'Position',[0.9250 0.1071 0.0476/2 0.8167]);
cyl = get(c,'YLim');
set(c,'YTick',cyl,'YTickLabel',{'Min','Max'});

linkaxes([ax1,ax2]);
end
