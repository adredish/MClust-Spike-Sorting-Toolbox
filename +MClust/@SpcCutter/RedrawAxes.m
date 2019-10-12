function RedrawAxes(self, ~, ~)

MCS = MClust.GetSettings();

% Something has changed in the control window, redraw as necessary...

if self.get_redrawStatus()
    
    % window for display
    if isempty(self.CC_displayWindow) || ~ishandle(self.CC_displayWindow)
        % create new drawing figure
        self.CC_displayWindow = ...
            figure('Name', 'Cluster Cutting Window',...
            'NumberTitle', 'off', ...
            'Tag', 'CHDrawingAxisWindow', ...
			'Position',MCS.CHDrawingAxisWindow_Pos);
        MCS.PlaceWindow(self.CC_displayWindow); % ADR 2013-12-12
    else
        % figure already exists -- select it
        figure(self.CC_displayWindow);
    end
        
    % get axes
    xFeat = self.Features{self.get_xAxis};
    yFeat = self.Features{self.get_yAxis};
    
    % get FD data
    xFD = xFeat.GetData();
    yFD = yFeat.GetData();
    
    clf;
    ax = axes('Parent', self.CC_displayWindow, ...
        'XLim', [min(xFD) max(xFD)], 'YLim', [min(yFD) max(yFD)]);
    hold on;
    
    % Recolor according to how many clusters there are
    AllClusters = self.getClusters();
    colors = hsv;
    ncd = 0; %number of clusters being displayed
    for iC = 2:length(AllClusters)  % count each cluster that isn't the all-points cluster
        if ~AllClusters{iC}.hide    % and isn't hidden
            ncd = ncd + 1;
        end
    end

    % Draw the clusters
    cncd = 1; %current number of clusters displayed
    for iC = 1:length(AllClusters)        
        if ~AllClusters{iC}.hide %for each cluster that isn't hidden
            if iC==1 %if this is the all-points cluster,
                AllClusters{iC}.color = [0, 0, 0]; %make it black
            else
                AllClusters{iC}.color = colors(round(size(colors,1)/ncd*cncd),:); %recolor
                cncd = cncd + 1;
            end
            AllClusters{iC}.PlotSelf(xFD, yFD, ax, xFeat, yFeat);  %plot it
        end
    end 
	    
    xlabel(xFeat.name,'interpreter','none');
    ylabel(yFeat.name,'interpreter','none');
    zoom on
end

end