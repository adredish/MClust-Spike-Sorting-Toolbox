function FindBestAxes(self)

% Finds the best pair of axes using L-ratio for one cell against the rest
% of the spikes or for two cells if two foci

axisSetX = get(self.xAxisLB, 'string'); nX = length(axisSetX);
axisSetY = get(self.yAxisLB, 'string'); nY = length(axisSetY);

bestLR = inf;
bestX = nan;
bestY = nan;

for iX = 1:nX
    if strcmp(axisSetX(iX),'_Time') % don't do time
        continue;
    end
    xFeat = self.Features{iX};
    FD(:,1) = xFeat.GetData();
    
    for iY = (iX+1):nY
        if strcmp(axisSetY(iY),'_Time') % don't do time
            continue;
        end
        yFeat = self.Features{iY};
        FD(:,2) = yFeat.GetData();
                        
        if ~isequal(self.whoHasComparison, self.whoHasFocus)
            LR = MClust.ClusterQuality.L_Ratio(FD, self.whoHasFocus.GetSpikes(), self.whoHasComparison.GetSpikes());
        else
            LR = MClust.ClusterQuality.L_Ratio(FD, self.whoHasFocus.GetSpikes());
        end            
        
        if LR.Lratio < bestLR
            bestLR = LR.Lratio;
            bestX = iX;
            bestY = iY;
        end
    end
end

if ~isnan(bestX) && ~isnan(bestY)
    self.set_xAxis(bestX);
    self.set_yAxis(bestY);
    self.RedrawAxes();
end