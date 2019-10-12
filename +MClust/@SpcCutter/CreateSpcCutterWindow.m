function CreateSpcCutterWindow(self)
% KKwikCutter.CreateCutterWindow

MCS = MClust.GetSettings();
MCD = MClust.GetData();

self.CreateCutterWindow@MClust.Cutter();  % call superclass to build initial window

%--------------------------------
% constants to make everything identical

uicHeight = self.uicHeight;
uicWidth  = self.uicWidth;
uicWidth0 = self.uicWidth0;
uicWidth1 = self.uicWidth1;
XLocs = self.XLocs;
dY = self.dY;
YLocs = self.YLocs;


% ----- Clusters
% These were made by superclass call, but need to re-make them for 
% hierarchical SPC
%-----------------------------------
disp('Turning the clusterpanel visibility off')
set(self.self.clusterPanel, 'Visible', 'off'); %hide old cluster panel 
set(self.uiScrollbar, 'Visible', 'off'); %and scrollbar
disp('it should be off...')

