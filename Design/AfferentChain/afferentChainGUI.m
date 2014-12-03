function afferentChainGUI
% GUI for plotting r* and psi* of afferent chains (one-to-one connections)
% consisting of fully expanded canonical models of Hopf bifurcation.
%
% Change parameter values by moving the sliders or by entering numbers
% in the text boxes. Min and max values for the sliders can also be
% changed.

hf = figure;
set(hf,'Name','Afferent Chain GUI','Visible','off','Toolbar','figure',...
  'Position',[100,100,1000,670]);
initialRun(hf)
movegui(hf,'center')
set(hf,'Visible','on')

% =========================================================================
function initialRun(hf)

global paramLim

% Determine font size
platform = computer;
vMATLAB = version('-release');
if strcmp(platform(1:3),'MAC') && strcmp(vMATLAB(1:4),'2011')
  fsize = 11;
else
  fsize = 10;
end

% Set initial values
numLayer = 2; % number of layers

param = [-.12       0; % alpha
            0 -150000; % beta1
            0       0; % beta2
            0       0]; % epsilon
force = .002;
slcLayer = 1; % selected layer
numOct = 1;
maxdf = 1;
layerString = [];
for nl = 1:numLayer
  layerString = [layerString 'Layer ' num2str(nl) '|'];
end
layerString = layerString(1:end-1);

paramLim = cell(numLayer,1);
for nl = 1:numLayer
  paramLim{nl} = [ [-1 1]*max(1,ceil(2*abs(param(1,nl))));
                    [-1 1]*max(1,ceil(2*abs(param(2,nl))));
                    [-1 0]*max(1,ceil(2*abs(param(3,nl))));
                    [ 0 2]*max(1,ceil(2*abs(param(4,nl)))) ];
end
forceLim = [0 1]*max(1,ceil(2*abs(force)));

% Create axes
clf
h.axr = axes('Parent',hf,'Units','normalized','Position',[.05 .55 .5 .4]);
h.axpsi = axes('Parent',hf,'Units','normalized','Position',[.05 .1 .5 .4]);

% Create a panel for layer control
panel1 = uipanel(hf,'Position',[.57 .4 .4 .55]);
uicontrol(panel1,'Style','text','String','Number of Layers:',...
  'FontSize',fsize,'Units','normalized','Position',[.02 .85 .25 .1],...
  'HorizontalAlignment','right');
h.numLayer = uicontrol(panel1,'Style','edit','String',numLayer,...
  'FontSize',fsize,'Units','normalized','Position',[.28 .9 .07 .07]);
uicontrol(panel1,'Style','text','String','Control Parameters for ',...
  'FontSize',fsize,'Units','normalized','Position',[.35 .85 .4 .1],...
  'HorizontalAlignment','right');
h.slcLayer = uicontrol(panel1,'Style','popup','String',layerString,...
  'Value',1,'FontSize',fsize,'Units','normalized',...
  'Position',[.75 .86 .25 .1]);

% Create a panel for parameter control
uipanel(panel1,'Position',[.02 .19 .96 .66]);
uicontrol(panel1,'Style','text','String','alpha',...
  'FontSize',fsize,'Units','normalized','Position',[.04 .75 .15 .05],...
  'HorizontalAlignment','left');
uicontrol(panel1,'Style','text','String','beta1',...
  'FontSize',fsize,'Units','normalized','Position',[.04 .6 .15 .05],...
  'HorizontalAlignment','left');
uicontrol(panel1,'Style','text','String','beta2',...
  'FontSize',fsize,'Units','normalized','Position',[.04 .45 .15 .05],...
  'HorizontalAlignment','left');
uicontrol(panel1,'Style','text','String','epsilon',...
  'FontSize',fsize,'Units','normalized','Position',[.04 .3 .15 .05],...
  'HorizontalAlignment','left');
uicontrol(panel1,'Style','text','String','forcing',...
  'FontSize',fsize,'Units','normalized','Position',[.04 .1 .15 .05],...
  'HorizontalAlignment','left');

h.slider(1) = uicontrol(panel1,'Style','slider','Value',param(1,slcLayer),...
  'Min',paramLim{slcLayer}(1,1),'Max',paramLim{slcLayer}(1,2),...
  'Units','normalized','Position',[.2 .75 .6 .05]);
h.slider(2) = uicontrol(panel1,'Style','slider','Value',param(2,slcLayer),...
  'Min',paramLim{slcLayer}(2,1),'Max',paramLim{slcLayer}(2,2),...
  'Units','normalized','Position',[.2 .6 .6 .05]);
h.slider(3) = uicontrol(panel1,'Style','slider','Value',param(3,slcLayer),...
  'Min',paramLim{slcLayer}(3,1),'Max',paramLim{slcLayer}(3,2),...
  'Units','normalized','Position',[.2 .45 .6 .05]);
h.slider(4) = uicontrol(panel1,'Style','slider','Value',param(4,slcLayer),...
  'Min',paramLim{slcLayer}(4,1),'Max',paramLim{slcLayer}(4,2),...
  'Units','normalized','Position',[.2 .3 .6 .05]);
h.slider(5) = uicontrol(panel1,'Style','slider','Value',force,...
  'Min',forceLim(1),'Max',forceLim(2),'Units','normalized',...
  'Position',[.2 .1 .6 .05]);

h.pval(1) = uicontrol(panel1,'Style','edit','String',param(1,slcLayer),...
  'FontSize',fsize,'Units','normalized','Position',[.82 .74 .15 .07]);
h.pval(2) = uicontrol(panel1,'Style','edit','String',param(2,slcLayer),...
  'FontSize',fsize,'Units','normalized','Position',[.82 .59 .15 .07]);
h.pval(3) = uicontrol(panel1,'Style','edit','String',param(3,slcLayer),...
  'FontSize',fsize,'Units','normalized','Position',[.82 .44 .15 .07]);
h.pval(4) = uicontrol(panel1,'Style','edit','String',param(4,slcLayer),...
  'FontSize',fsize,'Units','normalized','Position',[.82 .29 .15 .07]);
h.pval(5) = uicontrol(panel1,'Style','edit','String',force,...
  'FontSize',fsize,'Units','normalized','Position',[.82 .09 .15 .07]);

h.min(1) = uicontrol(panel1,'Style','edit',...
  'String',paramLim{slcLayer}(1,1),'FontSize',fsize-1,...
  'Units','normalized','Position',[.2 .7 .1 .05]);
h.min(2) = uicontrol(panel1,'Style','edit',...
  'String',paramLim{slcLayer}(2,1),'FontSize',fsize-1,...
  'Units','normalized','Position',[.2 .55 .1 .05]);
h.min(3) = uicontrol(panel1,'Style','edit',...
  'String',paramLim{slcLayer}(3,1),'FontSize',fsize-1,...
  'Units','normalized','Position',[.2 .4 .1 .05]);
h.min(4) = uicontrol(panel1,'Style','edit',...
  'String',paramLim{slcLayer}(4,1),'FontSize',fsize-1,...
  'Units','normalized','Position',[.2 .25 .1 .05]);
h.min(5) = uicontrol(panel1,'Style','edit','String',forceLim(1),...
  'FontSize',fsize,'Units','normalized','Position',[.2 .05 .1 .05]);

h.max(1) = uicontrol(panel1,'Style','edit',...
  'String',paramLim{slcLayer}(1,2),'FontSize',fsize-1,...
  'Units','normalized','Position',[.7 .7 .1 .05]);
h.max(2) = uicontrol(panel1,'Style','edit',...
  'String',paramLim{slcLayer}(2,2),'FontSize',fsize-1,...
  'Units','normalized','Position',[.7 .55 .1 .05]);
h.max(3) = uicontrol(panel1,'Style','edit',...
  'String',paramLim{slcLayer}(3,2),'FontSize',fsize-1,...
  'Units','normalized','Position',[.7 .4 .1 .05]);
h.max(4) = uicontrol(panel1,'Style','edit',...
  'String',paramLim{slcLayer}(4,2),'FontSize',fsize-1,...
  'Units','normalized','Position',[.7 .25 .1 .05]);
h.max(5) = uicontrol(panel1,'Style','edit','String',forceLim(2),...
  'FontSize',fsize-1,'Units','normalized','Position',[.7 .05 .1 .05]);

% Create a table for parameter display
pnames = {'alpha','beta1','beta2','epsilon'};
lnames = {};
for nl = 1:numLayer
  lnames{nl} = ['Layer ' num2str(nl)];
end
h.param = uitable(hf,'Data',param,'ColumnName',lnames,...
  'RowName',pnames,'FontSize',fsize,'Units','normalized',...
  'Position',[.57 .25 .4 .15],'ColumnEditable',false);

% Create a button group for frequency spacing type
h.bgFreq = uibuttongroup(hf,'Position',[.57 .15 .4 .1]);
uicontrol(h.bgFreq,'Style','text','String','Type of Frequency Spacing',...
  'FontSize',fsize,'Units','normalized','Position',[.02 .75 .9 .2],...
  'HorizontalAlignment','left');
h.log = uicontrol(h.bgFreq,'Style','radiobutton','FontSize',fsize,...
  'String','Log','Units','normalized','Position',[.05 .4 .2 .2]);
h.lin = uicontrol(h.bgFreq,'Style','radiobutton','FontSize',fsize,...
  'String','Lin','Units','normalized','Position',[.05 .1 .2 .2]);
set(h.bgFreq,'SelectedObject',h.log)
uicontrol(h.bgFreq,'Style','text',...
  'String','(Number of octaves each side:','FontSize',fsize,...
  'Units','normalized','Position',[.25 .4 .5 .2],...
  'HorizontalAlignment','right');
uicontrol(h.bgFreq,'Style','text','String',')','FontSize',fsize,...
  'Units','normalized','Position',[.9 .4 .1 .2],...
  'HorizontalAlignment','left');
uicontrol(h.bgFreq,'Style','text',...
  'String','(Max freq difference from input:','FontSize',fsize,...
  'Units','normalized','Position',[.25 .1 .5 .2],...
  'HorizontalAlignment','right');
uicontrol(h.bgFreq,'Style','text','String','Hz)','FontSize',fsize,...
  'Units','normalized','Position',[.9 .1 .1 .2],...
  'HorizontalAlignment','left');
h.numOct = uicontrol(h.bgFreq,'Style','edit','String',numOct,...
  'FontSize',fsize,'Units','normalized','Position',[.75 .4 .15 .3]);
h.maxdf = uicontrol(h.bgFreq,'Style','edit','String',maxdf,'Enable','off',...
  'FontSize',fsize,'Units','normalized','Position',[.75 .1 .15 .3]);

% Create a panel for misc controls
panel3 = uipanel(hf,'Position',[.57 .1 .4 .05]);
h.legend = uicontrol(panel3,'Style','checkbox','String','Show Legend',...
  'Value',true,'FontSize',fsize,'Units','normalized',...
  'Position',[.02 .1 .3 .8]);
h.spontAmp = uicontrol(panel3,'Style','checkbox',...
  'String','Use Spontaneous Amp for Unstable Fixed Pt',...
  'Value',false,'FontSize',fsize,'Units','normalized',...
  'Position',[.4 .1 .6 .8]);

% Set callback functions
set(h.slider(:),'Callback',{@slider_Callback,h})
set(h.pval(:),'Callback',{@pval_Callback,h})
set(h.min(:),'Callback',{@min_Callback,h})
set(h.max(:),'Callback',{@max_Callback,h})
set(h.numLayer,'Callback',{@numLayer_Callback,h})
set(h.slcLayer,'Callback',{@slcLayer_Callback,h})
set(h.bgFreq,'SelectionChangeFcn',{@bgFreq_Callback,h})
set(h.numOct,'Callback',{@drawAgain_Callback,h})
set(h.maxdf,'Callback',{@drawAgain_Callback,h})
set(h.legend,'Callback',{@legend_Callback,h})
set(h.spontAmp,'Callback',{@drawAgain_Callback,h})

drawAfferentChain(h)

% =========================================================================
function slider_Callback(source,eventdata,h)
ind = find(h.slider == source);
newValue = get(source,'Value');
if ind == 3 && newValue > 0
  oldValue = str2double(get(h.pval(ind),'String'));
  set(source,'Value',oldValue)
  warndlg('Beta2 cannot be negative')
  return
end
set(h.pval(ind),'String',num2str(newValue))

if ind < 5
  param = get(h.param,'Data');
  param(ind,get(h.slcLayer,'Value')) = newValue;
  set(h.param,'Data',param)
end
drawAfferentChain(h)

% =========================================================================
function pval_Callback(source,eventdata,h)
global paramLim
ind = find(h.pval == source);
slcLayer = get(h.slcLayer,'Value');
sliderMin = get(h.slider(ind),'Min');
sliderMax = get(h.slider(ind),'Max');
newValue = str2double(get(source,'String'));
if ind == 3 && newValue > 0
  oldValue = get(h.slider(ind),'Value');
  set(source,'String',oldValue)
  warndlg('Beta2 cannot be negative')
  return
end
  
if newValue < sliderMin
  set(h.slider(ind),'Min',newValue)
  set(h.min(ind),'String',num2str(newValue))
  paramLim{slcLayer}(ind,1) = newValue;
end
if newValue > sliderMax
  set(h.slider(ind),'Max',newValue)
  set(h.max(ind),'String',num2str(newValue))
  paramLim{slcLayer}(ind,2) = newValue;
end
set(h.slider(ind),'Value',newValue)

if ind < 5
  param = get(h.param,'Data');
  param(ind,slcLayer) = newValue;
  set(h.param,'Data',param)
end
drawAfferentChain(h)

% =========================================================================
function min_Callback(source,eventdata,h)
global paramLim
ind = find(h.min == source);
slcLayer = get(h.slcLayer,'Value');
newMin = str2double(get(source,'String'));
sliderValue = get(h.slider(ind),'Value');
currentMax = get(h.slider(ind),'Max');

if newMin >= currentMax
  newMax = newMin + .1;
  newValue = newMin;
  set(h.slider(ind),'Max',newMax)
  set(h.slider(ind),'Value',newValue)
  set(h.pval(ind),'String',num2str(newValue))
  set(h.max(ind),'String',num2str(newMax))
  paramLim{slcLayer}(ind,2) = newMax;
  if ind < 5
    param = get(h.param,'Data');
    param(ind,slcLayer) = newValue;
    set(h.param,'Data',param)
  end
  drawAfferentChain(h)
  
elseif newMin > sliderValue
  newValue = newMin;
  set(h.slider(ind),'Value',newValue)
  set(h.pval(ind),'String',num2str(newValue))
  if ind < 5
    param = get(h.param,'Data');
    param(ind,slcLayer) = newValue;
    set(h.param,'Data',param)
  end
  drawAfferentChain(h)
end

set(h.slider(ind),'Min',newMin)
paramLim{slcLayer}(ind,1) = newMin;

% =========================================================================
function max_Callback(source,eventdata,h)
global paramLim
ind = find(h.max == source);
slcLayer = get(h.slcLayer,'Value');
newMax = str2double(get(source,'String'));
sliderValue = get(h.slider(ind),'Value');
currentMin = get(h.slider(ind),'Min');

if newMax <= currentMin
  newMin = newMax - .1;
  newValue = newMax;
  set(h.slider(ind),'Min',newMin)
  set(h.slider(ind),'Value',newValue)
  set(h.pval(ind),'String',num2str(newValue))
  set(h.min(ind),'String',num2str(newMin))
  paramLim{slcLayer}(ind,1) = newMin;
  if ind < 5
    param = get(h.param,'Data');
    param(ind,get(h.slcLayer,'Value')) = newValue;
    set(h.param,'Data',param)
  end
  drawAfferentChain(h)
  
elseif newMax < sliderValue
  newValue = newMax;
  set(h.slider(ind),'Value',newValue)
  set(h.pval(ind),'String',num2str(newValue))
  if ind < 5
    param = get(h.param,'Data');
    param(ind,get(h.slcLayer,'Value')) = newValue;
    set(h.param,'Data',param)
  end
  drawAfferentChain(h)
end

set(h.slider(ind),'Max',newMax)
paramLim{slcLayer}(ind,2) = newMax;

% =========================================================================
function numLayer_Callback(source,eventdata,h)
global paramLim
newNum = str2double(get(h.numLayer,'String'));
param = get(h.param,'Data');
oldNum = size(param,2);
if mod(newNum,1) || newNum < 1
  set(h.numLayer,'String',oldNum)
  return
elseif newNum > oldNum
  param = [param param(:,end)*ones(1,newNum-oldNum)];
  set(h.param,'Data',param)
  for nl = oldNum+1:newNum
    paramLim{nl} = paramLim{nl-1};
  end
elseif newNum < oldNum
  param = param(:,1:newNum);
  set(h.param,'Data',param)
  paramLim = paramLim(1:newNum);
end

layerString = [];
for nl = 1:newNum
  layerString = [layerString 'Layer ' num2str(nl) '|'];
  lnames{nl} = ['Layer ' num2str(nl)];
end
layerString = layerString(1:end-1);
set(h.slcLayer,'String',layerString)
if get(h.slcLayer,'Value') > newNum
  set(h.slcLayer,'Value',newNum)
  slcLayer_Callback([],[],h)
end
set(h.param,'ColumnName',lnames)
drawAfferentChain(h)

% =========================================================================
function slcLayer_Callback(source,eventdata,h)
global paramLim
slcLayer = get(h.slcLayer,'Value');
param = get(h.param,'Data');
for ind = 1:4
  set(h.slider(ind),'Value',param(ind,slcLayer),...
    'Min',paramLim{slcLayer}(ind,1),'Max',paramLim{slcLayer}(ind,2))
  set(h.min(ind),'String',paramLim{slcLayer}(ind,1))
  set(h.max(ind),'String',paramLim{slcLayer}(ind,2))
  set(h.pval(ind),'String',param(ind,slcLayer))
end

% =========================================================================
function bgFreq_Callback(source,eventdata,h)
switch get(eventdata.NewValue,'String')
  case 'Log'
    set(h.numOct,'Enable','on')
    set(h.maxdf,'Enable','off')
  case 'Lin'
    set(h.numOct,'Enable','off')
    set(h.maxdf,'Enable','on')
end
drawAfferentChain(h)

% =========================================================================
function drawAgain_Callback(source,eventdata,h)
drawAfferentChain(h)

% =========================================================================
function legend_Callback(source,eventdata,h)
if get(source,'Value')
  drawAfferentChain(h)
else
  legend(h.axr,'off')
end

% =========================================================================
function drawAfferentChain(h)

% To deal w/ display issues
platform = computer;
vMATLAB = version('-release');
if strcmp(platform(1:3),'MAC') && strcmp(vMATLAB(1:4),'2011')
  msize = 7;
else
  msize = 10;
end

param = get(h.param,'Data');
a = param(1,:);
b1 = param(2,:);
b2 = param(3,:);
e = param(4,:);
F = str2double(get(h.pval(5),'String'));
freqSpac = get(get(h.bgFreq,'SelectedObject'),'String');

numOsc = 501; % number of oscillators in each layer
switch(freqSpac)
  case 'Log'
    numOct = str2double(get(h.numOct,'String'));
    fz = 2.^(linspace(-1,1,numOsc)*numOct);
    W = (fz-1)./fz*2*pi; % freq scaled omegas
  case 'Lin'
    maxdf = str2double(get(h.maxdf,'String'));
    fz = linspace(-1,1,numOsc)*maxdf;
    W = fz*2*pi;
end

nlayer = length(a);
maxnr = ones(size(a))*2;
if nlayer < 6; lcolor = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1];
else lcolor = jet(nlayer); end
rstar = cell(nlayer,1);
psistar = cell(nlayer,1);
hl = [];
rmax = zeros(nlayer,1);
switch(freqSpac)
  case 'Log'
    x = log2(fz);
  case 'Lin'
    x = fz;
end

for nl = 1:nlayer
  rstar{nl} = -ones(length(W),maxnr(nl));
  psistar{nl} = zeros(size(rstar{nl}));
  if nl == 1
    rin = ones(size(W))*F;
    psiin = zeros(size(W));
  else
    rin = rstar{nl-1};
    psiin = psistar{nl-1};
  end
  for nw = 1:length(W)
    if ~isnan(rin(nw))
      [rstarnw,psistarnw] = ...
        rstarfull(a(nl),b1(nl),b2(nl),e(nl),rin(nw),W(nw));
      nrstar = length(rstarnw); % # of r*s for W(nw)
      rstar{nl}(nw,1:nrstar) = rstarnw';
      psistar{nl}(nw,1:nrstar) = ...
        mod(real(psistarnw')+psiin(nw)+pi,2*pi)-pi;
    end
  end
  
  if get(h.spontAmp,'Value')
    ind1 = find(rstar{nl}(:,1) < 0); % indices w/o any stable fixed pt
    for ni = 1:length(ind1)
      r0 = spontAmp(a(nl),b1(nl),b2(nl),e(nl));
      rstar{nl}(ind1(ni),1:length(r0)) = r0;
      psistar{nl}(ind1(ni),1:length(r0)) = 0;
      W(ind1(ni)) = 0;
    end
  end
  ind2 = find(rstar{nl} < 0); % indices for unstable fixed pts
  rstar{nl}(ind2) = NaN;
  psistar{nl}(ind2) = NaN;
  rmax(nl) = max(rstar{nl}(:));
  
  for nr = 1:maxnr(nl)
    subplot(h.axr)
    hl(nl) = plot(x,rstar{nl}(:,nr),'.:','MarkerSize',msize,...
      'Color',lcolor(nl,:));
    hold on
    subplot(h.axpsi)
    plot(x,psistar{nl}(:,nr),'.:','MarkerSize',msize,...
      'Color',lcolor(nl,:));
    hold on
  end
end

subplot(h.axr)
hold off
grid on
if get(h.legend,'Value')
  lnames = get(h.param,'ColumnName');
  legend(hl,lnames,'Location','NorthWest')
end
set(gca,'XLim',[min(x) max(x)],'YLim',[0 max(rmax)*1.05])
ylabel('r*','FontSize',12)
title('Afferent Chains','FontSize',13)

subplot(h.axpsi)
hold off
grid on
switch(freqSpac)
  case 'Log'
    xlabel('log_2 ( f / f_0 )','FontSize',12)
  case 'Lin'
    xlabel('f - f_0 (Hz)','FontSize',12)
end
set(gca,'XLim',[min(x) max(x)],'Ylim',[-pi pi],...
  'YTick',[-pi,-pi/2,0,pi/2,pi],...
  'YTickLabel',{'-pi';'-pi/2';'0';'pi/2';'pi'})
ylabel('\psi*','FontSize',12)
drawnow
