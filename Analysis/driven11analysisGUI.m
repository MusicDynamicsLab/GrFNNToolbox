%% GREATE FIGURE WINDOW
function driven11AnalysisGUI
fh = figure;
set(fh,'Visible','off','Toolbar','figure','Color',[.8 .8 .8],...
  'Position',[100 100 1000 600]);
initialRun
movegui(fh,'center')
set(fh,'Visible','on')

%% INITIAL RUN
function initialRun
%% Create title and equation
axes('Units','normalized','Position',[.02 .79 .45 .19]);
lcolor = jet(5);
hc = zeros(5,1);
for nc = 1:5
  hc(nc) = plot(NaN,NaN,'-','Color',lcolor(6-nc,:),'LineWidth',2);
  hold on
end
hl = legend(hc,'Stable node','Stable spiral','Unstable node',...
  'Unstable spiral','Saddle point',...
  'Location','NorthEast');
lpos = get(hl,'Position');
lpos(1:2) = [.345 .605];
set(hl,'Position',lpos)
hold off
axis off
text('String','Single Oscillator Driven by a Sinusoid',...
  'Position',[0 1],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',13,'FontWeight','bold')
text('Interpreter','latex',...
	'String',['$$\frac{dz}{dt} = z\left(\alpha + \textrm{i}\omega + ' ...
  '\beta_1 |z|^2 + \frac{\epsilon\beta_2|z|^4}' ...
  '{1-\epsilon |z|^2}\right) + Fe^{\textrm{i}(\omega_0t+\theta_0)}$$'],...
	'Position',[.02 .25],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','bottom',...
	'FontSize',14)

%% Create Legend button
handles.legendButton = uicontrol('Parent',gcf,'Style','pushbutton',...
  'String','Legend','FontSize',10,'Units','normalized',...
  'Position',[.42 .93 .07 .05]);

%% Make text boxes for parameter input and Draw button
hp1 = uipanel('Title','Parameters','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Position',[.02 .6 .32 .2]);
handles.alpha = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[.025 .55 .2 .25]);
handles.beta1 = uicontrol('Parent',hp1,'Style','edit','String','-10',...
  'FontSize',11,'Units','normalized','Position',[.275 .55 .2 .25]);
handles.beta2 = uicontrol('Parent',hp1,'Style','edit','String','0',...
  'FontSize',11,'Units','normalized','Position',[.525 .55 .2 .25]);
handles.epsilon = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[.775 .55 .2 .25]);
handles.F = uicontrol('Parent',hp1,'Style','edit','String','.3',...
  'FontSize',11,'Units','normalized','Position',[.025 .05 .2 .25]);
handles.df = uicontrol('Parent',hp1,'Style','edit','String','.1',...
  'FontSize',11,'Units','normalized','Position',[.275 .05 .2 .25]);
handles.drawButton = uicontrol('Parent',hp1,'Style','pushbutton',...
  'String','Draw','FontSize',11,'Units','normalized',...
  'Position',[.65 .05 .325 .3]);

%% Labels for text boxes
axes('Parent',hp1,'Units','normalized','Position',[0 0 1 1])
axis off
text(.125,.85,'$$\alpha$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.375,.85,'$$\beta_1$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.625,.85,'$$\beta_2$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.875,.85,'$$\epsilon$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.125,.35,'$$F$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.375,.35,'$$\Omega/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)

%% Create axes
handles.ax1a = axes('Position',[.12 .07 .15 .45],'Units','normalized',...
  'XTick',[],'YTick',[],'box','on');
xlabel('{\itr}*','FontSize',11)
ylabel('{\itF}','FontSize',11)
text(-.5,.5,'Constant Frequency Difference (\Omega)','FontSize',11,...
  'FontWeight','bold','HorizontalAlignment','center',...
  'VerticalAlignment','middle','Units','normalized','Rotation',90)

handles.ax1b = axes('Position',[.34 .07 .15 .45],'Units','normalized',...
  'XTick',[],'YTick',[],'box','on');
xlabel('{\it\psi}*','FontSize',11)

handles.ax2a = axes('Position',[.57 .81 .41 .13],'Units','normalized',...
  'XTick',[],'YTick',[],'box','on');
ylabel('{\itr}*','FontSize',11)
title('Constant Forcing Amplitude ({\itF})','FontSize',11,...
  'FontWeight','bold')

handles.ax2b = axes('Position',[.57 .63 .41 .13],'Units','normalized',...
  'XTick',[],'YTick',[],'box','on');
%xlabel('\Omega/2\pi','FontSize',11)
ylabel('{\it\psi}*','FontSize',11)

handles.ax3 = axes('Position',[.57 .07 .41 .45],'Units','normalized',...
  'XTick',[],'YTick',[],'box','on');
xlabel('\Omega/2\pi','FontSize',11)
%ylabel('{\itF}','FontSize',11)
title('Resonance Region','FontSize',11,'FontWeight','bold')

%% Define callback functions
set(handles.drawButton,'Callback',{@drawAnalysis,handles})
set(handles.legendButton,'Callback',@showLegend)

%% DRAW ANALYSIS
function drawAnalysis(~,~,handles)
%% Get parameters
a = str2double(get(handles.alpha,'String'));
b1 = str2double(get(handles.beta1,'String'));
b2 = str2double(get(handles.beta2,'String'));
e = str2double(get(handles.epsilon,'String'));
F = str2double(get(handles.F,'String'));
W = str2double(get(handles.df,'String'))*2*pi;
lcolor = jet(5);

%% Prepare vectors of Forcing and Omega
if e < 0
  axes(handles.ax1a)
  cla
  set(handles.ax1a,'XTick',[],'YTick',[],'box','on')
  text(.5,.5,'{\it\epsilon} cannot be negative','FontSize',15,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Units','normalized')
  return
end
if F <= 0
  axes(handles.ax1a)
  cla
  set(handles.ax1a,'XTick',[],'YTick',[],'box','on')
  text(.5,.5,'{\itF} must be greater than 0','FontSize',15,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Units','normalized')
  return
end

FF = (.005:.005:1)*ceil(2.5*F)/2;
WW = (-1:.01:1)*ceil(2.5*abs(W+eps)/(2*pi))*pi;

%% Constant Omega plot
axes(handles.ax1a)
cla
set(gca,'XTick',[],'YTick',[],'box','on')
text(.5,.5,'Calculating...','FontSize',15,...
  'HorizontalAlignment','center','VerticalAlignment','middle',...
  'Units','normalized')
text(-.5,.5,'Constant Frequency Difference (\Omega)','FontSize',11,...
  'FontWeight','bold','HorizontalAlignment','center',...
  'VerticalAlignment','middle','Units','normalized','Rotation',90)
axes(handles.ax1b)
cla
set(gca,'XTick',[],'YTick',[],'box','on')
text(.5,.5,'Calculating...','FontSize',15,...
  'HorizontalAlignment','center','VerticalAlignment','middle',...
  'Units','normalized')
drawnow

rstar = cell(size(FF));
psistar = cell(size(FF));
stabtype = cell(size(FF));
rmax = 0;
for n = 1:length(FF)
  [rstar{n},psistar{n},~,stabtype{n}] = rStarDriven11(a,b1,b2,e,FF(n),W,1);
  stabtype{n} = ceil(stabtype{n});
  rmax = max([rmax max(rstar{n})]);
end

axes(handles.ax1a)
for n = 1:length(FF)
  for nr = 1:length(rstar{n})
    plot(rstar{n}(nr),FF(n),'.','MarkerSize',10,...
      'Color',lcolor(stabtype{n}(nr)+1,:))
    hold on
  end
end
if rmax == 0
  rmax = 1;
end
set(gca,'XLim',[0 rmax*1.1],'YLim',[min(FF) max(FF)])
plot(get(gca,'XLim'),[1 1]*F,'r--','LineWidth',2)
xlabel('{\itr}*','FontSize',11)
ylabel('{\itF}','FontSize',11)
text(-.5,.5,'Constant Frequency Difference (\Omega)','FontSize',11,...
  'FontWeight','bold','HorizontalAlignment','center',...
  'VerticalAlignment','middle','Units','normalized','Rotation',90)
grid on
hold off
set(gca,'XDir','Reverse')

axes(handles.ax1b)
for n = 1:length(FF)
  for nr = 1:length(psistar{n})
    plot(real(psistar{n}(nr)),FF(n),'.','MarkerSize',10,...
      'Color',lcolor(stabtype{n}(nr)+1,:))
    hold on
  end
end
set(gca,'XLim',[-1 1]*pi,'XTick',[-pi,-pi/2,0,pi/2,pi],...
  'XTickLabel',{'-pi';'-pi/2';'0';'pi/2';'pi'},'YLim',[min(FF) max(FF)])
plot(get(gca,'XLim'),[1 1]*F,'r--','LineWidth',2)
xlabel('{\it\psi}*','FontSize',11)
grid on
hold off
set(gca,'XDir','Reverse')

%% Constant Forcing plot
axes(handles.ax2a)
cla
set(gca,'XTick',[],'YTick',[],'box','on')
text(.5,.5,'Calculating...','FontSize',15,...
  'HorizontalAlignment','center','VerticalAlignment','middle',...
  'Units','normalized')
axes(handles.ax2b)
cla
set(gca,'XTick',[],'YTick',[],'box','on')
text(.5,.5,'Calculating...','FontSize',15,...
  'HorizontalAlignment','center','VerticalAlignment','middle',...
  'Units','normalized')
drawnow

rstar = cell(size(WW));
psistar = cell(size(WW));
stabtype = cell(size(WW));
rmax = 0;
for n = 1:length(WW)
  [rstar{n},psistar{n},~,stabtype{n}] = rStarDriven11(a,b1,b2,e,F,WW(n),1);
  stabtype{n} = ceil(stabtype{n});
  rmax = max([rmax max(rstar{n})]);
end

axes(handles.ax2a)
for n = 1:length(WW)
  for nr = 1:length(rstar{n})
    plot(WW(n)/(2*pi),rstar{n}(nr),'.','MarkerSize',10,...
      'Color',lcolor(stabtype{n}(nr)+1,:))
    hold on
  end
end
if rmax == 0
  rmax = 1;
end
set(gca,'XLim',[min(WW) max(WW)]./(2*pi),'YLim',[0 rmax*1.1])
plot([1 1]*W/(2*pi),get(gca,'YLim'),'b--','LineWidth',2)
ylabel('{\itr}*','FontSize',11)
title('Constant Forcing Amplitude ({\itF})','FontSize',11,...
  'FontWeight','bold')
grid on
hold off

axes(handles.ax2b)
for n = 1:length(WW)
  for nr = 1:length(rstar{n})
    plot(WW(n)/(2*pi),real(psistar{n}(nr)),'.','MarkerSize',10,...
      'Color',lcolor(stabtype{n}(nr)+1,:))
    hold on
  end
end
set(gca,'XLim',[min(WW) max(WW)]./(2*pi),'YLim',[-1 1]*pi,...
  'YTick',[-pi,-pi/2,0,pi/2,pi],...
  'YTickLabel',{'-pi';'-pi/2';'0';'pi/2';'pi'})
plot([1 1]*W/(2*pi),get(gca,'YLim'),'b--','LineWidth',2)
%xlabel('\Omega/2\pi','FontSize',11)
ylabel('{\it\psi}*','FontSize',11)
grid on
hold off

%% Draw Resonance Region
axes(handles.ax3)
cla
set(gca,'XTick',[],'YTick',[],'box','on')
text(.5,.5,'Calculating...','FontSize',15,...
  'HorizontalAlignment','center','VerticalAlignment','middle',...
  'Units','normalized')
drawnow
STABTYPE = zeros(length(FF),length(WW));
for nf = 1:length(FF)
  for nw = 1:length(WW)
    [~,~,stab,type] = rStarDriven11(a,b1,b2,e,FF(nf),WW(nw),1);
    if sum(stab) < 2
      STABTYPE(nf,nw) = max(type); % show stable type if any
    elseif sum(stab) == 2
      STABTYPE(nf,nw) = 5;
    else
      error('Three or more stable fixed points?')
    end
  end
end
imagesc(WW/(2*pi),FF,STABTYPE)
colormap([jet(5);.5 0 .5])
caxis([0 5])
hold on
plot(get(gca,'XLim'),[1 1]*F,'r--','LineWidth',2)
plot([1 1]*W/(2*pi),get(gca,'YLim'),'b--','LineWidth',2)
hold off
set(gca,'YDir','normal')
xlabel('\Omega/2\pi','FontSize',11)
%ylabel('{\itF}','FontSize',11)
title('Resonance Region','FontSize',11,'FontWeight','bold')
grid on

%% SHOW LEGEND
function showLegend(~,~)
hf = figure;
fpos = get(hf,'Position');

set(hf,'Position',[fpos(1) fpos(2) 350 200])
axes('Units','normalized','Position',[0 0 1 1]);
axis off
text('Interpreter','latex','String',{'$$z = re^{\textrm{i}\phi}$$,',...
  ['$$x = Fe^{\textrm{i}\theta}$$ '...
  'where $\theta = \omega_0t+\theta_0$,'],...
  '$$\psi = \phi - \theta$$, and',...
  '$$\Omega = \omega - \omega_0$$'},...
	'Position',[.05 .9],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',15)
text('Interpreter','latex','String',...
  {'*Purple in the Resonance Region indicates',...
  'the existence of two stable fixed points'},...
	'Position',[.05 .3],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',13)
