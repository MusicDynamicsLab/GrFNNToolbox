%% CREATE FIGURE WINDOW
function coupledkmAnalysisGUI
hf = figure;
set(hf,'Visible','off','Toolbar','figure','Color',[.8 .8 .8],...
  'Position',[100 100 1000 600]);
initialRun
movegui(hf,'center')
set(hf,'Visible','on')

%% INITIAL RUN
function initialRun
hf = gcf;
handles.hf = hf;

%% Create title and equation
axes('Parent',hf,'Units','normalized','Position',[.02 .79 .45 .19]);
lcolor = jet(5);
hc = zeros(5,1);
for nc = 1:5
  hc(nc) = plot(NaN,NaN,'-','Color',lcolor(6-nc,:),'LineWidth',2);
  hold on
end
hl = legend(hc,'Stable node','Stable spiral','Unstable node',...
  'Unstable spiral','Saddle point','Location','NorthEast');
lpos = get(hl,'Position');
lpos(1:2) = [.345 .56];
set(hl,'Position',lpos)
hold off
axis off
text('String','Two Oscillators with {\itk}:{\itm} Coupling ',...
  'Position',[.0 1],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',13,'FontWeight','bold')
text('Interpreter','latex',...
	'String',{['$$\frac{dz_1}{dt} = z_1\left(\alpha + ' ...
  '\textrm{i}\omega_1 + \beta_1 |z_1|^2 + ' ...
  '\frac{\epsilon\beta_2|z_1|^4}' ...
  '{1-\epsilon |z_1|^2}\right) + \epsilon^\frac{k+m-2}{2}' ...
  'cz_2^k \bar{z}_1^{m-1}$$'],...
  ['$$\frac{dz_2}{dt} = z_2\left(\alpha + \textrm{i}\omega_2 + ' ...
  '\beta_1 |z_2|^2 + ' ...
  '\frac{\epsilon\beta_2|z_2|^4}' ...
  '{1-\epsilon |z_2|^2}\right) + \epsilon^\frac{k+m-2}{2}' ...
  'cz_1^m \bar{z}_2^{k-1}$$']},...
	'Position',[.02 .7],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',12)

%% Create Legend button
handles.legendButton = uicontrol('Parent',hf,'Style','pushbutton',...
  'String','Legend','FontSize',10,'Units','normalized',...
  'Position',[.42 .93 .07 .05]);

%% Make text boxes for parameter input and Draw button
hp1 = uipanel('Parent',hf,'Title','Parameters','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Position',[.02 .56 .32 .2]);
nc = 5; % # of columns
rb = 15; % ratio of box length relative to (half) space
wb = 1/nc*rb/(rb+2); % normalized width of boxes
ws = 1/(nc*(rb+2)); % normalized width of (half) space
wu = wb + 2*ws;
handles.alpha = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[ws .55 wb .25]);
handles.beta1 = uicontrol('Parent',hp1,'Style','edit','String','-10',...
  'FontSize',11,'Units','normalized','Position',[ws+wu .55 wb .25]);
handles.beta2 = uicontrol('Parent',hp1,'Style','edit','String','-1',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*2 .55 wb .25]);
handles.epsilon = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*3 .55 wb .25]);
handles.df = uicontrol('Parent',hp1,'Style','edit','String','.1',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*4 .55 wb .25]);
handles.c = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[ws .05 wb .25]);
handles.k = uicontrol('Parent',hp1,'Style','edit','String','2',...
  'FontSize',11,'Units','normalized','Position',[ws+wu .05 wb .25]);
handles.m = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*2 .05 wb .25]);

handles.drawButton = uicontrol('Parent',hp1,'Style','pushbutton',...
  'String','Draw','FontSize',11,'Units','normalized',...
  'Position',[.65 .05 .325 .3]);

%% Labels for text boxes
axes('Parent',hp1,'Units','normalized','Position',[0 0 1 1])
axis off
text(wu/2,.85,'$$\alpha$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*3/2,.85,'$$\beta_1$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*5/2,.85,'$$\beta_2$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*7/2,.85,'$$\epsilon$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*9/2,.85,'$$\Omega/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu/2,.35,'$$c$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*3/2,.35,'$$k$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*5/2,.35,'$$m$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)

%% Create axes
handles.ax1a = axes('Parent',hf,'Units','normalized',...
  'Position',[.12 .07 .15 .45],'XTick',[],'YTick',[],'box','on');
xlabel('{\itr}*','FontSize',11)
ylabel('{\itc}','FontSize',11)
text(-.5,.5,'Constant Frequency Difference (\Omega)','FontSize',11,...
  'FontWeight','bold','HorizontalAlignment','center',...
  'VerticalAlignment','middle','Units','normalized','Rotation',90)

handles.ax1b = axes('Parent',hf,'Units','normalized',...
  'Position',[.34 .07 .15 .45],'XTick',[],'YTick',[],'box','on');
xlabel('{\it\psi}*','FontSize',11)

handles.ax2a = axes('Parent',hf,'Units','normalized',...
  'Position',[.57 .81 .41 .13],'XTick',[],'YTick',[],'box','on');
ylabel('{\itr}*','FontSize',11)
title('Constant Coupling Amplitude ({\itc})','FontSize',11,...
  'FontWeight','bold')

handles.ax2b = axes('Parent',hf,'Units','normalized',...
  'Position',[.57 .63 .41 .13],'XTick',[],'YTick',[],'box','on');
%xlabel('\Omega/2\pi','FontSize',11)
ylabel('{\it\psi}*','FontSize',11)

handles.ax3 = axes('Parent',hf,'Units','normalized',...
  'Position',[.57 .07 .41 .45],'XTick',[],'YTick',[],'box','on');
xlabel('\Omega/2\pi','FontSize',11)
%ylabel('{\itc}','FontSize',11)
title('Stability Region','FontSize',11,'FontWeight','bold')

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
W = str2double(get(handles.df,'String'))*2*pi;
c = str2double(get(handles.c,'String'));
k = str2double(get(handles.k,'String'));
m = str2double(get(handles.m,'String'));
lcolor = jet(5);

%% Prepare vectors of coupling strength and Omega
if e < 0
  axes(handles.ax1a)
  cla
  set(handles.ax1a,'XTick',[],'YTick',[],'box','on')
  text(.5,.5,'{\it\epsilon} cannot be negative','FontSize',15,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Units','normalized')
  return
end

if c <= 0
  axes(handles.ax1a)
  cla
  set(handles.ax1a,'XTick',[],'YTick',[],'box','on')
  text(.5,.5,'{\itc} must be greater than 0','FontSize',15,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Units','normalized')
  return
end

cc = (.01:.005:1)*ceil(2.5*c)/2;
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

r1star = cell(size(cc));
r2star = cell(size(cc));
psistar = cell(size(cc));
stabtype = cell(size(cc));
rmax = 0;
for n = 1:length(cc)
  [r1star{n},r2star{n},psistar{n},~,stabtype{n}] = ...
    rStarCoupledkm(a,b1,b2,e,cc(n),W,k,m,1);
  stabtype{n} = ceil(stabtype{n});
  rmax = max([rmax max(r1star{n})]);
end

axes(handles.ax1a)
for n = 1:length(cc)
  for nr = 1:length(r1star{n})
    hr = plot(r1star{n}(nr),cc(n),'.',...
      'Color',lcolor(stabtype{n}(nr)+1,:));
    if r1star{n}(nr) == r2star{n}(nr); set(hr,'MarkerSize',10);
      else set(hr,'MarkerSize',4); end
    hold on
  end
end
if rmax == 0
  rmax = 1;
end
plot([0 rmax*1.1],[1 1]*c,'r--','LineWidth',2)
set(gca,'XLim',[0 rmax*1.1],'YLim',[min(cc) max(cc)])
xlabel('{\itr}*','FontSize',11)
ylabel('{\itc}','FontSize',11)
text(-.5,.5,'Constant Frequency Difference (\Omega)','FontSize',11,...
  'FontWeight','bold','HorizontalAlignment','center',...
  'VerticalAlignment','middle','Units','normalized','Rotation',90)
grid on
hold off
set(gca,'XDir','Reverse')

axes(handles.ax1b)
for n = 1:length(cc)
  for nr = 1:length(psistar{n})
    hpsi = plot(real(psistar{n}(nr)),cc(n),'.',...
      'Color',lcolor(stabtype{n}(nr)+1,:));
    % thick if symmetric, thin if asymmetric
    if r1star{n}(nr) == r2star{n}(nr); set(hpsi,'MarkerSize',10);
      else set(hpsi,'MarkerSize',4); end
    hold on
  end
end
plot([-1 1]*pi,[1 1]*c,'r--','LineWidth',2)
set(gca,'XLim',[-1 1]*pi,'XTick',[-pi,-pi/2,0,pi/2,pi],...
  'XTickLabel',{'-pi';'-pi/2';'0';'pi/2';'pi'},'YLim',[min(cc) max(cc)])
xlabel('{\it\psi}*','FontSize',11)
grid on
hold off
set(gca,'XDir','Reverse')

%% Constant Coupling plot
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

r1star = cell(size(WW));
r2star = cell(size(WW));
psistar = cell(size(WW));
stabtype = cell(size(WW));
rmax = 0;

for n = 1:length(WW)
  [r1star{n},r2star{n},psistar{n},~,stabtype{n}] = ...
    rStarCoupledkm(a,b1,b2,e,c,WW(n),k,m,1);
  stabtype{n} = ceil(stabtype{n});
  rmax = max([rmax max(r1star{n})]);
end

axes(handles.ax2a)
for n = 1:length(WW)
  for nr = 1:length(r1star{n})
    hr = plot(WW(n)/(2*pi),r1star{n}(nr),'.',...
      'Color',lcolor(stabtype{n}(nr)+1,:));
    if r1star{n}(nr) == r2star{n}(nr); set(hr,'MarkerSize',10);
      else set(hr,'MarkerSize',4); end
    hold on
  end
end
if rmax == 0
  rmax = 1;
end
plot([1 1]*W/(2*pi),[0 rmax*1.1],'b--','LineWidth',2)
set(gca,'XLim',[min(WW) max(WW)]./(2*pi),'YLim',[0 rmax*1.1])
ylabel('{\itr}*','FontSize',11)
title('Constant Coupling Amplitude ({\itc})','FontSize',11,...
  'FontWeight','bold')
grid on
hold off

axes(handles.ax2b)
for n = 1:length(WW)
  for nr = 1:length(psistar{n})
    hpsi = plot(WW(n)/(2*pi),real(psistar{n}(nr)),'.',...
      'Color',lcolor(stabtype{n}(nr)+1,:));
    if r1star{n}(nr) == r2star{n}(nr); set(hpsi,'MarkerSize',10);
      else set(hpsi,'MarkerSize',4); end
    hold on
  end
end
plot([1 1]*W/(2*pi),[-1 1]*pi,'b--','LineWidth',2)
set(gca,'XLim',[min(WW) max(WW)]./(2*pi),'YLim',[-1 1]*pi,...
  'YTick',[-pi,-pi/2,0,pi/2,pi],...
  'YTickLabel',{'-pi';'-pi/2';'0';'pi/2';'pi'})
%xlabel('\Omega/2\pi','FontSize',11)
ylabel('{\it\psi}*','FontSize',11)
grid on
hold off

%% Draw Stability Region
axes(handles.ax3)
cla
set(gca,'XTick',[],'YTick',[],'box','on')
text(.5,.5,'Calculating...','FontSize',15,...
  'HorizontalAlignment','center','VerticalAlignment','middle',...
  'Units','normalized')
drawnow
STABTYPE = zeros(length(cc),length(WW));
for nc = 1:length(cc)
  for nw = 1:length(WW)
    [r1,r2,~,~,type] = rStarCoupledkm(a,b1,b2,e,cc(nc),WW(nw),k,m,1);
    typeSym = type(r1 == r2);
    ind = find(r1 ~= r2);
    rAsym = [r1(ind) r2(ind)];
    typeAsym = type(ind);
    for n = 1:length(ind)
      if typeAsym(n) > 0
        for nn = 1:length(ind)
          if sum(abs(rAsym(nn,:)-fliplr(rAsym(n,:)))) < 10^(-2)
            typeAsym(nn) = -1;
          end
        end
      end
    end
    typeAsym = typeAsym(typeAsym >= 0);
    type = [typeSym; typeAsym];
    stab = type >= 3;
    if isempty(type)
      STABTYPE(nc,nw) = -1;
    elseif sum(stab) < 2
      STABTYPE(nc,nw) = max(type); % show stable type if any
    elseif sum(stab) == 2
      STABTYPE(nc,nw) = 5;
    else
      error('Three or more stable fixed points?')
    end
  end
end
imagesc(WW/(2*pi),cc,STABTYPE)
colormap([1 1 1; jet(5); .5 0 .5])
      % purple for 2 stable pts, white for no fixed pt
caxis([-1 5])
hold on
plot(get(gca,'XLim'),[1 1]*c,'r--','LineWidth',2)
plot([1 1]*W/(2*pi),get(gca,'YLim'),'b--','LineWidth',2)
hold off
set(gca,'YDir','normal')
xlabel('\Omega/2\pi','FontSize',11)
%ylabel('{\itc}','FontSize',11)
title('Stability Region','FontSize',11,'FontWeight','bold')
grid on

%% SHOW LEGEND
function showLegend(~,~)
hf = figure;
fpos = get(hf,'Position');

set(hf,'Position',[fpos(1) fpos(2) 450 160])
axes('Units','normalized','Position',[0 0 1 1]);
axis off
text('Interpreter','latex','String',...
  {['$$z_i = r_ie^{\textrm{i}\phi_i}$$, ' ...
  '$$\psi = m\phi_1 - k\phi_2$$, and'],...
  '$$\Omega = m\omega_1 - k\omega_2$$'},...
	'Position',[.05 .9],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',15)
text('Interpreter','latex','String',...
  {'*Only symmetric solutions are shown ($$r^* \equiv r_1^* = r_2^*$$)',...
	'**Purple in the Stability Region indicates',...
  'the existence of two stable fixed points'},...
  'Position',[.05 .5],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',13)
