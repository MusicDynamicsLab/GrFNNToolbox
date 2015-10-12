%% CREATE FIGURE WINDOW
function drivenkmGUI
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
axes('Parent',hf,'Units','normalized','Position',[.025 .88 .5 .1]);
axis off
text('String',['Oscillator {\itk}:{\itm} Mode Locking ' ...
  'to a Sinusoid'],'Position',[0 1],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',13,'FontWeight','bold')
text('Interpreter','latex',...
	'String',['$$\frac{dz}{dt} = z\left(\alpha + \textrm{i}\omega + ' ...
  '(\beta_1+\textrm{i}\delta_1) |z|^2 + ' ...
  '\frac{\epsilon(\beta_2+\textrm{i}\delta_2)|z|^4}' ...
  '{1-\epsilon |z|^2}\right) + \epsilon^\frac{k+m-2}{2}' ...
  'x^k \bar{z}^{m-1}$$'],...
	'Position',[.02 .4],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',12.5)

%% Create Legend button
handles.legendButton = uicontrol('Parent',hf,'Style','pushbutton',...
  'String','Legend','FontSize',10,'Units','normalized',...
  'Position',[.455 .94 .07 .05]);

%% Make text boxes for parameter input
hp1 = uipanel('Parent',hf,'Title','Parameters','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Units','normalized',...
  'Position',[.025 .62 .5 .2]);
nc = 6; % # of columns
rb = 6; % ratio of box length relative to (half) space
wb = 1/nc*rb/(rb+2); % normalized width of boxes
ws = 1/(nc*(rb+2)); % normalized width of (half) space
wu = wb + 2*ws;
handles.alpha = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[ws .55 wb .25]);
handles.beta1 = uicontrol('Parent',hp1,'Style','edit','String','-1',...
  'FontSize',11,'Units','normalized','Position',[ws+wu .55 wb .25]);
handles.beta2 = uicontrol('Parent',hp1,'Style','edit','String','-1',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*2 .55 wb .25]);
handles.delta1 = uicontrol('Parent',hp1,'Style','edit','String','0',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*3 .55 wb .25]);
handles.delta2 = uicontrol('Parent',hp1,'Style','edit','String','0',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*4 .55 wb .25]);
handles.epsilon = uicontrol('Parent',hp1,'Style','edit','String','.8',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*5 .55 wb .25]);
handles.fz = uicontrol('Parent',hp1,'Style','edit','String','2.1',...
  'FontSize',11,'Units','normalized','Position',[ws .05 wb .25]);
handles.F = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[ws+wu .05 wb .25]);
handles.f0 = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*2 .05 wb .25]);
handles.theta0 = uicontrol('Parent',hp1,'Style','edit','String','0',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*3 .05 wb .25]);
handles.k = uicontrol('Parent',hp1,'Style','edit','String','2',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*4 .05 wb .25]);
handles.m = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[ws+wu*5 .05 wb .25]);

hp2 = uipanel('Parent',hf,'Title','Initial Values','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Units','normalized',...
  'Position',[.025 .48 .165 .12]);
handles.r0 = uicontrol('Parent',hp2,'Style','edit','String','.1',...
  'FontSize',11,'Units','normalized','Position',[.0625 .1 .375 .5]);
handles.phi0 = uicontrol('Parent',hp2,'Style','edit','String','0',...
  'FontSize',11,'Units','normalized','Position',[.5625 .1 .375 .5]);

hp3 = uipanel('Parent',hf,'Title','Duration','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Units','normalized',...
  'Position',[.19 .48 .085 .12]);
handles.duration = uicontrol('Parent',hp3,'Style','edit','String','10',...
  'FontSize',11,'Units','normalized','Position',[.125 .1 .75 .5]);

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
text(wu*7/2,.85,'$$\delta_1$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*9/2,.85,'$$\delta_2$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*11/2,.85,'$$\epsilon$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu/2,.35,'$$\omega/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalize','FontSize',15)
text(wu*3/2,.35,'$$F$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*5/2,.35,'$$\omega_0/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*7/2,.35,'$$\theta_0/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*9/2,.35,'$$k$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(wu*11/2,.35,'$$m$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)

axes('Parent',hp2,'Units','normalized','Position',[0 0 1 1])
axis off
text(.25,.7,'$$r_0$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.75,.7,'$$\phi_0/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalize','FontSize',15)

%% Create Run button
handles.runButton = uicontrol('Parent',hf,'Style','pushbutton',...
  'String','Run','FontSize',11,'Units','normalized',...
  'Position',[.38 .48 .145 .1]);

%% Create Time Plot panel
hp4 = uipanel('Parent',hf,'Title','Choose Time Plot','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Units','normalized',...
  'Position',[.025 .015 .95 .44]);
handles.plotType = uicontrol('Parent',hp4,'Style','popup',...
  'String',['Real & imaginary parts|Amplitude|Relative phase|' ...
  'Instantaneous frequency'],...
  'FontSize',11,'Units','normalized','Position',[.015 .88 .2 .1]);
handles.ax1 = axes('Parent',hp4,'Units','normalized',...
  'Position',[.075 .2 .9 .6],'XTick',[],'YTick',[],'box','on');
text(.5,.5,'Run simulation','FontSize',15,...
  'HorizontalAlignment','center','VerticalAlignment','middle',...
  'Units','normalized')

%% Create Phase Portrait panel
hp5 = uipanel('Parent',hf,'Title','Phase Portrait & Vector Field','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Units','normalized',...
  'Position',[.55 .45 .425 .54]);
handles.ax2 = axes('Parent',hp5,'Units','normalized',...
  'Position',[.2 0 .8 1],'XTick',[],'YTick',[],'box','on');
polar2(NaN,NaN,[0 1]);

%% Define callback functions
set(handles.runButton,'Callback',{@integrate,handles})
set(handles.plotType,'Callback',{@plotAxes1,handles})
set(handles.legendButton,'Callback',@showLegend)

%% INTEGRATE
function integrate(~,~,handles)
%% Get parameters
p.a = str2double(get(handles.alpha,'String'));
p.b1 = str2double(get(handles.beta1,'String'));
p.b2 = str2double(get(handles.beta2,'String'));
p.d1 = str2double(get(handles.delta1,'String'));
p.d2 = str2double(get(handles.delta2,'String'));
p.e = str2double(get(handles.epsilon,'String'));
p.fz = str2double(get(handles.fz,'String'));
p.k = str2double(get(handles.k,'String'));
p.m = str2double(get(handles.m,'String'));

s.F = str2double(get(handles.F,'String'));
s.f = str2double(get(handles.f0,'String'));
s.theta0 = str2double(get(handles.theta0,'String'))*2*pi;
s.endtime = str2double(get(handles.duration,'String'));

%% Integrate
fs = ceil(max([p.fz; s.f]))*20;
t = (0:1/fs:s.endtime*1)';
z0 = str2double(get(handles.r0,'String')) ...
  *exp(1i*2*pi*str2double(get(handles.phi0,'String')));

Z = ode4(@(t,z)zdot(t,z,p,s),t,z0);

s.x = s.F.*exp(1i*(2*pi*s.f*t+s.theta0));
M.Z = Z; M.t = t; M.p = p; M.s = s; M.fs = fs;
guidata(handles.hf,M)
plotAxes1([],[],handles)
plotAxes2([],[],handles)

%% DRAW TIME PLOT
function plotAxes1(~,~,handles)
axes(handles.ax1)
cla
M = guidata(handles.hf);
if ~isempty(M)
  Z = M.Z; t = M.t; p = M.p; s = M.s; fs = M.fs;
  val = get(handles.plotType,'Value');
  if val == 1 % real & imaginary parts
    plot(t,real(Z),'-',t,real(s.x),'-')
    hold on
    if isprop(gca,'ColorOrderIndex')
      set(gca,'ColorOrderIndex',1)
    end
    plot(t,imag(Z),':',t,imag(s.x),':')
    hold off
    ylabel('')
    ymax = max(max(abs(Z)),s.F);
    set(gca,'YLim',[-1 1]*ymax*1.1)
    hl = legend('Re({\itz})','Re({\itx})','Imag({\itz})','Imag({\itx})',...
      'Location','NorthEast','Orientation','horizontal');
  elseif val == 2 % amplitude
    plot(t,abs(Z))
    ylabel('{\itr}')
  elseif val == 3 % relative phase
    hrp = plot(t,angle(Z.^p.m.*conj(s.x).^p.k));
    set(hrp,'LineStyle','none','Marker','.','MarkerSize',5)
    set(gca,'YLim',pi*[-1 1],'YTick',[-pi,-pi/2,0,pi/2,pi],...
      'YTickLabel',{'-pi';'-pi/2';'0';'pi/2';'pi'})
    ylabel('{\it\psi}')
  elseif val == 4 % instantaneous frequency
    instFreq = angle(Z(2:end,:).*conj(Z(1:end-1,:)))*fs/(2*pi);
    plot(t(1:end-1),instFreq,'-',t([1 end-1]),[1 1]*s.f,'--')
    ylabel('Instantaneous frequency')
    df = abs(p.fz-s.f);
    if df > 0
      set(gca,'YLim',[min(s.f,p.fz)-df/2 max(s.f,p.fz)+df/2])
    end
    hold off
    hl = legend('{\itz}','{\itx}','Location','NorthEast',...
      'Orientation','horizontal');
  end
  xlabel('Time')
  set(gca,'XLim',[min(t) max(t)])
  grid on
  if exist('hl','var')
    lpos = get(hl,'Position');
    lpos(2) = .88;
    set(hl,'Position',lpos)
  end
else
  text(.5,.5,'Run simulation first','FontSize',15,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Units','normalized')
end

%% DRAW PHASE PORTRAIT
function plotAxes2(~,~,handles)
axes(handles.ax2)
cla
M = guidata(handles.hf);
Z = M.Z; p = M.p; s = M.s;

if p.e
  rmax = 1/sqrt(p.e)*.99;
else
  rmax = max(abs(Z))*1.2;
end

polar2(NaN,NaN,[0 rmax],'-');
hold on

%% Draw nullclines
W = (p.m*p.fz - p.k*s.f)*2*pi;
x = (-1:.01:1)*rmax;
y = (-1:.01:1)*rmax;
[X,Y] = meshgrid(x,y);
[PSI,R] = cart2pol(X,Y);
RDOT = p.a*R + p.b1*R.^3 + p.b2*p.e*R.^5./(1-p.e*R.^2) + ...
  p.e^((p.k+p.m-2)/2)*s.F^p.k*R.^(p.m-1).*cos(PSI);
PSIDOT = W + p.d1*R.^2 + p.d2*p.e*R.^4./(1-p.e*R.^2) - ...
  p.m*p.e^((p.k+p.m-2)/2)*s.F^p.k*R.^(p.m-2).*sin(PSI);
RDOT(find(R>rmax)) = NaN;
PSIDOT(find(R>rmax)) = NaN;
[~,hc(1)] = contour(x,y,RDOT,[0 0],'m--','LineWidth',1.5);
[~,hc(2)] = contour(x,y,PSIDOT,[0 0],'c--','LineWidth',1.5);

%% Draw vector field
ind1 = 1:10:length(x); % use sparser grid
X1 = X(ind1,ind1);
Y1 = Y(ind1,ind1);
R1 = R(ind1,ind1);
PSI1 = PSI(ind1,ind1);
RDOT1 = RDOT(ind1,ind1);
PSIDOT1 = PSIDOT(ind1,ind1);
ind = find(~isnan(RDOT1));

XDOT1 = RDOT1.*cos(PSI1) - R1.*PSIDOT1.*sin(PSI1);
YDOT1 = RDOT1.*sin(PSI1) + R1.*PSIDOT1.*cos(PSI1);

% Normalize and apply sigmoid
absXYDOT = sqrt(XDOT1.^2+YDOT1.^2);
unitXDOT = XDOT1./absXYDOT;
unitYDOT = YDOT1./absXYDOT;

XDOT1 = unitXDOT;
YDOT1 = unitYDOT;

quiver(X1(ind),Y1(ind),XDOT1(ind),YDOT1(ind),.5,'Color','g');

%% Draw trajectory
hz1 = polar2(angle(Z.^p.m.*conj(s.x).^p.k),abs(Z),[0 rmax],'-');
set(hz1,'LineWidth',2,'Color',lines(1))
hz2 = polar2(angle(Z(1).^p.m.*conj(s.x(1)).^p.k),abs(Z(1)),[0 rmax],'ko');
set(hz2,'MarkerSize',5)
hold off
hl = legend([hc,hz1,hz2],'r-nullcline','\psi-nullcline',...
  '{\itre}^{i{\it\psi}}','Initial point','Location','NorthWest');
lpos = get(hl,'Position');
lpos(1) = .02;
lpos(2) = .75;
set(hl,'Position',lpos)

%% SHOW LEGEND
function showLegend(~,~)
hf = figure;
fpos = get(hf,'Position');

set(hf,'Position',[fpos(1) fpos(2) 350 110])
axes('Units','normalized','Position',[0 0 1 1]);
axis off
text('Interpreter','latex','String',{'$$z = re^{\textrm{i}\phi}$$,',...
  ['$$x = Fe^{\textrm{i}\theta}$$ '...
  'where $\theta = \omega_0t+\theta_0$, and'],...
  '$$\psi = m\phi - k\theta$$'},...
	'Position',[.05 .9],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',15)

%% DIFFERENTIAL EQUATION
function dzdt = zdot(t,z,p,s)
x = s.F.*exp(1i*(2*pi*s.f*t+s.theta0));
dzdt = z.*(p.a + 1i*2*pi*p.fz + (p.b1+1i*p.d1)*abs(z).^2 + ...
  p.e*(p.b2+1i*p.d2)*abs(z).^4./(1-p.e*abs(z).^2)) + ...
  p.e^((p.k+p.m-2)/2)*x^p.k*conj(z)^(p.m-1);
