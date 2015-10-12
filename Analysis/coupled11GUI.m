%% CREATE FIGURE WINDOW
function coupled11GUI
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
text('String','Two Oscillators with Linear Coupling',...
  'Position',[0 1],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',13,'FontWeight','bold')
text('Interpreter','latex',...
	'String',{['$$\frac{dz_1}{dt} = z_1\left(\alpha + ' ...
  '\textrm{i}\omega_1 + (\beta_1+\textrm{i}\delta_1) |z_1|^2 + ' ...
  '\frac{\epsilon(\beta_2+\textrm{i}\delta_2)|z_1|^4}' ...
  '{1-\epsilon |z_1|^2}\right) + cz_2$$'],...
  ['$$\frac{dz_2}{dt} = z_2\left(\alpha + \textrm{i}\omega_2 + ' ...
  '(\beta_1+\textrm{i}\delta_1) |z_2|^2 + ' ...
  '\frac{\epsilon(\beta_2+\textrm{i}\delta_2)|z_2|^4}' ...
  '{1-\epsilon |z_2|^2}\right) + cz_1$$']},...
	'Position',[.02 .45],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',11)

%% Create Legend button
handles.legendButton = uicontrol('Parent',hf,'Style','pushbutton',...
  'String','Legend','FontSize',10,'Units','normalized',...
  'Position',[.455 .925 .07 .05]);

%% Make text boxes for parameter input
hp1 = uipanel('Parent',hf,'Title','Parameters & Duration','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Units','normalized',...
  'Position',[.025 .59 .5 .2]);
handles.alpha = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[.025 .55 .15 .25]);
handles.beta1 = uicontrol('Parent',hp1,'Style','edit','String','-10',...
  'FontSize',11,'Units','normalized','Position',[.225 .55 .15 .25]);
handles.beta2 = uicontrol('Parent',hp1,'Style','edit','String','0',...
  'FontSize',11,'Units','normalized','Position',[.425 .55 .15 .25]);
handles.delta1 = uicontrol('Parent',hp1,'Style','edit','String','0',...
  'FontSize',11,'Units','normalized','Position',[.625 .55 .15 .25]);
handles.delta2 = uicontrol('Parent',hp1,'Style','edit','String','0',...
  'FontSize',11,'Units','normalized','Position',[.825 .55 .15 .25]);
handles.epsilon = uicontrol('Parent',hp1,'Style','edit','String','0',...
  'FontSize',11,'Units','normalized','Position',[.025 .05 .15 .25]);
handles.fz1 = uicontrol('Parent',hp1,'Style','edit','String','1.1',...
  'FontSize',11,'Units','normalized','Position',[.225 .05 .15 .25]);
handles.fz2 = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[.425 .05 .15 .25]);
handles.c = uicontrol('Parent',hp1,'Style','edit','String','1',...
  'FontSize',11,'Units','normalized','Position',[.625 .05 .15 .25]);
handles.duration = uicontrol('Parent',hp1,'Style','edit','String','10',...
  'FontSize',11,'Units','normalized','Position',[.825 .05 .15 .25]);

hp2 = uipanel('Parent',hf,'Title','Initial Values','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Units','normalized',...
  'Position',[.025 .46 .4 .12]);
handles.r10 = uicontrol('Parent',hp2,'Style','edit','String','.6',...
  'FontSize',11,'Units','normalized','Position',[.0313 .1 .1875 .5]);
handles.phi10 = uicontrol('Parent',hp2,'Style','edit','String','.4',...
  'FontSize',11,'Units','normalized','Position',[.2813 .1 .1875 .5]);
handles.r20 = uicontrol('Parent',hp2,'Style','edit','String','.1',...
  'FontSize',11,'Units','normalized','Position',[.5313 .1 .1875 .5]);
handles.phi20 = uicontrol('Parent',hp2,'Style','edit','String','.7',...
  'FontSize',11,'Units','normalized','Position',[.7813 .1 .1875 .5]);

%% Labels for text boxes
axes('Parent',hp1,'Units','normalized','Position',[0 0 1 1])
axis off
text(.1,.85,'$$\alpha$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.3,.85,'$$\beta_1$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.5,.85,'$$\beta_2$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.7,.85,'$$\delta_1$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.9,.85,'$$\delta_2$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.1,.35,'$$\epsilon$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.3,.35,'$$\omega_1/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalize','FontSize',15)
text(.5,.35,'$$\omega_2/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.7,.35,'$$c$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.9,.35,'Duration','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',13)

axes('Parent',hp2,'Units','normalized','Position',[0 0 1 1])
axis off
text(.125,.7,'$$r_1$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.375,.7,'$$\phi_1/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalized','FontSize',15)
text(.625,.7,'$$r_2$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalize','FontSize',15)
text(.875,.7,'$$\phi_2/2\pi$$','Interpreter','Latex',...
  'HorizontalAlignment','center','VerticalAlignment','baseline',...
  'Units','normalize','FontSize',15)

%% Create Run button
handles.runButton = uicontrol('Parent',hf,'Style','pushbutton',...
  'String','Run','FontSize',11,'Units','normalized',...
  'Position',[.445 .46 .08 .1]);

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
hp5 = uipanel('Parent',hf,'Title','2D Phase Portrait','FontSize',11,...
  'BackgroundColor',[.8 .8 .8],'Units','normalized',...
  'Position',[.54 .45 .435 .54]);
handles.ax2 = axes('Parent',hp5,'Units','normalized',...
  'Position',[.2 0 .8 1],'XTick',[],'YTick',[],'box','on');
polar2(NaN,NaN,[0 1],'-');

%% Create 3D phase portrait button
handles.threeD = uicontrol('Parent',hp5,'Style','pushbutton',...
  'String','3D Portrait','FontSize',11,'Units','normalized',...
  'Position',[.01 .01 .2 .2]);

%% Define callback functions
set(handles.runButton,'Callback',{@integrate,handles})
set(handles.plotType,'Callback',{@plotAxes1,handles})
set(handles.threeD,'Callback',{@plot3D,handles})
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
p.fz1 = str2double(get(handles.fz1,'String'));
p.fz2 = str2double(get(handles.fz2,'String'));
p.c = str2double(get(handles.c,'String'));
endtime = str2double(get(handles.duration,'String'));
z10 = str2double(get(handles.r10,'String')) ...
  *exp(1i*2*pi*str2double(get(handles.phi10,'String')));
z20 = str2double(get(handles.r20,'String')) ...
  *exp(1i*2*pi*str2double(get(handles.phi20,'String')));

%% Integrate
fs = ceil(max(p.fz1,p.fz2))*20;
t = (0:1/fs:endtime*1)';

Z = ode4(@(t,z)zdot(z,p),t,[z10; z20]);

Z1 = Z(:,1);
Z2 = Z(:,2);
M.Z1 = Z1; M.Z2 = Z2;M.t = t; M.p = p; M.fs = fs;
guidata(handles.hf,M)
plotAxes1([],[],handles)
plotAxes2([],[],handles)

%% DRAW TIME PLOT
function plotAxes1(~,~,handles)
axes(handles.ax1)
cla
M = guidata(handles.hf);
if ~isempty(M)
  Z1 = M.Z1; Z2 = M.Z2; t = M.t; p = M.p; fs = M.fs;
  val = get(handles.plotType,'Value');
  if val == 1 % real & imaginary parts
    plot(t,real(Z1),t,real(Z2))
    hold on
    if isprop(gca,'ColorOrderIndex')
      set(gca,'ColorOrderIndex',1)
    end
    plot(t,imag(Z1),':',t,imag(Z2),':')
    hold off
    ylabel('')
    rmax = max(abs([Z1;Z2]));
    set(gca,'YLim',[-1 1]*rmax*1.1)
    hl = legend('Re({\itz}_1)','Re({\itz}_2)','Imag({\itz}_1)',...
      'Imag({\itz}_2)','Location','NorthEast','Orientation','horizontal');
  elseif val == 2 % amplitude
    plot(t,abs(Z1),t,abs(Z2))
    ylabel('{\itr}_i')
    hl = legend('{\itr}_1','{\itr}_2',...
      'Location','NorthEast','Orientation','horizontal');
  elseif val == 3 % relative phase
    hrp = plot(t,angle(Z1.*conj(Z2)));
    set(hrp,'LineStyle','none','Marker','.','MarkerSize',5)
    set(gca,'YLim',pi*[-1 1],'YTick',[-pi,-pi/2,0,pi/2,pi],...
      'YTickLabel',{'-pi';'-pi/2';'0';'pi/2';'pi'})
    ylabel('{\it\psi}')
  elseif val == 4 % instantaneous frequency
    instFreq1 = angle(Z1(2:end,:).*conj(Z1(1:end-1,:)))*fs/(2*pi);
    instFreq2 = angle(Z2(2:end,:).*conj(Z2(1:end-1,:)))*fs/(2*pi);
    plot(t(1:end-1),instFreq1,t(1:end-1),instFreq2)
    hold on
    if isprop(gca,'ColorOrderIndex')
      set(gca,'ColorOrderIndex',1)
    end
    plot(t([1 end-1]),[1 1]*p.fz1,'--',t([1 end-1]),[1 1]*p.fz2,'--')
    ylabel('Instantaneous frequency')
    df = abs(p.fz1-p.fz2);
    if df > 0
      set(gca,'YLim',[min(p.fz1,p.fz2)-df/2 max(p.fz1,p.fz2)+df/2])
    end
    hold off
    hl = legend('Instant freq {\itz}_1','Instant freq {\itz}_2',...
      'Natural freq {\itz}_1','Natural freq {\itz}_2',...
      'Location','NorthEast','Orientation','horizontal');
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
Z1 = M.Z1; Z2 = M.Z2; t = M.t; p = M.p; fs = M.fs;
if p.e
  rmax = 1/sqrt(p.e)*.99;
else
  rmax = max(abs([Z1;Z2]))*1.2;
end
lcolor = lines(2);
hz1 = polar2(angle(Z1.*conj(Z2)),abs(Z1),[0 rmax],'-');
set(hz1,'LineWidth',2,'Color',lcolor(1,:))
hold on
hz2 = polar2(angle(Z2.*conj(Z1)),abs(Z2),[0 rmax],'-');
set(hz2,'LineWidth',2,'Color',lcolor(2,:))
hz10 = polar2(angle(Z1(1).*conj(Z2(1))),abs(Z1(1)),[0 rmax],'ko');
hz20 = polar2(angle(Z2(1).*conj(Z1(1))),abs(Z2(1)),[0 rmax],'ko');
set([hz10 hz20],'MarkerSize',5)
hold off
hl = legend([hz1,hz2,hz10],'{\itr}_1{\ite}^{i{\it\psi}}',...
  '{\itr}_2{\ite}^{-i{\it\psi}}','Initial pts','Location','NorthWest');
lpos = get(hl,'Position');
lpos(1) = .01;
lpos(2) = .75;
set(hl,'Position',lpos)

%% DRAW 3D PHASE PORTRAIT
function plot3D(~,~,handles)
M = guidata(handles.hf);
if ~isempty(M)
  Z1 = M.Z1; Z2 = M.Z2; t = M.t; p = M.p; fs = M.fs;
  % Get (X,Y,Z) meshgrid
  if p.e
    xymax = 1/sqrt(p.e)*.99*sqrt(2);
  else
    xymax = max(abs([Z1;Z2]))*sqrt(2)*1.2;
  end
  zmax = xymax/2;
  x = (-1:.2:1)*xymax;
  y = (-1:.2:1)*xymax;
  z = (-1:.2:1)*zmax;
  [X,Y,Z] = meshgrid(x,y,z);

  % Get vector field
  R1 = (sqrt(X.^2+Y.^2)+Z)/sqrt(2);
  R2 = (sqrt(X.^2+Y.^2)-Z)/sqrt(2);
  cosPSI = X./sqrt(X.^2+Y.^2);
  sinPSI = Y./sqrt(X.^2+Y.^2);

  W = (p.fz1 - p.fz2)*2*pi;
  R1DOT = p.a*R1 + p.b1*R1.^3 + p.e*p.b2*R1.^5./(1-p.e*R1.^2)...
    + p.c*R2.*cosPSI;
  R2DOT = p.a*R2 + p.b1*R2.^3 + p.e*p.b2*R2.^5./(1-p.e*R2.^2)...
    + p.c*R1.*cosPSI;
  PSIDOT = W - p.c*(R1.^2+R2.^2)./(R1.*R2).*sinPSI;

  XDOT = ((R1DOT + R2DOT).*cosPSI - (R1 + R2).*sinPSI.*PSIDOT)/sqrt(2);
  YDOT = ((R1DOT + R2DOT).*sinPSI + (R1 + R2).*cosPSI.*PSIDOT)/sqrt(2);
  ZDOT = (R1DOT - R2DOT)/sqrt(2);
  
  absXYZDOT = sqrt(XDOT.^2+YDOT.^2+ZDOT.^2); % Normalize to unit vectors
  XDOT = XDOT./absXYZDOT;
  YDOT = YDOT./absXYZDOT;
  ZDOT = ZDOT./absXYZDOT;

  ind = intersect(find(R1 >= 0),find(R2 >= 0));
  ind = intersect(ind,find(sqrt(X.^2+Y.^2) <= xymax));
  if p.e
    ind = intersect(ind,find(R1 < 1/sqrt(p.e)));
    ind = intersect(ind,find(R2 < 1/sqrt(p.e)));
  end

  figure
  hq = quiver3(X(ind),Y(ind),Z(ind),XDOT(ind),YDOT(ind),ZDOT(ind),.7,...
    'Color','g');
  hold on

  % Draw fixed point from analysis
  [r1star,r2star,psistar] = rStarCoupled11(p.a,p.b1,p.b2,p.e,p.c,W,0);
  xstar = (r1star + r2star)/sqrt(2).*cos(psistar);
  ystar = (r1star + r2star)/sqrt(2).*sin(psistar);
  zstar = (r1star - r2star)/sqrt(2);
  plot3(xstar,ystar,zstar,'r.','MarkerSize',10)

  % Draw initial point
  z01 = Z1(1); z02 = Z2(1);
  psi0 = angle(z01*conj(z02));
  x0 = (abs(z01) + abs(z02))/sqrt(2).*cos(psi0);
  y0 = (abs(z01) + abs(z02))/sqrt(2).*sin(psi0);
  z0 = (abs(z01) - abs(z02))/sqrt(2);
  plot3(x0,y0,z0,'ko','MarkerSize',5)

  % Plot trajectories
  psi = angle(Z1.*conj(Z2));
  r1 = abs(Z1);
  r2 = abs(Z2);
  x = (r1+r2)/sqrt(2).*cos(psi);
  y = (r1+r2)/sqrt(2).*sin(psi);
  z = (r1-r2)/sqrt(2);
  plot3(x,y,z,'-','LineWidth',2,'Color',lines(1))
  
  % Draw axes
  plot3([0 zmax],[0 0],[0 zmax],'k-','LineWidth',2) % r1 axis
  plot3([0 zmax],[0 0],[0 -zmax],'k-','LineWidth',2) % r2 axis
  gridangle1 = (0:5)*pi/6;
  plot3([-1;1]*cos(gridangle1)*xymax,[-1;1]*sin(gridangle1)*xymax,...
    zeros(2,length(gridangle1)),'k:')
  text(-xymax*.05,0,0,'0')
  text(zmax,0,zmax*1.1,'{\itr}_1','FontSize',12,'FontWeight','bold')
  text(zmax,0,-zmax*1.1,'{\itr}_2','FontSize',12,'FontWeight','bold')
  text(xymax*1.2,0,0,'{\it\psi}','FontSize',12,'FontWeight','bold')
  for n = 0:11
    text(xymax*cos(n*pi/6)*1.1,xymax*sin(n*pi/6)*1.1,0,num2str(n*30),...
      'HorizontalAlignment','center','VerticalAlignment','middle')
  end

  patch([0 zmax xymax xymax zmax],[0 0 0 0 0],[0 zmax zmax -zmax -zmax],...
    [1 1 1]*.7,'FaceAlpha',.5,'EdgeColor','none') % plane w/ psi = 0
  gridangle2 = (0:.01:.99)*2*pi;
  patch(xymax*cos(gridangle2),xymax*sin(gridangle2),...
    zeros(size(gridangle2)),[1 1 1]*.7,'FaceAlpha',.5,'EdgeColor','none')
                                                  % mid plane (r1 = r2)
  axis equal
  axis off
  hold off
  set(gcf,'Color',[1 1 1])
else
  axes(handles.ax1)
  cla
  text(.5,.5,'Run simulation first','FontSize',15,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Units','normalized')
end

%% SHOW LEGEND
function showLegend(~,~)
hf = figure;
fpos = get(hf,'Position');

set(hf,'Position',[fpos(1) fpos(2) 300 80])
axes('Units','normalized','Position',[0 0 1 1]);
axis off
text('Interpreter','latex','String',...
  ['$$z_i = r_ie^{\textrm{i}\phi_i}$$ and ' ...
  '$$\psi = \phi_1 - \phi_2$$'],...
	'Position',[.05 .8],'Units','normalized',...
  'HorizontalAlignment','left','VerticalAlignment','top',...
	'FontSize',15)

%% DIFFERENTIAL EQUATION
function dzdt = zdot(z,p)
z1 = z(1);
z2 = z(2);
dz1dt = z1.*(p.a + 1i*2*pi*p.fz1 + (p.b1+1i*p.d1)*abs(z1).^2 + ...
  p.e*(p.b2+1i*p.d2)*abs(z1).^4./(1-p.e*abs(z1).^2)) + p.c*z2;
dz2dt = z2.*(p.a + 1i*2*pi*p.fz2 + (p.b1+1i*p.d1)*abs(z2).^2 + ...
  p.e*(p.b2+1i*p.d2)*abs(z2).^4./(1-p.e*abs(z2).^2)) + p.c*z1;
dzdt = [dz1dt; dz2dt];
