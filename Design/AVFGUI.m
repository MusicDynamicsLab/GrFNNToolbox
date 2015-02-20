function AVFGUI(alpha, beta1, beta2, epsilon, force)
% AVFGUI(alpha, beta1, beta2, epsilon, force)
%
% GUI for drawing amplitude vector field of the canonical version of
% Hoppensteadt model:
% drdt = r.*(alpha + beta1*r.^2 + beta2*epsilon*r.^4./(1-epsilon*r.^2)) + force;
%
% Change parameter values by moving the sliders or by entering numbers
% in the text boxes. Min and max values for the sliders can also be changed.
%
% 12/23/2011 JCK added markers for zero-crossings and local extremes

if nargin < 1
  alpha = -1;
end
if nargin < 2
  beta1 = 4;
end
if nargin < 3
  beta2 = -1;
end
if nargin < 4
  epsilon = 1;
end
if nargin < 5
  force = 0;
end

figure(998)
set(gcf,'Visible','off','Toolbar','figure','Color',[.8 .8 .8],...
  'Position',[100,100,1000,670]);
initialRun([],[],alpha,beta1,beta2,epsilon,force)
movegui(998,'center')
set(998,'Visible','on')

% =========================================================================
function initialRun(source,eventdata,alpha,beta1,beta2,epsilon,force)
clf
axes('Units','pixels','Position',[50 130 500 500]);

% Initial min and max for sliders
alphaLim = [-1 1]*max(1,ceil(2*abs(alpha)));
beta1Lim = [-1 1]*max(1,ceil(2*abs(beta1)));
beta2Lim = [-1 1]*max(1,ceil(2*abs(beta2)));
epsilonLim = [0 2]*max(1,ceil(2*abs(epsilon)));
forceLim = [0 1]*max(1,ceil(2*abs(force)));

% Set axis limit
rstarmax = max(rstar(alpha, beta1, beta2, epsilon, force));
rext = rextreme(alpha, beta1, beta2, epsilon);
rdotextmax = max(abs(drdt(rext, alpha, beta1, beta2, epsilon, force)));

if epsilon
  set(gca,'XLim',[-eps 1.1/sqrt(epsilon)])
elseif rstarmax
  set(gca,'XLim',[-eps 1.5*rstarmax])
elseif rext
  set(gca,'XLim',[-eps 1.5*rext])
else
  set(gca,'XLim',[-eps 1])
end
rLim = get(gca,'XLim');

if rdotextmax
  set(gca,'YLim',[-1 1]*rdotextmax*2)
else
  r = (0:.01:1)*rLim(2);
  rdot = drdt(r, alpha, beta1, beta2, epsilon, force);
  set(gca,'YLim',[-1 1]*max(abs(rdot)))
end
rdotLim = get(gca,'YLim');

% Construct UI components
slider1 = uicontrol('Style','slider','Value',alpha,...
  'Min',alphaLim(1),'Max',alphaLim(2),'Position',[600,600,250,25]);
slider2 = uicontrol('Style','slider','Value',beta1,...
  'Min',beta1Lim(1),'Max',beta1Lim(2),'Position',[600,500,250,25]);
slider3 = uicontrol('Style','slider','Value',beta2,...
  'Min',beta2Lim(1),'Max',beta2Lim(2),'Position',[600,400,250,25]);
slider4 = uicontrol('Style','slider','Value',epsilon,...
  'Min',epsilonLim(1),'Max',epsilonLim(2),'Position',[600,300,250,25]);
slider5 = uicontrol('Style','slider','Value',force,...
  'Min',forceLim(1),'Max',forceLim(2),'Position',[600,200,250,25]);

min1 = uicontrol('Style','edit','String',alphaLim(1),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[600,570,40,25]);
min2 = uicontrol('Style','edit','String',beta1Lim(1),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[600,470,40,25]);
min3 = uicontrol('Style','edit','String',beta2Lim(1),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[600,370,40,25]);
min4 = uicontrol('Style','edit','String',epsilonLim(1),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[600,270,40,25]);
min5 = uicontrol('Style','edit','String',forceLim(1),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[600,170,40,25]);

max1 = uicontrol('Style','edit','String',alphaLim(2),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[810,570,40,25]);
max2 = uicontrol('Style','edit','String',beta1Lim(2),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[810,470,40,25]);
max3 = uicontrol('Style','edit','String',beta2Lim(2),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[810,370,40,25]);
max4 = uicontrol('Style','edit','String',epsilonLim(2),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[810,270,40,25]);
max5 = uicontrol('Style','edit','String',forceLim(2),...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[810,170,40,25]);

text1 = uicontrol('Style','edit','String',alpha,...
  'FontSize',11,'Position',[900,600,70,25]);
text2 = uicontrol('Style','edit','String',beta1,...
  'FontSize',11,'Position',[900,500,70,25]);
text3 = uicontrol('Style','edit','String',beta2,...
  'FontSize',11,'Position',[900,400,70,25]);
text4 = uicontrol('Style','edit','String',epsilon,...
  'FontSize',11,'Position',[900,300,70,25]);
text5 = uicontrol('Style','edit','String',force,...
  'FontSize',11,'Position',[900,200,70,25]);
zeroX = uicontrol('Style','text','String',{'-'},'HorizontalAlignment','left',...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[160,50,600,25]);
extreme = uicontrol('Style','text','String',{'-'},'HorizontalAlignment','left',...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[160,20,600,25]);
rmax = uicontrol('Style','edit','String',rLim(2),...
  'FontSize',11,'Position',[600,100,70,25]);
rdotmax = uicontrol('Style','edit','String',max(abs(rdotLim)),...
  'FontSize',11,'Position',[700,100,70,25]);

startOver = uicontrol('Style', 'pushbutton', 'String', 'Start Over',...
  'FontSize',10,'Position',[600 40 100 30],...
  'Callback', {@initialRun,alpha,beta1,beta2,epsilon,force});
initialPosition = uicontrol('Style', 'pushbutton', 'String', 'Initial Position',...
  'FontSize',10,'Position',[720 40 100 30],...
  'Callback', {@initialPosition_Callback,alpha,beta1,beta2,epsilon,force});
holdControl = uicontrol('Style','checkbox','String','Hold On','Value',0,...
  'FontSize',10,'Position',[850 100 150 30],'BackgroundColor',[.8 .8 .8]);
fixedPtControl = uicontrol('Style','checkbox','String','Show fixed points','Value',1,...
  'FontSize',10,'Position',[850 70 150 30],'BackgroundColor',[.8 .8 .8]);
arrowControl = uicontrol('Style','checkbox','String','Show arrows','Value',1,...
  'FontSize',10,'Position',[850 40 150 30],'BackgroundColor',[.8 .8 .8]);

label1 = uicontrol('Style','text','String',{'Alpha'},...
  'BackgroundColor',[.8 .8 .8],'FontSize',11,'Position',[700,625,100,25]);
label2 = uicontrol('Style','text','String',{'Beta1'},...
  'BackgroundColor',[.8 .8 .8],'FontSize',11,'Position',[700,525,100,25]);
label3 = uicontrol('Style','text','String',{'Beta2'},...
  'BackgroundColor',[.8 .8 .8],'FontSize',11,'Position',[700,425,100,25]);
label4 = uicontrol('Style','text','String',{'Epsilon'},...
  'BackgroundColor',[.8 .8 .8],'FontSize',11,'Position',[700,325,100,25]);
label5 = uicontrol('Style','text','String',{'External Force'},...
  'BackgroundColor',[.8 .8 .8],'FontSize',11,'Position',[700,225,100,25]);
label6 = uicontrol('Style','text','String',{'Fixed Point:'},...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[50,50,100,25]);
label7 = uicontrol('Style','text','String',{'Local Extreme:'},...
  'BackgroundColor',[.8 .8 .8],'FontSize',10,'Position',[50,20,100,25]);
label8 = uicontrol('Style','text','String',{'r max'},...
  'BackgroundColor',[.8 .8 .8],'FontSize',11,'Position',[600,125,70,25]);
label9 = uicontrol('Style','text','String',{'|dr/dt| max'},...
  'BackgroundColor',[.8 .8 .8],'FontSize',11,'Position',[700,125,70,25]);

align([slider1,slider2,slider3,slider4,slider5,...
  label1,label2,label3,label4,label5],'Center','None')
set([label6,label7],'HorizontalAlignment','left')

% Create an array of handles
handles = [ slider1 min1 max1;
            slider2 min2 max2;
            slider3 min3 max3;
            slider4 min4 max4;
            slider5 min5 max5;
            text1   0    0;
            text2   0    0;
            text3   0    0;
            text4   0    0;
            text5   0    0;
            zeroX   0    0;
            extreme 0    0;
            holdControl fixedPtControl arrowControl;
            rmax rdotmax 0];

% Set callback functions
set(handles(1:5,1),'Callback',{@slider_Callback,handles})
set(handles(6:10,1),'Callback',{@text_Callback,handles})
set(handles(1:5,2),'Callback',{@min_Callback,handles})
set(handles(1:5,3),'Callback',{@max_Callback,handles})
set(handles(13,1),'Callback',{@holdControl_Callback,handles})
set(handles(13,2:3),'Callback',{@replot_Callback,handles})
set(handles(14,1:2),'Callback',{@resetLim_Callback,handles})

hold off

drawVectorField(alpha,beta1,beta2,epsilon,force,handles)

% =========================================================================
function slider_Callback(source,eventdata,handles)
[row,col] = find(handles == source);
set(handles(row+5,1),'String',num2str(get(source,'Value'),3))

parameterValues = get(handles(1:5,1),'Value');
drawVectorField(parameterValues{1},parameterValues{2},parameterValues{3},...
  parameterValues{4},parameterValues{5},handles)

% =========================================================================
function text_Callback(source,eventdata,handles)
[row,col] = find(handles == source);
sliderMin = get(handles(row-5,1),'Min');
sliderMax = get(handles(row-5,1),'Max');
textValue = str2double(get(source,'String'));

if textValue < sliderMin
  set(handles(row-5,1),'Min',textValue)
  set(handles(row-5,2),'String',num2str(textValue))
end

if textValue > sliderMax
  set(handles(row-5,1),'Max',textValue)
  set(handles(row-5,3),'String',num2str(textValue))
end

set(handles(row-5,1),'Value',str2double(get(source,'String')))

parameterValues = str2double(get(handles(6:10,1),'String'));
drawVectorField(parameterValues(1),parameterValues(2),parameterValues(3),...
  parameterValues(4),parameterValues(5),handles)

% =========================================================================
function min_Callback(source,eventdata,handles)
[row,col] = find(handles == source);
newMin = str2double(get(source,'String'));
sliderValue = get(handles(row,1),'Value');
currentMax = get(handles(row,1),'Max');

if newMin >= currentMax
  set(handles(row,1),'Value',newMin)
  set(handles(row+5,1),'String',num2str(newMin))
  set(handles(row,1),'Max',newMin+.1)
  set(handles(row,3),'String',num2str(newMin+.1))
  parameterValues = str2double(get(handles(6:10,1),'String'));
  drawVectorField(parameterValues(1),parameterValues(2),...
    parameterValues(3),parameterValues(4),parameterValues(5),handles)
elseif newMin > sliderValue
  set(handles(row,1),'Value',newMin)
  set(handles(row+5,1),'String',num2str(newMin))
  parameterValues = str2double(get(handles(6:10,1),'String'));
  drawVectorField(parameterValues(1),parameterValues(2),...
    parameterValues(3),parameterValues(4),parameterValues(5),handles)
end

set(handles(row,1),'Min',newMin)

% =========================================================================
function max_Callback(source,eventdata,handles)
[row,col] = find(handles == source);
newMax = str2double(get(source,'String'));
sliderValue = get(handles(row,1),'Value');
currentMin = get(handles(row,1),'Min');

if newMax <= currentMin
  set(handles(row,1),'Value',newMax)
  set(handles(row+5,1),'String',num2str(newMax))
  set(handles(row,1),'Min',newMax-.1)
  set(handles(row,2),'String',num2str(newMax-.1))
  parameterValues = str2double(get(handles(6:10,1),'String'));
  drawVectorField(parameterValues(1),parameterValues(2),...
    parameterValues(3),parameterValues(4),parameterValues(5),handles)
elseif newMax < sliderValue
  set(handles(row,1),'Value',newMax)
  set(handles(row+5,1),'String',num2str(newMax))
  parameterValues = str2double(get(handles(6:10,1),'String'));
  drawVectorField(parameterValues(1),parameterValues(2),...
    parameterValues(3),parameterValues(4),parameterValues(5),handles)
end

set(handles(row,1),'Max',newMax)

% =========================================================================
function initialPosition_Callback(source,eventdata,alpha,beta1,beta2,epsilon,force)
rstarmax = max(rstar(alpha, beta1, beta2, epsilon, force));
rext = rextreme(alpha, beta1, beta2, epsilon);
rdotextmax = max(abs(drdt(rext, alpha, beta1, beta2, epsilon, force)));

if epsilon
  set(gca,'XLim',[-eps 1.1/sqrt(epsilon)])
elseif rstarmax
  set(gca,'XLim',[-eps 1.5*rstarmax])
elseif rext
  set(gca,'XLim',[-eps 1.5*rext])
else
  set(gca,'XLim',[-eps 1])
end
rLim = get(gca,'XLim');

if rdotextmax
  set(gca,'YLim',[-1 1]*rdotextmax*2)
else
  r = (0:.01:1)*rLim(2);
  rdot = drdt(r, alpha, beta1, beta2, epsilon, force);
  set(gca,'YLim',[-1 1]*max(abs(rdot)))
end

% =========================================================================
function holdControl_Callback(source,eventdata,handles)
if get(source,'Value')
  hold on
else
  hold off
  parameterValues = str2double(get(handles(6:10,1),'String'));
  drawVectorField(parameterValues(1),parameterValues(2),...
    parameterValues(3),parameterValues(4),parameterValues(5),handles)
end

% =========================================================================
function replot_Callback(source,eventdata,handles)
parameterValues = str2double(get(handles(6:10,1),'String'));
drawVectorField(parameterValues(1),parameterValues(2),...
  parameterValues(3),parameterValues(4),parameterValues(5),handles)

% =========================================================================
function resetLim_Callback(source,eventdata,handles)
rmax = str2double(get(handles(14,1),'String'));
rdotmax = str2double(get(handles(14,2),'String'));
set(gca,'XLim',[-eps rmax],'YLim',[-1 1]*rdotmax)
parameterValues = str2double(get(handles(6:10,1),'String'));
drawVectorField(parameterValues(1),parameterValues(2),...
  parameterValues(3),parameterValues(4),parameterValues(5),handles)

% =========================================================================
function drawVectorField(alpha,beta1,beta2,epsilon,force,handles)

axisLimit = get(gca,{'XLim','YLim'});

realrstar = rstar(alpha, beta1, beta2, epsilon, force);
realrext = rextreme(alpha, beta1, beta2, epsilon);
rdotext = drdt(realrext, alpha, beta1, beta2, epsilon, force);

if epsilon
  r = (0:0.001:1)/sqrt(epsilon)*(1-eps);
else
  r = (0:0.001:1)*max(max(realrstar)*2,max(axisLimit{1})*1.5);
end
rdot = drdt(r, alpha, beta1, beta2, epsilon, force);

% Draw rdot curve
plot(r,rdot)
hold on
string1 = '';
for nr = 1:length(realrstar)
  string1 = [string1 '(' num2str(realrstar(nr)) ', 0)  '];
end

% Draw fixed points
if get(handles(13,2),'Value')
  for nr = 1:length(realrstar)
    h = plot(realrstar(nr),0,'ro','MarkerSize',7);
    if slope(realrstar(nr),alpha,beta1,beta2,epsilon) < 0
      set(h,'MarkerFaceColor','r')
    elseif drdt(realrstar(nr)-eps,alpha,beta1,beta2,epsilon,force) > 0 &&...
        drdt(realrstar(nr)+eps,alpha,beta1,beta2,epsilon,force) < 0
      set(h,'MarkerFaceColor','r')
    end
  end
end

% Draw arrows
if get(handles(13,3),'Value')
  if epsilon
    bound = unique([0; realrstar; 1/sqrt(epsilon)]);
  else
    bound = unique([0; realrstar; max(realrstar)*1.5]);
  end
  if length(bound) == 1
    bound = axisLimit{1};
  end
  midpt = bound(1:end-1)+diff(bound)/2;
  for nm = 1:length(midpt)
    rm = midpt(nm);
    drdtmid = drdt(rm, alpha, beta1, beta2, epsilon, force);
    if drdtmid < 0
      har = plot(rm,0,'k<');
    elseif drdtmid > 0
      har = plot(rm,0,'k>');
    elseif drdtmid == 0
      har = plot(rm,0,'ks');
    end
    set(har,'MarkerSize',7,'MarkerFaceColor','k')
  end
end

string2 = '';
for nr = 1:length(realrext)
  plot(realrext(nr),rdotext(nr),'bx','MarkerSize',7)
  string2 = [string2 '(' num2str(realrext(nr)) ', ' num2str(rdotext(nr)) ')  '];
end
plot([0 max([100, max(r)*2])],[0 0],'k-') % r-axis
plot([0 0],[-1 1]*max([100, max(rdot)*2]),'k-') % rdot-axis
if epsilon
  plot([1 1]/sqrt(epsilon),[-1 1]*max([100, max(rdot)*2]),'r--') % r bound
end
set(gca,'XLim',axisLimit{1},'YLim',axisLimit{2})
xlabel('{\itr}','FontSize',12);
ylabel('{\itdr}/{\itdt}','FontSize',12);
title('Amplitude Vector Field (\Omega = 0)','FontSize',14,'FontWeight','bold');
grid on
if ~get(handles(13,1),'Value')
  hold off
end
set(handles(11,1),'String',string1)
set(handles(12,1),'String',string2)

% =========================================================================
function realrstar = rstar(alpha, beta1, beta2, epsilon, force)
rstar = roots([epsilon*(beta2-beta1), 0, beta1-epsilon*alpha,...
  -epsilon*force, alpha, force]); % rstar (zero crossing)
rstar = unique(rstar);
realrstar = sort(rstar(intersect(find(imag(rstar)==0),find(real(rstar)>=-eps))));
                                  % only positive real solutions
if epsilon
  realrstar = realrstar(find(realrstar < 1/sqrt(epsilon)));
end

% =========================================================================
function realrext = rextreme(alpha, beta1, beta2, epsilon)
rext = sqrt(roots([3*epsilon^2*(beta1-beta2),...
  epsilon*(epsilon*alpha-6*beta1+5*beta2),...
  -2*epsilon*alpha+3*beta1, alpha])); % local extremes
realrext = sort(rext(intersect(find(imag(rext)==0),find(rext > 0))));
if epsilon
  realrext = realrext(find(realrext < 1/sqrt(epsilon)));
end
% =========================================================================
function rdot = drdt(r, alpha, beta1, beta2, epsilon, force)
rdot = alpha*r + beta1*r.^3 + beta2*epsilon*r.^5./(1-epsilon*r.^2) + force;

% =========================================================================
function drdotdr = slope(r, alpha, beta1, beta2, epsilon)
drdotdr = alpha + 3*beta1*r.^2 + ...
  (5*epsilon*beta2*r.^4-3*epsilon^2*beta2*r.^6)./((1-epsilon*r.^2).^2);
