% Plots r* and psi* for afferent chain w/ one-to-one connections
a =  [   0   0.01  0.01];
b1 = [-100 -10    -.1];
b2 = [  -1  -1    -1];
e =  [   0.0025   0.107   0.350];
typeLabel = {'Layer 1 (critical LC)','Layer 2 (supercritical LC)',...
  'Layer 3 (subcritical DLC)'};

%F = (.1:.1:1)*1;
F = dB2Pa(80);

W = (-10:.001:10)*2*pi*.5;

figure(101)
nlayer = length(a);
maxnr = [1 1 2];
lcolor = {'b','g','r'};
rstar = cell(nlayer,1);
psistar = cell(nlayer,1);
rmax = zeros(nlayer,1);

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
  ind = find(rstar{nl} < 0); % indices for unstable fixed pts
  rstar{nl}(ind) = NaN;
  psistar{nl}(ind) = NaN;
  rmax(nl) = max(rstar{nl}(:));
  for nr = 1:maxnr(nl)
    subplot(2,1,1)
    h(nl) = plot(W/(2*pi),rstar{nl}(:,nr),'.--','MarkerSize',5,...
      'Color',lcolor{nl});
    hold on
    subplot(2,1,2)
    plot(W/(2*pi),psistar{nl}(:,nr),'.--','MarkerSize',5,...
      'Color',lcolor{nl});
    hold on
  end
end

subplot(2,1,1)
hold off
grid on
legend(h,typeLabel,'Location','NorthWest')
set(gca,'YLim',[0 max(rmax)*1.05])
%set(gca,'XTick',(-.5:.1:.5))
xlabel('\Omega/2\pi (Hz)','FontSize',12)
ylabel('r*','FontSize',12)
title(['Afferent Chain (F = ' num2str(F) ')'],'FontSize',13)

subplot(2,1,2)
hold off
grid on
set(gca,'Ylim',[-pi pi],'YTick',[-pi,-pi/2,0,pi/2,pi],...
    'YTickLabel',{'-pi';'-pi/2';'0';'pi/2';'pi'})
%set(gca,'XTick',(-.5:.1:.5))
xlabel('\Omega/2\pi (Hz)','FontSize',12)
ylabel('\psi*','FontSize',12)