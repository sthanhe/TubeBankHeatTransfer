fname='extLit.mat';


%% Load data
s=load(fname);

% load('intMeas.mat')

%Create response variable y and regressor matrix X: only data from Grewal
idx=strcmp(s.tab.Author,'Eder');
tabEder=s.tab(idx,:);
pisEder=s.pis(idx,:);


s=load('intMeas.mat');
tabOwn=s.htc;
pisOwn=s.pis;


%%
C3=1.63634185340022;
C4=0.616116450106231;

pi5=linspace(0,30,1000)';
pi10=linspace(0,2,1000)';


pisMeanEder=pisEder(1,:);
pisMeanEder{:,:}=mean(pisEder{:,:},1);

pisMeanOwn=pisOwn(1,:);
pisMeanOwn{:,:}=mean(pisOwn{:,:},1);


d_cf5=@(pi5,pis) 1-tanh((pi5./pis.pi6).^C3.*...
    (pi5./pis.pi10).^C4.*...
    (1-pis.pi9).^(C3+C4));

d_cf10=@(pi10,pis) 1-tanh((pis.pi5./pis.pi6).^C3.*...
    (pis.pi5./pi10).^C4.*...
    (1-pis.pi9).^(C3+C4));


%% Mean factor values
f1=@(pis) (pis.pi5./pis.pi6).^C3;
f2=@(pis) (pis.pi5./pis.pi10).^C4;
f3=@(pis) (1-pis.pi9).^(C3+C4);


authors={'Own','Eder'};
names=[{'Source'},compose('pi%d',[6,9,10]),compose('f%d',1:3)];
tab=table('Size',[length(authors),length(names)],...
    'VariableTypes',[{'string'},repmat({'double'},1,length(names)-1)],...
    'VariableNames',names);

tab.Source=authors';
for i=authors
    switch i{1}
        case 'Own'
            pis=pisMeanOwn;
        case 'Eder'
            pis=pisMeanEder;
    end

    idx=strcmp(tab.Source,i);

    tab.pi6(idx)=pis.pi6;
    tab.pi9(idx)=pis.pi9;
    tab.pi10(idx)=pis.pi10;

    tab.f1(idx)=mean(f1(pis));
    tab.f2(idx)=mean(f2(pis));
    tab.f3(idx)=mean(f3(pis));
end




%%
figidx=8;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


mkrSize=27;

plot(ax,pi5,[d_cf5(pi5,pisMeanOwn),d_cf5(pi5,pisMeanEder)]);

scatter(ax,pisOwn.pi5,d_cf5(pisOwn.pi5,pisOwn),mkrSize,colors(1,:),'o');

scatter(ax,pisEder.pi5,d_cf5(pisEder.pi5,pisEder),mkrSize,colors(2,:),'x');


legItems=repmat(line(ax,'Visible','off'),2,1);
legItems(1)=plot(ax,NaN,NaN,'Color',colors(1,:),...
    'Marker','o','MarkerSize',sqrt(mkrSize));
legItems(2)=plot(ax,NaN,NaN,'Color',colors(2,:),...
    'Marker','x','MarkerSize',sqrt(mkrSize));

hold(ax,'off');


lgd=legend(ax,legItems,{'Test rig','Eder et al.'},'Location','northeast');

xlabel(ax,subsz('\pi_5 (-)',6));
ylabel(ax,subsz('d_{cf} (-)',6));


%Text size
ax.FontSize=7;
lgd.FontSize=7;


%Export figure for manuscript
fname=['Figures',filesep,'Figure',num2str(figidx)];

t.Units='centimeters';
t.OuterPosition=[0,0,9,9];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
t.OuterPosition=[0,0,9,9];      %1 column
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.eps']);




%%
% C3=1.63634185340022;
% C4=0.616116450106231;
% 
% f1=@(pis) (pis.pi5./pis.pi6).^C3;
% f2=@(pis) (pis.pi5./pis.pi10).^C4;
% f3=@(pis) (1-pis.pi9).^(C3+C4);
% 
% d_cf=@(pis) 1-tanh(f1(pis).*f2(pis).*f3(pis));
% 
% 
% names={'Author','Factor','Min','Max'};
% tab=table('Size',[8,length(names)],...
%     'VariableTypes',[repmat({'string'},1,2),repmat({'double'},1,2)],...
%     'VariableNames',names);
% 
% tab.Author=repmat({'Eder';'Own'},height(tab)/2,1);
% tab.Factor=sort(repmat([compose('f%d',1:3)';{'d_cf'}],2,1));
% 
% for i=1:height(tab)
%     val=eval([tab.Factor{i},'(pis',tab.Author{i},')']);
%     tab.Min(i)=min(val);
%     tab.Max(i)=max(val(~isinf(val)));
% end

% tab.Min=arrayfun(@(i) ...
%     min(eval([tab.Factor{i},'(pis',tab.Author{i},')'])),...
%     1:height(tab));











%%
% tabEder.Nu_cf=arrayfun(@(i) ...
%     pisEder.pi1(i)-pisEder.pi1(tabEder.w==tabEder.w(i) & tabEder.w_p==0),...
%     1:height(tabEder))';



%%
% fig=figure(8);
% clf(fig);
% ax=gca();
% colors=ax.ColorOrder;
% hold(ax,'on');
% 
% 
% % s=scatter(ax,tab.Nu_cf,tab.Nu_cfExt);
% s=scatter(ax,tabEder.Nu_cf./pisEder.pi1,tabEder.Nu_cfExt./tabEder.Nu_mixExt);
% 
% eq=linspace(ax.XLim(1),ax.XLim(2),100);
% 
% plot(ax,eq,eq,'Color','k');
% % plot(ax,eq,eq.*1.2,'Color','k','LineStyle','--');
% % plot(ax,eq,eq./1.2,'Color','k','LineStyle','--');
% 
% hold(ax,'off');


% ax.YLim=[-1,7].*10^-3;
% ax.XLim=ax.YLim;