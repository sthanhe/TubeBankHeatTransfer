%% Set data directories
dirFigures='Figures';       %Path to directory where figures should be stored
fname='primData.mat';


%% Load data
load(fname);


%% Bed temperature sensor analysis
names={'Strategy','mode','T1','T3','DeltaT','P_el','Tsurf'};
strat={'Pel=const.';'Tsurf-Tbed=const.'};
Tbed=table('Size',[length(strat),length(names)],...
    'VariableTypes',[{'string','logical'},repmat({'double'},1,length(names)-2)],...
    'VariableNames',names);


Tbed.Strategy=strat;
Tbed.mode(2)=true;

Tbed.T1=arrayfun(@(tf) mean(prim.T1(prim.mode==tf)),Tbed.mode);
Tbed.T3=arrayfun(@(tf) mean(prim.T3(prim.mode==tf)),Tbed.mode);
Tbed.Tsurf=arrayfun(@(tf) mean(prim.Tsurf(prim.mode==tf)),Tbed.mode);
Tbed.P_el=arrayfun(@(tf) mean(prim.P_el(prim.mode==tf)),Tbed.mode);
Tbed.DeltaT=Tbed.T3-Tbed.T1;


%Set up figure
fig=figure(912);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot
plot(ax,Tbed.mode,Tbed.T3,'Color',colors(2,:));
plot(ax,Tbed.mode,Tbed.T1,'Color',colors(1,:));

hold(ax,'off');


%Format axes and add legend
ax.XTick=[0,1];
grid(ax,'on');

xlabel(ax,'Mode');
ylabel(ax,'Temperature (K)');

legend(ax,{'T_3','T_1'},'Location','north');


%Size figure for repository
t.Units='centimeters';
t.OuterPosition=[0,0,17,8.5];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;


%Add arrows
x=[0,0];
y=[Tbed.T1(1),Tbed.T3(1)];
figaux.arrow(ax,x,y,'doublearrow');
text(ax,x(1),mean(y),compose(' \\DeltaT = %.2f K',Tbed.DeltaT(1)));

x=[1,1];
y=[Tbed.T1(2),Tbed.T3(2)];
figaux.arrow(ax,x,y,'doublearrow');
text(ax,x(1),mean(y),compose('\\DeltaT = %.2f K ',Tbed.DeltaT(2)),...
    'HorizontalAlignment','right');


%Export figure
exportgraphics(fig,[dirFigures,filesep,'T1T3.tiff'],'Resolution',600);




