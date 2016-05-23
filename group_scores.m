%% Calculate group scores and plot for all group games

load experiment
load peoeva

% Find trial numbers
Ind50  = find(gamestatus(:,1)==0 & gamestatus(:,2)==50); Ind50  = Ind50(~ismember(Ind50,missing));
Ind25  = find(gamestatus(:,1)==0 & gamestatus(:,2)==25); Ind25  = Ind25(~ismember(Ind25,missing));
Ind5   = find(gamestatus(:,1)==0 & gamestatus(:,2)==5);  Ind5   = Ind5(~ismember(Ind5,missing));
Ind    = sort(vertcat(Ind50,Ind25,Ind5));
FTG_5  = find(gamestatus(:,1)==1 & gamestatus(:,3)==5);  FTG_5  = FTG_5(~ismember(FTG_5,missing));
LTG_5  = find(gamestatus(:,1)==2 & gamestatus(:,3)==5);  LTG_5  = LTG_5(~ismember(LTG_5,missing));
MV_5   = find(gamestatus(:,1)==3 & gamestatus(:,3)==5);  MV_5   = MV_5(~ismember(MV_5,missing));
FTG_25 = find(gamestatus(:,1)==1 & gamestatus(:,3)==25); FTG_25 = FTG_25(~ismember(FTG_25,missing));
LTG_25 = find(gamestatus(:,1)==2 & gamestatus(:,3)==25); LTG_25 = LTG_25(~ismember(LTG_25,missing));
MV_25  = find(gamestatus(:,1)==3 & gamestatus(:,3)==25); MV_25  = MV_25(~ismember(MV_25,missing));
Group  = sort(vertcat(FTG_5,LTG_5,MV_5,FTG_25,LTG_25,MV_25));

Ind(:,2:4) = gamestatus(Ind(:,1),:);
Group(:,2:4) = gamestatus(Group(:,1),:);
FTG = Group(Group(:,2) == 1,:);
LTG = Group(Group(:,2) == 2,:);
MV = Group(Group(:,2) == 3,:);

ftg5sc = peoeva(:,FTG_5);
ftg25sc = peoeva(:,FTG_25);
ltg5sc = peoeva(:,LTG_5);
ltg25sc = peoeva(:,LTG_25);
mv5sc = peoeva(:,MV_5);
mv25sc = peoeva(:,MV_25);

% Calculate scores
ftg5sc = calc_scores(ftg5sc,FTG_5,gameinfo);
ftg25sc = calc_scores(ftg25sc,FTG_25,gameinfo);
ltg5sc = calc_scores(ltg5sc,LTG_5,gameinfo);
ltg25sc = calc_scores(ltg25sc,LTG_25,gameinfo);
mv5sc = calc_scores(mv5sc,MV_5,gameinfo);
mv25sc = calc_scores(mv25sc,MV_25,gameinfo);
allscores = calc_scores(peoeva,1:160,gameinfo);
allscores(:,missing) = 0;

% Sum scores
groupsc = horzcat(ftg5sc,ftg25sc,ltg5sc,ltg25sc,mv5sc,mv25sc);
group_tot = sum(groupsc,2);
ftg5_tot = sum(ftg5sc,2);
ftg25_tot = sum(ftg25sc,2);
ltg5_tot = sum(ltg5sc,2);
ltg25_tot = sum(ltg25sc,2);
mv5_tot = sum(mv5sc,2);
mv25_tot = sum(mv25sc,2);

% Sort by individual score
indsc = allscores(:,Ind(:,1));
[~,indSort] = sort(sum(indsc,2));
group_tot_sort = group_tot(indSort);

% Plot scores
figure()
plot(group_tot_sort,'o');
hold on
plot(ftg5_tot(indSort)+100,'o');
plot(ftg25_tot(indSort)+100,'o');
plot(ltg5_tot(indSort)+300,'o');
plot(ltg25_tot(indSort)+250,'o');
plot(mv5_tot(indSort)+80,'o');
plot(mv25_tot(indSort)+150,'o');

%% Calculate group scores and plot, only using group games after trial 44

load experiment
load peoeva


% Find trial numbers
Ind50  = find(gamestatus(:,1)==0 & gamestatus(:,2)==50); Ind50  = Ind50(~ismember(Ind50,missing));
Ind25  = find(gamestatus(:,1)==0 & gamestatus(:,2)==25); Ind25  = Ind25(~ismember(Ind25,missing));
Ind5   = find(gamestatus(:,1)==0 & gamestatus(:,2)==5);  Ind5   = Ind5(~ismember(Ind5,missing));
Ind    = sort(vertcat(Ind50,Ind25,Ind5));
FTG_5  = find(gamestatus(:,1)==1 & gamestatus(:,3)==5);  FTG_5  = FTG_5(~ismember(FTG_5,missing));
LTG_5  = find(gamestatus(:,1)==2 & gamestatus(:,3)==5);  LTG_5  = LTG_5(~ismember(LTG_5,missing));
MV_5   = find(gamestatus(:,1)==3 & gamestatus(:,3)==5);  MV_5   = MV_5(~ismember(MV_5,missing));
FTG_25 = find(gamestatus(:,1)==1 & gamestatus(:,3)==25); FTG_25 = FTG_25(~ismember(FTG_25,missing));
LTG_25 = find(gamestatus(:,1)==2 & gamestatus(:,3)==25); LTG_25 = LTG_25(~ismember(LTG_25,missing));
MV_25  = find(gamestatus(:,1)==3 & gamestatus(:,3)==25); MV_25  = MV_25(~ismember(MV_25,missing));
Group  = sort(vertcat(FTG_5,LTG_5,MV_5,FTG_25,LTG_25,MV_25));

Ind(:,2:4) = gamestatus(Ind(:,1),:);
Group(:,2:4) = gamestatus(Group(:,1),:);
FTG = Group(Group(:,2) == 1,:);
LTG = Group(Group(:,2) == 2,:);
MV = Group(Group(:,2) == 3,:);

FTG_5 = FTG_5(FTG_5 > 44);
FTG_25 = FTG_25(FTG_25 > 44);
LTG_5 = LTG_5(LTG_5 > 44);
LTG_25 = LTG_25(LTG_25 > 44);
MV_5 = MV_5(MV_5 > 44);
MV_25 = MV_25(MV_25 > 44);

ftg5sc = peoeva(:,FTG_5);
ftg25sc = peoeva(:,FTG_25);
ltg5sc = peoeva(:,LTG_5);
ltg25sc = peoeva(:,LTG_25);
mv5sc = peoeva(:,MV_5);
mv25sc = peoeva(:,MV_25);

% Calculate scores
ftg5sc = calc_scores(ftg5sc,FTG_5,gameinfo);
ftg25sc = calc_scores(ftg25sc,FTG_25,gameinfo);
ltg5sc = calc_scores(ltg5sc,LTG_5,gameinfo);
ltg25sc = calc_scores(ltg25sc,LTG_25,gameinfo);
mv5sc = calc_scores(mv5sc,MV_5,gameinfo);
mv25sc = calc_scores(mv25sc,MV_25,gameinfo);
allscores = calc_scores(peoeva,1:160,gameinfo);
allscores(:,missing) = 0;

% Sum scores
groupsc = horzcat(ftg5sc,ftg25sc,ltg5sc,ltg25sc,mv5sc,mv25sc);
group_tot = sum(groupsc,2);
ftg5_tot = sum(ftg5sc,2);
ftg25_tot = sum(ftg25sc,2);
ltg5_tot = sum(ltg5sc,2);
ltg25_tot = sum(ltg25sc,2);
mv5_tot = sum(mv5sc,2);
mv25_tot = sum(mv25sc,2);

% Sort by individual score
indsc = allscores(:,Ind(:,1));
[~,indSort] = sort(sum(indsc,2));
group_tot_sort = group_tot(indSort);

% Plot scores
figure('position',[0 0 700 525])
plot(group_tot_sort,'o','MarkerSize',6);
hold on
blank = plot(1:50,1:50,'o','color',[1 1 1]);
set(blank,'visible','off')
ax = gca; ax.ColorOrderIndex = 2; 
plot(ftg5_tot(indSort)+230,'o','MarkerSize',6);
h = plot(ftg25_tot(indSort)+230,'o','MarkerSize',6);
hcolor = get(h,'color');
plot(ltg5_tot(indSort)+230,'o','MarkerSize',6);
plot(ltg25_tot(indSort)+230,'o','MarkerSize',6);
plot(mv5_tot(indSort)+230,'o','MarkerSize',6);
plot(mv25_tot(indSort)+230,'o','MarkerSize',6);
set(gca,'ticklength',[0.005 0.005])
overlaps = find(ftg25_tot(indSort)==38); % some of the points coincide, make them more obvious
plot(overlaps,38+230,'o','MarkerSize',4,'color',hcolor)

ylim([250 440])
set(gca,'ytick',[250:20:330 370:20:440])
set(gca,'yticklabel',{'20','40','60','80','100','370','390','410','430'})
set(gca,'fontsize',16)
xlabel('Participant (ordered by individual score)','fontsize',18)
ylabel('Group score','fontsize',18)

columnlegend(4,{'All groups','','FTG-5','FTG-25','LTG-5','LTG-25','MV-5','MV-25'},'fontsize',16,'location','best','boxon')

% savefig('groupscores')
