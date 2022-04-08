%%% 20211026 plot eye trace of one trial, also plot the eye fixation result
%%% of all valid 38 subjects

%% ===  plot xy eye trace during reading
load('Z:\Lexical\Analyse_data\20190608_b54e\Targ60_TW1000\epoch_WrdOff.mat')
trl_id = 3;
tpw = epochdata.trialinfo(trl_id,7); %trig_pretarg_wrdoff (meg)
meg_range = tpw-500:tpw+500;
eye_range = meg_range+(Event.SentenceOnTrigger(trl_id,2)-Event.SentenceOnTrigger(trl_id,3));
% just get xy from .asci
eyefile = 'Z:\Lexical\RawData\EyeLink_data\b54e.asc';
fid = fopen(eyefile);
% reading file line by line -- loop
BlockEnd = 0;
ContinuousEyeData = [];
ii =  1;
while ~BlockEnd && ii<100000
    strline = fgetl(fid);
    isevent = 0;
    %%% if this line is continuous xy data
    content = sscanf(strline, '%f %f %f %f %f');
    if size(content,1) == 5
        ContinuousEyeData = [ContinuousEyeData; content(1:3)'];
    end
    
    if contains(strline, 'end of block')
        BlockEnd = 1;
    end
    ii = ii+1;
end
fclose(fid);
eyeid = find(ContinuousEyeData(:,1)==eye_range(1)):find(ContinuousEyeData(:,1)==eye_range(end));
xy = ContinuousEyeData(eyeid,2:3);
sacdur = 37; %unit in timepoint
%%% plot eye trace from difference sources
%%%% pretarget fixation:[1809604 1809808]
%%%% saccade to target: [1809809 1809845]
%%%% target fixation:   [1809846 1810060]
%%%% trig_pretarg_wrdoff:time 0, datapoint 500
%%%% trig_pretarg_wrdoff:time 0.37, datapoint 537
%% plot
% converted xy into visual angels aligned with the coords of target word
xy(:,1) = xy(:,1)-xy(500,1);
xy(:,2) = xy(:,2)-xy(500,2);
xy = xy./34.59;% convert pixels into visual degrees
%%
figtitle = 'Eye Trace';
h = figure('Name',figtitle,'color',[1 1 1],'position',[100 100 420 220]);
% x
subplot(2,1,1);
plot(xy(:,1),'k','LineWidth',1)
set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
set(gca,'box','off','LineWidth',1)
xlim([0 1000]);
set(gca,'XTick',0:100:1000);
set(gca,'XTickLabel',{'-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.2',...
    '0.3','0.4','0.5'},'FontWeight','normal','FontSize',7,'FontName','Arial');
xlabel('Time (s)','FontSize',7,'FontWeight','normal','FontName','Arial');
hold on;
yrange = get(gca,'YLim');
plot([500 500],yrange,':k','LineWidth',1)
plot([500+sacdur 500+sacdur],yrange,':k','LineWidth',1)
ylabel('X coords distance from target centre (visual degree)','FontSize',7,'FontWeight','normal','FontName','Arial');
% y
subplot(2,1,2);
plot(xy(:,2),'k','LineWidth',1)
set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
set(gca,'box','off','LineWidth',1)
xlim([0 1000])
set(gca,'XTick',0:100:1000);
set(gca,'XTickLabel',{'-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.2',...
    '0.3','0.4','0.5'},'FontWeight','normal','FontSize',7,'FontName','Arial');
xlabel('Time (s)','FontSize',7,'FontWeight','normal','FontName','Arial');
hold on;
ylim([-2 2])
yrange = get(gca,'YLim');
plot([500 500],yrange,':k','LineWidth',1)
plot([500+sacdur 500+sacdur],yrange,':k','LineWidth',1)
ylabel('Y coords distance from target centre (visual degree)','FontSize',7,'FontWeight','normal','FontName','Arial');
set(gcf, 'renderer', 'painters')
PPath.FigPath = 'Z:\Saccade_AlphaPhase\figs\';
saveas(h,[PPath.FigPath figtitle]);
saveas(h,[PPath.FigPath figtitle],'svg');


%% plot eye movement data
%%%% first fixation duration and gaze duration for pre-target and target in
%%%% low and high condition
load('Z:\Lexical\Results\Behavioral\CmbDatasets\CmbDatasets_BehaData'); %% load in 'BehaData'
% remove eye data that are from invalid subs
load('Z:\Saccade_AlphaPhase\Results\AllSubs'); %all valid subs
load('Z:\Lexical\Analyse_data\SubInfo.mat')
% [SubInfo.subjects.Targ60 subs]; %manually find the missing/invalid subid
s = 28; %id of the bad sub in total 39 sub index

%% ----- violin plot
addpath(genpath('Z:\Saccade_AlphaPhase\Analyse_codes\Violinplot\'));
colmat = [0 114 189;217 83 25]./255;
figtitle = 'Fixation Duration';
h = figure('Name',figtitle,'color',[1 1 1],'position',[300 300 320 230]);
tmp = BehaData.FirstFix;
%%% remove subjects without 2 datasets
tmp.Pre(tmp.Pre(:,1)==0,:) = [];
tmp.Tag(tmp.Tag(:,1)==0,:) = [];
%%% remove the invalid sub 28
tmp.Pre(s,:) = [];
tmp.Tag(s,:) = [];
nsub = length(tmp.Tag);
group = [cellstr(repmat('PT-l',nsub,1));cellstr(repmat('PT-h',nsub,1));cellstr(repmat('T-l',nsub,1));cellstr(repmat('T-h',nsub,1))];
grouporder={'PT-l','PT-h','T-l','T-h'};
vdata = [tmp.Pre(:,1);tmp.Pre(:,2);tmp.Tag(:,1);tmp.Tag(:,2)];
vp = violinplot(vdata, group,'GroupOrder',grouporder);
vp(1).ViolinColor = colmat(1,:);
vp(2).ViolinColor = colmat(2,:);
vp(3).ViolinColor = colmat(1,:);
vp(4).ViolinColor = colmat(2,:);
vp(1).ShowMean = 1; vp(2).ShowMean = 1;
vp(3).ShowMean = 1; vp(4).ShowMean = 1;
vp(1,1).MeanPlot.LineWidth = 2;
vp(1,2).MeanPlot.LineWidth = 2;
vp(1,3).MeanPlot.LineWidth = 2;
vp(1,4).MeanPlot.LineWidth = 2;
vp(1,1).BoxWidth = 0.01; vp(1,2).BoxWidth = 0.01;
vp(1,3).BoxWidth = 0.01; vp(1,4).BoxWidth = 0.01;

%%% ====== set labels and other stuff
set(gca,'LineWidth',1,'FontSize',7,'FontWeight','normal','FontName','Arial')
subtitle = 'First fixation duration(ms)';
ylabel(subtitle,'FontSize',7,'FontWeight','normal','FontName','Arial')

%%% save out
set(gcf, 'renderer', 'painters')
PPath.FigPath = 'Z:\Saccade_AlphaPhase\figs\';
saveas(h,[PPath.FigPath figtitle]);
saveas(h,[PPath.FigPath figtitle],'svg');

%%% statistics
[~,p,~,stat] = ttest(tmp.Pre(:,1),tmp.Pre(:,2));
stat.pvalue = p; %[t=0.054, p=0.958, df = 37]
Pre_stat = stat;
[~,p,~,stat] = ttest(tmp.Tag(:,1),tmp.Tag(:,2));
stat.pvalue = p; %[t=7.099, p=2.095e-8, df = 37]
Tag_stat = stat;

%%% saving out
BehaData = [];
BehaData.subs = SubInfo.subjects.Targ60([1:s-1 s+1:end]);
BehaData.FirstFix.Pre = tmp.Pre;
BehaData.FirstFix.Tag = tmp.Tag;
BehaData.FirstFix.Pre_stat = Pre_stat;
BehaData.FirstFix.Tag_stat = Tag_stat;
save('Z:\Saccade_AlphaPhase\Results\BehaData.mat','BehaData')

%%% copy the data into first_fixation.csv then run paired_Ttest.R
n = length(BehaData.subs);
subid = kron(ones(4,1),[1:n]');
data = [BehaData.FirstFix.Pre(:,1);BehaData.FirstFix.Pre(:,2);...
    BehaData.FirstFix.Tag(:,1);BehaData.FirstFix.Tag(:,2)];
frqid = [ones(n,1);2*ones(n,1);ones(n,1);2*ones(n,1)];
locid = [ones(2*n,1);2*ones(2*n,1)];

[subid data frqid locid];




%% ----- bar plot
addpath(genpath('Z:\Saccade_AlphaPhase\Analyse_codes\plot_funs\'));
load('Z:\Lexical\Results\Behavioral\CmbDatasets\CmbDatasets_BehaData'); %% load in 'BehaData'
load('Z:\Saccade_AlphaPhase\Results\AllSubs'); %all valid subs
load('Z:\Lexical\Analyse_data\SubInfo.mat')
s = 28; %id of the bad sub in total 39 sub index
tmp = BehaData.FirstFix;
tmp.Pre(tmp.Pre(:,1)==0,:) = [];
tmp.Tag(tmp.Tag(:,1)==0,:) = [];
tmp.Pre(s,:) = [];
tmp.Tag(s,:) = [];
y = [mean(tmp.Pre); mean(tmp.Tag)];         % random y values (2 groups of 2 parameters)
errY = [std(tmp.Pre); std(tmp.Tag)]./sqrt(size(tmp.Pre,1));
% plot
figtitle = 'Fixation Duration';
h = figure('Name',figtitle,'color',[1 1 1],'position',[300 300 300 210]);
a = barwitherr(cat(3,-errY,errY), y);    % Plot with errorbars
% set labels and other stuff
set(gca,'box','off','LineWidth',1)
set(gca,'LineWidth',1,'FontSize',7,'FontWeight','normal','FontName','Arial')
subtitle = 'First fixation duration(ms)';
ylabel(subtitle,'FontSize',7,'FontWeight','normal','FontName','Arial')
set(gca,'XTickLabel',{'Pre-target','Target'},'FontSize',7,'FontWeight',...
    'normal','FontName','Arial')
legendflex(a,[{'Low frequency target'},{'High frequency target'}],'anchor',...
    {'ne','ne'}, 'buffer', [5 -5],'FontSize',7,'FontWeight','normal',...
    'FontName','Arial','xscale',0.8,'box', 'off');
ylim([180 260])
%%% save out
set(gcf, 'renderer', 'painters')
PPath.FigPath = 'Z:\Saccade_AlphaPhase\figs\';
saveas(h,[PPath.FigPath figtitle],'svg');




























