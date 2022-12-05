%% 20210825 clean up scripts 
% run ITC for all 39 subs, pre-target and target words, aligned with saccade_on (WrdOff) and fixation_on (WrdOn)
% the range of fixation durations for epochs is [0.08 1]s defined in the pre-processing

function ITC_WrdFrq(sid)
tic
server = 1;
%%%% running setting
runconds = {'WrdOff','WrdOn'};
EpochType = {'PreTarg','Targ'};
targid = [-1, 0];
conds = {'low','high','diff'};
%%% time window parameters
timwin = [-0.5 0.5]; %%% select epochs aligned with WrdOff and WrdOn
stat_win = [-0.2 0];
%%% parameter for hilbert filter
freqrange = 4:1:30;
freqwidth_factor = 0.3; % freq widh: ±0.3*fcenter(Hz)

%%% basic settingup
if server
    rootdir = '/rds/projects/2018/jenseno-reading/';
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
    ft_defaults
    addpath(genpath('/rds/projects/2018/jenseno-reading/Analyse_codes/helper_functions/'))
else
    %%% local
    rootdir = 'Z:\';
    addpath(genpath('Z:\Saccade_AlphaPhase\codes'));
    addpath Z:\fieldtrip-20200220\
    ft_defaults
    %%% for colormap
    cmap = colormap(cbrewer('div','RdBu',32));
    cmap = colormap(flipud(cmap));
end

%%% paths
PPath.Data = [rootdir 'Lexical' filesep 'Analyse_data' filesep];%data are under the folder of lexical 
PPath.FigPath = [rootdir 'Saccade_AlphaPhase' filesep 'Results' filesep 'ITC' filesep];
%%% get the subinformation
load([PPath.Data 'SubInfo.mat']);

if server == 1    
    %%%%%%==================== epoch alignment loop
    for rrr = 1:length(runconds)
        RunCond = runconds{rrr};
        
        %%%%%%===================== subject loop
        for sss = sid
            fprintf(['********** s' num2str(sid) ' ******** \n\n']);
            
            %%% set output
            ITC = [];
            ITC.conds = conds;
            ITC.sub = SubInfo.subjects.Targ60{sss};
            
            %%%%%%%%%%========= getting the combined dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%=== JEP 60
            PPath.SaveData = [PPath.Data SubInfo.subjects.JEP60{sss} filesep 'JEP60_TW1000' filesep];
            load([PPath.SaveData 'epoch_' RunCond]);
            eval(['epochdata_jep = epoch_' RunCond ';']);
            %%%=== targ 60
            PPath.SaveData = [PPath.Data SubInfo.subjects.Targ60{sss} filesep 'Targ60_TW1000' filesep];
            load([PPath.SaveData 'epoch_' RunCond]);
            eval(['epochdata = epoch_' RunCond ';']);
            %%%=== combining 2 datasets together
            epochdata.trialinfo = [epochdata.trialinfo; epochdata_jep.trialinfo];
            epochdata.trial = [epochdata.trial epochdata_jep.trial];
            epochdata.time = [epochdata.time epochdata_jep.time];
            clear epochdata_jep
            clear epoch_*
            
            %%% get first fixation trials
            validduration = find(epochdata.trialinfo(:,10) == 1); % FirstPassFix
            cfg = [];
            cfg.trials = validduration;
            cfg.channel = 'MEGGRAD';
            epochdata = ft_selectdata(cfg, epochdata);
            
            for mmm = 1:length(EpochType)
                %%% get equal trial number across all conditions
                %%%% seperate trials based on target-freq
                idx_low = find(epochdata.trialinfo(:,3) == targid(mmm) & epochdata.trialinfo(:,12) == 1); %%[loc2targ SentenceCondition]
                idx_high = find(epochdata.trialinfo(:,3) == targid(mmm) & epochdata.trialinfo(:,12) == 2);
                %%% make sure same number of trials for different conditions
                nlow = length(idx_low);
                nhigh = length(idx_high);
                if nlow > nhigh
                    idx_low = idx_low(randperm(nlow));
                    idx_low = idx_low(1:nhigh);
                elseif nlow < nhigh
                    idx_high= idx_high(randperm(nhigh));
                    idx_high = idx_high(1:nlow);
                end
                cfg         = [];
                cfg.latency = timwin;
                cfg.trials  = idx_low;
                data_low    = ft_selectdata(cfg, epochdata);
                cfg.trials  = idx_high;
                data_high   = ft_selectdata(cfg, epochdata);
                
                %%% get trl numbers in each condition
                for cid = 1:length(conds)-1
                    eval(['ITC.' EpochType{mmm} '_TrlInfo{1,cid} = data_' conds{cid} '.trialinfo;']);
                    eval(['ITC.' EpochType{mmm} '_RT(1,cid) = nanmean(data_' conds{cid} '.trialinfo(:,8));']);
                    eval(['ITC.' EpochType{mmm} '_TrlIdx_equ{1,cid} = idx_' conds{cid} ';']);
                end
                
                %%%%%%%%%%%=================== Run the analysis===============%%%%%%%%%%%%
                for ccc = 1:length(conds)-1
                    %%% setting outputs
                    itc        = [];
                    itc.label  = data_low.label;
                    itc.freq   = freqrange;
                    itc.time   = data_low.time{1,1};
                    itc.dimord = 'chan_freq_time';
                    %%% get phases from hilbert
                    nfreq = length(freqrange);
                    ntrl = length(data_low.trial);
                    freq_hilbert = zeros(ntrl,length(itc.label),nfreq,length(itc.time));
                    for ifreq = 1:nfreq
                        fc = freqrange(ifreq);
                        hw = fc*freqwidth_factor; %freq_halfwidth
                        %%%% get hilbert filert data
                        cfg = [];
                        cfg.bpfilter   = 'yes';
                        cfg.bpfreq     = [fc-hw fc+hw];
                        cfg.hilbert    = 'complex';
                        cfg.keeptrials = 'yes';
                        eval(['fltdata = ft_preprocessing(cfg,data_' conds{ccc} ');']); %% angles/phases
                        %%% get hilbert analytic signal from each trial
                        for ttt = 1:ntrl
                            freq_hilbert(ttt,:,ifreq,:) = fltdata.trial{1,ttt};
                        end
                    end
                    clear fltdata 
                    %%%===== calculating ITC in the fieldtrip method 
                    % compute inter-trial phase coherence (itpc)
                    itc.itpc      = freq_hilbert./abs(freq_hilbert);         % divide by amplitude
                    itc.itpc      = sum(itc.itpc,1);   % sum angles
                    itc.itpc      = abs(itc.itpc)/ntrl;   % take the absolute value and normalize
                    itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
                    itc.powspctrm = itc.itpc;
                    itc = rmfield(itc,'itpc');
                    %%% combine planar sensors(as if it was a powspctrm)
                    itc = ft_combineplanar([],itc);
                    %%% saving out
                    eval(['ITC.' EpochType{mmm} '_' conds{ccc} '= itc;']);
                end
                clear data_*
            end
            save([PPath.FigPath RunCond '_ITC_' num2str(sid)],'ITC','-v7.3');
            toc
            fprintf(['********** s' num2str(sid) ': done ******** \n\n']);
        end
    end
    close all
else
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run locally get the whole structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for rrr = 1:length(runconds)
        RunCond = runconds{rrr};
        %%% check!subjects have more than 50 trials in either cond
        %%% remove subs that has too few trials (sub28 only has 40 trials
        %%% per condition)
        load([PPath.FigPath RunCond '_ITC_1']); % ITC
        ITC_all = [];
        ss = 0;  
        %%% only run group analysis with 39 subjects who has both datasets
        for ps = 1:39 
            fprintf(['*** s' num2str(ps) ' *** \n']);
            load([PPath.FigPath RunCond '_ITC_' num2str(ps)]);
            nntrl = size(ITC.PreTarg_TrlInfo{1, 1},1);
            if nntrl >= 50 % nntrl trls for each epoch under each cond
                ss = ss + 1;
                ITC_all.subs{ss,1} = ITC.sub;
                for fd = 1:length(EpochType)
                    eval(['ITC_all.' EpochType{fd} '_TrlIdx_equ(ss,:) = ITC.' EpochType{fd} '_TrlIdx_equ;']);
                    eval(['ITC_all.' EpochType{fd} '_RT(ss,:) = ITC.' EpochType{fd} '_RT;']);
                    eval(['ITC_all.' EpochType{fd} '_TrlInfo(ss,:) = ITC.' EpochType{fd} '_TrlInfo;']);
                    eval(['ITC_all.' EpochType{fd} '_low{ss} = ITC.' EpochType{fd} '_low;']);
                    eval(['ITC_all.' EpochType{fd} '_high{ss} = ITC.' EpochType{fd} '_high;']);
                end
            else
                fprintf(['*** s ' num2str(ps) ': only ' num2str(nntrl) ' trls **** \n\n']);
            end
        end
        %%% in total, 38 subs
        clear ITC
        ITC = ITC_all;
        clear ITC_all
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Group statistics-----ITC %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for mmm = 1:length(EpochType)
            %%%%% pre-define the output
            eval(['Data_low = ITC.' EpochType{mmm} '_low;']);
            eval(['Data_high = ITC.' EpochType{mmm} '_high;']);
               
            %%%%%%%%% permutation statistics %%%%%%%%%
            load([rootdir filesep 'neuromag306cmb_neighb.mat'])           
            cfg                  = [];
            cfg.method           = 'montecarlo';
            cfg.statistic        = 'ft_statfun_depsamplesT';
            cfg.parameter        = 'powspctrm';
            cfg.latency          = stat_win;
            cfg.frequency        = [9 13];
            cfg.channel          = 'all';
            cfg.avgoverfreq      = 'yes';
            cfg.avgovertime      = 'yes';
            cfg.avgoverchan      = 'no';
            cfg.correctm         = 'cluster';
            cfg.clusterstatistic = 'maxsum';
            cfg.clusteralpha     = 0.05;
            cfg.tail             = 0;
            cfg.clustertail      = 0;
            cfg.alpha            = 0.05;
            cfg.correcttail      = 'prob';
            cfg.minnbchan        = 1;
            cfg.neighbours       = neighbours;
            cfg.numrandomization = 5000;
            nsub = length(Data_low);
            design = [1:nsub 1:nsub; ones(1,nsub) ones(1,nsub)*2];
            cfg.design           = design;
            cfg.uvar             = 1; % row in design indicating subjects, repeated measure
            cfg.ivar             = 2; % row in design indicating condition for contrast
            stat = ft_freqstatistics(cfg, Data_low{:}, Data_high{:});
            any(any(any(stat.mask)))
            %%% save stat to the ITC struct
            eval(['ITC.' EpochType{mmm} '_stat = stat']);
            
            %%%============== plot topography figure
            if any(stat.mask)
                %%% get the ITC_diff for the group
                grandavg_low = ft_freqgrandaverage(cfg, Data_low{:});
                grandavg_high = ft_freqgrandaverage(cfg, Data_high{:});
                cfg           = [];
                cfg.parameter = 'powspctrm';
                cfg.operation = 'subtract';
                itc_diff      = ft_math(cfg,grandavg_low, grandavg_high);
                
                %%% topography for the sig sensors
                cfg                  = [];
                % only positive cluster is sig
                cfg.highlightchannel = stat.label(stat.posclusterslabelmat == 1);
                cfg.parameter        = 'powspctrm';
                cfg.layout           = 'neuromag306cmb.lay';
                cfg.maskstyle        = 'saturation';
                cfg.highlight        = 'on'; %'numbers';
                cfg.contournum       = 0;
                cfg.colormap         = cmap;
                cfg.zlim             = [-.03 .03];
                cfg.xlim             = stat_win;
                cfg.ylim             = [9 13];
                cfg.gridscale        = 360;
                cfg.marker           = 'no'; %'numbers';
                cfg.comment          = 'no';
                cfg.style            = 'straight';
                figname              = [RunCond '_Grp_ITC_diff_' EpochType{mmm} '_sig_1nb'];
                h                    = figure('name',figname,'color',[1 1 1],'position',[0 0 440 300]);
                ft_topoplotTFR(cfg,itc_diff);
                title(['\Delta ITC aligned to saccade on to target'],'FontWeight','normal','FontSize',7,'FontName','Arial');
                set(gcf, 'renderer', 'painters')
                saveas(h,[PPath.FigPath figname]);
                saveas(h,[PPath.FigPath figname],'svg');
            end
        end
        save([PPath.FigPath RunCond '_ITC'],'ITC','-v7.3');
        % rename
        eval([RunCond '_ITC = ITC']);
        clear ITC
    end
end

%% save out the sig sens of WrdOff_PreTarg
SigSens = [];
SigSens.label = ITC.PreTarg_stat.label(ITC.PreTarg_stat.mask);
SigSens.mask = ITC.PreTarg_stat.mask;
save([PPath.FigPath 'SigSens_WrdOffPreTarg'],'SigSens');

%% get averaged ITC-diff over WrdOff_PreTarg sigsens for WrdOff_PreTarg and WrdOn_Targ
runconds = {'WrdOff','WrdOn'};
EpochType = {'PreTarg','Targ'};
plot_0ms = {'saccade on','fixation on'};
xline = -0.2;
for m = 1:length(EpochType)
    % get itc-diff of all sens
    eval(['Data_low =' runconds{m} '_ITC.' EpochType{m} '_low;']);
    eval(['Data_high =' runconds{m} '_ITC.' EpochType{m} '_high;']);
    grandavg_low = ft_freqgrandaverage(cfg, Data_low{:});
    grandavg_high = ft_freqgrandaverage(cfg, Data_high{:});
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'subtract';
    itc_diff      = ft_math(cfg,grandavg_low, grandavg_high);
    % averaged itc-diff over sig sensors
    itc_avg = itc_diff.powspctrm(WrdOff_ITC.PreTarg_stat.mask,:,:);
    itc_avg = squeeze(mean(itc_avg,1));
    % plot
    figname = ['ITCdiff_' runconds{m} '_' EpochType{m}];
    h = figure('Name',figname,'color',[1 1 1],'position',[0 0 220 150]);
    pcolor(itc_diff.time, itc_diff.freq, itc_avg);
    shading interp;
    %%% for colormap
    cmap = colormap(cbrewer('div','RdBu',32));
    cmap = colormap(flipud(cmap));
    colorbar('FontWeight','normal','FontSize',7,'FontName','Arial')
    caxis([-0.03 0.03])
    xlim([-0.4 0.4])
    % title
    title({['Averaged \Delta ITC'],[plot_0ms{m} ' to target']},'FontWeight','normal','FontSize',7,'FontName','Arial');
    % plot some lines
    hold on;
    plot([-0.5 0.5],[9 9],':k','LineWidth',1)
    plot([-0.5 0.5],[13 13],':k','LineWidth',1)
    plot([xline xline],[4 30],':k','LineWidth',1)
    plot([0 0],[4 30],':k','LineWidth',1)
    set(gca,'XTick',-0.4:0.2:0.4);
    set(gca,'XTickLabel',{'-0.4','-0.2',plot_0ms{m},'0.2','0.4'},'FontWeight','normal','FontSize',7,'FontName','Arial')
    set(gca,'YTick',5:5:30);
    set(gca,'YTickLabel',{'5','','','20','25','30'},'FontWeight','normal','FontSize',7,'FontName','Arial')
    text([-0.45 -0.47],[9 13],[{'9'},{'13'}],'FontWeight','normal','FontSize',7,'FontName','Arial')
    ylabel('Frequency (Hz)','FontWeight','normal','FontSize',7,'FontName','Arial')
    xlabel('Time (s)','FontWeight','normal','FontSize',7,'FontName','Arial')
    text(0.8,6,{'\Delta Inter-trial coherence (r^2)','(low-high lexical freq target)'},'FontWeight','normal','FontSize',7,'FontName','Arial','Rotation',90)
    text(-0.18,23,{['pre-target'],['interval']},'FontWeight','normal','FontSize',7,'FontName','Arial')
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figname]);
    saveas(h,[PPath.FigPath figname],'svg');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%====Supplementary Fig S3: plot the raw PLI for low/high condition ===%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figname = 'Raw PLI_pretarget_high' ;
% itc_avg = grandavg_high.powspctrm(SigSens.mask,:,:);
% itc_avg = squeeze(mean(itc_avg,1));
figname = 'Raw PLI_pretarget_low' ;
itc_avg = grandavg_low.powspctrm(SigSens.mask,:,:);
itc_avg = squeeze(mean(itc_avg,1));
% plot
h = figure('Name',figname,'color',[1 1 1]);
pcolor(grandavg_high.time, grandavg_high.freq, itc_avg);
shading interp;
cmap = colormap(cbrewer('div','RdBu',32));
cmap = colormap(flipud(cmap));
%%% for colormap
colorbar('FontWeight','normal','FontSize',12,'FontName','Arial')
caxis([0.1 0.4])
xlim([-0.4 0.4])
% title
title('Averaged raw PLI, saccade on to high freq target','FontWeight','normal','FontSize',14,'FontName','Arial');
% plot some lines
hold on;
plot([-0.5 0.5],[9 9],':k','LineWidth',1)
plot([-0.5 0.5],[13 13],':k','LineWidth',1)
plot([xline xline],[4 30],':k','LineWidth',1)
plot([0 0],[4 30],':k','LineWidth',1)
set(gca,'XTick',-0.4:0.2:0.4);
set(gca,'XTickLabel',{'-0.4','-0.2',plot_0ms{m},'0.2','0.4'},'FontWeight','normal','FontSize',12,'FontName','Arial')
set(gca,'YTick',5:5:30);
set(gca,'YTickLabel',{'5','','','20','25','30'},'FontWeight','normal','FontSize',12,'FontName','Arial')
text([-0.43 -0.445],[9 13],[{'9'},{'13'}],'FontWeight','normal','FontSize',12,'FontName','Arial')
ylabel('Frequency (Hz)','FontWeight','normal','FontSize',12,'FontName','Arial')
xlabel('Time (s)','FontWeight','normal','FontSize',12,'FontName','Arial')
saveas(h,[PPath.FigPath figname]);
% zoom in
figname = [figname '_zoomin'];
h = figure('Name',figname,'color',[1 1 1]);
pcolor(grandavg_high.time, grandavg_high.freq, itc_avg);
shading interp;
cmap = colormap(cbrewer('div','RdBu',32));
cmap = colormap(flipud(cmap));
%%% for colormap
colorbar('FontWeight','normal','FontSize',12,'FontName','Arial')
xlim([-0.2 0])
caxis([0.15 0.2])
ylim([8 30])
set(gca,'XTick',[-0.2 0]);
set(gca,'XTickLabel',{'-0.2','saccade on'},'FontWeight','normal','FontSize',12,'FontName','Arial')
set(gca,'YTick',5:5:30);
set(gca,'YTickLabel',{'','','','20','25','30'},'FontWeight','normal','FontSize',12,'FontName','Arial')
text([-0.21 -0.21 -0.214],[8 9 13],[{'8'},{'9'},{'13'}],'FontWeight','normal','FontSize',12,'FontName','Arial')
ylabel('Frequency (Hz)','FontWeight','normal','FontSize',12,'FontName','Arial')
xlabel('Time (s)','FontWeight','normal','FontSize',12,'FontName','Arial')
hold on
plot([-0.2 0],[9 9],':k','LineWidth',1)
plot([-0.2 0],[13 13],':k','LineWidth',1)
saveas(h,[PPath.FigPath figname]);











