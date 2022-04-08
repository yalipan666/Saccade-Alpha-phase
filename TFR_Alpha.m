%%%% 20210830 run TFR (2-30Hz) for all sensors using Hilbert
% to show raw alpha power (without BL correction) over all sensors during
% pre-target fixations, aligned with fix_on(wrd_on)
% get rid of sub28 due to few trials

function TFR_Alpha(sid)
tic
server = 1;
%%%% running setting
runconds = {'WrdOff','WrdOn'};
EpochType = {'PreTarg','Targ'};
targid = [-1, 0];
minfix = 0; %% minimum fixation duration (80ms threshold in the pre-processing)
conds = {'low','high','diff'};
%%% time window parameters
timwin = [-0.5 0.5]; %%% {[-0.5 0];[0 0.5]}; select epochs aligned with WrdOff and WrdOn
stat_win = [-0.2 0];
%%% parameter for hilbert filter
freqrange = 4:1:30;
freqwidth_factor = 0.3; % freq widh: ±0.3*fcenter(Hz)

%%% basic settingup
if server
    rootdir = '/rds/projects/2018/jenseno-reading/';
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
    ft_defaults
else
    %%% local
    rootdir = 'Z:\';
    addpath(genpath('Z:\Saccade_AlphaPhase\Analyse_codes'));
    addpath Z:\fieldtrip-20200220\
    ft_defaults
    %%% for colormap
    cmap = colormap(cbrewer('div','RdBu',32));
    cmap = colormap(flipud(cmap));
    close
end

%%% paths
PPath.Data = [rootdir 'Lexical' filesep 'Analyse_data' filesep];%data are under the folder of lexical
PPath.FigPath = [rootdir 'Saccade_AlphaPhase' filesep 'Results' filesep 'TFR' filesep];
%%% get the subinformation
load([PPath.Data 'SubInfo.mat']);

if server == 1
    %%%%%%==================== epoch alignment loop
    for rrr = 1:length(runconds)
        RunCond = runconds{rrr};
        
        %%%%%%===================== subject loop
        for sss = sid
            fprintf(['********** s' num2str(sss) ' ******** \n\n']);
            
            %%% set output
            TFR = [];
            TFR.conds = conds;
            TFR.sub = SubInfo.subjects.Targ60{sss};
            
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
                cfg         = [];
                cfg.latency = timwin;
                cfg.trials  = idx_low;
                data_low    = ft_selectdata(cfg, epochdata);
                cfg.trials  = idx_high;
                data_high   = ft_selectdata(cfg, epochdata);
                
                %%% get trl numbers in each condition
                for cid = 1:length(conds)-1
                    eval(['TFR.' EpochType{mmm} '_TrlInfo{1,cid} = data_' conds{cid} '.trialinfo;']);
                    eval(['TFR.' EpochType{mmm} '_RT(1,cid) = nanmean(data_' conds{cid} '.trialinfo(:,8));']);
                    eval(['TFR.' EpochType{mmm} '_TrlIdx_equ{1,cid} = idx_' conds{cid} ';']);
                end
                
                %%%%%%%%%%%=================== Run the analysis===============%%%%%%%%%%%%
                for ccc = 1:length(conds)-1
                    %%% setting outputs
                    pow        = [];
                    pow.label  = data_low.label;
                    pow.freq   = freqrange;
                    pow.time   = data_low.time{1,1};
                    pow.dimord = 'chan_freq_time';
                    %%% get power from hilbert
                    nfreq = length(freqrange);
                    eval(['ntrl = length(data_' conds{ccc} '.trial);']);
                    pow.powspctrm  = zeros(length(pow.label),nfreq,length(pow.time));
                    for ifreq = 1:nfreq
                        %%% pre-allocate output
                        tempow = zeros(ntrl,length(pow.label),length(pow.time));
                        %%% setting the filter parameters
                        fc = freqrange(ifreq);
                        hw = fc*freqwidth_factor; %freq_halfwidth
                        %%%% get amplitude
                        cfg = [];
                        cfg.bpfilter   = 'yes';
                        cfg.bpfreq     = [fc-hw fc+hw];
                        cfg.hilbert    = 'abs';
                        cfg.keeptrials = 'yes';
                        eval(['fltdata = ft_preprocessing(cfg,data_' conds{ccc} ');']); %% amplitude
                        %%% get power from each trial
                        for ttt = 1:ntrl
                            tempow(ttt,:,:) = fltdata.trial{1,ttt}.^2;
                        end
                        %%% get the avergaed power across trls
                        pow.powspctrm(:,ifreq,:) = mean(tempow,1);
                    end
                    clear fltdata
                    %%% combine planar sensors(as if it was a powspctrm)
                    pow = ft_combineplanar([],pow);
                    %%% saving out
                    eval(['TFR.' EpochType{mmm} '_' conds{ccc} '= pow;']);
                end
                clear data_*
            end
            save([PPath.FigPath RunCond '_TFR_' num2str(sid)],'TFR','-v7.3');
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
        %%% check!subjects have more than 30 trials in either cond
        %%% remove subs that has too few trials (sub28 only has 40 trials
        %%% per condition)
        %%% get the group data cell
        load([PPath.FigPath RunCond '_TFR_1']); % TFR
        ss = 0;
        for ps = 1:39
            fprintf(['*** s' num2str(ps) ' *** \n']);
            load([PPath.FigPath RunCond '_TFR_' num2str(ps)]);
            nntrl = size(TFR.PreTarg_TrlInfo{1, 1},1);
            if nntrl >= 50 % 50 nntrl for each epoch under each cond
                ss = ss + 1;
                TFR_all.subs{ss,1} = TFR.sub;
                for fd = 1:length(EpochType)
                    eval(['TFR_all.' EpochType{fd} '_TrlIdx_noequ(ss,:) = TFR.' EpochType{fd} '_TrlIdx_equ;']);
                    eval(['TFR_all.' EpochType{fd} '_RT(ss,:) = TFR.' EpochType{fd} '_RT;']);
                    eval(['TFR_all.' EpochType{fd} '_TrlInfo(ss,:) = TFR.' EpochType{fd} '_TrlInfo;']);
                    eval(['TFR_all.' EpochType{fd} '_low{ss} = TFR.' EpochType{fd} '_low;']);
                    eval(['TFR_all.' EpochType{fd} '_high{ss} = TFR.' EpochType{fd} '_high;']);
                end
            else
                fprintf(['*** s ' num2str(ps) ': only ' num2str(nntrl) ' trls **** \n\n']);
            end
        end
        clear TFR
        TFR = TFR_all;
        clear TFR_all
        save([PPath.FigPath RunCond '_TFR' ],'TFR','-v7.3');
        % rename
        eval([RunCond '_TFR = TFR;']);
        clear TFR
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%====== plotting raw alpha power for PreTarg_SacOn ======= %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%get averaged TFR over subs
    Data_low = WrdOff_TFR.PreTarg_low;
    Data_high = WrdOff_TFR.PreTarg_high;
    % raw data, without removing baseline
    cfg = [];
    cfg.parameter = 'powspctrm';
    plot_low = ft_freqgrandaverage(cfg,Data_low{:});
    plot_high = ft_freqgrandaverage(cfg,Data_high{:});
    plot_all = plot_low;
    plot_all.powspctrm = (plot_low.powspctrm + plot_high.powspctrm)./2;
    %%% transfer the unit from (T/m)^2 to (fT/m)^2
    plot_all.powspctrm = plot_all.powspctrm.*10e24;
    
% %     %% plotting TFR for epochs across all conditions for pre-target wrd_on
% %     %%% for colormap
% %     cmap = colormap(cbrewer('div','RdBu',32));
% %     cmap = colormap(flipud(cmap));
% %     
% %     %%%%%%%%%  plot figures %%%%%%%%%
% %     figname = 'TFR_raw_WrdOn_PreTarg';
% %     h = figure('name',figname,'color',[1 1 1],'position',[0 0 220 150]);
% %     cfg = [];
% %     cfg.baseline = 'no';
% %     cfg.showlabels = 'no';
% %     cfg.comment = 'no';
% %     cfg.layout = 'neuromag306cmb.lay';
% %     cfg.colormap = cmap;
% %     cfg.dataname = figname;
% %     cfg.xlim = [-0.4 0.4]; %% only plot ±0.4s to avoid the 'ugly' edge artefacts
% %     % cfg.zlim = [-0.4 0.4];
% %     colorbar;
% %     ft_multiplotTFR(cfg, plot_all);
% %     title('Raw power during pre-target interval','FontSize',7);
% %     set(gcf, 'renderer', 'painters')
% %     saveas(h,[PPath.FigPath figname]);
    

    %%averaged TFR for all chans
    figname = 'TFR_raw_WrdOff_PreTarg_OverAllSens';
    h = figure('Name',figname,'color',[1 1 1],'position',[0 0 250 150]);
    pcolor(plot_all.time, plot_all.freq, squeeze(mean(plot_all.powspctrm,1)));
    shading interp;
    %%% for colormap
    cmap = colormap(cbrewer('div','RdBu',32));
    cmap = colormap(flipud(cmap));
    title('Averaged power over all planar sensors','FontSize',7,'FontName','Arial')
    % plot some lines
    hold on;
    plot([-0.5 0.5],[9 9],':k','LineWidth',1)
    plot([-0.5 0.5],[13 13],':k','LineWidth',1)
    plot([-0.2 -0.2],[4 30],':k','LineWidth',1)
    plot([0 0],[4 30],':k','LineWidth',1)
    xlim([-0.4 0.4])
    set(gca,'XTick',-0.4:0.2:0.4);
    set(gca,'XTickLabel',{'-0.4','-0.2','saccade on','0.2','0.4'},'FontWeight','normal','FontSize',7,'FontName','Arial')
    set(gca,'YTick',5:5:30);
    set(gca,'YTickLabel',{'5','','','20','25','30'},'FontWeight','normal','FontSize',7,'FontName','Arial')
    text([-0.45 -0.48],[9 13],[{'9'},{'13'}],'FontWeight','normal','FontSize',7,'FontName','Arial')
    ylabel('Frequency (Hz)','FontWeight','normal','FontSize',7,'FontName','Arial')
    xlabel('Time (s)','FontWeight','normal','FontSize',7,'FontName','Arial')
    caxis([0 250])
    colorbar
    text(-0.22,23,{['pre-target'],['interval']},'FontWeight','normal','FontSize',7,'FontName','Arial')
    text(0.72,7,'Raw power (fT/m)^2','FontWeight','normal','FontSize',7,'FontName','Arial','Rotation',90)
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figname]);
    saveas(h,[PPath.FigPath figname],'svg');
    
    %%plot topography figure
    occip_left = {'MEG1642+1643';'MEG1712+1713';'MEG1722+1723';'MEG1732+1733';...
        'MEG1742+1743';'MEG1912+1913';'MEG1922+1923';'MEG1932+1933';'MEG1942+1943'};
    occip_right = {'MEG2312+2313';'MEG2322+2323';'MEG2332+2333';'MEG2342+2343';...
        'MEG2432+2433';'MEG2512+2513';'MEG2522+2523';'MEG2532+2533';'MEG2542+2543'};
    cfg            = [];
    cfg.parameter  ='powspctrm';
    cfg.layout     = 'neuromag306cmb.lay';
    cfg.maskstyle  = 'saturation';
    cfg.contournum = 0;
    cfg.colormap   = cmap;
    cfg.zlim       = [0 600];
    cfg.xlim       = [-0.2 0];
    cfg.ylim       = [9 13];
    cfg.gridscale  = 360;
    cfg.marker     = 'no'; %'numbers';
    cfg.comment    = 'no';
    cfg.style      = 'straight';
    cfg.highlightchannel = [occip_left;occip_right]; %OccipSens;
    cfg.highlight        = 'on';
    figname        = 'TFR_raw_WrdOff_PreTarg_topo';
    h              = figure('name',figname,'color',[1 1 1],'position',[0 0 220 150]);
    ft_topoplotTFR(cfg,plot_all);
    title({['Topography for alpha power'],['during the pre-target interval']},'FontWeight','normal','FontSize',7,'FontName','Arial');
    colorbar
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figname]);
    saveas(h,[PPath.FigPath figname],'svg');
        
    %%plot lateralized alpha power
    figname = 'TFR_WrdOff_PreTarg_alpha_LI';
    % select left/right alpha power
    cfg = [];
    cfg.channel = occip_left;
    alpha_left = ft_selectdata(cfg,plot_all);
    cfg.channel = occip_right;
    alpha_right = ft_selectdata(cfg,plot_all);
    % replace powspctrm with lateralization index
    plot_all.powspctrm = (alpha_right.powspctrm-alpha_left.powspctrm)./(alpha_right.powspctrm+alpha_left.powspctrm);
    h = figure('Name',figname,'color',[1 1 1],'position',[0 0 250 150]);
    pcolor(plot_all.time, plot_all.freq, squeeze(mean(plot_all.powspctrm,1)));
    shading interp;
    caxis([-0.07 0.07])
    %%% for colormap
    cmap = colormap(cbrewer('div','RdBu',32));
    cmap = colormap(flipud(cmap));
    title('Averaged power over all planar sensors','FontSize',7,'FontName','Arial')
    % plot some lines
    hold on;
    plot([-0.5 0.5],[9 9],':k','LineWidth',1)
    plot([-0.5 0.5],[13 13],':k','LineWidth',1)
    plot([-0.2 -0.2],[4 30],':k','LineWidth',1)
    plot([0 0],[4 30],':k','LineWidth',1)
    xlim([-0.4 0.4])
    set(gca,'XTick',-0.4:0.2:0.4);
    set(gca,'XTickLabel',{'-0.4','-0.2','saccade on','0.2','0.4'},'FontWeight','normal','FontSize',7,'FontName','Arial')
    set(gca,'YTick',5:5:30);
    set(gca,'YTickLabel',{'5','','','20','25','30'},'FontWeight','normal','FontSize',7,'FontName','Arial')
    text([-0.45 -0.48],[9 13],[{'9'},{'13'}],'FontWeight','normal','FontSize',7,'FontName','Arial')
    ylabel('Frequency (Hz)','FontWeight','normal','FontSize',7,'FontName','Arial')
    xlabel('Time (s)','FontWeight','normal','FontSize',7,'FontName','Arial')
    caxis([0 250])
    colorbar
    text(-0.22,23,{['pre-target'],['interval']},'FontWeight','normal','FontSize',7,'FontName','Arial')
    text(0.72,7,'Raw power (fT/m)^2','FontWeight','normal','FontSize',7,'FontName','Arial','Rotation',90)
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figname]);
    saveas(h,[PPath.FigPath figname],'svg');
    
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%====== plotting difference alpha power for PreTarg_SacOn ======= %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%plot power-diff for WrdOff_pre-target interval over sig chans
    % get the sig sens
    load('Z:\Saccade_AlphaPhase\Results\ITC\SigSens_WrdOffPreTarg') %SigSens 
    % get power for both conditions
    Data_low = WrdOff_TFR.PreTarg_low;
    Data_high = WrdOff_TFR.PreTarg_high;
    grandavg_low = ft_freqgrandaverage(cfg, Data_low{:});
    grandavg_high = ft_freqgrandaverage(cfg, Data_high{:});
    % get power-diff
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'subtract';
    tfr_diff      = ft_math(cfg,grandavg_low, grandavg_high);
    % get index
    freq_id = dsearchn([grandavg_low.freq]',[9:13]');
    time_id = dsearchn([grandavg_low.time]',[-0.2:0.001:0]');
    chan_id = SigSens.mask;
    % get averaged tfr_diff over sig sens
    tfr_diff_avg = squeeze(mean(tfr_diff.powspctrm(chan_id,:,:),1));
    tfr_diff_avg = tfr_diff_avg.*10e24;%transfer unit from (T/m)^2 to (fT/m)^2
    %% plot
    figname = 'TFR_diff_WrdOff_PreTarg';
    h = figure('Name',figname,'color',[1 1 1],'position',[0 0 250 150]);
    pcolor(tfr_diff.time, tfr_diff.freq,tfr_diff_avg);
    shading interp;
    % for colormap
    cmap = colormap(cbrewer('div','RdBu',32));
    cmap = colormap(flipud(cmap));
    title({['Averaged \Delta power over sig sensors'],['saccade on to target']},'FontWeight','normal','FontSize',7,'FontName','Arial')
    % plot some lines
    hold on;
    plot([-0.5 0.5],[9 9],':k','LineWidth',1)
    plot([-0.5 0.5],[13 13],':k','LineWidth',1)
    plot([-0.2 -0.2],[4 30],':k','LineWidth',1)
    plot([0 0],[4 30],':k','LineWidth',1)
    xlim([-0.4 0.4])
    set(gca,'XTick',-0.4:0.2:0.4);
    set(gca,'XTickLabel',{'-0.4','-0.2','saccade on','0.2','0.4'},'FontWeight','normal','FontSize',7,'FontName','Arial')
    set(gca,'YTick',5:5:30);
    set(gca,'YTickLabel',{'5','','','20','25','30'},'FontWeight','normal','FontSize',7,'FontName','Arial')
    text([-0.45 -0.48],[9 13],[{'9'},{'13'}],'FontWeight','normal','FontSize',7,'FontName','Arial')
    ylabel('Frequency (Hz)','FontWeight','normal','FontSize',7,'FontName','Arial')
    xlabel('Time (s)','FontWeight','normal','FontSize',7,'FontName','Arial')
    caxis([-30 30])
    colorbar
    text(-0.19,23,{['pre-target'],['interval']},'FontWeight','normal','FontSize',7,'FontName','Arial')
    text(0.66,6,{'\Delta power(T/cm)^2','(low-high lexical target)'},'FontWeight','normal','FontSize',7,'FontName','Arial','Rotation',90)
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figname]);
    saveas(h,[PPath.FigPath figname],'svg');
    
    
    %% simple ttest for alpha power for WrdOff_pre-target interval over sig chans
    nsub = length(WrdOff_TFR.PreTarg_low);
    alpha_pow = nan(nsub,2);
    for s = 1:nsub
        chan_id = cellfun(@(x) find(strcmp(x,WrdOff_TFR.PreTarg_low{1,s}.label)),SigSens.label,'Uni',true);
        tmp_low = WrdOff_TFR.PreTarg_low{1,s}.powspctrm(chan_id,freq_id,time_id);
        tmp_low = mean(mean(mean(tmp_low,1),2),3);
        tmp_high = WrdOff_TFR.PreTarg_high{1,s}.powspctrm(chan_id,freq_id,time_id);
        tmp_high = mean(mean(mean(tmp_high,1),2),3);
        alpha_pow(s,:) = [tmp_low tmp_high];
    end
    alpha_pow = alpha_pow.*10e24;%transfer unit from (T/m)^2 to (fT/m)^2
    [~,p,~,stat] = ttest(alpha_pow(:,1),alpha_pow(:,2)); % [t=-1.189 p=0.242]
    stat.p = p
    % save data
    AlphaPow_WrdOff_PreTarg.alpha_pow = alpha_pow;
    AlphaPow_WrdOff_PreTarg.stat = stat;
    save([PPath.FigPath 'AlphaPow_WrdOff_PreTarg'],'AlphaPow_WrdOff_PreTarg','-v7.3');
    %% % plot
    colmat = [0 114 189;217 83 25]./255;
    figtitle = 'AlphaPow_WrdOff_PreTarg_violin';
    h = figure('Name',figtitle,'color',[1 1 1],'position',[0 0 220 150]);
    % violin
    group = [cellstr(repmat('Low',nsub,1)); cellstr(repmat('High',nsub,1))];
    grouporder={'Low','High'};
    vdata = [alpha_pow(:,1); alpha_pow(:,2)];
    vp = violinplot(vdata, group,'GroupOrder',grouporder);
    vp(1).ViolinColor = colmat(1,:);
    vp(2).ViolinColor = colmat(2,:);
    vp(1).ShowMean = 1; vp(2).ShowMean = 1;
    vp(1,1).MeanPlot.LineWidth = 1.5;
    vp(1,2).MeanPlot.LineWidth = 1.5;
    vp(1,1).BoxWidth = 0.01; vp(1,2).BoxWidth = 0.01;
    %%% plot the line linking each subject
    hold on;
    x1 = vp(1,1).ScatterPlot.XData;
    y1 = vp(1,1).ScatterPlot.YData;
    x2 = vp(1,2).ScatterPlot.XData;
    y2 = vp(1,2).ScatterPlot.YData;
    plot([x1; x2],[y1; y2],'Color',[.8 .8 .8],'linewidth',0.3)    
    plot([1 2],[1900 1900],'k','LineWidth',1)
    text(1.5,2000,'n.s.','FontWeight','normal','FontSize',7,'FontName','Arial')
    ylabel('Power (T/cm)^2','FontSize',7,'FontWeight','normal','FontName','Arial');
    xlabel('Target lexical frequency','FontSize',7,'FontWeight','normal','FontName','Arial');
    title('Averaged alpha power for pre-target','FontSize',7,'FontWeight','normal','FontName','Arial');
    set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
    set(gca,'box','off','LineWidth',1)
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figtitle]);
    saveas(h,[PPath.FigPath figtitle],'svg');
    
    %%%% fow R
    subid = [[1:38]';[1:38]'];
    wrdfrq = [ones(38,1); 2*ones(38,1)];
    power = [alpha_pow(:,1);alpha_pow(:,2)];
    [subid wrdfrq power];
    
end





























