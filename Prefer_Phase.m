%% 20220906 control analyses required by reviewers
% get the actual preferred phase for each conditions
% use exactly the same trials as in the main PLV analysis
% only for the sig planar sensors (not combination)
% only for alpha band (9-13Hz), [-0.2 0]s time window

function Prefer_Phase(sid)
tic
server = 1;
%%%% running setting
runconds = {'WrdOff'};
EpochType = {'PreTarg'};
targid = -1;
conds = {'low','high'};
%%% time window parameters
presac_win = [-0.2 0]; %%% select epochs aligned with WrdOff
%%% parameter for hilbert filter
freqrange = 9:1:13;
freqwidth_factor = 0.3; % freq widh: ±0.3*fcenter(Hz)

%%% basic settingup
if server
    rootdir = '/rds/projects/2018/jenseno-reading/';
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
    ft_defaults
else
    %%% local
    rootdir = 'Z:\';
    addpath Z:\fieldtrip-20200220\
    ft_defaults
end

%%% paths
PPath.Data = [rootdir 'Lexical' filesep 'Analyse_data' filesep];%data are under the folder of lexical
PPath.FigPath = [rootdir 'Saccade_AlphaPhase' filesep 'Results' filesep 'ITC' filesep];

%%% get the sig sensors
load([PPath.FigPath 'SigSens_WrdOffPreTarg.mat'])
sens1 = cellfun(@(x) x(1:7),SigSens.label,'Uni',false);
sens2 = cellfun(@(x) x([1:3 9:12]),SigSens.label,'Uni',false);
Sens = cell(2*length(sens1),1);
Sens(1:2:end) = sens1;
Sens(2:2:end) = sens2;

%%% get the trl id
load([PPath.FigPath 'WrdOff_ITC.mat'])
Subs = ITC.subs;
TrlID = ITC.PreTarg_TrlIdx_equ;
clear ITC


%% start to get the preferred phase
if server == 1
    for rrr = 1:length(runconds)
        RunCond = runconds{rrr};
        %%% subject loop
        for sss = sid
            fprintf(['********** s' num2str(sid) ' ******** \n\n']);
            %%% set output
            ITC = [];
            ITC.conds = conds;
            ITC.sub = Subs{sss};
            %%%%%%%%%%========= getting the combined dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%=== JEP 60
            PPath.SaveData = [PPath.Data ITC.sub '_JEP' filesep 'JEP60_TW1000' filesep];
            load([PPath.SaveData 'epoch_' RunCond]);
            eval(['epochdata_jep = epoch_' RunCond ';']);
            %%%=== targ 60
            PPath.SaveData = [PPath.Data ITC.sub filesep 'Targ60_TW1000' filesep];
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
            cfg.channel = Sens;%'MEGGRAD';
            epochdata = ft_selectdata(cfg, epochdata);
            
            for mmm = 1:length(EpochType)
                %%% get equal trial number across all conditions
                cfg         = [];
                cfg.trials  = TrlID{sss,1};
                data_low    = ft_selectdata(cfg, epochdata);
                cfg.trials  = TrlID{sss,2};
                data_high   = ft_selectdata(cfg, epochdata);
                
                %%%%%%%%%%%=================== Run the analysis===============%%%%%%%%%%%%
                for ccc = 1:length(conds)
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
                    %%% get preferred
                    itc.phase = squeeze(angle(itc.itpc));
                    % get ITPC
                    itc.itpc      = abs(itc.itpc)/ntrl;   % take the absolute value and normalize
                    itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
                    itc.powspctrm = itc.itpc;
                    itc = rmfield(itc,'itpc');
                    % only select the presac_win
                    cfg         = [];
                    cfg.latency = presac_win;
                    itc    = ft_selectdata(cfg, itc);
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
        for ps = 1:38
            fprintf(['*** s' num2str(ps) ' *** \n']);
            load([PPath.FigPath RunCond '_ITC_' num2str(ps)]);
            ss = ss + 1;
            ITC_all.subs{ss,1} = ITC.sub;
            for fd = 1:length(EpochType)
                eval(['ITC_all.' EpochType{fd} '_low{ss} = ITC.' EpochType{fd} '_low;']);
                eval(['ITC_all.' EpochType{fd} '_high{ss} = ITC.' EpochType{fd} '_high;']);
            end
        end
        %%% in total, 38 subs
        clear ITC
        ITC_PreferPhase = ITC_all;
        clear ITC_all
        save([PPath.FigPath 'ITC_PreferPhase'],'ITC_PreferPhase','-v7.3');
    end
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%====Suppleentary Fig S4: preferred phase for PLI ===%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%plot figs for all sub
    fid = 2; %10Hz
    tim = ITC_PreferPhase.PreTarg_low{1,1}.time; 
    sen = ITC_PreferPhase.PreTarg_low{1,1}.label;
    sub = 1:length(ITC_PreferPhase.subs);
    for c = 5 %loop over sensors, 5--MEG2032
        plv_low = zeros(length(sub),length(tim));
        plv_high = zeros(length(sub),length(tim));
        pha_low = zeros(length(sub),length(tim));
        pha_high = zeros(length(sub),length(tim));
        for s = sub
            plv_low(s,:) = squeeze(ITC_PreferPhase.PreTarg_low{1,s}.powspctrm(c,fid,:));
            plv_high(s,:) = squeeze(ITC_PreferPhase.PreTarg_high{1,s}.powspctrm(c,fid,:));
            pha_low(s,:) = squeeze(ITC_PreferPhase.PreTarg_low{1,s}.phase(c,fid,:));
            pha_high(s,:) = squeeze(ITC_PreferPhase.PreTarg_high{1,s}.phase(c,fid,:));
        end
        % plot
        figname = ['PreferPhae-Sensor-' sen{c} '-10Hz'];
        figure('name',figname,'color',[1 1 1]);
        subplot(2,2,1)
        pcolor(tim,sub,plv_low)
        shading interp
        title('PLV: low wordfreq')
        caxis([0 0.2])
        xlabel('Time (s)')
        ylabel('Subject ID')
        colorbar
        
        subplot(2,2,2)
        pcolor(tim,sub,pha_low)
        shading interp
        title('Preferred phase: low wordfreq')
        xlabel('Time (s)')
        ylabel('Subject ID')
        colorbar
        
        subplot(2,2,3)
        pcolor(tim,sub,plv_high)
        shading interp
        title('PLV: high wordfreq')
        caxis([0 0.2])
        xlabel('Time (s)')
        ylabel('Subject ID')
        colorbar
       
        subplot(2,2,4)
        pcolor(tim,sub,pha_high)
        shading interp
        title('Preferred phase: high wordfreq')
        xlabel('Time (s)')
        ylabel('Subject ID')
        colorbar
        saveas(gcf,[PPath.FigPath figname]);
    end
    
    
% plot preferred phases at MEG2032 at 10Hz at -0.1s over all subs
sen = 5;
tim = 101;
frq = 2;
pha_low = cellfun(@(x) x.phase(sen,frq,tim),ITC_PreferPhase.PreTarg_low,'uni', true);
pha_high = cellfun(@(x) x.phase(sen,frq,tim),ITC_PreferPhase.PreTarg_high,'uni', true);
figname = 'PreferPhae-MEG032-10Hz--100ms';
figure('name',figname,'color',[1 1 1]);
subplot(1,2,1)        
edges = -4:0.5:4;
h = histogram(pha_low,edges,'Normalization','probability');
hold on;
h = histogram(pha_high,edges,'Normalization','probability');
legend('low','high')
title('Actual preferred phases','fontsize',12) 
xlabel('Alpha phase (radian)','fontsize',11)
ylabel('Proportion','fontsize',11)
subplot(1,2,2)        
h = histogram(pha_low-pha_high,edges,'Normalization','probability');
title('Difference of preferred phases (low-high)','fontsize',12) 
xlabel('Alpha phase (radian)','fontsize',11)
ylabel('Proportion','fontsize',11)
set(h,'FaceColor',[.7 .7 .7])
saveas(gcf,[PPath.FigPath figname]);
 
% test weather the distribution of the phase_diff is a uniform distribution
% using Kolmogorov-Smirnov test
pd = makedist('uniform');
pd.Lower = min(h.Values);
pd.Upper = max(h.Values);
[~,p] = kstest(h.Values,'cdf',pd) %p=0.34
      
end


















