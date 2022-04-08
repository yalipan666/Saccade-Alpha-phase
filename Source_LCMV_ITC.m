%%% doing source modeling for the inter-trial-coherence using LCMV (Linearly Constrained Minimum Variance)
%%% https://www.fieldtriptoolbox.org/workshop/aarhus/beamformingerf/#eeg-source-analysis

%% 20210903
% 1. Just do source modeling for the 38 valid subs that were used
%in the ITC_alpha analysis (criteria is based on trial numbers)
% 2. Using the same trls from ITC sensor level analysis for the source level
% 3. Create beamformer based on alpha 9-13Hz, [-0.2 0]s only
% 4. play around the lamda

function Source_LCMV_ITC(sid)
%%% local
%%% basic settingups
server = 1;
if server
    saccade_path = '/rds/projects/2018/jenseno-reading/Saccade_AlphaPhase/';
    lexical_path = '/rds/projects/2018/jenseno-reading/Lexical/';
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
    ft_defaults
else
    %%% local
    saccade_path = 'Z:\Saccade_AlphaPhase\';
    lexical_path = 'Z:\Lexical\';
    addpath(genpath([saccade_path 'codes']));
    addpath Z:\fieldtrip-20200220\
    ft_defaults
end

%%% setting path
PPath.FigPath = [saccade_path 'Results' filesep 'source' filesep];
PPath.HMPath = [lexical_path 'data sharing' filesep 'HeadModel' filesep];
PPath.MEGPath = [lexical_path 'RawData' filesep 'MEG_data' filesep];
PPath.Data = [lexical_path 'Analyse_data' filesep];%data are under the folder of lexical
figname = '13Hz02s_fourier';

%%% get sub information
load([saccade_path 'Results' filesep 'ITC' filesep 'WrdOff_ITC.mat'])
sub_all = ITC.subs;
TrlIdx_equ = ITC.PreTarg_TrlIdx_equ; % [low high]
clear ITC
nsub = length(sub_all);
RunCond = 'WrdOff';
%%% setup outputs
Source.Low = cell(nsub,1);
Source.High = cell(nsub,1);

%% %================= subject loop
for sss = sid
    disp(['******* analyzing sub ' num2str(sss) '******'])
    tic
    sub = sub_all{sss,1};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% preparing headmodel %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([PPath.HMPath 'Hdm_MRIaligned_' sub]); %%hdm,mri_aligned
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% loading combined dataset %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%=== JEP 60
    PPath.SaveData = [PPath.Data sub '_JEP' filesep 'JEP60_TW1000' filesep];
    load([PPath.SaveData 'epoch_' RunCond]);
    eval(['epochdata_jep = epoch_' RunCond ';']);
    %%%=== targ 60
    PPath.SaveData = [PPath.Data sub filesep 'Targ60_TW1000' filesep];
    load([PPath.SaveData 'epoch_' RunCond]);
    eval(['epochdata = epoch_' RunCond ';']);
    %%%=== combining 2 datasets together
    epochdata.trialinfo = [epochdata.trialinfo; epochdata_jep.trialinfo];
    epochdata.trial = [epochdata.trial epochdata_jep.trial];
    epochdata.time = [epochdata.time epochdata_jep.time];
    clear epochdata_jep
    clear epoch_*
    
    %%% get first fixation trials (pre_target + target + post_targ).Because
    %%% TrlIdx_equ is according to epochdata_all
    validduration = find(epochdata.trialinfo(:,10) == 1); % FirstPassFix
    cfg = [];
    cfg.trials = validduration;
    cfg.channel = 'MEGGRAD';
    epochdata = ft_selectdata(cfg, epochdata);
    %%% get the hdr.grad
    megfile = [PPath.MEGPath sub filesep sub(3:8) filesep sub(end-3:end) '-1.fif'];
    cfg            = [];
    cfg.dataset    = megfile;
    hdr            = ft_read_header(cfg.dataset);
    epochdata.grad = hdr.grad;
    %%% select low and high condition using trl_id from ITC_alpha sensor
    %%% level analysis, where the trl_id will select pre-target
    %%% automatically
    cfg         = [];
    cfg.trials  = TrlIdx_equ{sss,1};
    epochdata_low = ft_selectdata(cfg, epochdata);
    cfg.trials  = TrlIdx_equ{sss,2};
    epochdata_high = ft_selectdata(cfg, epochdata);
    %%% select pre-target for epoch_all
    cfg        = [];
    cfg.trials = epochdata.trialinfo(:,3)==-1;
    epochdata  = ft_selectdata(cfg, epochdata);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% computing the leadfield %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% dir for stabdard MNI spac,where the template source models are
    [~, ftdir] = ft_version;
    templatedir = fullfile(ftdir, 'template', 'sourcemodel');
    template = load(fullfile(templatedir, 'standard_sourcemodel3d5mm')); % 5mm spacing grid
    %%% inverse-warp the template grid to subject specific coordinates
    cfg           = [];
    cfg.warpmni   = 'yes';
    cfg.template  = template.sourcemodel;
    cfg.nonlinear = 'yes'; % use non-linear normalization
    cfg.mri       = mri_aligned;
    sourcemodel   = ft_prepare_sourcemodel(cfg);
    % prepare leadfields
    cfg = [];
    cfg.channel = epochdata.label;
    cfg.grad = epochdata.grad; % sensor position
    cfg.sourcemodel = sourcemodel; %source model from ft
    cfg.headmodel = hdm; % volume conductor head model from first step
    cfg.reducerank = 2; %'no', or number(default = 3 for EEG, 2 for MEG)
    leadfield = ft_prepare_leadfield(cfg);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% get LCMV filters %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% bandpass filtering for creating the beamformer
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [9 13];
    epochdata = ft_preprocessing(cfg, epochdata);
    % calculate the common filter using all data
    cfg                  = [];
    cfg.covariance       = 'yes';
    cfg.covariancewindow = [-0.2 0]; %the whole time window
    avg                  = ft_timelockanalysis(cfg,epochdata); %use overall data to create a common filter
    cfg                  = [];
    cfg.method           = 'lcmv';
    cfg.sourcemodel      = leadfield;
    cfg.headmodel        = hdm;
    cfg.lcmv.keepfilter  = 'yes'; %keep filter to apply it to diff cond data later
    cfg.lcmv.fixedori    = 'yes';
    cfg.lambda           = '5%'; %smoothing parameter for covariance matrix regularization
    sourceavg            = ft_sourceanalysis(cfg, avg);
    sourceavg.pos        = template.sourcemodel.pos;
    sourceavg.dim        = template.sourcemodel.dim;
    
    % second and third call to ft_sourceanalysis now applying the
    % precomputed filters to epochs from diff condition
    cfg                    = [];
    cfg.method             = 'lcmv';
    cfg.sourcemodel        = leadfield;
    cfg.headmodel          = hdm;
    cfg.sourcemodel.filter = sourceavg.avg.filter;
    ITC_low                = ft_sourceanalysis(cfg, epochdata_low); %%% get the fake pow struct
    ITC_high               = ft_sourceanalysis(cfg, epochdata_high);
    ITC_low.pos            = template.sourcemodel.pos;
    ITC_low.dim            = template.sourcemodel.dim;
    ITC_high.pos           = template.sourcemodel.pos;
    ITC_high.dim           = template.sourcemodel.dim;
    
    %%% applying the filter
    W = sourceavg.avg.filter(sourceavg.inside); % only get the non-empty filter
    %%%% get hilbert filert data with the 1s long epoch, avoding filtering
    %%%% artefacts
    cfg = [];
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [9 13];
    cfg.hilbert    = 'complex';
    cfg.keeptrials = 'yes';
    epochdata_low  = ft_preprocessing(cfg,epochdata_low);
    epochdata_high = ft_preprocessing(cfg,epochdata_high);
    cfg = [];
    cfg.latency = [-0.2 0];
    epochdata_low = ft_selectdata(cfg, epochdata_low);
    epochdata_high = ft_selectdata(cfg, epochdata_high);
    ntrl = length(epochdata_low.trial);
    ngrid = length(W);
    ntim = length(epochdata_low.time{1,1});
    FourierSpct_low = zeros(ntrl, ngrid, ntim);
    FourierSpct_high = zeros(ntrl, ngrid, ntim);
    for tt = 1:ntrl
        FourierSpct_low(tt,:,:) = cell2mat(cellfun(@(x) x*epochdata_low.trial{tt}, W, 'Uni', false));
        FourierSpct_high(tt,:,:) = cell2mat(cellfun(@(x) x*epochdata_high.trial{tt}, W, 'Uni', false));
    end
    
    %%% compute inter-trial phase coherence (itpc)
    itc = FourierSpct_low./abs(FourierSpct_low); % divide by amplitude
    itc = sum(itc,1); % sum angles
    itc = abs(itc)/ntrl; % take the absolute value and normalize
    ITC_low.avg.pow(sourceavg.inside) = mean(squeeze(itc),2);
    itc = FourierSpct_high./abs(FourierSpct_high);         % divide by amplitude
    itc = sum(itc,1);   % sum angles
    itc = abs(itc)/ntrl;   % take the absolute value and normalize
    ITC_high.avg.pow(sourceavg.inside) = mean(squeeze(itc),2);
    save([PPath.FigPath 'Source_' num2str(sss)],'ITC_low','ITC_high','-v7.3')
    
    %%% plot fig for each sub
    % get diff
    cfg = [];
    cfg.parameter     = 'avg.pow';
    cfg.operation     = 'x1-x2';
    source_wrdfrq_avg = ft_math(cfg, ITC_low, ITC_high);
    % get template_mri
    [~, ftdir]            = ft_version;
    templatedir           = fullfile(ftdir, 'template', 'headmodel');
    template_mri          = ft_read_mri(fullfile(templatedir, 'standard_mri.mat'));
    template_mri.coordsys = 'mni'; % we know it's in MNI space
    % interpolate
    cfg               = [];
    cfg.parameter     = 'pow';
    cfg.interpmethod  = 'nearest';
    cfg.coordsys      = 'mni';
    source_wrdfrq_avg = ft_sourceinterpolate(cfg, source_wrdfrq_avg,template_mri);
    % ortho plot
    cfg               = [];
    cfg.method        = 'ortho';
    cfg.funparameter  = 'pow';
    cfg.funcolormap   = 'hot';
    cfg.funcolorlim    = [0 max(source_wrdfrq_avg.pow)];
    h = ft_sourceplot(cfg, source_wrdfrq_avg);
    saveas(h,[PPath.FigPath figname '_S' num2str(sss)],'bmp');
    close all
    disp(['******* Done sub-' num2str(sss) '******'])
    toc
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%============= compute grand average =============%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if server == 0
    for sss = 1:nsub
        load([PPath.FigPath 'Source_' num2str(sss)])
        Source.Low{sss} = ITC_low;
        Source.High{sss} = ITC_high;
    end
    clear ITC_*
    
    %%%%%%%%===== Source Statistics =====%%%%%%%%%%%
    %%% run statistics over subjects %https://www.fieldtriptoolbox.org/example/source_statistics/
    cfg                  = [];
    cfg.dim              = Source.Low{1}.dim;
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.parameter        = 'avg.pow';
    cfg.correctm         = 'cluster';
    cfg.clusterstatistic = 'maxsum';
    cfg.clusteralpha     = 0.05;
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.05;
    cfg.correcttail      = 'prob';
    cfg.numrandomization = 1000;
    nsubj = numel(Source.Low);
    cfg.design = [1:nsub 1:nsub; ones(1,nsub) ones(1,nsub)*2];
    cfg.uvar             = 1; % row in design indicating subjects, repeated measure
    cfg.ivar             = 2; % row in design indicating condition for contrast
    stat = ft_sourcestatistics(cfg, Source.Low{:}, Source.High{:});
    any(stat.mask)
    
    %%%%=== get grand average source model and plot figures
    cfg                    = [];
    cfg.keepindividual     = 'no';
    cfg.parameter          = 'avg.pow';
    source_low_avg         = ft_sourcegrandaverage(cfg, Source.Low{:});
    source_high_avg        = ft_sourcegrandaverage(cfg, Source.High{:});
    source_low_avg.dimord  = 'pos';
    source_high_avg.dimord = 'pos';
    clear Source
    %get diff
    cfg = [];
    cfg.parameter     = 'pow';
    cfg.operation     = 'x1-x2';
    source_wrdfrq_avg = ft_math(cfg, source_low_avg, source_high_avg);
    % get template_mri
    [~, ftdir]            = ft_version;
    templatedir           = fullfile(ftdir, 'template', 'headmodel');
    template_mri          = ft_read_mri(fullfile(templatedir, 'standard_mri.mat'));
    template_mri.coordsys = 'mni'; % we know it's in MNI space
    % interpolate
    cfg               = [];
    cfg.parameter     = 'pow';
    cfg.interpmethod  = 'nearest';
    cfg.coordsys      = 'mni';
    source_wrdfrq_avg = ft_sourceinterpolate(cfg, source_wrdfrq_avg,template_mri);
    % ortho plot
    cfg               = [];
    cfg.method        = 'ortho';
    cfg.funparameter  = 'pow';
    cfg.maskparameter = 'pow';
    cfg.funcolormap    = 'hot';
    cfg.funcolorlim    = [0 0.03];
    h = ft_sourceplot(cfg, source_wrdfrq_avg);
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figname]);
    saveas(h,[PPath.FigPath figname '.svg']);
    % surface plot
    cfg = [];
    cfg.nonlinear     = 'no';
    sourceWrdfrqIntNorm = ft_volumenormalise(cfg, source_wrdfrq_avg);
    cfg = [];
    cfg.method         = 'surface';
    cfg.funparameter   = 'pow';
    cfg.maskparameter  = cfg.funparameter;
    cfg.funcolormap    = 'jet';
    cfg.opacitymap     = 'rampup';
    cfg.projmethod     = 'nearest';
    cfg.surfdownsample = 10;
    cfg.funcolorlim    = [0 0.02];
    cfg.surffile       = 'surface_white_both.mat';
    h = ft_sourceplot(cfg, sourceWrdfrqIntNorm);
    saveas(h,[PPath.FigPath figname '_surface']);
end

%% ========= simple t-test for each source-grid, no cluster-correction
% get t-vlues for the whole source grid and plot it
PPath.FigPath = 'Z:\Saccade_AlphaPhase\Results\source_P\';
nsub = 38;
for sss = 1:nsub
    load([PPath.FigPath 'Source_' num2str(sss)])
    Source.Low{sss} = ITC_low;
    Source.High{sss} = ITC_high;
end
% get the dummy struct to plot source in MNI space
cfg                    = [];
cfg.keepindividual     = 'no';
cfg.parameter          = 'avg.pow';
source_low_avg         = ft_sourcegrandaverage(cfg, Source.Low{:});
source_high_avg        = ft_sourcegrandaverage(cfg, Source.High{:});
source_low_avg.dimord  = 'pos';
source_high_avg.dimord = 'pos';
cfg = [];
cfg.parameter     = 'pow';
cfg.operation     = 'x1-x2';
source_wrdfrq_avg = ft_math(cfg, source_low_avg , source_high_avg);
% put the t-values into the dummy struct
p_values = source_wrdfrq_avg.pow;
for gg = 1:length(source_wrdfrq_avg.pow) %loop over grids
    if ~isnan(source_wrdfrq_avg.pow(gg))
        data_low = cellfun(@(x) x.avg.pow(gg),Source.Low,'UniformOutput',false);
        data_high = cellfun(@(x) x.avg.pow(gg),Source.High,'UniformOutput',false);
        [~,p,~,stats] = ttest(cell2mat(data_low'),cell2mat(data_high'));
        source_wrdfrq_avg.pow(gg) = stats.tstat;
        p_values(gg) = 1-p;
    end
end
% get template_mri
[ftver, ftdir]        = ft_version;
templatedir           = fullfile(ftdir, 'template', 'headmodel');
template_mri          = ft_read_mri(fullfile(templatedir, 'standard_mri.mat'));
template_mri.coordsys = 'mni'; % we know it's in MNI space

%% plot t-values
figname = 'Tvalues_NOcorrect';
% interpolate
cfg               = [];
cfg.parameter     = 'pow';
cfg.interpmethod  = 'nearest';
cfg.coordsys      = 'mni';
source_wrdfrq_avg = ft_sourceinterpolate(cfg, source_wrdfrq_avg,template_mri);
% ortho plot
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap    = 'hot';
cfg.funcolorlim    = [2.03 4]; %t=2.026 for two-tailed ttest with df=37
h = ft_sourceplot(cfg, source_wrdfrq_avg);
saveas(h,[PPath.FigPath figname]);
% surface plot
cfg = [];
cfg.nonlinear     = 'no';
sourceWrdfrqIntNorm = ft_volumenormalise(cfg, source_wrdfrq_avg);
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolormap    = 'jet';
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surfdownsample = 10;
cfg.funcolorlim    = [2.03 2.5];
cfg.surffile       = 'surface_white_both.mat';
h = ft_sourceplot(cfg, sourceWrdfrqIntNorm);
saveas(h,[PPath.FigPath figname '_surface']);

%% plot p-values
figname = 'Pvalues_NOcorrect';
cfg = [];
cfg.parameter     = 'pow';
cfg.operation     = 'x1-x2';
source_wrdfrq_p = ft_math(cfg, source_low_avg , source_high_avg);
source_wrdfrq_p.pow = p_values;
% interpolate
cfg               = [];
cfg.parameter     = 'pow';
cfg.interpmethod  = 'nearest';
cfg.coordsys      = 'mni';
source_wrdfrq_p = ft_sourceinterpolate(cfg, source_wrdfrq_p,template_mri);
% ortho plot
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap    = 'hot';
cfg.funcolorlim    = [0.95 1];
h = ft_sourceplot(cfg, source_wrdfrq_p);
saveas(h,[PPath.FigPath figname]);
% surface plot
cfg = [];
cfg.nonlinear     = 'no';
sourceWrdfrqIntNorm = ft_volumenormalise(cfg, source_wrdfrq_p);
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolormap    = 'jet';
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surfdownsample = 10;
cfg.funcolorlim    = [0.9 1];
cfg.surffile       = 'surface_white_both.mat';
h = ft_sourceplot(cfg, sourceWrdfrqIntNorm);
saveas(h,[PPath.FigPath figname '_surface']);







