%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading Experiment with Rapid Invisible Frequency Tagging
% 20210216 Yai Pan (yalipan666@gmail.com)
% presentation:
% § one line of sentence on the screen
% § using 'OpenOffscreenwindow' to avoid any delays relatd to screen flip during
%   rapid frequency mode
% § applying a gaussian window to smooth the edges of flickering patch and make it 'invisible'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function RIFT_reading_demo

cfg.debugmode = 1; % debug mode, no Propix

% Input dialog
prompt = {'Subject code:', ...
    'Subject Number:', ...
    'Screen width (cm): ', ...
    'Screen height (cm): ', ...
    'Screen distance (cm): '};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'A','2','70.6','39.5','145'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

%Physical screen parameters (cm) (input)
cfg.width=str2double(answer{3});   %projection screen width in cm
cfg.height=str2double(answer{4});  %projection screen height in cm
cfg.dist=str2double(answer{5});    %distance from subject eye to screen in cm

%%%%%=====PTB screen Initialization=====%%%%%%%
AssertOpenGL;
PsychDefaultSetup(2);    % call some default settings for setting up Psychtoolbox

%Open screen
screens = Screen('Screens'); % Get the screen numbers
cfg.screenNumber = max(screens); %select screen
cfg.ScrBgc = [0.5 0.5 0.5];% backgroud color
%%%  Get the size of the on screen window and set resolution
sc_info = Screen('Resolution', cfg.screenNumber);
resx = sc_info.width;
resy = sc_info.height;
cfg.resx = resx;
cfg.resy = resy;

%%%%%======= open PTB window
[window] = PsychImaging('OpenWindow', cfg.screenNumber, cfg.ScrBgc,[0 0 800 800]);
cfg.window = window;
%%% for the rapid drawing in the offscreen window
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');

% Query the maximum priority level
topPriorityLevel = MaxPriority(window);
HideCursor;
Priority(topPriorityLevel);


% enable alpha blending, also required for using offscreen window
Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% Flip to clear
Screen('Flip', window);

%%%%%=====Stimuli settings=====%%%%%%
%screen setup
cfg.WordSpace = 0.35; % unit in visual angle, equal space between each word;
WordSpace = usrDeg2Pix(cfg.WordSpace,cfg);
cfg.WordStart = 4*cfg.WordSpace;%6*cfg.WordSpace; %% unit in visual angle, the start point of sentence
WordStart = usrDeg2Pix(cfg.WordStart,cfg);
cfg.TextStyle = 1;                  %0=normal,1=bold,2=italic,4=underline,8=outline,32=condense,64=extend.
cfg.TextFont = 'Courier New';
cfg.TextSize = 22;%22; %22--0.3556 visdeg -> 12.3pix; 24--0.3951
cfg.TextColor = [0 0 0];

%Calculate the stimulus centers
xPos = resx/4;
yPos = resy/4;
pos_1 = round([xPos yPos]); % left upper
pos_2 = round([3*(xPos) yPos]); % right-upper
pos_3 = round([xPos 3*(yPos)]); % left-lower
pos_4 = round([3*(xPos) 3*(yPos)]); % right-lower
qcenters =[pos_1 ; pos_2 ; pos_3 ; pos_4]; % centers for each quandrant

%for every quandrant, calculate rects (useful for text placement)
q_rects(1,:)=[0 0 qcenters(1,:)*2];
q_rects(2,:)=[resx/2 0 resx qcenters(2,2)*2];
q_rects(3,:)=[0 resy/2 resx/2 resy];
q_rects(4,:)=[resx/2 resy/2 resx resy];

%%%%====Frequency Tagging Timecourse===%%%%%
f1 = 60; % frequency for RIFT
cfg.FreqMat = f1;
cfg.WaveShape = 'sin'; %Shape of the waveform 'sin' for sinusoidal, 'square' for square wave

%Get the maximal amount of frames to calculate the timecourse for
ifi = Screen('GetFlipInterval', window); % Query the frame duration
flicker_duration = 1000; % unit in second
max_trialframes = round((flicker_duration)/ifi);%max frames per trial

%Here we define the phase timecourse of the frequencies. To enable trial
%averaging in the time domain, phase will be 0 at t=0 (zeroframe), stimulus (figure)
%onset
frame_mult=12; %every refresh is 12 frames

%Frequency timecourse parameters
cfg.patch_amplitude = 0.5;
cfg.patch_startPhase = 0;
cfg.f_offset = 0;

%initialize the table
freqTable=NaN(length(cfg.FreqMat),(max_trialframes*frame_mult));

for f = 1:length(cfg.FreqMat)
    patch_frequency = cfg.FreqMat(f);
    patch_angFreq = 2 * pi * patch_frequency;
    start_time = (1)*-1;
    frametime = start_time:ifi/frame_mult:(max_trialframes*frame_mult)*(ifi/frame_mult)+start_time;
    frametime = frametime(1:max_trialframes*frame_mult);
    %sinusoidal
    freqTable(f,:)= cfg.patch_amplitude * sin(patch_angFreq * frametime + cfg.patch_startPhase) + cfg.patch_amplitude + cfg.f_offset;
end

%%%%%%%===Propixx initialization=====%%%%%%%%%%
% Setup Propixx 1440 Hz
if  ~cfg.debugmode  
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram', 5); % 2 for 480, 5 for 1440 Hz, 0 for normal
    Datapixx('RegWrRd');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%===================== EXPERIMENTAL LOOP =======================%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up Text parameters
Screen('TextFont', window, cfg.TextFont);
Screen('TextSize', window, cfg.TextSize);
Screen('TextStyle',window, cfg.TextStyle);

%%%set frequency table
cur_freqTable = freqTable;

%%%DRAW STIMULI
%%% get the sentence and its x-coordinate in this trial
Words = [{'welcome'},{'to'},{'this'},{'study'},{'.'}]; %% words in the present sentences
NumWord = length(Words);
TxtXcoord = zeros(1,NumWord); %% the x-coordinate of the first word
WordWid = zeros(1,NumWord);
WordHit = zeros(1,NumWord);
for www = 1:NumWord
    textmp = Words{www};
    WordWid(www) = RectWidth(Screen('TextBounds',window,textmp));
    WordHit(www) = RectHeight(Screen('TextBounds',window,textmp));
    if www ~= 1
        TxtXcoord(www) = TxtXcoord(www-1)+WordWid(www-1)+WordSpace;
    end
end
TxtXcoord = TxtXcoord + WordStart;
MostHeight = mode(WordHit);
% Compute bounding box of textstring:
TextureBox1 = q_rects(1,:);

%%%%% draw sentence in the offscreen window
%%all words in offscreen 1
woff1 = Screen('OpenOffscreenwindow', window, [cfg.ScrBgc 0],TextureBox1);
Screen('TextFont', woff1, cfg.TextFont);
Screen('TextSize', woff1, cfg.TextSize);
Screen('TextStyle',woff1, cfg.TextStyle);
for www = 1:NumWord
    Screen('DrawText', woff1, Words{www},TxtXcoord(www),0.5*TextureBox1(4)+0.25*MostHeight,[1 1 1], [], 1);
end
%%gaussian mask of the small flickering patch underlying target word 1 in offscreen 2
id_wrd = 3; % index for the flickering target word -- 'this'
bcgheight = round(TxtXcoord(id_wrd+1)-TxtXcoord(id_wrd));
msy = round(bcgheight)/2;
woff2 = Screen('OpenOffscreenwindow', window, [cfg.ScrBgc 0],TextureBox1);
rectwidth = round(TxtXcoord(id_wrd+1)-TxtXcoord(id_wrd));
rectheight = bcgheight;
x_start = round(TxtXcoord(id_wrd)-0.5*WordSpace);
y_start = round(0.5*TextureBox1(4)-rectheight/2);
x_coords = [x_start (x_start+rectwidth)];
y_coords = [y_start (y_start+rectheight)];
Screen('FillRect', woff2,[1 1 1],[x_coords(1) y_coords(1) x_coords(end) y_coords(end)]);
% We create a Luminance+Alpha matrix for use as transparency mask:
% Layer 1 (Luminance) is filled with luminance value 'gray' of the background.
ms = round(1.5*WordSpace + WordWid(id_wrd))/2;
transLayer = 2;
%%%% square mask
[x,y] = meshgrid(-ms:ms, -msy:msy);
maskblob = uint8(ones(2*msy+1, 2*ms+1, transLayer) * 128);
% Layer 2 (Transparency aka Alpha) is filled with gaussian transparency mask.
xsd = ms/1.2; %bigger than 1.2, more concentrated blob, less smoothing area
ysd = msy/1.2;
maskblob(:,:,transLayer) = uint8(round(255 - exp(-((x/xsd).^2)-((y/ysd).^2))*255));
% Build a single transparency mask texture in offscreen 3
woff3 = Screen('MakeTexture', window, maskblob);
wordcenter_x = (x_coords(1)+x_coords(end))/2;
wordcenter_y = (y_coords(1)+y_coords(end))/2;
rect4offscreen = [wordcenter_x-ms  wordcenter_y-msy  wordcenter_x+ms  wordcenter_y+msy];

%%%% introduction
KbReleaseWait;  %%%  make sure no key unreleased in debug
%welcome screen at first trial
for q = 1 :4
    %for some weird ass reason, this only works with supressed output arguments [~,~,~]
    [~,~,~] = DrawFormattedText(window, 'Press  Any  Key  To  END trial !', 'center', 'center',cfg.TextColor,[],[],[],[],[],q_rects(q,:));
end
vbl = Screen('Flip', window);
KbWait; %% waiting for key pressing

%%%%%========= frames loops in each trial
j = 1; %% index of the frames to change the frequency table in each frame
fff = 1;
while 1  %% word frames during one trial
    %%% put sentences onscreen
    for q = 1:4 %for all quadrants
        %%% flicking bkg of word1
        colortmp2 = [cur_freqTable(1,(((j-1)*12)+q)), cur_freqTable(1,(((j-1)*12)+q+4)),  cur_freqTable(1,(((j-1)*12)+q+8))];
        Screen('DrawTexture', window, woff2,[], q_rects(q,:),0, [], 1, colortmp2);
        Screen('DrawTexture', window, woff3,[], rect4offscreen+[q_rects(q,1) q_rects(q,2) q_rects(q,1) q_rects(q,2)]);
        %%% draw all words
        colortmp1 = cfg.TextColor;
        Screen('DrawTexture', window, woff1,[], q_rects(q,:),0, [], 1, colortmp1);
    end
    %%% flip the frame
    [vbl] = Screen('Flip', window, vbl + 0.5 * ifi);
 
    %CHECK RESPONSES
    if KbCheck
        cleanup(cfg);
        break;
    end
    j = j + 1;
    fff = fff + 1;
end

%%% close offscreens
Screen('Close', woff1)
Screen('Close', woff2)
Screen('Close', woff3)

ShowCursor;
%set propixx to normal state
if ~cfg.debugmode || cfg  .DataPixxOnly
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');
    Datapixx('close');
end
toc
%Close screen
Screen('CloseAll');
end

function [] = cleanup(cfg)
%Return propixx to normal state
if  ~cfg.debugmode
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');
    Datapixx('close');
end

%lower priority
if ~cfg.debugmode
    Priority(0);
end

%close screen
Screen('CloseAll');
%ListenChar(0);
ShowCursor;
%throw warning due to prematurely aborted experiment56
warning('Experiment aborted');
end

function pix = usrDeg2Pix(degree,cfg)
ProjectorResolution_y = 0.5*cfg.resy; %% in rapid mode
pix = tan(degree/2/180*pi)*cfg.dist/cfg.height*2*ProjectorResolution_y;
end