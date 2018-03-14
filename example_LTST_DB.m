function [ann,ant,ast,chn,num,com] = example_LTST_DB(varargin)
  %% Example usage of get_coupled_QRTwaves
  %% Do not run it as is! It is a long-long process...
  %% Just take as a template.
  %%
  %% This example annotates (using ecgpuwave) and processes the
  %% annotations of the Long-Term ST Database from PhysioNet (file
  %% RECORDS must be present and contain the list of the LTST DB
  %% signal files.
  
  annotate = 0;

  i = 1;
  while(i <= length(varargin))
    if(strcmpi("annotate",varargin{i}))
      annotate = 1;
    end
    i = i + 1;
  end
  
  prefix = "ltstdb/";
  ## rec = importdata("RECORDS");
  rec = textread("RECORDS","%s");%% textread is more stable in Octave
  more off;

  if annotate
    fprintf("Annotation in progress:\n");
    for i = 1:length(rec)
      fprintf("\t %s (%d/%d)\n",rec{i},i,length(rec));
      ecgpuwave([prefix rec{i}],'ann',[],[],[],1);
    end
  end

  ann = cell(1,length(rec));
  ant = cell(1,length(rec));
  ast = cell(1,length(rec));
  chn = cell(1,length(rec));
  num = cell(1,length(rec));
  com = cell(1,length(rec));
  %% Global RR threshold, used for getting consecutive heartbeats
  RRthr = 1.5;% in sec

  fprintf("Reading annotations:\n");
  for i = 1:length(rec)
    fprintf("\t %s (%d/%d)\n",rec{i},i,length(rec));
    %% Read annotation
    [ann{i},ant{i},ast{i},chn{i},num{i},com{i}] = rdann([prefix rec{i}],'ann');
    %% Read the sampling frequency
    siginfo = wfdbdesc([prefix rec{i}]);
    Fs = double(siginfo(1).SamplingFrequency);%% force double-type
    %% Determine the location of the crucial waves
    ri = find(ant{i} == 'N'); % R-wave
    qi = find((ant{i} == '(') & (num{i} == 1)); % Q-wave
    ti = find(ant{i} == 't'); % T-wave peak
    ti2 = find((ant{i} == ')') & (num{i} == 2));% T-wave end
    %% Remove uncoupled waves
    fprintf("Remove uncoupled waves...\n");
    [qi, ri, ti, ti2] = get_coupled_QRTwaves(ann{i},qi,ri,ti,ti2);
    %% Form first RR's
    RR0 = diff(ann{i}(ri)/Fs);
    %% Get consecutive heart beats by applying global RR threshold
    I = find(RR0 < RRthr);
    %% Get filtered RR
    RR = RR0(I);
    %% Get QTstart
    Qwaves = ann{i}(qi); Twaves = ann{i}(ti2);
    QTstart = (Twaves(I) - Qwaves(I))/Fs;
    %% Get QTend (NOTE: RRthr automatically assures the next beat is
    %% within the reach)
    QTend = (Twaves(I+1) - Qwaves(I+1))/Fs;
    %% Output
    dlmwrite(['data/' rec{i} '.txt'],[QTstart RR QTend],' ');
  end
  
end



