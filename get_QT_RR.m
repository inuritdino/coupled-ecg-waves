function Xout = get_QT_RR(ann,ant,num,Fs,outfn)
%% 

  %% Options
  RRthr = 1.5; % sec, basic filtering for consecutive beats

  %%
  if nargin < 5
    outfn = 0;
    if nargin < 4
      Fs = 1.0;
    end
  end

  %% Find the location indices of the waves
  ri = find(ant == 'N'); % R-wave
  qi = find((ant == '(') & (num == 1)); % Q-wave
  ti = find(ant == 't'); % T-wave peak
  te = find((ant == ')') & (num == 2));% T-wave end

  %% remove uncoupled waves
  [qi, ri, ti, te, ii, ni] = get_coupled_QRTwaves(ann,Fs,qi,ri,ti,te);

  %% Form QT and RR intervals
  RR = diff(ann(ri)/Fs);
  %% RR0 = diff(ann(ri)/Fs);
  %% I = find(RR0 < RRthr);% very basic filtering
  %% RR = RR0(I);
  %% NOTE: RRthr automatically assures the next beat is within the
  %% reach
  Qwaves = ann(qi);
  Twaves = ann(te);
  QT = (Twaves - Qwaves)./Fs;
  QTstart = QT(1:end-1);
  QTend = QT(2:end);

  Xout = [QTstart RR QTend].*1000;% make in msec
  %% Output
  if outfn ~= 0
    dlmwrite(outfn,Xout,' ');
  
end
