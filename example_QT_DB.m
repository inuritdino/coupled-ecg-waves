%%% Simple demo script (NOTE: one must have WFDB Toolbox of PhysioNet)

%% Prefix of the PhysioNet database
prefix = 'qtdb/';

%% Take a sample and read the annotations (already prepared in the database)
record = 'sel100';
fprintf("Reading the sample %s from QT DB...",record);
[ann,ant,ast,chn,num,com] = rdann([prefix record],'pu0');

%% Read the sampling requency and, possibly, the signal
[sig,Fs,tm] = rdsamp([prefix record]);
fprintf("Done.\n");

%% Determine the location of the crucial waves
ri0 = find(ant == 'N'); % R-wave
qi0 = find((ant == '(') & (num == 1)); % Q-wave
ti0 = find(ant == 't'); % T-wave peak
te0 = find((ant == ')') & (num == 2));% T-wave end

%% Remove uncoupled waves
fprintf("Removing uncoupled waves...\n");
[qi, ri, ti, te, ii, ni] = get_coupled_QRTwaves(ann,Fs,qi0,ri0,ti0,te0);

fprintf("Number of the interruptions: %d\n",length(ni));
fprintf("Number of beats in the interruptions: %d\n",unique(ni));
fprintf("Check this visually on the plot...\n");

%% Plot the signal along with the detection
plot(tm,sig(:,1));hold on;
plot(tm(ann(ri)),sig(ann(ri),1),'or');% R-waves
plot(tm(ann(qi)),sig(ann(qi),1),'og');% Q-waves
plot(tm(ann(ti)),sig(ann(ti),1),'om');% T-wave peaks
plot(tm(ann(te)),sig(ann(te),1),'ok');% T-wave ends
plot(tm(ann(ii)),sig(ann(ii),1),'xk','linewidth',2);% Interruption points
legend('Signal','R-wave','Q-wave','T-wave peak','T-wave end',...
       'Last continuous R-wave')
%% text(tm(ann(ii))+0.1,sig(ann(ii),1),strsplit(sprintf("%d Interrupt Beats  ",ni),"  ")(1:end-1));

%% Iterate over the interruption points
for n = 1:length(ii)
  xlim([tm(ann(ii(n)))-5 tm(ann(ii(n)))+5]);
  pause(3.0);
end




