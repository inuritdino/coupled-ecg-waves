function [qi, ri, ti, te, ii] = get_coupled_QRTwaves(ann,qi,ri,ti,te,sig,Fs,tm)
  %% Remove uncoupled Q-, R-, and T-waves from WFDB annotations.
  %%  [qi, ri, ti, te] = GET_COUPLED_WAVES(ann, qi, ri, ti, te, sig,
  %%                                       Fs, tm)
  %% INPUT:
  %%  ann - annotation sample number
  %%  qi - indices of Q-waves (e.g. as qi = find((ant == '(') & (num == 1));)
  %%  ri - indices of R-waves (e.g. as ri = find(ant == 'N');)
  %%  ti - indices of T-wave peaks (e.g. as ti = find(ant == 't');)
  %%  te - indeces of T-wave ends (e.g. as ti2 = find((ant == ')') & (num == 2));)
  %%  NOTE: ann, ant and num are assumed to be corresponding outputs of
  %%  WFDB's RDANN as in [ann,ant,ast,chn,num,com] = rdann(...);
  %%  sig (optional) - signal name in WFDB notation for debugging
  %%  purposes only (Takes too long for 24-hour signal as it calls rdsamp(...))
  %%  If Fs and tm are supplied sig, Fs, and tm are the output from
  %%  WFDB RDSAMP as [sig, Fs, tm] = rdsamp(...);
  %%
  %% OUTPUT:
  %%  The similar indices (qi, ri, ti, te) but filtered to be
  %%  corresponding to each other (coupled/synchronised).
    
  debug = 1;
  if(nargin < 7)
    Fs = [];
    tm = [];
    if(nargin < 6)
      debug = 0;
    end
  end
  more off; %% for printing the progress...
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% VER 3 (straightforward implementation, non-Google way)
  %% Find timings (in num of samples) of the annotations
  tri = double(ann(ri));
  tqi = double(ann(qi));
  tti = double(ann(ti));
  tte = double(ann(te));

  %% Thresholds for Q-waves
  al = 0.3;
  bet = 0.01;
  %% Pointers to be conditionally moved
  i = 1;% R pointer
  j = 1;% Q pointer
  m = 1;% T-peak pointer
  n = 1;% T-end pointer
  k = 1;% coupled counter
  p = 1;% interruption counter
  print_count = 1;
  %% Global upper limits for pointers
  N = length(ri);
  M = length(qi);
  P = length(ti);
  O = length(te);
  %% Output synchronized/coupled indices
  qidx = zeros(1,length(qi));
  ridx = zeros(1,length(ri));
  tidx = zeros(1,length(ti));
  tedx = zeros(1,length(te));
  %% Continuous segment flag
  continuous_flag = false;
  %% Interruption point indices
  iidx = zeros(1,length(ri));
  %% Irregular beat flag
  irregular = false;
  %% Main loop until we reach the end of the signal
  while((i <= N) && (j < M) && (m < P) && (n < O))
    %% debug
    ## if k > 5
    ##   return;
    %% end
    %% Small conditioning on the 1st and last beats
    if (i == 1)
      iprev = i+1;
      inext = i+1;
    elseif (i == N)
      iprev = i-1;
      inext = i-1;
    else
      iprev = i-1;
      inext = i+1;
    end
    %% Regular beat...
    %% if( ( tqi(j) <= tri(i) + bet*(abs(tri(inext)-tri(i))) ) && ...
    %% 	( tqi(j) >= tri(i) -  al*(abs(tri(i)-tri(iprev))) )  && ...
    %% 	( tqi(j+1) <= tri(i) + bet*(abs(tri(inext)-tri(i))) ) && ...
    %% 	( tqi(j+1) >= tri(i) -  al*(abs(tri(i)-tri(iprev))) )  && ...
    %%     ((tti(m) > tri(i)) && (tti(m) < tri(inext))) && ...
    %%     ((tti(m+1) > tri(i)) && (tti(m+1) < tri(inext))) && ...
    %%     ((tte(n) > tri(i)) && (tte(n) < tri(inext))) && ...
    %% 	((tte(n+1) > tri(i)) && (tte(n+1) < tri(inext))) )
    if( ( tqi(j) <= tri(i) + bet*(abs(tri(inext)-tri(i))) ) && ...
	( tqi(j) >= tri(i) -  al*(abs(tri(i)-tri(iprev))) )  && ...
        ((tti(m) > tri(i)) && (tti(m) < tri(inext))) && ...
        ((tte(n) > tri(i)) && (tte(n) < tri(inext))) )
      qidx(k) = j;
      ridx(k) = i;
      tidx(k) = m;
      %% fprintf("Tp: %.0f\n",tti(m));
      tedx(k) = n;
      if ~continuous_flag, continuous_flag = true; end
      %% fprintf("Regular beat\n");
      if irregular
	%% fprintf("Reg: Qi%.0f,Ri%.0f,Tpi%.0f,Tei%.0f,R(i+1)%.0f\n",tri(i),...
	%% 	tqi(j),tti(m),tte(n),tri(i+1));
      end
      irregular = false;
      k = k + 1;
    else%% Irregular beat...
      %% fprintf("Irreg: Qi%.0f,Ri%.0f,Tpi%.0f,Tei%.0f,R(i+1)%.0f\n",tqi(j),tri(i),tti(m),tte(n),tri(i+1));
      %% fprintf("\tQi <= Ri+bet*(...): %.0f <= %.0f\n",tqi(j),tri(i)+bet*(abs(tri(inext)-tri(i))));
      %% fprintf("\tQi >= Ri-al*(...): %.0f >= %.0f\n",tqi(j),tri(i)-al*(abs(tri(i)-tri(iprev))));
      %% fprintf("\tTpi > Ri & Tpi < R(i+1): %.0f > %.0f & %.0f < %.0f\n",tti(m),tri(i),tti(m),tri(inext));
      %% fprintf("\tTei > Ri & Tei < R(i+1): %.0f > %.0f & %.0f < %.0f\n",tte(n),tri(i),tte(n),tri(inext));
      if irregular
	%% fprintf("Warning: two irregular beats\n");
	%% return;
      end
      if ~irregular
	irregular = true;
      end
      %% Designate the end of continuous segment, i.e. interruption
      %% point
      if continuous_flag
	%% fprintf("%.0f continuous stop\n",tri(ri(ridx(k-1))));
	iidx(p) = k-1; p = p + 1;
	continuous_flag = false;
      end
      %%Move R-, T- and Q-pointers forward depending on conditions,
      %%until the coupled waves are found
      while( ...
	  ( (i < N) && (j <= M) && (m <= P) && (n <= O) ) && ...
	  ( ~((tqi(j) <= tri(i) + bet*(abs(tri(inext)-tri(i)))) && ...
	      (tqi(j) >= tri(i) - al*(abs(tri(i)-tri(iprev)))) && ...
	      (tti(m) > tri(i) && tti(m) < tri(inext)) && ...
              (tte(n) > tri(i) && tte(n) < tri(inext)) ) ) ...
	)
	%% Missed Q-wave ==> skip the R-wave (the beat)
	if( tqi(j) > tri(i) + bet*(tri(inext)-tri(i)) )
	  i = i + 1;
	  inext = i+1;
	  iprev = i-1;
	end
	%% Missed R-wave ==> move Q-pointer forward
	if( tqi(j) < tri(i) -  al*(tri(i)-tri(iprev)) )
	  j = j + 1;
	end
	%% Missed T-wave ==> skip the R-wave (the heartbeat)
	if( (tti(m) >= tri(inext)) && (tte(n) >= tri(inext)) )
	  i = i + 1;
	  inext = i+1;
	  iprev = i-1;
	elseif( tte(n) >= tri(inext) ) %% only T-end is missing
	  i = i + 1;
	  inext = i+1;
	  iprev = i-1;
	  m = m + 1;
	elseif( tti(m) >= tri(inext) ) %% only T-peak is missing
	  i = i + 1;
	  inext = i+1;
	  iprev = i-1;
	  n = n + 1;
	end
	%% Missed R-wave ==> move T-peak pointer forward
	if( tti(m) <= tri(i) )
	  m = m + 1;
	end
	%% Missed R-wave ==> move T-end pointer forward
	if( tte(n) <= tri(i) )
	  n = n + 1;
	end
      end
      %% If aligned all waves and not reached the end of signal, store
      %% the index of the aligned waves
      if( (i < N) && (j <= M) && (m <= P) && (n <= O) )
	qidx(k) = j;
	ridx(k) = i;
	tidx(k) = m;
	tedx(k) = n;
	if ~continuous_flag, continuous_flag = true; end
	k = k + 1;
      end
    end
    %% Move all indices forward
    i = i + 1;
    j = j + 1;
    m = m + 1;
    n = n + 1;

    if( (i/N)*100 >= 10*print_count)
      fprintf("%.0f%%  ",(i/N)*100);
      print_count = print_count + 1;
    end
  end
  fprintf("Done.\n");
  %% Remove positions that were not occupied
  qidx = qidx(qidx ~= 0);
  ridx = ridx(ridx ~= 0);
  tidx = tidx(tidx ~= 0);
  tedx = tedx(tedx ~= 0);
  iidx = iidx(iidx ~= 0);
  %% Update the wave indices
  qi = qi(qidx);
  ri = ri(ridx);
  ti = ti(tidx);
  te = te(tedx);
  ii = ri(iidx);% interruption points

  %% Some debugging
  if debug
    %% Reading sample takes a lot of time...
    if isempty(Fs)
      fprintf("Reading signal...");
      [sig,Fs,tm] = rdsamp(sig);
      fprintf("Done.\n");
    end
    plot(tm,sig(:,1));
    hold on;
    plot(tm(ann(ri)),sig(ann(ri),1),'or');
    plot(tm(ann(qi)),sig(ann(qi),1),'og');
    plot(tm(ann(ti)),sig(ann(ti),1),'om');
    plot(tm(ann(te)),sig(ann(te),1),'ok');
    hold off;
    drawnow();
    n_frames = 10;
    frame_delay = 3.0;% in sec
    rnd_frame = randi(floor(max(tm(ann(ri)))/30),1,n_frames);
    for l = rnd_frame
      xlim([(l-1)*30 l*30]);
      id = find((tm(ann(ri)) > (l-1)*30) & (tm(ann(ri)) < l*30));
      for p = 1:length(id)
	text(tm(ann(ri(id(p)))),sig(ann(ri(id(p))),1),num2str(p));
	text(tm(ann(qi(id(p)))),sig(ann(qi(id(p))),1),num2str(p));
	text(tm(ann(ti(id(p)))),sig(ann(ti(id(p))),1),num2str(p));
	text(tm(ann(te(id(p)))),sig(ann(te(id(p))),1),num2str(p));
      end
      title(['Random frame ' num2str(find(rnd_frame == l)) '/' num2str(n_frames)]);
      pause(frame_delay);
    end
  end



  %% VER 2
  %% %% Remove unordered Q and R in the beginning
  %% while (ri(1) < qi(1))
  %%   ri = ri(2:end);
  %% end
  %% %% Form the difference
  %% maxlen = min(length(qi),length(ri));
  %% ri = ri(1:maxlen);
  %% qi = qi(1:maxlen);
  %% D = ri - qi;
  %%
  %% %% Form the diff of D
  %% Dd = diff(D);
  %%
  %% %% Ignore places where Dd is restored in consecutive positions:
  %% %% e.g. Dd = 0 0 -n n 0 0 these shifts indicate detection errors,
  %% %% not missing annotations
  %% Ipos = find(Dd > 0);
  %% Ineg = find(Dd < 0);
  %% J = find(diff(abs(Dd)) == 0);
  %% J = [J J+1];%% we need two consecutive numbers to ignore
  %% %% Drop unpaired R- and Q-waves
  %% ri(setdiff(Ineg,J)) = [];
  %% qi(setdiff(Ipos,J)) = [];

  
end
