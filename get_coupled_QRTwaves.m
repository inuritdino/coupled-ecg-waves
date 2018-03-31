function [qi, ri, ti, te, ii, nuri] = get_coupled_QRTwaves(ann,fs,qi,ri,ti,te)
  %% Remove uncoupled Q-, R-, and T-waves from WFDB annotations. 
  %% USAGE:
  %% [qi, ri, ti, te, ii, nuri] = GET_COUPLED_WAVES(ann, fs, qi, ri,
  %% ti, te)
  %% INPUT:
  %% ann - annotation sample number
  %% fs - sampling frequency of the signal
  %% qi - indices of Q-waves (e.g. as qi = find((ant == '(') & (num == %% 1));)
  %% ri - indices of R-waves (e.g. as ri = find(ant == 'N');)
  %% ti - indices of T-wave peaks (e.g. as ti = find(ant == 't');)
  %% te - indices of T-wave ends (e.g. as te = find((ant == ')') & (num == 2));)
  %% NOTE: ann, ant and num are assumed to be corresponding outputs of
  %% WFDB's RDANN as in [ann,ant,ast,chn,num,com] = rdann(...);
  %%
  %% OUTPUT: The similar indices (qi, ri, ti, te) but filtered to be
  %% corresponding to each other (coupled/synchronised). Additionally,
  %% the interruption point indices are returned (ii). (ii) are the
  %% indices of the R-waves where the continuous segment interruptions
  %% occur. (nuri) is the number of hear beats that are in the
  %% interruption regions (i.e. (nuri-1) R-waves).
    
  more off; %% for printing the progress...
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Find timings (in num of samples) of the annotations
  tri = double(ann(ri))./fs;
  tqi = double(ann(qi))./fs;
  tti = double(ann(ti))./fs;
  tte = double(ann(te))./fs;

  %% Thresholds for Q-waves
  al = 0.3;
  bet = 0.01;
  %% Thresholds for RR checks
  RR_low_trsh = 0.3;
  RR_up_trsh = 1.5;
  %% Pointers to be conditionally moved
  i = 1;% R pointer
  j = 1;% Q pointer
  m = 1;% T-peak pointer
  n = 1;% T-end pointer
  k = 1;% coupled counter
  p = 1;% interruption point counter
  print_count = 1;
  %% Flags
  %% Double Q/T-wave flags (2+ Q/T-waves determined for the same beat)
  double_wave = false;
  %% Continuous flag
  continuous = true;
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
  iidx = zeros(1,length(ri));% interruption points
  nuri = zeros(1,length(ri));% num of R-waves in interruptions
  %% Main loop until we reach the end of the signal
  while((i < N) && (j <= M) && (m <= P) && (n <= O))
    %% Small conditioning on the 1st and last beats
    if (i == 1)
      iprev = i+1;
      inext = i+1;
    else
      iprev = i-1;
      inext = i+1;
    end
    %% Regular beat...
    if( ( tqi(j) <= tri(i) + bet*(abs(tri(inext)-tri(i))) ) && ...
    	( tqi(j) >= tri(i) -  al*(abs(tri(i)-tri(iprev))) )  && ...
    	( (tti(m) > tri(i)) && (tti(m) < tri(inext)) ) && ...
        ( (tte(n) > tri(i)) && (tte(n) < tri(inext)) ) )
      qidx(k) = j;
      ridx(k) = i;
      tidx(k) = m;
      tedx(k) = n;
      k = k + 1;
      if ~continuous, continuous = true; end
    else%% Irregular beat...
      %%Move R-, T- and Q-pointers forward depending on conditions,
      %%until the coupled waves are found
      double_wave = false;
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
	  if( i > 1 )
	    if( (tri(i) - tri(ridx(k-1))) > RR_low_trsh && ...
		(tri(i) - tri(ridx(k-1))) < RR_up_trsh )
	      %% Double-Q beat...
	      double_wave = true;
	    end
	  end
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
	  if( i > 1 )
	    if( (tri(i) - tri(ridx(k-1))) > RR_low_trsh && ...
		(tri(i) - tri(ridx(k-1))) < RR_up_trsh )
	      %% Double T-wave...
	      double_wave = true;
	    end
	  end
	  m = m + 1;
	end
	%% Missed R-wave ==> move T-end pointer forward
	if( tte(n) <= tri(i) )
	  if( i > 1 )
	    if( (tri(i) - tri(ridx(k-1))) > RR_low_trsh && ...
		(tri(i) - tri(ridx(k-1))) < RR_up_trsh )
	      %% Double T-wave...
	      double_wave = true;
	    end
	  end
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
	k = k + 1;
	if ~double_wave & continuous
	  continuous = false;
	  iidx(p) = ridx(k - 2);% last regular/continuous beat
	  nuri(p) = i - ridx(k - 2);
	  p = p + 1;
	end
	## double_wave = false;% reset to default
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
  nuri = nuri(nuri ~= 0);
  %% Update the wave indices
  ii = ri(iidx);
  qi = qi(qidx);
  ri = ri(ridx);
  ti = ti(tidx);
  te = te(tedx);

  %% Some debugging
  %% if debug
  %%   %% Reading sample takes a lot of time...
  %%   if isempty(Fs)
  %%     fprintf("Reading signal...");
  %%     [sig,Fs,tm] = rdsamp(sig);
  %%     fprintf("Done.\n");
  %%   end
  %%   plot(tm,sig(:,1));
  %%   hold on;
  %%   plot(tm(ann(ri)),sig(ann(ri),1),'or');
  %%   plot(tm(ann(qi)),sig(ann(qi),1),'og');
  %%   plot(tm(ann(ti)),sig(ann(ti),1),'om');
  %%   plot(tm(ann(te)),sig(ann(te),1),'ok');
  %%   hold off;
  %%   drawnow();
  %%   n_frames = 10;
  %%   frame_delay = 3.0;% in sec
  %%   rnd_frame = randi(floor(max(tm(ann(ri)))/30),1,n_frames);
  %%   for l = rnd_frame
  %%     xlim([(l-1)*30 l*30]);
  %%     id = find((tm(ann(ri)) > (l-1)*30) & (tm(ann(ri)) < l*30));
  %%     for p = 1:length(id)
  %% 	text(tm(ann(ri(id(p)))),sig(ann(ri(id(p))),1),num2str(p));
  %% 	text(tm(ann(qi(id(p)))),sig(ann(qi(id(p))),1),num2str(p));
  %% 	text(tm(ann(ti(id(p)))),sig(ann(ti(id(p))),1),num2str(p));
  %% 	text(tm(ann(te(id(p)))),sig(ann(te(id(p))),1),num2str(p));
  %%     end
  %%     title(['Random frame ' num2str(find(rnd_frame == l)) '/' num2str(n_frames)]);
  %%     pause(frame_delay);
  %%   end
  %% end



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
