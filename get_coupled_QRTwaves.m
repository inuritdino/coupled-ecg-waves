function [qi, ri, ti, te, ii, nuri] = get_coupled_QRTwaves(ann,fs,qi,ri,ti,te,varargin)
  %% Remove uncoupled Q-, R-, and T-waves from WFDB annotations.
  %%
  %% USAGE:
  %%
  %%    [qi, ri, ti, te, ii, nuri] =
  %%         GET_COUPLED_QRTWAVES(ann, fs, qi, ri, ti, te, varargin)
  %%
  %% INPUT:
  %%
  %% ann -- annotation sample number
  %% fs -- sampling frequency of the signal
  %% qi -- indices of Q-waves (e.g. as qi = find((ant == '(') & (num == %% 1));)
  %% ri -- indices of R-waves (e.g. as ri = find(ant == 'N');)
  %% ti -- indices of T-wave peaks (e.g. as ti = find(ant == 't');)
  %% te -- indices of T-wave ends (e.g. as te = find((ant == ')') & (num == 2));)
  %% varargin -- see below.
  %%
  %% NOTE: ann, ant and num are assumed to be corresponding outputs of
  %% WFDB's RDANN as in [ann,ant,ast,chn,num,com] = rdann(...);
  %%
  %% OUTPUT: The similar indices (qi, ri, ti, te) but filtered to be
  %% corresponding to each other (coupled/synchronised). Additionally,
  %% the interruption point indices are returned (ii). (ii) are the
  %% indices of the R-waves where the continuous segment interruptions
  %% occur (the last R of the continuous region). (nuri) is the number 
  %% of heart beats that are in the interruption regions (i.e. (nuri-1) 
  %% R-waves).
  %%
  %% NOTE: the interruption points are determined accounting for
  %% double Q/T-waves, that is such annotations that show 2 or more
  %% Q/T-waves for the single R-wave (annotation problem). Namely, no
  %% interruption occurs if the double wave is detected. The double
  %% waves are detected by using RR intervals subject to
  %% thresholds. The thresholds are specified by the VARARGIN:
  %%
  %% VARARGIN:
  %%
  %% get_coupled_QRTwaves(...,'RRlow',Value,...) -- the lower boundary
  %% for RR intervals is set to Value (default: 0.3 sec).
  %%
  %% get_coupled_QRTwaves(...,'RRup',Value,...) -- the upper boundary
  %% for RR intervals is set to Value (default: 1.5 sec).
  %%
  %% Data processing note: 'RRlow' and 'RRup' options impose a
  %% processing/cleaning step to the data and, hence, must be taken
  %% with caution. All other procedures to find the coupled waves are
  %% (semi)qualitative not requiring thresholds, except finding the
  %% contiguous regions. Thus, contiguous region procedure can be
  %% performed by external tools depending on application.
  %%
  %% Other VARARGIN options:
  %%
  %% ...,'MinInterruptBeats',Value,... -- min number of the
  %%     interruption segment beats considered an interruption
  %%     (default: 1)
  %%
  %% ...,'ReturnCells',... -- turns on the cell arrays in the return
  %%     of the function: each cell holds the wave indices
  %%     corresponding to each contiguous region of heart beats.
  %%
  %% ...,'MinWavesACell',Value,... -- sets the min number of waves per
  %%     each contiguous region to return in cells (only when ReturnCells
  %%     is true; default: 1).
    
  more off; %% for printing the progress...

  %% Thresholds for Q-waves alignment to R
  al = 0.3;
  bet = 0.01;
  %% Thresholds for RR checks
  RR_low_trsh = 0.3;
  RR_up_trsh = 1.5;
  %% Min number of the interruption region beats to be considered as
  %% the interruption
  min_interrupt_beats = 1;
  %% Return cell arrays
  return_cell = false;
  %% Min number of waves in a cell (only if return_cell=true)
  min_waves_a_cell = 1;

  vi = 1;
  while( vi <= length(varargin))
    if( strcmpi('RRlow',varargin{vi}) )
      RR_low_trsh = varargin{vi+1};
      vi = vi + 1;
    elseif( strcmpi('RRup',varargin{vi}) )
      RR_up_trsh = varargin{vi+1};
      vi = vi + 1;
    elseif( strcmpi('MinInterruptBeats',varargin{vi}) )
      min_interrupt_beats = varargin{vi+1};
      vi = vi + 1;
    elseif( strcmpi('ReturnCells',varargin{vi}) )
      return_cell = true;
    elseif( strcmpi('MinWavesACell',varargin{vi}) )
      min_waves_a_cell = varargin{vi+1};
      vi = vi + 1;
    else
      fprintf("Warning: VARARGIN argument is unknown.\n");
    end
    vi = vi + 1;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Find timings of the annotations
  tri = double(ann(ri))./fs;
  tqi = double(ann(qi))./fs;
  tti = double(ann(ti))./fs;
  tte = double(ann(te))./fs;

  %% Pointers to be conditionally moved
  i = 1;% R pointer
  j = 1;% Q pointer
  m = 1;% T-peak pointer
  n = 1;% T-end pointer
  k = 1;% coupled counter
  p = 1;% interruption point counter
  print_count = 1;
  %% Double Q/T-wave flag (2+ Q/T-waves determined for the same beat)
  double_wave = false;
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
	%% Missed R-wave ==> move Q-pointer forward
	elseif( tqi(j) < tri(i) -  al*(tri(i)-tri(iprev)) )
	  if( k > 1 )
	    if( (tri(i) - tri(ridx(k-1))) > RR_low_trsh && ...
		(tri(i) - tri(ridx(k-1))) < RR_up_trsh )
	      %% Double-Q beat...
	      double_wave = true;
	    end
	  end
	  j = j + 1;
	%% Missed T-wave ==> skip the R-wave (the heartbeat)
	elseif( (tti(m) >= tri(inext)) && (tte(n) >= tri(inext)) )
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
	%% Missed R-wave ==> move T-peak pointer forward
	elseif( tti(m) <= tri(i) )
	  if( k > 1 )
	    if( (tri(i) - tri(ridx(k-1))) > RR_low_trsh && ...
		(tri(i) - tri(ridx(k-1))) < RR_up_trsh )
	      %% Double T-wave...
	      double_wave = true;
	    end
	  end
	  m = m + 1;
	%% Missed R-wave ==> move T-end pointer forward
	elseif( tte(n) <= tri(i) )
	  if( k > 1 )
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
	%% Get the interruption point
	if ~double_wave & (k > 2)
	  if( (i - ridx(k-2)) >= min_interrupt_beats )
	    iidx(p) = ridx(k - 2);% last continuous beat
	    nuri(p) = i - ridx(k - 2);
	    p = p + 1;
	  end
	end
      end
    end
    %% Move all indices forward
    i = i + 1;
    j = j + 1;
    m = m + 1;
    n = n + 1;

    %% Progress printing...
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

  if return_cell
    %% Sort the indices into a cell array
    ric = cell(1,length(ii)+1);
    qic = cell(1,length(ii)+1);
    tic = cell(1,length(ii)+1);
    tec = cell(1,length(ii)+1);
    for n = 1:(length(ii)+1)
      if n == 1
	I0 = 1;
      else
	I0 = find(ri == ii(n-1),1) + 1;
      end
      if n == (length(ii)+1)
	I1 = length(ri);
      else
	I1 = find(ri == ii(n),1);
      end
      if( (I1 - I0 + 1) >= min_waves_a_cell )
	ric{n} = ri(I0:I1);
	qic{n} = qi(I0:I1);
	tic{n} = ti(I0:I1);
	tec{n} = te(I0:I1);
      end
    end
    %% Form the output
    ri = ric(~cellfun(@isempty,ric));
    qi = qic(~cellfun(@isempty,qic));
    ti = tic(~cellfun(@isempty,tic));
    te = tec(~cellfun(@isempty,tec));
  end
  
end
