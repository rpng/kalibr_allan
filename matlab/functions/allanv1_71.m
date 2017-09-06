function [retval, s, errorb] = allanv1_71(data,tau,name)

% Compute the Allan deviation for a set of time-domain frequency data
% [RETVAL, S, ERRORB] = ALLAN(DATA,TAU,NAME)
% DATA should be a struct and have the following fields:
%  DATA.freq    The frequency measurements in Hz
%  DATA.rate or DATA.time
%               The sampling rate in Hertz (DATA.rate) or a timestamp for
%               each measurement in seconds (DATA.time). Computation is
%               shorter when the rate is known, but if the rate is
%               inaccurate, then the Allan plot will be skewed.
%               DATA.rate is used if both fields are present.
%               If DATA.rate == 0, then the timestamps are used.
% TAU is an array of tau values for computing Allan deviation.
% NAME is a label that is added to the plot titles.
%
% RETVAL is the array of Allan deviation values at each TAU.
% S is an optional output of other statistical measures of the data (mean, std, etc).
% ERRORB is an optional output containing the error estimates for a 1-sigma
%   confidence interval. Error bars are plotted as vertical lines at each point
%   (MATLAB-style error bars cannot be used due to a known bug in MATLAB R14SP3).
%
% Example:
%
% To compute the Allan deviation for the data in the variable "lt":
% >> lt
% lt = 
%     freq: [1x86400 double]
%     rate: 0.50
%
% Use:
%
% >> ad = allan(lt,[1 10 100],'lt data');
%
% The Allan deviation will be computed and plotted at tau = 1,10,100 seconds.
%  1-sigma confidence intervals will be indicated by vertical lines.
%
% Notes:
%  No pre-processing of the data is performed. 
%  For rate-based data, AD is computed only for tau values greater than the
%   minimum time between samples and less than the half the total time. For
%   time-stamped data, only tau values greater than the maximum gap between
%   samples and less than half the total time are used.
%  The calculation for rate-based data is *much* faster than for time-stamp
%   data. You may wish to run the rate-based calculation first, then
%   compare with time-stamp-based. Often the differences are insignificant.
%  To plot the "tau bins", uncomment the code at the beginning of the
%   "plot" section (last section of code, search for "TAUBIN").
%  This function has been validated using the test data from NBS Monograph
%   140 and the 1000-point test data set given by Riley [1].
%   If you have other validation results, please let me know!
%
% For more information, see:
% [1] W. J. Riley, "Addendum to a test suite for the calculation of time domain
%  frequency stability," presented at IEEE Frequency Control Symposium,
%  1996.
% Available on the web:
%  http://www.ieee-uffc.org/freqcontrol/paper1ht.html
%
%
% M.A. Hopcroft
%      hopcroft at mems stanford edu
%
%
% MH APR2008
% v1.61  improve error handling, plotting
%        fix bug in regular data calc for high-rate data
%        fix bug in timestamp data calc for large starting gap
%         (thanks to C. B. Ruiz for identifying these bugs)
%        uses timestamps for DATA.rate=0
%        progress indicator for large timestamp data processing

versionstr = 'allan v1.71';
% FCz OCT2008
% v1.71 'lookfor' gives now useful comments; script and sample data are
%       availabie on www.nbi.dk/~czerwin/files/allan.zip
% v1.7  Improve program performance by mainly predefining matrices outside
%       of loops (avoiding memory allocation within loops); no changes to manual
%
% MH JUN2007
% v1.54 Improve data plotting and optional bin plotting
% MH FEB2007
% v1.5  use difference from median for plotting
%       added MAD calculation for outlier detection
% MH JAN2007
% v1.48 plotting typos fixes
% MH DEC2006
% v1.46 hack to plot error bars
% v1.44 further validation (Riley 1000-pt)
%       plot mean and std
% MH NOV2006
% v1.42 typo fix comments
% v1.4  fix irregular rate algorithm
%       irregular algorithm rejects tau less than max gap in time data
%       validate both algorithms using test data from NBS Monograph 140
% v1.3  fix time calc if data.time not present
%       add error bars (not possible due to bug in MATLAB R14SP3)
%       remove offset calculation
% v1.24 improve feedback
% MH SEP2006
% v1.22 updated comments
% v1.2  errors and warnings
% v1.1  handle irregular interval data

% Suggestion for further improvement: The limiting performance factor for the timestamp
% case is the slow performance of 'find()'. Here might logical indexing be
% an appropriate solution.
% Pre-allocation of f, fa and fs does not improve performance
% significantly.


% defaults
if nargin < 3, name=''; end
if nargin < 2, tau=[1 10]; end %one could in principle construct a nice TAU matrix which incorporates knowledge about rate, halftime, etc.

fprintf(1,'allan: %s\n\n',versionstr);

%% Basic statistical tests on the data set
s.numpoints=length(data.freq);
s.max=max(data.freq);
s.min=min(data.freq);
s.mean=mean(data.freq);
s.median=median(data.freq);
if isfield(data,'time')
    s.linear=polyfit(data.time,data.freq,1);
elseif isfield(data,'rate')
    s.linear=polyfit(1/data.rate:1/data.rate:length(data.freq)/data.rate,data.freq,1);
else
    error('Either "time" or "rate" must be present in DATA. See "help allan" for details. [err1]');
end
s.std=std(data.freq);

if nargout < 2, disp(s); end


dfreq=data.freq;
% scale to median for plotting
mfreq=dfreq-s.median;

%pre-defining sm and sme here to prevent pre-allocation within loop
sm=zeros(1, length(tau)); sme=zeros(1, length(tau));

% Screen for outliers using 5x Median Absolute deviation (MAD) criteria
MAD = median(abs(dfreq-s.median)/0.6745);

%%%%
% Two cases - regular or irregular data

%% Regular - Is there a regular interval between measurements?
if isfield(data,'rate') && data.rate > 0 % if there is a regular interval
    fprintf(1, 'allan: regular data (rate = %g Hz)\n',data.rate);
    
    tmstep = 1/data.rate;
    
    if isfield(data,'time')
        % adjust time to remove any starting gap
        dtime=data.time-data.time(1)+mean(diff(data.time));        
        if (data.rate - 1/mean(diff(dtime))) > 1e-6
            fprintf(1,'allan: WARNING: data.rate (%f Hz) does not match recorded data rate (%f Hz)\n',data.rate,1/mean(diff(dtime)));
        end
        fprintf(1,'allan: Length of time record: %f sec\n',dtime(end)-dtime(1));
        % find halfway point
        halftime = fix((dtime(end)-dtime(1))/2);
    else
        halftime = fix(tmstep*length(data.freq)/2);
    end
    fprintf(1, 'allan: max. tau value: %g sec. (time/2)\n',halftime);

    % time axis data
    dtime=tmstep:tmstep:length(dfreq)*tmstep;
    
    
    % truncate tau to appropriate values
    tau = tau(find(tau >= tmstep & tau <= halftime));
    % size(tau)
    fac = round(tau/tmstep);
    
    fprintf(1,'allan: calculating Allan deviation...\n');    

    % calculate the Allan deviation for each value of tau
    k=0; tic;
    for i = tau
        k=k+1;
        %% Instead of defining a matrix which allocated lateron, it is
        %% pre-allocated here and used at #2#
        %f=[];
        f= zeros( fac(k), floor(length(data.freq)/fac(k)) );
        
        % truncate frequency set to an even multiple of this tau value
        freq=dfreq(1:end-rem(length(dfreq),fac(k)));
        length(freq);
        % group the measurements into columns by tau 
        for j=1:fac(k)
            f(j,:)=freq(j:fac(k):end);  % #2#
        end
        % save the binning points for plotting
        fs(k,1:length(freq)/fac(k))=fac(k):fac(k):length(freq); fval{k}=mean(f,1);
        % average in each "tau group"
        fa=mean(f,1);
        % first finite difference
        fd=diff(fa);
        % calculate Allan deviation for this tau
        m=length(fa);
        sm(k)=sqrt(0.5/(m-1)*(sum(fd.^2)));

        % estimate error bars
        sme(k)=sm(k)/sqrt(m);
        
        fprintf(1,'%d ',i);
        
    end
    fprintf(1,'\n'); toc
    % plot the frequency results
    %figure
    %plot(dtime,mfreq,'.b');

    % string for plot title
    name=[name ' (' num2str(data.rate) ' Hz)'];
        
    
%% Irregular data, no fixed interval    
elseif isfield(data,'time')
    % the interval between measurements is irregular
    %  so we must group the data by time
    
    %the two following matrices don't have to be generated, +2 because of
    %beginning and end bin
    time=zeros(1, 2*(length(data.time)+2)); freq=zeros(1, 2*(length(data.freq)+2));
    
    % adjust time to remove any starting gap or zero
    dtime=data.time-data.time(1)+mean(diff(data.time));
    
    fprintf(1, 'allan: irregular data (no rate).\n');
    fprintf(1,'allan: Length of time record: %f sec\n',dtime(end)-dtime(1));
    fprintf(1, '       Average rate: %g Hz (%g sec/measurement); Max. gap: %g sec at position %d\n',...
        1/mean(diff(dtime)),mean(diff(dtime)),max(diff(dtime)),find(diff(dtime)==max(diff(dtime))));
    

    % truncate tau to appropriate values
    % find halfway point
    halftime = fix((dtime(end)-dtime(1))/2);
    tau = tau(find(tau >= max(diff(dtime)) & tau <= halftime));
    if isempty(tau)
        error('allan: ERROR: no appropriate tau values (> %g s, < %g s)\n',max(diff(dtime)),fix(dtime(end)/2));
    end
        
    
    fprintf(1,'allan: calculating Allan deviation...');

    k=0; tic;
    for i = tau
        fprintf(1,'\n%d ',i);
        
        k=k+1; fa=[]; 
        f=[];
        jk=0; km=0;
        
        % truncate data set to an even multiple of this tau value
        freq=dfreq(find(dtime <= dtime(end)-rem(dtime(end),i)));
        time=dtime(find(dtime <= dtime(end)-rem(dtime(end),i)));
        
        % break up the data into groups of tau length
        while i*km < time(end)
            km=km+1;
                        
            % progress bar
            if rem(km,100)==0, fprintf(1,'.'); end
            if rem(km,1000)==0, fprintf(1,'%g/%g\n',km,round(time(end)/i)); end

            f = freq(find(i*(km-1) < time & time <= i*km));
            if ~isnan(any(f)) && any(f) ~= 0
                fa(km)=mean(f);
            else
                fa(km)=0;
            end
            
            % save the binning points for plotting
            if find(time <= i*km) > 0
                fs(k,km)=max(time(find(time <= i*km)));
            else
                fs(k,km)=0;
            end
            fval{k}=fa;
            
        end

        % first finite difference of the averaged results
        fd=diff(fa);
        % calculate Allan deviation for this tau
        m=length(fa);
        sm(k)=sqrt(0.5/(m-1)*(sum(fd.^2)));
        
        % estimate error bars
        sme(k)=sm(k)/sqrt(m);
        

    end
    fprintf(1,'\n'); toc
    
    % string for plot title
    name=[name ' (timestamp)'];

%     % plot the frequency results
%     if ~isempty(time) && ~isempty(freq)
%         figure
%         plot(dtime,mfreq,'.b');
%         name=['(time) ' name];
%     else
%         fprintf(1,'allan: ERROR: no appropriate tau values (> %g s)\n',max(diff(data.time)));
%     end

else
    error('allan: WARNING: no DATA.rate or DATA.time! Type "help allan" for more information. [err2]');
end


%%%%%%%%
%% plotting

% plot the frequency data, centered on median
figure
plot(dtime,mfreq,'.b');
hold on;


    %% Optional plot    TAUBIN
    % plot the time divisions
%     [rfs,cfs]=size(fs);
%     colororder=get(gca,'ColorOrder');
%     axis tight; ap=axis; kc=2;
%     for j=1:rfs
%         kc=kc+1; if rem(kc,length(colororder))==1, kc=2; end
%         plot the tau division boundaries
%         for b=1:max(find(fs(j,:)));
%             plot([fs(j,b) fs(j,b)],[ap(3)*1.1 ap(4)*1.1],'-','Color',colororder(kc,:));
%             if b == 1
%                 plot([dtime(1) fs(j,b)],[fval{j}(b)-s.median fval{j}(b)-s.median],'-','Color',colororder(kc,:),'LineWidth',4);
%             else
%                 plot([fs(j,b-1) fs(j,b)],[fval{j}(b)-s.median fval{j}(b)-s.median],'-','Color',colororder(kc,:),'LineWidth',4);
%             end
%         end
%     end
%     axis auto
    %% End optional plot
    

fx = xlim;
% plot([fx(1) fx(2)],[s.median s.median],'-k');
plot([fx(1) fx(2)],[0 0],':k');
% plot([fx(1) fx(2)],[s.mean+s.std s.mean+s.std],'-r');
% plot([fx(1) fx(2)],[s.mean-s.std s.mean-s.std],'-r');
% plot([fx(1) fx(2)],[3*MAD 3*MAD],'--r');
% plot([fx(1) fx(2)],[-3*MAD -3*MAD],'--r');
plot([fx(1) fx(2)],[5*MAD 5*MAD],'-r');
plot([fx(1) fx(2)],[-5*MAD -5*MAD],'-r');
title(['Data: ' name],'FontSize',16,'FontName','Arial');
%set(get(gca,'Title'),'Interpreter','none');
xlabel('time (sec)','FontSize',14,'FontName','Arial');
ylabel('Data - Median [frequency]','FontSize',14,'FontName','Arial');
set(gca,'FontSize',14,'FontName','Arial');

if ne(isempty(sm), 1)
    figure
    plotlinewidth=2;
    
    
    % plot with loglog scale
    loglog(tau, sm, '.-b', 'LineWidth', plotlinewidth, 'MarkerSize', 24);
    errorbar(tau,sm,sme,'.-b', 'Markersize', 18); set(gca,'XScale','log'); set(gca,'YScale','log');
    % plot with y log scale
    % semilogy(tau, sm, '.-b', 'LineWidth', plotlinewidth, 'MarkerSize', 24);
    % errorbar(tau,sm,sme,'.-b', 'Markersize', 18); set(gca,'YScale','log');
    % plot with x log scale
    % semilogx(tau,sm,'.-b','LineWidth',plotlinewidth,'MarkerSize',24);
    % errorbar(tau,sm,sme,'.-b', 'Markersize', 18); set(gca,'XScale','log');
    
    % There's still a bug that screws up the error bars on a loglog and semilog plots.
    % When this is fixed, set(...) becomes unnecessary
     
    grid on;
    title(['Allan Deviation: ' name],'FontSize',16,'FontName','Arial');
    %set(get(gca,'Title'),'Interpreter','none');
    xlabel('\tau (sec)','FontSize',14,'FontName','Arial');
    ylabel('\sigma_y(\tau) (Hz)','FontSize',14,'FontName','Arial');
    set(gca,'FontSize',14,'FontName','Arial');
    % expand the x axis a little bit so that the errors bars look nice
    adax = axis;
    axis([adax(1)*0.9 adax(2)*1.1 adax(3) adax(4)]);
else
    fprintf(1,'allan: WARNING: no values calculated. Check that TAU > 1/DATA.rate\n');
    fprintf(1,'Type "help allan" for more information.\n\n');
end
        
retval = sm;
errorb = sme;

return
