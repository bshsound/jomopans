function [s]=create_jomopans_s_struct_v201905(spec_path,instr_settings,proc_settings)
% Program originally developed for data analysis in BIAS/LIFE/ENV/SE/841 
% by Leif Persson 2014.
% It is a small program for power spectrum estimation for measurement data.
% Jomopans are not responsible for the results and the use of this script.
%
% [s]=create_jomopans_s_struct_v2(spec_path,instr_settings,proc_settings)
%   instr_settings require information on 
%         gain = 0 dB
%         voltage correction=1, (i.e. what 1 V means on the normalised wav-file), 
%         channel = 1
%         type of instrument = e.g. Soundtrap or DSG or rtsys
%         sensitivity = a matrix in the form [freq sens].
%
%   proc_settings require information on if windowing will be applied
%       proc_setting.wind = 0 or 1. If 1 = hanning window
%       proc_settings.overlap = 0 (however, this does not work so skip it
%       proc_settings.zeropadding = 0, integer number adds zeros with value of
%       integer number * length of vector x.
%
%  Following scripts are needed in able to run the script:
%     "find_files_in_folder" 
%     "prcentile" 
%
% output:
%         s = struct for
%         SPL_1smean_time_domain'
%         SPL_1smean_third_octave_bands = All octave bands from 10 Hz - 20 kHz 
%         ters_center = third octave band center frequency
%         ters_upper_lim = third octave band upper limit frequency
%         ters_lowerlim = third octave band lower limit frequency
%         PSD_values (1s average)
%         Frequency_scale
%         file_size
%         clipping_statistic_max
%         clipping_statistic_min
%         number_of_files
%         file_time
%         file_time_datenum
%         time_datenum_second
%         Number_of_files_in_folder
% %
% Emilia Lalander 20190521: Change script for JOMOPANS data
%-------------------------------------------
% Processing settings:
wind=proc_settings.wind; % if a Hann window is used change to wind=1, wind=0 mean no window, i.e. boxcar window;
noverlap=proc_settings.noverlap; % if overlap is used set noverlap to the number of overlapping samples
max_level=proc_settings.max_level;%0.9; % maximum level on the sensors in output wavefile, used for testing of clipping
min_level=proc_settings.min_level;%-0.9; % minimum level on the sensors in output wavefile, used for testing of clipping
%spec_out=1; % flag for if spec_out=1, the all 1sec spectras will be saved in the struct 
clipping_test=1; % flag for test (0) means no test, if clipping test is done change to (1)
zero_padding=1;
qq=0; % Uses all files in the folder

% Instrumet settings
instrument=instr_settings.instrument;
gain=instr_settings.gain;
volt_correction=instr_settings.volt_correction;
channel=instr_settings.channel; 
sensitivity_matrix=instr_settings.sensitivity;


%-------------------------------------------------------------------------------------------------
disp('Monthly data specified')

if strcmp(instrument(1:5),'Hydro')
    nrFilesInDir=dir(fullfile(spec_path,'/*.mat')); %count the files
    load([spec_path '\' nrFilesInDir(1).name]); % load a test data fil
    x=data;
    fs=instr_settings.fs;
else
    nrFilesInDir=dir(fullfile(spec_path,'/*.wav')); %count the files
    [x,fs]=audioread([spec_path nrFilesInDir(1).name]); % load a test data fil
end

loopsize_normal=max(size(nrFilesInDir));
tot_loopsize=loopsize_normal;

find_files_in_folder

% Quality control ---------------------------------------------------
%x  =  x(~isnan(x(:,1))); % Remove missing observations indicated by NaN's.

if isempty(x)
   error('NotEnoughData',...
         'Input samples in first file has no valid data.');
end
%--------------------------------------------------------------------
if strcmp(instrument,'Soundtrap') 
    x=x(3*fs:end,channel);
elseif strcmp(instrument,'DSG-LS1') 
    x=x(1*fs:end,channel);
else
    x=x(:,channel);
end

N=fs; % Number of samples in 1 sec of data
x=x-mean(x); %removing the DC offset, this is the global mean of one file
% [f_0,~]=fftspectrum(x(1:fs),N,fs,wind,noverlap); % first spectra for specify the bins
wind_size=1;

xx=x(1:N);
[f1,~]=spec_test(x(1:N),N,fs,wind,noverlap); % first spectra for specify the bins

% For lowest frequencies -> zeropadding *1
xx_0=[xx;zeros(length(xx)*zero_padding,1)];
N_0=length(xx_0);
[f_0,~]=spec_test(xx_0,N_0,fs,wind,noverlap);

% For lowest frequencies -> zeropadding *4
xx_00=[xx;zeros(length(xx)*zero_padding*4,1)];
N_00=length(xx_00);
[f_00,~]=spec_test(xx_00,N_00,fs,wind,noverlap);

sensitivityI=load_sensitivity(sensitivity_matrix,f1);
sensitivityI_0=load_sensitivity(sensitivity_matrix,f_0);
sensitivityI_00=load_sensitivity(sensitivity_matrix,f_00);

% Define tersband
n=-20:13;
ters_center_tot=1000*10.^(n/10);
ters_upperlim=ters_center_tot*10^(1/20); % Tersband upper limit
ters_lowerlim=ters_center_tot*10^(-1/20); % Tersband upper limit

nr_frek_tot=max(size(ters_center_tot));
low_freq_pad=100;   %Frequency to zeropadd 4 times N
middle_freq_pad=200;    %Frequency to zeropadd 1 times N

for fx=1:nr_frek_tot
    if ters_center_tot(fx)< middle_freq_pad  %Above 300 Hz --> no zeropadding
        if ters_center_tot(fx)< low_freq_pad
            octave_bands_idx_00{fx}=find(f_00>ters_lowerlim(fx) & f_00<=ters_upperlim(fx)); % define the bins tersband fx
            if length(sensitivityI) > 1
                mSens(fx)=-mean(sensitivityI_00(octave_bands_idx_00{fx}));
            else
                mSens(fx)=sensitivityI;
            end
        else
            octave_bands_idx_0{fx}=find(f_0>ters_lowerlim(fx) & f_0<=ters_upperlim(fx)); % define the bins tersband fx
            if length(sensitivityI) > 1
                mSens(fx)=-mean(sensitivityI_0(octave_bands_idx_0{fx}));
            else
                mSens(fx)=sensitivityI;
            end
        end
    else
        octave_bands_idx{fx}=find(f1>ters_lowerlim(fx) & f1<=ters_upperlim(fx)); % define the bins tersband fx        
        if length(sensitivityI) > 1
            mSens(fx)=-mean(sensitivityI(octave_bands_idx{fx}));
        else
            mSens(fx)=sensitivityI;
        end
    end
end

% Predefine variables
% ptot=[];
% ptot2=[];
nr_loop1=fix(max(size(x))/N); % total file length in number of 1s block

test_file_size=zeros(tot_loopsize,1);
if clipping_test==1
    clipping_max=zeros(nr_loop1*tot_loopsize,1);
    clipping_min=zeros(nr_loop1*tot_loopsize,1);
else
    clipping_max=0;
    clipping_min=0;
end
value_count=1;
empty=0; 
nan_full=0;
%-------------------------------------------------------------------------
%------------------ loop over all specified BIAS files -------------------
%-------------------------------------------------------------------------
%0-matrix
SPL_1smean_time_domain=zeros(1,tot_loopsize*nr_loop1);    %All bands
psd_all=zeros(nr_frek_tot,tot_loopsize*nr_loop1-2);    %All bands
SPL_1smean_third_octave_bands=zeros(nr_frek_tot,tot_loopsize*nr_loop1-2);    %All bands

file_time=zeros(6,tot_loopsize);
file_time_datenum=zeros(1,tot_loopsize*nr_loop1);
time_datenum_second=zeros(1,tot_loopsize*nr_loop1);
for loop_aaa=1:tot_loopsize % start the loop
%profile on
    disp([ 'File ' num2str(loop_aaa) ' of totally ' num2str(tot_loopsize)])
    if strcmp(instrument(1:5),'Hydro')
        load([spec_path '\' nrFilesInDir(loop_aaa).name]); % load a test data fil
        x=data;
        fs=instr_settings.fs;
    else
        [x,N]=audioread([spec_path '\' nrFilesInDir(loop_aaa).name]); % load the data monthly file 
    end
    %[spec_path nrFilesInDir(1).name]
    if isempty(x)==1
        error([spec_path  '\' nrFilesInDir(loop_aaa).name ' is empty'])
        return
    end

    if strcmp(instrument,'Soundtrap') 
        x=x(3*fs:end,channel);
    elseif strcmp(instrument,'DSG-LS1') 
        x=x(1*fs:end,channel);
    else
        x=x(:,channel); 
    end
    x=x*volt_correction;


    % extract the time from filename1
    if strcmp(instrument(1:5),'Hydro')
        filename=nrFilesInDir(1).name; %count the files
    else
        info=audioinfo([spec_path nrFilesInDir(loop_aaa).name]); %read the filename
        filename=info.Filename; %filename
    end


    date1=get_datenum_from_filename(filename,date_format);  
    %Date_format is defined in "find_files_in_folders.m"
    year=str2double(datestr(date1,'yyyy'));
    mon=str2double(datestr(date1,'mm'));
    day=str2double(datestr(date1,'dd'));
    sec=str2double(datestr(date1,'SS'));
    minute=str2double(datestr(date1,'MM'));
    hour=str2double(datestr(date1,'HH'));

    file_time(1,loop_aaa)=year;
    file_time(2,loop_aaa)=mon;
    file_time(3,loop_aaa)=day;
    file_time(4,loop_aaa)=hour;
    file_time(5,loop_aaa)=minute;
    file_time(6,loop_aaa)=sec;

    file_time_datenum=date1;
% Simple Quality control ---------------------------------------------------

    if isempty(x)
        warning('create_jomopans_s_struct:NotEnoughData',...
             'Input sample X is empty and has no valid data.');
         empty=1;
    else    
        x  =  x(~isnan(x)); % Remove missing observations indicated by NaN's.
        if isempty(x)
           warning('create_jomopans_s_struct:NotEnoughData',...
                 'Input sample X has no valid data (all NaN''s).');
             nan_full=1;
        end
    end

    if nan_full==1 || empty==1 
        error(['create_jomopans_s_struct:Terminated remove file ' nrFilesInDir(loop_aaa).name])
    end


    ind1=1:(N*wind_size); % 1 second of data
    %setting up the variables
    
    nr_loop_new=fix(max(size(x))/N); % total file length in number of 1s block


    % if file size is smaller than normal
    if loop_aaa>1
        if nr_loop_new<nr_loop1
            %making files larger
            disp(['Filesize of file ' num2str(loop_aaa) ' is smaller'])
            nr_loop1=nr_loop_new;
        elseif nr_loop_new>nr_loop1
            disp(['Filesize of file ' num2str(loop_aaa) ' is larger'])
            nr_loop1=nr_loop_new;
%         elseif nr_loop_new==nr_loop1
%     %         test_file=0;
        end    
        test_file_size(loop_aaa)=nr_loop1;  %saving the file sizes
    else
        test_file_size(loop_aaa)=nr_loop_new;  %saving the file sizes
    end


    x=x-mean(x); %remove the DC offset, this is the global mean


    for loop1=1:nr_loop1
        % Loop per sekund för varje fil
        
        if value_count>nr_loop1*tot_loopsize
            disp('Warning: the file size is increasing with 1s window')
        end
        try
        xx=x(ind1);
        catch
            keyboard
        end
        xx=xx-mean(xx); %removal of possible local offset caused by non-stationarities
        
        %Stop when reaching the end
        if ind1(end)<=max(size(x))
        %-------- testing for clipping --------------------------------
            if clipping_test==1
                [aw_max,~]=find(xx>=max_level);
                [aw_min,~]=find(xx<=min_level);
                if isempty(aw_max)
                    clipping_max(value_count)=0;
                else
                    clipping_max(value_count)=100*max(size(aw_max))/max(size(xx));
                end

                if isempty(aw_min)
                clipping_min(value_count)=0;
                else
                clipping_min(value_count)=100*max(size(aw_min))/max(size(xx));
                end
            end
            %------- end of testing for clipping --------------------------
            

            %Frequency domain calculation
            psd_sum_octaveband=zeros(nr_frek_tot,1);
            SPL_level=zeros(nr_frek_tot,1);
            
            f_0idxStart=find(ters_center_tot<low_freq_pad);
            f_idxStart=find(ters_center_tot<middle_freq_pad);
            for fx=1:nr_frek_tot
                if ters_center_tot(fx) < low_freq_pad
                    if fx==1
                        xx=x(ind1);
                        xx=xx-mean(xx); %removal of possible local offset caused by non-stationarities
                        xx=[xx;zeros(length(xx)*zero_padding*4,1)];
                        N=length(xx);
                        [~,psdx]=spec_test(xx,N,fs,wind,noverlap);
                    end
                    octave_bands_idxx=octave_bands_idx_00{fx};
                elseif ters_center_tot(fx) < middle_freq_pad
                    if fx == f_0idxStart(end)+1
                        xx=x(ind1);
                        xx=xx-mean(xx); %removal of possible local offset caused by non-stationarities
                        xx=[xx;zeros(length(xx)*zero_padding,1)];
                        N=length(xx);
                        [~,psdx]=spec_test(xx,N,fs,wind,noverlap);
                    end
                    octave_bands_idxx=octave_bands_idx_0{fx};
                else
                    if fx == f_idxStart(end)+1
                        xx=x(ind1);
                        xx=xx-mean(xx); %removal of possible local offset caused by non-stationarities
                        N=length(xx);
                        [~,psdx]=spec_test(xx,N,fs,wind,noverlap);
                    end
                    octave_bands_idxx=octave_bands_idx{fx};
                end
                %Total P in tersband
                %p_t=sqrt((psdx.').^2);  %abs(psdx)
                %p_ters_tot(fx)=(sum(abs(psdx(octave_bands_idxx)).^2));

                psd_sum_octaveband(fx)=sum(psdx(octave_bands_idxx)); 
                SPL_level(fx)=mSens(fx)-gain+10*log10(psd_sum_octaveband(fx)); 
            end
            %Time domain calculation
            Lp= 10 * log10(1/(wind_size*fs)*sum(xx.^2))-mean(sensitivityI);
            SPL_1smean_time_domain(value_count)=Lp';    %All bands
            
            psd_all(:,value_count)=psd_sum_octaveband';    %All bands
            SPL_1smean_third_octave_bands(:,value_count)=SPL_level';    %All bands
            %Value_count uppdateras hela tiden så för varje ny timme
            %fortsätter value_count att ticka uppåt...
            time_datenum_second(value_count)=date1+(value_count-1)/(24*3600);


        %--------------- wrong in data ------------------------------------------
        else
            disp('Warning size of datafile reduced')
            clipping_min(value_count)=0;
            clipping_max(value_count)=0;
            SPL_1smean_time_domain(value_count)=0;    %All bands          
            psd_all(:,value_count)=0;    %All bands
            SPL_1smean_third_octave_bands(:,value_count)=0;    %All bands
            time_datenum_second(value_count)=date1+(value_count-1)/(24*3600);
        end
        ind1=ind1+(N*wind_size);
        value_count=value_count+1;
    end % end of one file

    if qq==0
        q_tot=tot_loopsize;
    end
%         profile viewer
end

s = struct('SPL_1smean_time_domain',SPL_1smean_time_domain,...
    'SPL_1smean_third_octave_bands',SPL_1smean_third_octave_bands,...
    'ters_center',ters_center_tot','ters_upper_lim',ters_upperlim',...
    'ters_lowerlim',ters_lowerlim','PSD_values',psd_all,'Frequency_scale',f_0,...
    'file_size',test_file_size,'clipping_statistic_max',clipping_max',...
    'clipping_statistic_min',clipping_min','number_of_files',tot_loopsize,...
    'file_time',file_time,'file_time_datenum',file_time_datenum,...
    'time_datenum_second',time_datenum_second,'Number_of_files_in_folder',q_tot);%,...
%     'SPL_1smean_third_octave_bandsPwelch',SPL_1smean_third_octave_bandsPwel);

disp([' Number of lines in output header ' num2str(round(sum(test_file_size)./20))])
end

function [f,pxx]=spec_test(x,N,fsamp,wind,noverlap)

% minimum version of PSD periodogram Leif Persson 20140607 for BIAS/LIFE/ENV/SE/841 
% with time domain windowing 20140622, only Hann window
if wind==1
    window=hann(N);
     if noverlap==0
         k=1.633;
     else
         k=1.0;
     end
else
    window=ones(N,1);
    k=1.0; 
end

xx=x.*window;

xdft=k*fft(xx);	%Computes the discrete Fourier Transform
xdft=xdft(1:N/2+1);	%Use only one-side of the spektrum

 
psdx=(1/(fsamp*N)).*abs(xdft).^2;	% Divide by number of samples and 
                                    % square it to get PSD.
 
psdx(2:end-1)=2*psdx(2:end-1);	% Due to loss of energy when only looking 
                                % at one-sided spectrum, AC component have 
                                % to be multiplied by 2. PSDX(1)= DC-component

pxx=psdx;
 
f=0:fsamp/length(x):fsamp/2;
end

function sensitivity_2=load_sensitivity(sensitivity,f)
sensitivity_2=zeros(size(f));
size_sens=size(sensitivity(:,:),1);
if size_sens > 1
    if size_sens <= 2
        sensitivity=sensitivity';
    end
    sens_freq=sensitivity(:,1);
    maxfreq=min([max(sens_freq) max(f)]);
    minfreq=max([min(sens_freq) min(f)]);
    if max(max(f)) < max(sens_freq)
        warning('Number of octave bands will be less than specified due to limited sampling frequency')
    end
    freq_idx=find(f<=maxfreq & f >= minfreq);

    %Mixtrar med känsligheten
    sens_f_idx=find(sens_freq<=maxfreq & sens_freq >= minfreq);
    sens_int=interp1(sens_freq(sens_f_idx)',sensitivity(sens_f_idx,2)',f(freq_idx));
    sensitivity_2(1:freq_idx(1)-1)=sensitivity(1,2);
    sensitivity_2(freq_idx)=sens_int;
    sensitivity_2(freq_idx(end)+1:length(f))=sensitivity(end,2);
    % Slutmixtrat
else
    sensitivity_2=sensitivity;
end
end

function filedatenum = get_datenum_from_filename(filename,timeFormat)
%
% function filedatenum = get_datenum_from_filename(filename,timeFormat)
%
% strips off possible (and impossible stuff from filename) 
% to retrieve the Matlabtime from the filename
%

% test how many index the time code is in the files: 
i0 = 1; 
%not_found_start_index = 1; 
gg = []; 

while isempty(gg)
    try
        gg = datenum(filename(i0:end),timeFormat);   
    catch err      
    end
    i0 = i0 + 1;
end
i0 = i0-1; 
filedatenum = datenum(filename(i0:end),timeFormat); 
end

