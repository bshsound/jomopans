function [gaps,gap_in,file_time,size_mbytes,size_min]=bias_gaps20180920(spec_path,instrument)

% call in Matlab (only use this script with Matlab version 2007b or later)
% [gaps,gap_in,file_time,size_mbytes]=bias_gaps20140710;
%
%
% Simple script for BIAS/LIFE/ENV/SE/841 pre-processing data analysis. 
% bias_gaps.m is small program for preprocessing of BIAS/LIFE/ENV/SE/841 data files
% and is a script for use by the members in the BIAS/LIFE/ENV/SE/841 project.
% BIAS/LIFE/ENV/SE/841 are not responsible for the results and the use of this script.
% Only the users are responsible for the use of this script and the results 
% in a appropriate way in the BIAS/LIFE/ENV/SE/841 project.
%
% output:
% gaps = gaps between files in seconds
% gap_in = gaps in samples within file in samples
% file_time = starting time for each file
% size_mbytes = size in Mb of each file
% dsg_files = 1 (DSG files) or 2 (SM2M files)
%
% Leif Persson revised for BIAS/LIFE/ENV/SE/841 20140305 version 1.0
% Leif Persson modified for matching other scripts in BIAS/LIFE/ENV/SE/841 20140710
% version 2
% Leif Persson BIAS/LIFE/ENV/SE/841, error in testing sm2m files corrected version
% bias_gaps20140711
% Leif Persson BIAS/LIFE/ENV/SE/841 20141203 can handle recovered DSG files and original DSG files start
% with 0000.DSG... if data_3=3;
% Emilia 20181003 modified and removed a lot of code concerning DSG since they do not exist any more...
%-------------------------------------------
% data_s=2;    %Data type dsg files nonrecovered or recovered (1) or SM2M (2) (3) DSG odd/original files, i.e. starting with 0000.DSG...;
% Set the direcory name

% checking for filetype SM2M, DSG, soundtrap

% if data_s==3
% test4=dir(fullfile(spec_path,'/0*.wav')); %count the original DSG files
% else
% test1=dir(fullfile(spec_path,'/DSG*.wav')); %count the recovered DSG files
% test2=dir(fullfile(spec_path,'/cDSG*.wav')); %count the unrecovered DSG files
% test3=dir(fullfile(spec_path,'/*.wav')); %count the SM2M files
% end
% 

% if data_s==1 % checking foor DSG files
% 
% if isempty(test1)==1 && isempty(test2)==1 
%     warning([' No DSG or cDSG files in Path ' spec_path])
%     if isempty(test3)==0
%     error([' SM2M files in Path ' spec_path ' Please change filetype to (1)'])
%     else
%         return
%     end
% end
% 
% elseif data_s==2 % checking for SM2M files   
%     if isempty(test3)==1
%     warning([' No SM2M files in Path ' spec_path])
%     end
%     if isempty(test1)==0 || isempty(test2)==0
%     error([' DSG files in Path ' spec_path ' Please change filetype to (2)'])      
%     end
% 
% elseif data_s==3 % checking for SM2M files   
%     if isempty(test4)==1
%     warning([' No original DSG files in Path ' spec_path])
%     end
% else
%     error(['Wrong data file specified: DSG(1) or SM2M(2) or original DSG(3) files'])
%     return    
% end



find_files_in_folder

loopsize=max(size(aaa));

if isempty(aaa)==1
    error([' Path ' spec_path ' are without wav files'])
    return
end

disp(['Total numbers of wav files: ' num2str(loopsize)  ' in ' spec_path])

size_min=zeros(1,loopsize);
gaps=zeros(1,loopsize); %gaps between files
gap_in=zeros(1,loopsize); %gaps within files
size_mbytes=gaps;
file_time=zeros(6,loopsize);
gap_old=0;



for loop_aaa=1:loopsize
    size_mbytes(loop_aaa)=aaa(loop_aaa).bytes/1000000;
    info=audioinfo([spec_path aaa(loop_aaa).name]);
    size_min(loop_aaa)=info.TotalSamples/info.SampleRate/60; %size of file in minutes
    filename=info.Filename; %filename
    filename_length=length(filename);

    date1=get_datenum_from_filename(filename,date_format);
    len_in_sec1=info.TotalSamples/info.SampleRate;
    len_in_sec2=info.Duration;
    gap_in(1,loop_aaa)=(len_in_sec1-len_in_sec2)*info.SampleRate;
    gap=size_min(1,loop_aaa)+str2double(datestr(date1,'SS'))/60;

    if gap<60
       gaps(1,loop_aaa)=60-gap;
    elseif gap==60
       gaps(1,loop_aaa)=0;
    elseif gap>60
       if gap_old<60
           gaps(1,loop_aaa)=120-gap-gap_old;
       else
           gaps(1,loop_aaa)=0;
       end
    elseif gap2>120 && gap2<=180   
       if gap2_old<120
           gaps2(1,loop_aaa2)=180-gap2-gap_old2;
       else
           gaps2(1,loop_aaa2)=0;
       end
    end
    file_time(1,loop_aaa)=str2double(datestr(date1,'yyyy'));
    file_time(2,loop_aaa)=str2double(datestr(date1,'mm'));
    file_time(3,loop_aaa)=str2double(datestr(date1,'dd'));
    file_time(4,loop_aaa)=str2double(datestr(date1,'HH'));
    file_time(5,loop_aaa)=str2double(datestr(date1,'MM'));
    file_time(6,loop_aaa)=str2double(datestr(date1,'SS'));
    gap_old=gap;
end

% if dsg_files<2
disp(' Analysis of one file type finished')

tsc=decyear(file_time(1:6,:)');


sc1=1:max(size(size_mbytes));
%    figure(1),plot(sc1,size_mbytes,'ro--');
figure(1),plot(tsc(1:loopsize),size_mbytes,'ro--');
ylabel('Filesize in Mbyte','FontSize',14),xlabel('Decimal year','FontSize',14)
title(['Filesizes, instrument: ' instrument])

%counting the time between files
t_gaps=zeros(1,max(size(tsc)));
t_gaps(1)=0;
for tcheck=2:max(size(tsc))
    t_gaps(tcheck)=24*(datenum(file_time(:,tcheck)')-datenum(file_time(:,tcheck-1)')); % in hour
%    t_gaps(tcheck)=tsc(tcheck)-tsc(tcheck-1);
end


%figure(2),plot(tsc(1:loopsize),t_gaps,'ro--')
%xlabel('Decimal year','FontSize',14),ylabel('Gaps between files in hour','FontSize',14)

figure(2),plot(tsc(1:loopsize),gaps,'ro')
xlabel('Decimal year','FontSize',14),ylabel('Gaps between files in minutes','FontSize',14)
title(['Gaps between recorded BIAS files, instrument: ' instrument])

figure(3),plot(tsc(1:loopsize),gap_in,'ro--')
xlabel('Decimal year','FontSize',14),ylabel('Gaps within files in minutes','FontSize',14)
title(['Gaps within recorded BIAS files, instrument: ' instrument])

[n_files,b_files]=hist(size_mbytes);
figure(5)
bar(b_files,n_files)
xlabel('Filesize Mbyte','FontSize',14),ylabel('Number of files','FontSize',14)
[n_gaps,b_gaps]=hist(gaps);
title(['Filesizes, instrument: ' instrument])

figure(6)
bar(b_gaps,n_gaps)
xlabel('Gaps between recorded BIAS files [min]','FontSize',14),ylabel('Number of files','FontSize',14)
title(['Gaps between recorded BIAS files [min], instrument: ' instrument])





function [decyr] = decyear(date)

% function [decyr] = decyear(date)  
%
% This Matlab function converts dates in year,month,day into 
% year+fraction of yr. Takes into account leap years.
% date is a matrix with dates as yr, mn, dy , etc 
% It returns a n-vector (decyr) with results
%
% If you need days after first of january of each year do the following
%        days = (decyr - date(:,cy))*365;
%  where cy is the column in matrix date corresponding to years
%
%  ------------------                         
% Last modification 6/95, A.Allmann
% 
% transformation including seconds added, Samuel Neukomm, 3/12/04

%disp('This is /src/decyear.m');

mday= [0,31,59,90,120,151,181,212,243,273,304,334]';%cumulative days in one year
mdayl=[0,31,60,91,121,152,182,213,244,274,305,335]'; %leapyear
l=1:length(date(:,1));
if length(date(1,:))==3 % yr,mon,day!!!
    decyr = mday(date(l,2))/365+(date(l,3)-1)/365 + floor(date(l,1));
elseif length(date(1,:))==5 % yr,mon,day,hr,min !!!
    decyr = mday(date(l,2))/365+((date(l,3)-1) + date(l,4)/24 + date(l,5)/1440)/365 + floor(date(l,1));
elseif length(date(1,:))==6 % yr,mon,day,hr,min,sec !!!
    decyr = mday(date(l,2))/365+((date(l,3)-1) + date(l,4)/24 + date(l,5)/1440+date(l,6)/86400)/365 + floor(date(l,1));
end
%
% test for leap year 
%
leapy = rem(fix(date(:,1)),4) == 0 & rem(fix(date(:,1)),100) ~= 0 | rem(fix(date(:,1)),400) == 0;
ones = find(leapy);
if length(ones) >= 1
    if length(date(1,:))==3
        decyr(ones) = (mdayl(date(ones,2))+(date(ones,3)-1))/366 + floor(date(leapy,1));
    elseif length(date(1,:))==5
        decyr(ones)=(mdayl(date(ones,2))+(date(ones,3)-1)+date(ones,4)/24+date(ones,5)/1440)/366 + floor(date(ones,1));
    elseif length(date(1,:))==6
        decyr(ones)=(mdayl(date(ones,2))+(date(ones,3)-1)+date(ones,4)/24+date(ones,5)/1440+date(ones,6)/86400)/366 + floor(date(ones,1));
    end
end


