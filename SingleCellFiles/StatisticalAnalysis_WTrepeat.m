%% README
%  File for analysis of WT data (WTrepeat) for Gould, Domijan et al.
%  "Coordination of robust single cell rhythms in the Arabidopsis circadian 
%  clock via spatial waves of gene expression" (BioRxiv)
%  doi https://doi.org/10.1101/208900s
% 
%  by Mirela Domijan (U. of Liverpool)
%%  
name='WTrepeat'; %Experiment to run. 
%% add in relevant paths to the data
addpath([pwd '/Data_singlecell/'])
addpath([pwd '/Data_singlecell/' name '_final_coordinates'])
% For some plots , install Red Blue Colormap package from MATHWORKS from:
% https://uk.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap?focused=5138232&tab=function
% and place it in the Data_singlecell folder. 
addpath([pwd '/redblue/'])
%% is there a folder where to save all the data?
if exist(['Doneanalysis/' name  '/summary'],'file')~=7
    mkdir(['Doneanalysis/' name  '/summary']);
end
%% is there a folder where to save all the figures?
if exist(['Doneanalysis/' name  '/figures'],'file')~=7
    mkdir(['Doneanalysis/' name  '/figures']);
end
%% pull out all the period estimates and the data
[aa, texts]=xlsread('Data_singlecell//WTrepeat_final_coordinates/interpolated/Individual results_1h.xlsx');
aa(1:4, :)=[];
texts(1:4, :)=[];
rowstorem=isnan(aa(:,1));
[rowstorem2]=find(rowstorem==1);
aa(rowstorem2,:)=[];
texts(rowstorem2,:)=[];
count=[];
for j=1: size(texts,1)-1
   if isequal(texts(j,2),texts(j+1,2))~=1
       count=[count j];
   end      
end
clear j
separators=[];
posi=13:20;
for j=1:length(posi)
    xs=csvread(['Data_singlecell/WTrepeat_final_coordinates/interpolated/interpolateddata_WTrepeat_pos' num2str(posi(j)) '_spacing1h.csv']);
    separators=[separators size(xs,2)-1];
end
clear j rowstorem rowstorem2  xs ind x
ss=cumsum(separators); 
separates=ss;
pp=reshape(aa(:,12),3,length(aa(:,12))/3)';
ppgof=reshape(aa(:,34),3,length(aa(:,34))/3)';
for i=1:length(texts(:,11))
    if strcmp(texts(i,11), 'IGNORED')==1
        ignor(i)=1;
    else
        ignor(i)=0;
    end
end
ppignore=reshape(ignor,3,length(ignor)/3)';
separates=[0 separates];
for j=1:length(posi)
    k=posi(j);
    num1=separates(j)+1;
    num2=separates(j+1);
    pos{k}.periods=pp(num1:num2,:);
    pos{k}.periodsgofs=ppgof(num1:num2,:);
    pos{k}.periodsignore=ppignore(num1:num2,:);  
end
clearvars -except pos name
%% determine which positions you are looking at
if strcmp(name,'WT')==1
    posi=1:6;
    sect={ {'Root tip', ''},{'root', '(sect 1)'}, {'root', '(sect 2)'},{'hypocotyl', '(sect 1)'},{'hypocotyl', '(sect 2)'},{'cotyledon', ''}};
elseif strcmp(name,'WTrepeat')==1
    posi=13:20;    
    sect={{'Root tip', ''},{'Root', '(sect 1)'}, {'Root', '(sect 2)'},{'Root', '(sect 3)'}, {'Root', '(sect 4)'} ,{'Hypocotyl', '(sect 1)'},{'Hypocotyl', '(sect 2)'},{'Cotyledon', ''}}; 
end
%% now which positions you will look at 
start_pos=1;
numpos=8;
%% for each file you pull out, write down the positions; the YFP levels and the periods. 
for k= start_pos:numpos
%    %two types of periods: (i) periods(:,1)== FFT unscaled; (ii) periods(:,2)== FFT unscaled RAE; (iii) periods(:,9)== LS acceptance 10%
     i=posi(k); 
     if strcmp(name, 'WTrepeat')==1
         name2=['pos',num2str(i)];
         if i==14 || i==19
             posfile=dir([pwd , '/Data_singlecell/' name '_final_coordinates/',['*_pos_' name2 '*tracked*.csv']]);
         else
             posfile=dir([pwd , '/Data_singlecell/' name '_final_coordinates/',[ name2 '*tracked*.csv']]);
         end
         fid= fopen( posfile.name,'r');
         C=textscan(fid, repmat('%s',1,10), 'delimiter',',', 'CollectOutput',true);
         if i==3
             C{1}(676, :)=[];
         end
         if i==4
             triple=str2double(C{1}(3:end,1:3));
             LL=logical(diff(str2double(C{1}(3:end,7))));
         else
             triple=str2double(C{1}(4:end,1:3));
             LL=logical(diff(str2double(C{1}(4:end,7))));
         end
         fa=find(LL==1);
 % Note on Labels:
 % triple is the x,y,z triple of the coordinates. 
 % LL keeps clear of when the timing changes 
 % fa identifies where the new time point begins
 % res reshapes the positions so every three columns (i.e. x,y,z) you have
 % a new time point i.e, dimensions are no. of coordinates* time points. 
 % num_cells is the number of cells
 % cell_labels is the labels of the cells 
if i==16
     LL(1800:end)=[];
     triple(1801:end,:)=[];
     fa(end)=[];
 end
if numel(unique(diff([0; fa])))==1 && (size(LL,1)-fa(end))== (fa(1)-1)
    res=reshape(triple', fa(1)*3, (length(fa)+1));
    % verify that everything is pulled out correctly and that you have correct cell IDs;
    num_cells=size(res,1)/3;
    cell_labels=str2double(C{1}(4:end,8));
    if i==16
        cell_labels(1801:end)=[];
    end
    %now make sure they are always increasing and OK (quick method):    
    diffy=diff(cell_labels)-1;
    diffywo = diffy(diffy~=0);
    %Note: diffy  keeps track of differences in cell parents (i.e. if
    %same time point, cell parents/labels should be increasing by one)
    %when new time point starts diffywo records how many cells there
    %are recorded in each time point.
    allok(i)=isequal(abs(diffywo),ones(size(res,2)-1,1)*num_cells);
    % allok tells us that thre are right number of cells within each time point 
    data.pos=res';
    pos{i}.x=data.pos(:,1:3:end);
    pos{i}.y=data.pos(:,2:3:end); 
    pos{i}.z=data.pos(:,3:3:end);
end
    % clear everything not needed:
    clear diffy diffywo res  fa fid LL cell_labels triple
    % now put in the YFP files  
    YFPfile=dir([pwd , '/Data_singlecell/'  name '_final_coordinates/',['*mean*' name2 '*YFP_*.csv']]); 
    fid= fopen( YFPfile.name,'r');
    clear C
    C=textscan(fid, repmat('%s',1,num_cells+1), 'delimiter',',', 'CollectOutput',true);
    if i==4
         pos{i}.YFP=str2double(C{1}(3:end,2:end));
         cell_labels2=str2double(C{1}(2,2:end));
    else    
    pos{i}.YFP=str2double(C{1}(4:end,2:end));
    cell_labels2=str2double(C{1}(3,2:end));
    end
    % verify that everything is pulled out correctly and make sure the cells are in increasing order
    % repeat of above just wiht the new YFP file. 
    diffy2=diff(cell_labels2)-1;
    diffywo2 = diffy2(diffy2~=0);
    allok2(i)=isempty(diffywo2);
    end
    clear cell_labels2  fid YFPfile diffy2 diffywo2 num_cells
% time points 
if strcmp(name, 'WT')==1
    Times=xlsread('timestamp_WT.xlsx', 'C8:DR15');
    pos{i}.time= Times(i,:);
elseif strcmp(name, 'WTrepeat')==1
    Times=xlsread('timestamp_WTrepeat.xlsx', 'B8:AG27');
    pos{i}.time= Times(i,2:end);
end
clear data FFTrecs periods
end
clear i k  txt2 posfile bb a Times C name2 name3 
%% removal of doubles in overlapping sections: (these are calculated in Analysis_WTrepeat in dontkeep variable)
rem2=[  5     7    22    25    26    27    28    29    30    31    32    33    34    35  36    40    49    50    51    52    58];
rem3=[ 1     22   ];
rem4=[  3     4     7     9    12    16    18    19    23    24    25    28    33    35    36    37    38    39    40    41    52    53    54    55    57];
rem5=[ 1     2     3     4     6     8     9    11    12    24    25    26    27    28   29    32    36    37    38    45    46    55    65    66    67];
for  kk= start_pos:numpos
    k=posi(kk);
    if kk==2
        pos{k}.periods(rem2,:)=[];
        pos{k}.periodsgofs(rem2,:)=[];
        pos{k}.periodsignore(rem2,:)=[];
        pos{k}.x(:, rem2)=[];
        pos{k}.y(:,rem2)=[];
        pos{k}.z(:, rem2)=[];
        pos{k}.YFP(:,rem2)=[];
    elseif kk==3
        pos{k}.periods(rem3,:)=[];
        pos{k}.periodsgofs(rem3,:)=[];
        pos{k}.periodsignore(rem3,:)=[];
        pos{k}.x(:, rem3)=[];
        pos{k}.y(:,rem3)=[];
        pos{k}.z(:, rem3)=[];
        pos{k}.YFP(:,rem3)=[];
    elseif kk==4 
        pos{k}.periods(rem4,:)=[];
        pos{k}.periodsgofs(rem4,:)=[];
        pos{k}.periodsignore(rem4,:)=[];
        pos{k}.x(:, rem4)=[];
        pos{k}.y(:,rem4)=[];
        pos{k}.z(:, rem4)=[];
        pos{k}.YFP(:,rem4)=[];   
    elseif kk==5
        pos{k}.periods(rem5,:)=[];
        pos{k}.periodsgofs(rem5,:)=[];
        pos{k}.periodsignore(rem5,:)=[];
        pos{k}.x(:, rem5)=[];
        pos{k}.y(:,rem5)=[];
        pos{k}.z(:, rem5)=[];
        pos{k}.YFP(:,rem5)=[];
    end  
end
%% Amplitude
counter=NaN(1,numpos);counter2=NaN(1,numpos);counter3=NaN(1,numpos);counter4=NaN(1,numpos);
counter6=NaN(1,numpos);counter7=NaN(1,numpos);counter5=NaN(1,numpos);
mpkdis=5;
mpkdis2=5;
if exist(['Doneanalysis/' name  '/summary/dist' num2str(mpkdis2) '.mat'], 'file')~=2
    mkdir(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/']);
    gof_cutoff=0.9; 
    gof_cutoff2=1;
    for k= start_pos:numpos
        i=posi(k);
        Amps_for=[]; Amps_back=[];
        Amps=[];
        mkdir(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i)]);
        counter(k)=0; %counter for the number of cells that will be tested
        counter2(k)=0; %counter for the cells that have more than 2 peaks
        counter3(k)=0; %counter for the cells that have more than 2 peaks and where the number max and min correctly determined
        counter4(k)=0; %counter for the cells that have more than 3 peaks and where the number max and min correctly determined
        counter5(k)=0; %counter for the cells that have more than 2 peaks and where the number max and min correctly determined, with trough no=peak no with trough before peak
        counter7(k)=0;  %as counter 5 but trough after peak
        counter6(k)=0; %counter for the cells that have more than 2 peaks and where the number max and min correctly determined, with trough no> peak no       
        if i==16 || i==17 || i==18 || i==19 || i==20
            xs=size(pos{i}.time,2)-1;
        else
        xs=size(pos{i}.time,2);
        end
        for j= 1:size(pos{i}.YFP,2)
        if  pos{i}.periodsgofs(j,1)<gof_cutoff && pos{i}.periodsgofs(j,2)<gof_cutoff && pos{i}.periodsgofs(j,3)<gof_cutoff2  && abs(pos{i}.periods(j,1)- pos{i}.periods(j,3))<2.5 && abs(pos{i}.periods(j,2)- pos{i}.periods(j,3))<2.5 && abs(pos{i}.periods(j,2)- pos{i}.periods(j,1))<2.5 && sum(pos{i}.periodsignore(j,:))==0
            counter(k)=counter(k)+1;
            mpkdis=5;
            if exist(['Doneanalysis/' name  '/summary/dist' num2str(mpkdis2)],'file')~=7
                mkdir(['Doneanalysis/' name  '/summary/' 'dist' num2str(mpkdis2)]);
            end
            % finding the peaks strategies:
            % detrending out the base value and then the following settings
             MinPeakDist=mpkdis2*(pos{i}.time(2)-pos{i}.time(1));
             MinPkHeight=0; 
             MinPkProm=20;
             MinPkPromTr=2;
            % Find the peaks:
            clear dtdata p trend
            p = polyfit(pos{i}.time(1:xs)',pos{i}.YFP(1:xs,j),1);  % find trend as polynomial of degree 'd'
            trend = polyval(p,pos{i}.time(1:xs));
            dtdata = pos{i}.YFP(1:xs,j) - trend';  
            [a,b]= findpeaks(dtdata, pos{i}.time(1:xs),  'MinPeakDistance',MinPeakDist, 'MinPeakHeight', MinPkHeight, 'MinPeakProminence',MinPkProm); 
            [a2,b2]= findpeaks(-dtdata, pos{i}.time(1:xs), 'MinPeakDistance',MinPeakDist, 'MinPeakHeight',MinPkHeight, 'MinPeakProminence',MinPkPromTr);
            pos{i}.index(j)=200;  %this keeps track of cells that will have enough amplitudes for analysis adn tells you how many amplitudes
            if length(b)>2
                pos{i}.index(j)=100;
                counter2(k)=counter2(k)+1;
                if abs(length(b)-length(b2))<2
                    counter3(k)=counter3(k)+1;
                    pos{i}.index(j)=max([length(b) length(b2)]);
                    if length(b2)== length(b)
                        if b2(1)<b(1)
                       counter5(k)=counter5(k)+1; 
                        elseif b2(1)>b(1)
                            counter7(k)=counter7(k)+1; 
                        end
                    elseif length(b2)== length(b)+1
                         counter6(k)=counter6(k)+1; 
                    end
                    if length(b)>3
                     counter4(k)=counter4(k)+1;
                    end
                end         
            end
            % convert b's and b2's into time point values:
            clear bk  bk2  pkTime troughTime
            pkTime=NaN(1,length(b)); troughTime=NaN(1,length(b));
            for s=1: length(b)
                pkTime(s)=find(pos{i}.time==b(s));
            end
            clear s
            for s=1: length(b2)
                troughTime(s)=find(pos{i}.time==b2(s));
            end
            %now replace 
            pkTrueTimes=b; 
            troughTrueTime=b2;          
            if isempty(b)==1
                L1=0;
            else
                L1=length(pkTime);
            end
            if isempty(b2)==1
                L2=0;
            else
                L2=length(troughTime);
            end
            clear b b2 s   
            mm= min([L1 L2]);
            if mm>1 
                haystack=sort([troughTime pkTime], 'ascend');
                clear u l  u2 ul2
                [l, u]=ismember(troughTime, haystack);
                tooManyTroughs=sum(ismember(diff(u), 1) ); 
                [l2, u2]=ismember(pkTime, haystack);
                tooManyPks=sum(ismember(diff(u2), 1) ); 
                clear haystack l u l2 u2
                if tooManyTroughs==0 && tooManyPks==0 && abs(length(tooManyTroughs)-length(tooManyPks))<2
                % need to measure amplitude by two ways trough to peak and peak to trough. 
                    if length(pkTime)==length(troughTime) && length(pkTime)>1% b is peak and b2 is trough
                        if pkTime(1)<troughTime(1) 
                            diffs_forward=abs(pos{i}.YFP(troughTime(1:end-1),j)-pos{i}.YFP(pkTime(2:end),j));
                            diffs_back=abs(pos{i}.YFP(troughTime,j)-pos{i}.YFP(pkTime,j));
                            if length(diffs_forward)>2 
                                pos{i}.val(j)=mean(diffs_forward);
                                pos{i}.WCV_for(j)=std(diffs_forward)/mean(diffs_forward);
                            else
                                pos{i}.val(j)=NaN;
                                pos{i}.WCV_for(j)=NaN;
                            end
                            if length(diffs_back)>2 
                                pos{i}.val2(j)=mean(diffs_back);
                                pos{i}.WCV_back(j)=std(diffs_back)/mean(diffs_back);
                            else
                                pos{i}.val2(j)=NaN;
                                pos{i}.WCV_back(j)=NaN;
                            end
                            Amps_for=[Amps_for diffs_forward' ];
                            Amps_back=[Amps_back diffs_back' ];
                            if  length(diffs_forward)>2 
                                Amps=[Amps diffs_forward'];
                            else
                                if length(diffs_back)>2 
                                    Amps=[Amps diffs_back'];
                                end
                            end 
                            pos{i}.Ampsfor(j, 1:2)=diffs_forward([ 1 end]); 
                            pos{i}.Ampsback(j, 1:2)=diffs_back([ 1 end]); 
                            if length(pkTrueTimes)>3 
                       pos{i}.ppmeans(j)=mean(diff(pkTrueTimes)); %mean period per cell
                       pos{i}.CVp(j)=std(diff(pkTrueTimes))/mean(diff(pkTrueTimes));
                       pos{i}.pkTrueTimesN{j}=pkTrueTimes; 
                   else
                       pos{i}.ppmeans(j)=NaN; %mean period per cell
                       pos{i}.CVp(j)=NaN;
                       pos{i}.pkTrueTimesN{j}=NaN; 
                            end
                         elseif pkTime(1)>troughTime(1)
                            diffs_forward=abs(pos{i}.YFP(troughTime,j)-pos{i}.YFP(pkTime,j));
                            diffs_back=abs(pos{i}.YFP(troughTime(2:end),j)-pos{i}.YFP(pkTime(1:end-1),j));
                            if length(diffs_forward)>2 
                                pos{i}.val(j)=mean(diffs_forward);
                                pos{i}.WCV_for(j)=std(diffs_forward)/mean(diffs_forward);
                            else
                                pos{i}.val(j)=NaN;
                                pos{i}.WCV_for(j)=NaN;
                            end
                            if length(diffs_back)>2 
                                 pos{i}.val2(j)=mean(diffs_back);
                                 pos{i}.WCV_back(j)=std(diffs_back)/mean(diffs_back);
                            else
                                pos{i}.val2(j)=NaN;
                                pos{i}.WCV_back(j)=NaN;
                            end
                            Amps_for=[Amps_for diffs_forward' ];
                            Amps_back=[Amps_back diffs_back' ];
                            if  length(diffs_forward)>2 
                                Amps=[Amps diffs_forward'];
                            else
                                if length(diffs_back)>2 
                                    Amps=[Amps diffs_back'];
                                end
                            end
                            pos{i}.Ampsfor(j, 1:2)=diffs_forward([ 1 end]); 
                            pos{i}.Ampsback(j, 1:2)=diffs_back([ 1 end]); 
                            if length(pkTrueTimes)>3
                                pos{i}.ppmeans(j)=mean(diff(pkTrueTimes)); %mean period per cell
                                pos{i}.CVp(j)=std(diff(pkTrueTimes))/mean(diff(pkTrueTimes));
                                pos{i}.pkTrueTimesN{j}=pkTrueTimes; 
                            else
                                pos{i}.ppmeans(j)=NaN; %mean period per cell
                                pos{i}.CVp(j)=NaN;
                                pos{i}.pkTrueTimesN{j}=NaN; 
                            end                                           
                        end
            elseif length(pkTime)>length(troughTime) && length(troughTime)==length(pkTime)-1 && length(pkTime)>1
                if pkTime(1)<troughTime(1)
                    diffs_forward=abs(pos{i}.YFP(troughTime,j)-pos{i}.YFP(pkTime(2:end),j));
                    diffs_back=abs(pos{i}.YFP(troughTime,j)-pos{i}.YFP(pkTime(1:end-1),j));
                    if length(diffs_forward)>2 
                        pos{i}.val(j)=mean(diffs_forward);
                        pos{i}.WCV_for(j)=std(diffs_forward)/mean(diffs_forward);
                    else
                        pos{i}.val(j)=NaN;
                        pos{i}.WCV_for(j)=NaN;
                    end
                    if length(diffs_back)>2 
                        pos{i}.val2(j)=mean(diffs_back);
                        pos{i}.WCV_back(j)=std(diffs_back)/mean(diffs_back);
                    else
                        pos{i}.val2(j)=NaN;
                        pos{i}.WCV_back(j)=NaN;
                    end
                    Amps_for=[Amps_for diffs_forward' ];
                    Amps_back=[Amps_back diffs_back' ]; 
                    if  length(diffs_forward)>2 
                        Amps=[Amps diffs_forward'];
                    else
                        if length(diffs_back)>2 
                            Amps=[Amps diffs_back'];
                        end
                    end
                    pos{i}.Ampsfor(j, 1:2)=diffs_forward([ 1 end]); 
                    pos{i}.Ampsback(j, 1:2)=diffs_back([ 1 end]);      
                    if length(pkTrueTimes)>3 
                        pos{i}.ppmeans(j)=mean(diff(pkTrueTimes)); %mean period per cell
                        pos{i}.CVp(j)=std(diff(pkTrueTimes))/mean(diff(pkTrueTimes));
                        pos{i}.pkTrueTimesN{j}=pkTrueTimes; 
                    else
                        pos{i}.ppmeans(j)=NaN; %mean period per cell
                        pos{i}.CVp(j)=NaN;
                        pos{i}.pkTrueTimesN{j}=NaN; 
                    end   
                 else
                    pos{i}.val(j)=NaN;
                    pos{i}.val2(j)=NaN; 
                    Amps_for=[Amps_for NaN ];
                    Amps_back=[Amps_back NaN ];
                    Amps=[Amps NaN];
                    pos{i}.Ampsfor(j, 1:2)=NaN; 
                    pos{i}.Ampsback(j, 1:2)=NaN; 
                    pos{i}.WCV_for(j)=NaN;
                    pos{i}.WCV_back(j)=NaN;
                    pos{i}.ppmeans(j)=NaN; %mean period per cell
                    pos{i}.CVp(j)=NaN;
                    pos{i}.pkTrueTimesN{j}=NaN;          
                end
            elseif length(pkTime)<length(troughTime) && length(pkTime)==length(troughTime)-1
                if pkTime(1)>troughTime(1)
                    diffs_forward=abs(pos{i}.YFP(troughTime(1:end-1),j)-pos{i}.YFP(pkTime,j));
                    diffs_back=abs(pos{i}.YFP(troughTime(2:end),j)-pos{i}.YFP(pkTime,j));
                    if length(diffs_forward)>2 
                        pos{i}.val(j)=mean(diffs_forward);
                        pos{i}.WCV_for(j)=std(diffs_forward)/mean(diffs_forward);
                    else
                      pos{i}.val(j)=NaN;
                      pos{i}.WCV_for(j)=NaN;
                    end
                    if length(diffs_back)>2 
                        pos{i}.val2(j)=mean(diffs_back);
                        pos{i}.WCV_back(j)=std(diffs_back)/mean(diffs_back);
                    else
                        pos{i}.val2(j)=NaN;
                        pos{i}.WCV_back(j)=NaN;
                    end 
                    Amps_for=[Amps_for diffs_forward' ];
                    Amps_back=[Amps_back diffs_back' ];
                    if  length(diffs_forward)>2 
                        Amps=[Amps diffs_forward'];
                    else
                        if length(diffs_back)>2 
                            Amps=[Amps diffs_back'];
                        end
                    end
                    pos{i}.Ampsfor(j, 1:2)=diffs_forward([ 1 end]); 
                    pos{i}.Ampsback(j, 1:2)=diffs_back([ 1 end]);            
                   if length(pkTrueTimes)>3
                       pos{i}.ppmeans(j)=mean(diff(pkTrueTimes)); %mean period per cell
                       pos{i}.CVp(j)=std(diff(pkTrueTimes))/mean(diff(pkTrueTimes));
                       pos{i}.pkTrueTimesN{j}=pkTrueTimes; 
                   else
                       pos{i}.ppmeans(j)=NaN; %mean period per cell
                       pos{i}.CVp(j)=NaN;
                       pos{i}.pkTrueTimesN{j}=NaN; 
                   end   
                else
                    pos{i}.val(j)=NaN;
                    pos{i}.val2(j)=NaN; 
                    pos{i}.WCV_for(j)=NaN;
                    pos{i}.WCV_back(j)=NaN;
                    Amps_for=[Amps_for NaN ];
                    Amps_back=[Amps_back NaN ];
                    Amps=[Amps NaN];
                    pos{i}.Ampsfor(j, 1:2)=NaN; 
                    pos{i}.Ampsback(j, 1:2)=NaN;  
                    pos{i}.ppmeans(j)=NaN; %mean period per cell
                    pos{i}.CVp(j)=NaN;
                    pos{i}.pkTrueTimesN{j}=NaN;   
                        end
                    end
                clear pkTimes troughTimes a a2  pkTrueTimes
                else
                pos{i}.val(j)=NaN;
                pos{i}.val2(j)=NaN; 
                pos{i}.WCV_for(j)=NaN;
                pos{i}.WCV_back(j)=NaN;
                Amps_for=[Amps_for NaN ];
                Amps_back=[Amps_back NaN ];
                Amps=[Amps NaN];
                pos{i}.Ampsfor(j, 1:2)=NaN; 
                pos{i}.Ampsback(j, 1:2)=NaN; 
                pos{i}.ppmeans(j)=NaN; %mean period per cell
                pos{i}.CVp(j)=NaN;
                pos{i}.pkTrueTimesN{j}=NaN; 
                end   
            else
                pos{i}.val(j)=NaN;
                pos{i}.val2(j)=NaN; 
                pos{i}.WCV_for(j)=NaN;
                pos{i}.WCV_back(j)=NaN;
                Amps_for=[Amps_for NaN ];
                Amps_back=[Amps_back NaN ]; 
                Amps=[Amps NaN];
                pos{i}.Ampsfor(j, 1:2)=NaN; 
                pos{i}.Ampsback(j, 1:2)=NaN; 
                pos{i}.ppmeans(j)=NaN; %mean period per cell
                pos{i}.CVp(j)=NaN;
                pos{i}.pkTrueTimesN{j}=NaN; 
            end
            figure(1)
            if exist(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/good/'],'file')~=7
                mkdir(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/good/']);
            end
            if exist(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/notenough/'],'file')~=7
                mkdir(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/notenough/']);
            end
            if sum(isnan(dtdata))==xs
                set(gcf,'DefaultLineLineWidth',3)
                plot(pos{i}.time(1:xs),pos{i}.YFP(1:xs,j),'o-', 'Linewidth',3)
                xlabel('Time (h)')
                ylabel('Levels of undetrendable time series')  
                print(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/notenough/Cell_' num2str(j) '.png'],'-dpng','-r0')
            else
                set(gcf,'DefaultLineLineWidth',3)
                findpeaks(dtdata, pos{i}.time(1:xs),  'MinPeakDistance',mpkdis2*(pos{i}.time(2)-pos{i}.time(1)), 'MinPeakHeight', 0, 'MinPeakProminence',MinPkProm)
                [a2,b2]=findpeaks(-dtdata, pos{i}.time(1:xs),   'MinPeakDistance',mpkdis2*(pos{i}.time(2)-pos{i}.time(1)), 'MinPeakHeight', 0, 'MinPeakProminence',MinPkPromTr);
                title(['Cell=' num2str(j) '  Pos=' num2str(i)])
                hold on
                plot( b2, -a2, 'ro')
                xlabel('Time (h)')
                ylabel('Levels of detrended time series')  
                if isnan(pos{i}.val(j))==1 && isnan(pos{i}.val2(j))==1
                    print(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/notenough/Cell_' num2str(j) '.png'],'-dpng','-r0')    
                else
                    print(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/good/Cell_' num2str(j) '.png'],'-dpng','-r0')
                end
            end
            clf
        else
            pos{i}.val(j)=NaN;
            pos{i}.val2(j)=NaN; 
            Amps_for=[Amps_for NaN ];
            Amps_back=[Amps_back NaN ];
            pos{i}.Ampsfor(j, 1:2)=NaN; 
            pos{i}.Ampsback(j, 1:2)=NaN; 
            pos{i}.WCV_for(j)=NaN;
            pos{i}.WCV_back(j)=NaN;
            pos{i}.ppmeans(j)=NaN; %mean period per cell
            pos{i}.CVp(j)=NaN;
            pos{i}.pkTrueTimesN{j}=NaN; 
            if exist(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/notperiodic/'],'file')~=7
                mkdir(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/notperiodic/']);
            end
            set(gcf,'DefaultLineLineWidth',3)
            plot(pos{i}.time(1:xs),pos{i}.YFP(1:xs,j),'o-', 'Linewidth',3)
            xlabel('Time (h)')
            ylabel('Levels of undetrendable time series')  
            print(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/notperiodic/Cell_' num2str(j) '.png'],'-dpng','-r0')
        end %% this is when you have a period restriction
        clf
        end %this is so you go through ALL THE CELLS IN THE POSITION
        Amps_forward{i}=Amps_for;
        Amps_backward{i}= Amps_back; 
        Amps_mix{i}=Amps;
    end
    save(['Doneanalysis/' name  '/summary/dist' num2str(mpkdis2) '.mat'], 'pos', 'counter', 'counter2', 'counter3','counter4', 'Amps_forward', 'Amps_backward', 'Amps_mix');
else    
    a=load([ 'Doneanalysis/' name  '/summary/dist' num2str(mpkdis2) '.mat']);
    pos=a.pos;
    counter=a.counter;
    counter2=a.counter2;
    counter3=a.counter3;
    Amps_forward= a.Amps_forward;
    Amps_backward = a.Amps_backward;
    Amps_mix=a.Amps_mix;
    clear a;
end
%% split pos13 into two sets: 
%pos13 missing two indexes at end: 
pos{13}.index(369:370)=200;  %program only takes into account up to point where rhythmic consideration passed
pos{10}=pos{13}; %save pos13 into pos10 and now split it up:
pos{13}=[];
xfin=find((pos{10}.y(1,:))<max(max(pos{10}.y(1,:)))-179.97 ); % y cut-off for root tip comparabel to WT. 
pos{13}.periods=pos{10}.periods(xfin, :); 
pos{13}.periodsgofs=pos{10}.periodsgofs(xfin, :);
pos{13}.periodsignore=pos{10}.periodsignore(xfin, :);
pos{13}.x=pos{10}.x(:, xfin);
pos{13}.y=pos{10}.y(:, xfin);
pos{13}.z=pos{10}.z(:, xfin);
pos{13}.YFP=pos{10}.YFP(:, xfin);
pos{13}.time=pos{10}.time;
pos{13}.val=pos{10}.val( xfin); 
pos{13}.WCV_for=pos{10}.WCV_for(xfin);
pos{13}.val2=pos{10}.val2(xfin);
pos{13}.WCV_back=pos{10}.WCV_back(xfin);
pos{13}.Ampsfor=pos{10}.Ampsfor(xfin, :);
pos{13}.Ampsback=pos{10}.Ampsback(xfin, :);
pos{13}.ppmeans=pos{10}.ppmeans(xfin);
pos{13}.CVp=pos{10}.CVp(xfin);
pos{13}.index=pos{10}.index(xfin);
xfintot=1:size(pos{10}.y,2);
xfintot(xfin)=[];
pos{12}.periods=pos{10}.periods(xfintot, :); 
pos{12}.periodsgofs=pos{10}.periodsgofs(xfintot, :);
pos{12}.periodsignore=pos{10}.periodsignore(xfintot, :);
pos{12}.x=pos{10}.x(:,xfintot);
pos{12}.y=pos{10}.y(:, xfintot);
pos{12}.z=pos{10}.z(:, xfintot);
pos{12}.YFP=pos{10}.YFP(:, xfintot);
pos{12}.time=pos{10}.time;
pos{12}.val=pos{10}.val(xfintot); 
pos{12}.WCV_for=pos{10}.WCV_for(xfintot);
pos{12}.val2=pos{10}.val2(xfintot);
pos{12}.WCV_back=pos{10}.WCV_back(xfintot);
pos{12}.Ampsfor=pos{10}.Ampsfor(xfintot, :);
pos{12}.Ampsback=pos{10}.Ampsback(xfintot, :);
pos{12}.ppmeans=pos{10}.ppmeans(xfintot);
pos{12}.CVp=pos{10}.CVp(xfintot);
pos{12}.pkTrueTimesN=pos{10}.pkTrueTimesN(xfintot);
pos{12}.index=pos{10}.index(xfintot);
alltimes=[];
for j=2:numpos % pos12 and 13 have same time points
    i=posi(j);   
    alltimes=[alltimes pos{i}.time];
end
alltimes=sort(alltimes, 'ascend');
%now interpolate only over the needed region
posi= 12:20; 
numpos=length(posi);
for k=1:numpos
    i=posi(k);
    if  i==17 || i==18 || i==19 || i==20
        for j=1: size(pos{i}.YFP,2)
            pos{i}.YFPinterp(:,j)=interp1(pos{i}.time(1:end-1), pos{i}.YFP(:,j), alltimes);%when outside the range, gives nan
        end
    else
        for j=1: size(pos{i}.YFP,2)
            pos{i}.YFPinterp(:,j)=interp1(pos{i}.time, pos{i}.YFP(:,j), alltimes);%when outside the range, gives nan
        end  
    end
end
%%
pos{1}=pos{15}; 
pos{2}=pos{16}; 
pos{3}=pos{17}; 
pos{15}.periods= [pos{15}.periods; pos{16}.periods; pos{17}.periods]; 
pos{15}.periodsgofs= [pos{15}.periodsgofs; pos{16}.periodsgofs;  pos{17}.periodsgofs];
pos{15}.periodsignore= [pos{15}.periodsignore; pos{16}.periodsignore; pos{17}.periodsignore] ;
pos{15}.x= [pos{15}.x(1:end-1, :)  pos{16}.x  pos{17}.x] ;
pos{15}.y= [pos{15}.y(1:end-1, :) pos{16}.y  pos{17}.y];
pos{15}.z= [pos{15}.z(1:end-1, :) pos{16}.z  pos{17}.z];
pos{15}.YFP= [pos{15}.YFP(1:end-1, :) pos{16}.YFP(1:end-1, :)  pos{17}.YFP] ;
pos{15}.time= [pos{15}.time(1:end-1) pos{16}.time pos{17}.time] ;
pos{15}.val= [pos{15}.val pos{16}.val pos{17}.val] ; 
pos{15}.WCV_for= [pos{15}.WCV_for pos{16}.WCV_for pos{17}.WCV_for];
pos{15}.val2= [pos{15}.val2 pos{16}.val2 pos{17}.val2];
pos{15}.WCV_back= [pos{15}.WCV_back pos{16}.WCV_back pos{17}.WCV_back];
pos{15}.Ampsfor= [pos{15}.Ampsfor; pos{16}.Ampsfor; pos{17}.Ampsfor];
pos{15}.Ampsback= [pos{15}.Ampsback; pos{16}.Ampsback; pos{17}.Ampsback];
pos{15}.ppmeans= [pos{15}.ppmeans pos{16}.ppmeans pos{17}.ppmeans];
pos{15}.CVp= [pos{15}.CVp pos{16}.CVp pos{17}.CVp];
pos{15}.pkTrueTimesN= [pos{15}.pkTrueTimesN pos{16}.pkTrueTimesN pos{17}.pkTrueTimesN];
pos{15}.index= [pos{15}.index pos{16}.index pos{17}.index];
pos{15}.YFPinterp= [ pos{15}.YFPinterp pos{16}.YFPinterp pos{17}.YFPinterp]; 
pos{15}.YFPinterp(1:3,:)=nan; 
pos{15}.YFPinterp(207:217,:)=nan; 
%%
numpos=7;
counter(1:7)=zeros(1,7);
posi=[ 12:15 18:20]; 
%% redo counter again: 
gof_cutoff=0.9; 
gof_cutoff2=1;
for k=start_pos:numpos
    i=posi(k);
    pos{i}.indy=[];
    pos{i}.indygt24=[];
    for j= 1:size(pos{i}.YFP,2)
      if  pos{i}.periodsgofs(j,1)<gof_cutoff && pos{i}.periodsgofs(j,2)<gof_cutoff && pos{i}.periodsgofs(j,3)<gof_cutoff2  && abs(pos{i}.periods(j,1)- pos{i}.periods(j,3))<2.5 && abs(pos{i}.periods(j,2)- pos{i}.periods(j,3))<2.5 && abs(pos{i}.periods(j,2)- pos{i}.periods(j,1))<2.5 && sum(pos{i}.periodsignore(j,:))==0%<0.9  & pos{i}.periods(j,5)<=35 & pos{i}.periods(j,5)>=18 & pos{i}.periods(j,9)>0 & abs(pos{i}.periods(j,5)- pos{i}.periods(j,1))<=5% pos{i}.periods(j,2)<0.5
           counter(k)=counter(k)+1;
           pos{i}.indy=[pos{i}.indy j];
           if pos{i}.periods(j,1)>24
               pos{i}.indygt24=[pos{i}.indygt24 j];
           end
       end
    end
 end
%% Figure 1 source data 1. (Table data for the percentage rhythmic cells)
% TableOfRhythmics is  the table shown in Figure 1 source data 1.
nonrhyth2=zeros(7,6); 
xb=zeros(1,6);xa=zeros(1,6);
for k=start_pos:numpos
    i=posi(k); 
    vv1=(pos{i}.periodsgofs(:,1)<0.9& pos{i}.periodsignore(:,1)==0); %+(pos{i}.periodsignore(:,1)==0);
    vv2=(pos{i}.periodsgofs(:,2)<0.9 & pos{i}.periodsignore(:,2)==0); %+(pos{i}.periodsignore(:,2)==0);
    vv3=(pos{i}.periodsgofs(:,3)<1 & pos{i}.periodsignore(:,3)==0); %+(pos{i}.periodsignore(:,3)==0);
    vv4=(abs(pos{i}.periods(:,1)- pos{i}.periods(:,3))<2.5 & abs(pos{i}.periods(:,2)- pos{i}.periods(:,3))<2.5 & abs(pos{i}.periods(:,1)- pos{i}.periods(:,2))<2.5); 
    vv5=(pos{i}.periodsgofs(:,1)<0.9& pos{i}.periodsignore(:,1)==0  & pos{i}.periodsgofs(:,2)<0.9 & pos{i}.periodsignore(:,2)==0 & pos{i}.periodsgofs(:,3)<1 & pos{i}.periodsignore(:,3)==0);
    pos{i}.siz=size(pos{i}.periods(vv5==1 & vv4==1 ), 1);
    nonrhyth2(k, 1)=size(pos{i}.periods(vv1==1), 1)/size(pos{i}.periods,1); %RAE cut-off at 0.6
    nonrhyth2(k, 2)=size(pos{i}.periods(vv2==1), 1)/size(pos{i}.periods,1);
    nonrhyth2(k, 3)=size(pos{i}.periods(vv3==1), 1)/size(pos{i}.periods,1);
    nonrhyth2(k,4)=size(pos{i}.periods(vv4==1), 1)/size(pos{i}.periods,1);
    nonrhyth2(k,5)=size(pos{i}.periods(vv5==1), 1)/size(pos{i}.periods,1);
    nonrhyth2(k,6)=size(pos{i}.periods(vv5==1 & vv4==1 ), 1)/size(pos{i}.periods,1);
    xb(i)=size(pos{i}.periods(vv5==1 & vv4==1), 1);
    xa(i)=size(pos{i}.periods,1);
end
% following is the table of rhythmic values (as appearing in the paper):
  TableOfRhythmics=nonrhyth2(end:-1:1,[6 1:5]);
%%
% mean of   amplitudes per cell, just taking into account that they have 3 or more  
nnforw=NaN(1,numpos);nnback=NaN(1,numpos);
for k=start_pos:numpos
    i=posi(k);
    nnforw(k)=nanmean(pos{i}.val(pos{i}.index<100));
    %make a meta- forward/backward mix called nnforw
    if isnan(nnforw(k))==1 %update with bac value if forward is NaN
        nnforw(k)=nanmean(pos{i}.val2(pos{i}.index<100));
    end
    nnback(k)=nanmean(pos{i}.val2(pos{i}.index<100));
end
%%  amplitude for mean trace: need to interpolate it
val=NaN(1, numpos);val2=NaN(1, numpos);
for k=start_pos:numpos 
    i=posi(k);
    clear r c 
    [r,c]=find(~isnan(pos{i}.val));
    pos{i}.meanlevsX=mean(pos{i}.YFPinterp(:,c),2);
    clear ks
    ks = ~isnan(pos{i}.meanlevsX);
    pos{i}.meanlevs= pos{i}.meanlevsX(ks);
    pos{i}.meanlevstime= alltimes(ks);
end
for j= start_pos:numpos
    i=posi(j);
    clear p trend dtdata
    p = polyfit(pos{i}.meanlevstime',pos{i}.meanlevs,1);  % find trend as polynomial of degree 'd'
    trend = polyval(p,pos{i}.meanlevstime);
    dtdata = pos{i}.meanlevs - trend';  
    MinPeakDist=mpkdis2*(pos{i}.meanlevstime(2)-pos{i}.meanlevstime(1));
    MinPeakDist2=3*(pos{i}.meanlevstime(2)-pos{i}.meanlevstime(1));
    if i~=13
        [~,b]=findpeaks(dtdata,pos{i}.meanlevstime,  'MinPeakDistance',8,  'MinPeakProminence',20,'MinPeakHeight', 0);
        [~,b2]=findpeaks(-dtdata,pos{i}.meanlevstime, 'MinPeakDistance',8, 'MinPeakProminence',20); 
    else
        [~,b]=findpeaks(dtdata,pos{i}.meanlevstime,  'MinPeakDistance',10,  'MinPeakProminence',0,'MinPeakHeight', 0);
        [~,b2]=findpeaks(-dtdata,pos{i}.meanlevstime, 'MinPeakDistance',10, 'MinPeakProminence',10); 
    end
    figure(1)
    findpeaks(dtdata, pos{i}.meanlevstime, 'MinPeakDistance',MinPeakDist,'MinPeakProminence',10);
    hold on
    figure(2)
    findpeaks(-dtdata, pos{i}.meanlevstime, 'MinPeakDistance',mpkdis,'MinPeakProminence',10);
    figure(3)
    for s=1: length(b)
       pkTime(s)=find(pos{i}.meanlevstime==b(s));
    end
    clear s
    for s=1: length(b2)
       troughTime(s)=find(pos{i}.meanlevstime==b2(s));
    end
    clear s
    b=pkTime;
    b2=troughTime;
    clf
    set(gcf,'DefaultLineLineWidth',3)
    if i~=13
       findpeaks(dtdata,pos{i}.meanlevstime,  'MinPeakDistance',8,  'MinPeakProminence',20,'MinPeakHeight', 0);
       [ac2,bc2]=findpeaks(-dtdata,pos{i}.meanlevstime, 'MinPeakDistance',8, 'MinPeakProminence',20); %troughsbar
       title([' Mean of Pos=' num2str(i)])
       hold on
       plot( bc2, -ac2, 'ro')
    else
       findpeaks(dtdata,pos{i}.meanlevstime,  'MinPeakDistance',10,  'MinPeakProminence',0,'MinPeakHeight', 0);
       [ac2,bc2]=findpeaks(-dtdata,pos{i}.meanlevstime, 'MinPeakDistance',10, 'MinPeakProminence',10); %troughsbar
       title([' Mean of Pos=' num2str(i)])
       hold on
       plot( bc2, -ac2, 'ro')
    end
    xlabel('Time (h)')
    ylabel('Levels of detrended mean time series')  
    print(['Doneanalysis/' name  '/figures/pos' num2str(i)  '__rootlumpmean.png'],'-dpng','-r0')    
    clf
    clear pkTime troughTime
    % need to measure amplitude by two ways trough to peak and peak to trough. 
    if length(b)==length(b2) % b is peak and b2 is trough
    if b(1)<b2(1) 
    diffs_forward=abs(pos{i}.meanlevs(b2(1:end-1))-pos{i}.meanlevs(b(2:end)));
    diffs_back=abs(pos{i}.meanlevs(b2)-pos{i}.meanlevs(b));
    diffs_forward=[diffs_forward; diffs_back]; 
    if length(diffs_forward)>2
    val(i)=mean(diffs_forward);
    else
     val(i)=NaN;   
    end
    if length(diffs_back)>2
    val2(i)=mean(diffs_back);
     else
     val2(i)=NaN;   
    end
     pos{i}.AmpsforMean( 1:2)=diffs_forward([ 1 end]); 
    pos{i}.AmpsbackMean( 1:2)=diffs_back([ 1 end]); 
    
    elseif b(1)>b2(1)
    diffs_forward=abs(pos{i}.meanlevs(b2)-pos{i}.meanlevs(b));
    diffs_back=abs(pos{i}.meanlevs(b2(2:end))-pos{i}.meanlevs(b(1:end-1)));
    diffs_forward=[diffs_forward; diffs_back]; 
    if length(diffs_forward)>2
    val(i)=mean(diffs_forward);
    else
     val(i)=NaN;   
    end
    if length(diffs_back)>2
    val2(i)=mean(diffs_back);
    else
    val2(i)=NaN;   
    end
    pos{i}.AmpsforMean( 1:2)=diffs_forward([ 1 end]); 
    pos{i}.AmpsbackMean( 1:2)=diffs_back([ 1 end]); 
    end    
    elseif length(b)>length(b2)
      if b(1)<b2(1)
      diffs_forward=abs(pos{i}.meanlevs(b2)-pos{i}.meanlevs(b(2:end)));
      diffs_back=abs(pos{i}.meanlevs(b2)-pos{i}.meanlevs(b(1:end-1)));
      diffs_forward=[diffs_forward; diffs_back]; 
     if length(diffs_forward)>2
    val(i)=mean(diffs_forward);
    else
     val(i)=NaN;   
    end
    if length(diffs_back)>2
    val2(i)=mean(diffs_back);
     else
     val2(i)=NaN;   
    end 
     pos{i}.AmpsforMean( 1:2)=diffs_forward([ 1 end]); 
    pos{i}.AmpsbackMean( 1:2)=diffs_back([ 1 end]); 
      end
    elseif length(b)<length(b2)
      if b(1)>b2(1)
       diffs_forward=abs(pos{i}.meanlevs(b2(1:end-1))-pos{i}.meanlevs(b));
       diffs_back=abs(pos{i}.meanlevs(b2(2:end))-pos{i}.meanlevs(b));
       diffs_forward=[diffs_forward; diffs_back]; 
      if length(diffs_forward)>2
          val(i)=mean(diffs_forward);
      else
          val(i)=NaN;
      end
        if length(diffs_back)>2
            val2(i)=mean(diffs_back);
        else
            val2(i)=NaN;   
        end
     pos{i}.AmpsforMean( 1:2)=diffs_forward([ 1 end]); 
    pos{i}.AmpsbackMean( 1:2)=diffs_back([ 1 end]); 
      end
      
    end
    if isnan(val(i))~=1
    valmix(i)=val(i);
else
    if isnan(val2(i))~=1
     valmix(i)=val2(i);
    else
      valmix(i)=NaN;  
    end
    end
    clear a b a2 b2 
end
%%
% BAR PLOT OF Amplitude of mean trace vs  Mean of cell amplitudes
% note that amplitude of any trace is a mean of diffs_forward and
% diffs_backward
btwn=NaN(1,numpos); witn=NaN(1,numpos);
s=1;
for k=start_pos:numpos
    i=posi(k);
    btwn(s)=nanstd(pos{i}.ppmeans)/nanmean(pos{i}.ppmeans);
    witn(s)=nanmean(pos{i}.CVp);
    s=s+1;
end
%% Figure 1 Supp3 A
YFPint{1}=pos{12}.YFPinterp;
YFPint{2}=[pos{13}.YFPinterp pos{14}.YFPinterp pos{15}.YFPinterp]; 
YFPint{3}=[ pos{18}.YFPinterp pos{19}.YFPinterp];
YFPint{4}=[pos{20}.YFPinterp];
YFPintot=[YFPint{1} YFPint{2} YFPint{3} YFPint{4}];
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1600 3600]/300)
maxylim=1400;
for j=1:4
    set(gca, 'Fontsize', 7)
    i=posi(j);
    subplot(4,1,5-j)
    hold on 
    y = [0 0 maxylim maxylim];
    x = [36 48 48 36];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [60 72 72 60];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [84 96 96 84];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none') 
    x = [108 120 120 108];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [132 144 144 132];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [156 168 168 156];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [180 192 192 180];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')   
    x = [204 216 216 204];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')   
    ylabel('Mean Fluor. Intensity', 'Fontsize', 7) 
    ylim([0 maxylim])  
    xlim([72 168])  
    clear x y 
    if j==1 || j==4
        y=nanmean(YFPint{j}, 2);
    elseif j==2
        y=nanmean(YFPint{j}, 2);
        y([1:3 207:217])=nan; 
    elseif j==3
        y=nanmean(YFPint{j}, 2);
        y([1:5 209:217])=nan; 
    end  
    x=alltimes;
    if j==1
        plot(x,y, 'Color',[0.6,0.6,0.6])
        set(gca, 'LineWidth', 1.2)
        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
     elseif j==2
        plot(x,y, 'k')
        set(gca, 'LineWidth', 1.2)
        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
     elseif j==3
        plot(x,y, 'b')
        set(gca, 'LineWidth', 1.2)
        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
     elseif j==4
        plot(x,y, 'r')
        set(gca, 'LineWidth', 1.2)
        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
    end
    hold on
    ax = gca;
    if j==1
        ax.XTick = 0:24:220;
    else
       set(gca,'xtick',[])
    end    
    if j==1
        xlabel('Time (h)', 'Fontsize', 7) 
    end
    hold on
    set(gca, 'Fontsize', 7)
    set(gca, 'LineWidth', 1.2)
    set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
end
if exist('Figure1Supp3A.pdf', 'file')~=2
    print('-Painters',kkq, 'Figure1Supp3A','-dpdf','-r300')
    movefile('Figure1Supp3A.pdf', 'Doneanalysis/WTrepeat/figures')
end
%%  Figure 1 Supp3 B
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1600 3600]/300)
maxylim=7000;
for j=1:4
    set(gca, 'Fontsize', 7)
    i=posi(j);
    subplot(4,1,5-j)
    hold on 
    y = [0 0 maxylim maxylim];
    x = [36 48 48 36];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [60 72 72 60];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [84 96 96 84];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none') 
    x = [108 120 120 108];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [132 144 144 132];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [156 168 168 156];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
    x = [180 192 192 180];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')   
    x = [204 216 216 204];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')   
    ylabel('Mean Fluor.', 'Fontsize', 7) 
    ylim([0 6000])  
    xlim([72 168])  
    clear x y 
    hold on
    if j==1
        cmap=colormap(gray(size(pos{13}.YFP,2)));
        for k= 1:size(pos{13}.YFP,2)
            plot(pos{13}.time, pos{13}.YFP(:,k), 'Color', cmap(randi([1 size(pos{13}.YFP,2)]), :)*diag([1 1 1]), 'Linewidth', 0.05);
            hold on 
            set(gca, 'LineWidth', 1.2)
            set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
        end
        ylim([0 4300]);
    elseif j==2
         clear concat
         cmap=colormap(gray(size(pos{14}.YFP,2)));
         for k= 1:size(pos{14}.YFP,2)  
             clear cmap
             cmap=colormap(gray(size(pos{14}.YFP,2)));
             plot(pos{14}.time, pos{14}.YFP(:,k), 'Color', cmap(randi([1 size(pos{14}.YFP,2)]), :)*diag([1 1 1]), 'Linewidth', 0.05);
             hold on 
             set(gca, 'LineWidth', 1.2)
             set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
         end   
         for k= 1:size(pos{1}.YFP,2)  
                    clear cmap
                     cmap=colormap(gray(size(pos{15}.YFP,2)));
                     plot(pos{1}.time, pos{1}.YFP(:,k), 'Color', cmap(randi([1 size(pos{1}.YFP,2)]), :)*diag([1 1 1]), 'Linewidth', 0.05);
                     hold on 
                     set(gca, 'LineWidth', 1.2)
                     set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
          end
         for k= 1:size(pos{2}.YFP,2)  
                    clear cmap
                     cmap=colormap(gray(size(pos{2}.YFP,2)));
                     plot(pos{2}.time, pos{2}.YFP(:,k), 'Color', cmap(randi([1 size(pos{2}.YFP,2)]), :)*diag([1 1 1]), 'Linewidth', 0.05);
                  hold on 
               set(gca, 'LineWidth', 1.2)
            set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
                 end
         for k= 1:size(pos{3}.YFP,2)  
         clear cmap
                     cmap=colormap(gray(size(pos{3}.YFP,2)));
                  %  if mod(k,1)==0
                     plot(pos{3}.time(1:end-1), pos{3}.YFP(:,k), 'Color', cmap(randi([1 size(pos{3}.YFP,2)]), :)*diag([1 1 1]), 'Linewidth', 0.05);
                  hold on 
                set(gca, 'LineWidth', 1.2)
            set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
                 %   end
                end
         ylim([0 2800]);
    elseif j==3
         concat=[pos{18}.YFP pos{19}.YFP];
         lev=size(concat,2)*2;
         cc=redblue(lev);
         cmap= cc(1:lev/2,:);
         for k= 1:size(pos{18}.YFP,2)
             plot(pos{18}.time(1:end-1), pos{18}.YFP(:,k), 'Color',cmap(randi([1 lev/2]), :)*diag([1 1 1]), 'Linewidth', 0.05);
             hold on 
             set(gca, 'LineWidth', 1.2)
             set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
         end
         for k= 1:size(pos{19}.YFP,2)
             plot(pos{19}.time(1:end-1), pos{19}.YFP(:,k), 'Color', cmap(randi([1 lev/2]), :)*diag([1 1 1]), 'Linewidth', 0.05);
             hold on 
             set(gca, 'LineWidth', 1.2)
             set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
          end
         ylim([0 4000]);
    elseif j==4
         lev=size(pos{20}.YFP,2)*2;
         clear cc
         cc=redblue(lev);
         cmap= cc(end:-1:lev/2+1,:);
         for k= 1:size(pos{20}.YFP,2)
             plot(pos{20}.time(1:end-1), pos{20}.YFP(:,k), 'Color', cmap(randi([1 lev/2]), :)*diag([1 1 1]), 'Linewidth', 0.05);
             hold on 
             set(gca, 'LineWidth', 1.2)
             set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
          end  
                  ylim([0 2500]);
    end
    xlim([72 168])     
    ax = gca;  
    if j==1
        ax.XTick=0:24:220;
    else
       set(gca,'xtick',[])
    end
    if j==1
        xlabel('Time (h)', 'Fontsize', 7) 
    end
    hold on
    set(gca, 'Fontsize', 7)
end
if exist('Figure1Supp3B.pdf', 'file')~=2
    print('-Painters',kkq, 'Figure1Supp3B','-dpdf','-r300')
    movefile('Figure1Supp3B.pdf', 'Doneanalysis/WTrepeat/figures')
end
close all
%% Figure 2 Supp1A
Bigcombo=NaN([ numpos size((pos{13}.val(pos{13}.index<100))') ]); 
for k=1:numpos
    i=posi(k); 
    Bigcombo(k,1:size((pos{i}.val(pos{i}.index<100))'),1)=(pos{i}.val(pos{i}.index<100));
    Bigcombo(k,size((pos{i}.val(pos{i}.index<100))',1)+1:end)=nan;  
end
Bigcombo=Bigcombo';
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  3000 1800]/300)
axes('FontSize',7);
bh=boxplot(Bigcombo, 'symbol','');
set(bh, 'Linewidth', 1)
secto={ 'Root tip'; 'Root (up from tip)'; 'Root (sect 1)'; 'Root (sect 2, 3, 4)';'Hypocotyl (sect 1)'; 'Hypocotyl (sect 2)'; 'Cotyledon'};
ylabel('Amplitude levels')
set(gca, 'XTick', 1:numpos, 'XTickLabel', secto,'FontSize',7);
hold on
for nums=1:numpos
    line(nums-0.35:0.7:nums+0.35, [valmix(posi(nums)) valmix(posi(nums))], 'Linewidth', 2, 'Color', 'g')
end
ylim([0 2000])
box off
xlim([0.5 numpos+0.5])
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
ax=gca;
ax.XTickLabelRotation=25;
xlim([0.5 numpos+0.5])
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure2Supp1A.pdf', 'file')~=2
    print('-Painters',kkq, 'Figure2Supp1A','-dpdf','-r300')
    movefile('Figure2Supp1A.pdf', 'Doneanalysis/WTrepeat/figures')
end
%%
nnforw=NaN(1,6);
nnback=NaN(1,6);
nnforwstd=NaN(1,6);
nnbackstd=NaN(1,6);
for k=start_pos:numpos
    i=posi(k);
    nnforw(k)=nanmean(pos{i}.val(pos{i}.index<100));
    nnback(k)=nanmean(pos{i}.val2(pos{i}.index<100));
end
for k=start_pos:numpos
    i=posi(k);
    nnforwstd(k)=nanstd(pos{i}.val(pos{i}.index<100));
    nnbackstd(k)=nanstd(pos{i}.val2(pos{i}.index<100));
end
s=1;
for k=start_pos:numpos
    i=posi(k);
    btwn(k)=nanstd(pos{i}.ppmeans)/nanmean(pos{i}.ppmeans);
    witn(k)=nanmean(pos{i}.CVp);
    s=s+1;
end
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  3000 1200]/300)
h=bar([btwn' witn'],'BarWidth', 1);
set(gca, 'FontSize', 7)
set(gca, 'XTick', 1:numpos, 'XTickLabel', cell(numpos,1));
set(h(1),'facecolor','k') 
set(h(2),'facecolor','w') 
h(2).LineWidth=1;
h=legend({' Between-cell', ' Within-cell'}, 'Position', [0.75,0.9,0.01,0.01], 'box', 'off');
box off
ylabel({' Period variability (CV)'})
xlim([0.5 numpos+0.5])
set(gca, 'YTick', 0:0.05:2, 'YTickLabel', 0:0.05:2);
secto={ 'Root tip'; 'Root (up from tip)'; 'Root (sect 1)'; 'Root (sect 2, 3, 4)';'Hypocotyl (sect 1)'; 'Hypocotyl (sect 2)'; 'Cotyledon'};
set(gca, 'XTick', 1:numpos, 'XTickLabel', secto);
set(gca, 'Fontsize', 7)
box off
ax=gca; 
ax.XTickLabelRotation=25;
xlim([0.5 numpos+0.5])
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure2Supp1B.pdf', 'file')~=2
    print('-Painters',kkq, 'Figure2Supp1B','-dpdf','-r300')
    movefile('Figure2Supp1B.pdf', 'Doneanalysis/WTrepeat/figures')
end
 %%  Figure 1 supp6 B 


for k= start_pos:numpos
    i=posi(k); 
    % isnan(pos{i}.val)~=1 is just to make sure only taking itno account
    % cells where more than 2 amplitudes. call it AA, mainly based on
    % forward but when missing then on backward
    ccount=0;
    for j=1:size(pos{i}.YFP,2)
    if isnan(pos{i}.val(j))~=1
        pos{i}.AA(j, :)=pos{i}.Ampsfor(j, :);
    else
         if isnan(pos{i}.val2(j))~=1
             pos{i}.AA(j, :)=pos{i}.Ampsback(j, :);
         else
              pos{i}.AA(j, :)=[NaN NaN];
              ccount=ccount+1;
         end   
    end
    end 
    perc{i}=(size(pos{i}.YFP,2)-ccount)/size(pos{i}.YFP,2);
    ccount_all{i}=ccount;
end

posnew{1}.AA=pos{12}.AA;
posnew{2}.AA=[pos{13}.AA;pos{14}.AA; pos{15}.AA];
posnew{3}.AA=[pos{18}.AA;pos{19}.AA];
posnew{4}.AA=pos{20}.AA;
ccounts{1}=ccount_all{12};
ccounts{2}=ccount_all{13}+ccount_all{14}+ccount_all{15};
ccounts{3}=ccount_all{18}+ccount_all{19};
ccounts{4}=ccount_all{20}; 
sectnew={ {'Root tip'},{'Root (rest)'}, {'Hypocotyl'},{'Cotyledon'}};

kkq=figure;
axes('FontSize',7);
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  5400 3200]/300)
for i=1:4
    subplot(1,4,i)
    bh=boxplot(posnew{i}.AA, 'symbol','');
    p = prctile(posnew{i}.AA,[9 91]);
    h = flipud(findobj(gca,'Tag','Upper Whisker'));
    for j=1:length(h)
        ydata = get(h(j),'YData');
        ydata(2) = p(2,j);
        set(h(j),'YData',ydata);
    end
    h = flipud(findobj(gca,'Tag','Upper Adjacent Value'));
    for j=1:length(h)
        ydata = get(h(j),'YData');
        ydata(:) = p(2,j);
        set(h(j),'YData',ydata);
    end
    h = flipud(findobj(gca,'Tag','Lower Whisker'));
    for j=1:length(h)
        ydata = get(h(j),'YData');
        ydata(1) = p(1,j);
        set(h(j),'YData',ydata);
    end
    h = flipud(findobj(gca,'Tag','Lower Adjacent Value'));
    for j=1:length(h)
        ydata = get(h(j),'YData');
        ydata(:) = p(1,j);
        set(h(j),'YData',ydata);
    end
    secto={'first amplitude', 'last amplitude'};
    set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
    set(gca,'FontSize',7);
    ax=gca;
    ax.XTickLabelRotation=25;
    ax.XTickLabel = secto;
    if i==1     
        ylim([0 1300])
        set(gca, 'LineWidth', 1.2)
        ylabel('Amplitude levels')
    elseif i==2
        ylim([0 1200])
        set(gca, 'LineWidth', 1.2)
    elseif i==3
        ylim([0 1400])
        set(gca, 'LineWidth', 1.2)
    elseif i==4
        ylim([0 1000])
        set(gca, 'LineWidth', 1.2)
    end    
    set(gca, 'Fontsize',7);
    title(sectnew{i})
end
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
if exist('Figure1Supp6B.pdf', 'file')~=2
    print('-Painters',kkq, 'Figure1Supp6B','-dpdf','-r300')
    movefile('Figure1Supp6B.pdf', 'Doneanalysis/WTrepeat/figures')
end
for i=1:4 
    NumberNotUsed(i)=sum(isnan(posnew{i}.AA(:,1)));
    All(i)=size(posnew{i}.AA,1);
end
Percentages=(All-NumberNotUsed)./All*100; %percentage of all cells used in the analysis in the last boxplot. 
close all



