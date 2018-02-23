%% README
%  File for statistical analysis of WT data for Gould, Domijan et al.
%  "Coordination of robust single cell rhythms in the Arabidopsis circadian 
%  clock via spatial waves of gene expression" (BioRxiv)
%  doi https://doi.org/10.1101/208900s
% 
% file is supplemented by: Analysis_WT.m
%  by Mirela Domijan (U. of Liverpool)
%% which experiment we are looking at 
name='WT'; %user chosen
%% first add in useful paths
addpath([pwd '/Data_singlecell/'])
addpath([pwd '/Data_singlecell/' name '_final_coordinates'])
% For some plots , install Red Blue Colormap package from MATHWORKS from:
% https://uk.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap?focused=5138232&tab=function
% and place it in the Data_singlecell folder. 
addpath([pwd '/redblue/'])
%% is there a folder where to save all the processed data?
if exist(['Doneanalysis/' name  '/summary'],'file')~=7
    mkdir(['Doneanalysis/' name  '/summary']);
end
%% is there a folder where to save all the figures?
if exist(['Doneanalysis/' name  '/figures'],'file')~=7
    mkdir(['Doneanalysis/' name  '/figures']);
end
addpath(genpath([pwd '/Doneanalysis/']))
%% pull out all the period estimates and the data
[aa, texts]=xlsread(['Data_singlecell/' name '_final_coordinates/interpolated/Individual results.xlsx']);
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

posi=1:6;
separators=NaN(1,length(posi));
for j=1:length(posi)
    xs=csvread(['Data_singlecell/' name '_final_coordinates/interpolated/interpolateddata_WT_pos' num2str(posi(j)) '.csv']);
    separators(j)= size(xs,2)-1;
end
clear j rowstorem rowstorem2  xs ind x
ss=cumsum(separators); 
separates=ss;
pp=reshape(aa(:,12),3,length(aa(:,12))/3)';
ppgof=reshape(aa(:,34),3,length(aa(:,34))/3)';
ignor=NaN(1,length(texts(:,11)));
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
    sect={ {'Root tip'},{'Lower root'}, {'Upper root'},{'Lower hypocotyl'},{'Upper hypocotyl'},{'Cotyledon'}};
end
%% now which positions you will look at 
start_pos=1;
numpos=6;
%% for each file you pull out, write down the positions; the YFP levels and the periods. 
for k=start_pos:numpos
%two types of periods: (i) periods(:,1)== FFT unscaled; (ii) periods(:,2)== FFT unscaled RAE; (iii) periods(:,9)== LS acceptance 10%
    i=posi(k); 
    if strcmp(name, 'WT')==1
        name2=['pos',num2str(i)];
        posfile=dir([pwd , '/Data_singlecell/' name '_final_coordinates/',[name2 '*tracked*.csv']]);
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

    if numel(unique(diff([0; fa])))==1 && (size(LL,1)-fa(end))== (fa(1)-1) 
        res=reshape(triple', fa(1)*3, (length(fa)+1));
        % verify that everything is pulled out correctly and that you have correct cell IDs;
        num_cells=size(res,1)/3;
        cell_labels=str2double(C{1}(4:end,8));
        %now make sure they are always increasing and OK (quick method):    
        diffy=diff(cell_labels)-1;
        diffywo = diffy(diffy~=0);
        %Note: diffy  keeps track of differences in cell parents (i.e. if
        %same time point, cell parents/labels should be increasing by one)
        %when new time point starts diffywo records how many cells there
        %are recorded in each time point.
        allok=isequal(abs(diffywo),ones(size(res,2)-1,1)*num_cells);
        % allok tells us that thre are right number of cells within each time point 
        data.pos=res';
        pos{i}.x=data.pos(:,1:3:end);
        pos{i}.y=data.pos(:,2:3:end); 
        pos{i}.z=data.pos(:,3:3:end);
    end
    
    % clear everything not needed:
    clear diffy diffywo res  fa fid LL cell_labels triple
    % now put in the YFP files  
    if i==4
    YFPfile=dir([pwd , '/Data_singlecell/'  name '_final_coordinates/',['mean*' name2 '*YFP_median*.csv']]);
    else
     YFPfile=dir([pwd ,'/Data_singlecell/'  name '_final_coordinates/',['*_' name2 '*YFP_median*.csv']]);
    end
    
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
    allok2=isempty(diffywo2);
    end
    clear cell_labels2 allok allok2 fid YFPfile diffy2 diffywo2 num_cells
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
clear i k  posfile bb a Times C name2 name3 
%% removal of doubles in overlapping sections: (these are calculated in Analysis_WT.m as 'dontkeep' variables)
rem2=[21    25    27    30   37   102];
rem3=[ 1     2     3     5     8    13    14    16    20    30    33    34    37    40    41    47    67    68]; 
rem5=[3     6    40    41    53];
% now remove them: 
for k=1:numpos
    if k==2
        pos{k}.periods(rem2,:)=[];
        pos{k}.periodsgofs(rem2,:)=[];
        pos{k}.periodsignore(rem2,:)=[];
        pos{k}.x(:, rem2)=[];
        pos{k}.y(:,rem2)=[];
        pos{k}.z(:, rem2)=[];
        pos{k}.YFP(:,rem2)=[];    
    elseif k==3
        pos{k}.periods(rem3,:)=[];
        pos{k}.periodsgofs(rem3,:)=[];
        pos{k}.periodsignore(rem3,:)=[];
        pos{k}.x(:, rem3)=[];
        pos{k}.y(:,rem3)=[];
        pos{k}.z(:, rem3)=[];
        pos{k}.YFP(:,rem3)=[];     
    elseif k==5   
        pos{k}.periods(rem5,:)=[];
        pos{k}.periodsgofs(rem5,:)=[];
        pos{k}.periodsignore(rem5,:)=[];
        pos{k}.x(:, rem5)=[];
        pos{k}.y(:,rem5)=[];
        pos{k}.z(:, rem5)=[];
        pos{k}.YFP(:,rem5)=[];
    end 
end
%% Figure 1 source data 1. (Table data for the percentage rhythmic cells)
% TableOfRhythmics is  the table shown in Figure 1 source data 1.
gof_cutoff=0.9;
gof_cutoff2=1;
nonrhyth2=zeros(7,6); 
for k=start_pos:numpos
    i=posi(k); 
    vv1=(pos{i}.periodsgofs(:,1)<0.9& pos{i}.periodsignore(:,1)==0); %+(pos{i}.periodsignore(:,1)==0);
    vv2=(pos{i}.periodsgofs(:,2)<0.9 & pos{i}.periodsignore(:,2)==0); %+(pos{i}.periodsignore(:,2)==0);
    vv3=(pos{i}.periodsgofs(:,3)<1 & pos{i}.periodsignore(:,3)==0); %+(pos{i}.periodsignore(:,3)==0);
    vv4=(abs(pos{i}.periods(:,1)- pos{i}.periods(:,3))<2.5 & abs(pos{i}.periods(:,2)- pos{i}.periods(:,3))<2.5 & abs(pos{i}.periods(:,1)- pos{i}.periods(:,2))<2.5); 
    vv5=(pos{i}.periodsgofs(:,1)<0.9& pos{i}.periodsignore(:,1)==0  & pos{i}.periodsgofs(:,2)<0.9 & pos{i}.periodsignore(:,2)==0 & pos{i}.periodsgofs(:,3)<1 & pos{i}.periodsignore(:,3)==0);
    pos{i}.siz=size(pos{i}.periods(vv5==1 & vv4==1 ), 1);
    nonrhyth2(k, 1)=size(pos{i}.periods(vv1==1), 1)/size(pos{i}.periods,1); %RAE cut-off 
    nonrhyth2(k, 2)=size(pos{i}.periods(vv2==1), 1)/size(pos{i}.periods,1);
    nonrhyth2(k, 3)=size(pos{i}.periods(vv3==1), 1)/size(pos{i}.periods,1);
    nonrhyth2(k,4)=size(pos{i}.periods(vv4==1), 1)/size(pos{i}.periods,1);
    nonrhyth2(k,5)=size(pos{i}.periods(vv5==1), 1)/size(pos{i}.periods,1);
    nonrhyth2(k,6)=size(pos{i}.periods(vv5==1 & vv4==1 ), 1)/size(pos{i}.periods,1);
end
% following is the table of rhythmic values (as appearing in the paper):
  TableOfRhythmics=nonrhyth2(end:-1:1,[6 1:5]);
%% Amplitude
  counter=NaN(1,numpos);counter2=NaN(1,numpos);counter3=NaN(1,numpos);counter4=NaN(1,numpos);
  mpkdis2=15;
if exist(['Doneanalysis/' name  '/summary/dist' num2str(mpkdis2) '.mat'], 'file')~=2
    mkdir(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/']);
    for k=  start_pos:numpos
        i=posi(k);
        Amps_for=[]; 
        Amps_back=[];     
        counter(k)=0; %counter for the number of cells that will be tested
        counter2(k)=0; %counter for the cells that have more than 2 peaks
        counter3(k)=0; %counter for the cells that have more than 2 peaks and where the number max and min correctly determined
        counter4(k)=0; %counter for the cells that have more than 2 peaks and where the number max and min correctly determined and peaks nad troughs come at right time
        xs=size(pos{k}.time,2);
        pos{i}.indy=[];
        for j= 1:size(pos{i}.YFP,2)
            if  pos{i}.periodsgofs(j,1)<gof_cutoff && pos{i}.periodsgofs(j,2)<gof_cutoff && pos{i}.periodsgofs(j,3)<gof_cutoff2  && abs(pos{i}.periods(j,1)- pos{i}.periods(j,3))<2.5 && abs(pos{i}.periods(j,2)- pos{i}.periods(j,3))<2.5 && abs(pos{i}.periods(j,2)- pos{i}.periods(j,1))<2.5 && sum(pos{i}.periodsignore(j,:))==0%<0.9  & pos{i}.periods(j,5)<=35 & pos{i}.periods(j,5)>=18 & pos{i}.periods(j,9)>0 & abs(pos{i}.periods(j,5)- pos{i}.periods(j,1))<=5% pos{i}.periods(j,2)<0.5
                counter(k)=counter(k)+1;
                pos{i}.indy=[pos{i}.indy j];
                mpkdis=15; 
                if exist(['Doneanalysis/' name  '/summary/dist' num2str(mpkdis2)],'file')~=7
                    mkdir(['Doneanalysis/' name  '/summary/' 'dist' num2str(mpkdis2)]);
                end
                % finding the peaks strategies:
                % detrending out the base value and then the following settings
                MinPeakDist=mpkdis2*(pos{i}.time(2)-pos{i}.time(1));
                MinPeakDist2=mpkdis2*(pos{i}.time(2)-pos{i}.time(1));
                MinPkHeight=-200; 
                MinPkProm=20;
                MinPkPromTr=1;
               % Find the peaks:
                clear dtdata p trend
                p = polyfit(pos{i}.time(1:xs)',pos{i}.YFP(1:xs,j),1);  % find trend as polynomial of degree 'd'
                trend = polyval(p,pos{i}.time(1:xs));
                dtdata = pos{i}.YFP(1:xs,j) - trend';  
                dtdata= smooth(pos{i}.time(1:xs),dtdata,0.05,'rloess');
                [a,b]= findpeaks(dtdata, pos{i}.time(1:xs),  'MinPeakDistance',16, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkProm); 
                [a2,b2]= findpeaks(-dtdata, pos{i}.time(1:xs), 'MinPeakDistance',16, 'MinPeakHeight',-50, 'MinPeakProminence',MinPkPromTr);
                pos{i}.index(j)=200;  %this keeps track of cells that will have enough amplitudes for analysis an tells you how many amplitudes
                if length(b)>2
                    pos{i}.index(j)=100;
                    counter2(k)=counter2(k)+1;
                    if abs(length(b)-length(b2))<2
                        counter3(k)=counter3(k)+1;
                        pos{i}.index(j)=max([length(b) length(b2)]);   
                    elseif abs(length(b)-length(b2))>=2  && b2(1)<b(1) && b2(2)<b(1) && i==1
                    % in root tip some cells start oscillating later, so
                    % this is to process them:
                        b2(1)=[];
                        if abs(length(b)-length(b2))<2
                            counter3(k)=counter3(k)+1;
                            pos{i}.index(j)=max([length(b) length(b2)]);
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
                clear b b2 s        
                mm= min([length(pkTime) length(troughTime)]);
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
                                indTtake=min(find(pos{i}.time(troughTime(1:end-1))>72));
                                pos{i}.Ampsmid(j)=diffs_forward(indTtake);  
                                pos{i}.Ampsfor(j, 1:2)=diffs_forward([ 1 end]); 
                                diffs_forward=[diffs_forward; diffs_back]; 
                                pos{i}.val(j)=mean(diffs_forward);
                                pos{i}.val2(j)=mean(diffs_back);   
                                Amps_for=[Amps_for diffs_forward' ];
                                Amps_back=[Amps_back diffs_back' ];   
                                pos{i}.Ampsback(j, 1:2)=diffs_back([ 1 end]); 
                                pos{i}.WCV_for(j)=std(diffs_forward)/mean(diffs_forward);
                                pos{i}.WCV_back(j)=std(diffs_back)/mean(diffs_back);
                                pos{i}.ppmeans(j)=mean(diff(pkTrueTimes)); %mean period per cell
                                pos{i}.CVp(j)=std(diff(pkTrueTimes))/mean(diff(pkTrueTimes));
                                pos{i}.pkTrueTimesN{j}=pkTrueTimes; 
                                counter4(k)=counter4(k)+1;
                            elseif pkTime(1)>troughTime(1)
                                diffs_forward=abs(pos{i}.YFP(troughTime,j)-pos{i}.YFP(pkTime,j));
                                diffs_back=abs(pos{i}.YFP(troughTime(2:end),j)-pos{i}.YFP(pkTime(1:end-1),j));
                                indTtake=min(find(pos{i}.time(troughTime)>72));
                                pos{i}.Ampsmid(j)=diffs_forward(indTtake);
                                pos{i}.Ampsfor(j, 1:2)=diffs_forward([ 1 end]);
                                diffs_forward=[diffs_forward; diffs_back]; 
                                pos{i}.val(j)=mean(diffs_forward);
                                pos{i}.val2(j)=mean(diffs_back);   
                                Amps_for=[Amps_for diffs_forward' ];
                                Amps_back=[Amps_back diffs_back' ];
                                pos{i}.Ampsback(j, 1:2)=diffs_back([ 1 end]); 
                                pos{i}.WCV_for(j)=std(diffs_forward)/mean(diffs_forward);
                                pos{i}.WCV_back(j)=std(diffs_back)/mean(diffs_back);
                                pos{i}.ppmeans(j)=mean(diff(pkTrueTimes)); %mean period per cell
                                pos{i}.CVp(j)=std(diff(pkTrueTimes))/mean(diff(pkTrueTimes));    
                                pos{i}.pkTrueTimesN{j}=pkTrueTimes; 
                                counter4(k)=counter4(k)+1;       
                            end
                          elseif length(pkTime)>length(troughTime) && length(troughTime)==length(pkTime)-1 && length(pkTime)>1
                            if pkTime(1)<troughTime(1)
                                diffs_forward=abs(pos{i}.YFP(troughTime,j)-pos{i}.YFP(pkTime(2:end),j));
                                diffs_back=abs(pos{i}.YFP(troughTime,j)-pos{i}.YFP(pkTime(1:end-1),j));
                                indTtake=min(find(pos{i}.time(troughTime)>72));
                                pos{i}.Ampsmid(j)=diffs_forward(indTtake);  
                                pos{i}.Ampsfor(j, 1:2)=diffs_forward([ 1 end]); 
                                diffs_forward=[diffs_forward; diffs_back]; 
                                pos{i}.val(j)=mean(diffs_forward);
                                pos{i}.val2(j)=mean(diffs_back);     
                                Amps_for=[Amps_for diffs_forward'];
                                Amps_back=[Amps_back diffs_back']; 
                                pos{i}.Ampsback(j, 1:2)=diffs_back([ 1 end]); 
                                pos{i}.WCV_for(j)=std(diffs_forward)/mean(diffs_forward);
                                pos{i}.WCV_back(j)=std(diffs_back)/mean(diffs_back); 
                                pos{i}.ppmeans(j)=mean(diff(pkTrueTimes)); %mean period per cell
                                pos{i}.CVp(j)=std(diff(pkTrueTimes))/mean(diff(pkTrueTimes));
                                pos{i}.pkTrueTimesN{j}=pkTrueTimes;
                                counter4(k)=counter4(k)+1;  
                            else
                                pos{i}.val(j)=NaN;
                                pos{i}.val2(j)=NaN; 
                                Amps_for=[Amps_for NaN ];
                                Amps_back=[Amps_back NaN ];
                                pos{i}.Ampsfor(j, 1:2)=NaN; 
                                pos{i}.Ampsmid(j)=NaN;
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
                                indTtake=min(find(pos{i}.time(troughTime(1:end-1))>72));
                                pos{i}.Ampsmid(j)=diffs_forward(indTtake);
                                pos{i}.Ampsfor(j, 1:2)=diffs_forward([ 1 end]); 
                                diffs_forward=[diffs_forward; diffs_back]; 
                                pos{i}.val(j)=mean(diffs_forward);
                                pos{i}.val2(j)=mean(diffs_back);
                                Amps_for=[Amps_for diffs_forward' ];
                                Amps_back=[Amps_back diffs_back' ];
                                pos{i}.Ampsback(j, 1:2)=diffs_back([ 1 end]); 
                                pos{i}.WCV_for(j)=std(diffs_forward)/mean(diffs_forward);
                                pos{i}.WCV_back(j)=std(diffs_back)/mean(diffs_back);
                                pos{i}.ppmeans(j)=mean(diff(pkTrueTimes)); %mean period per cell
                                pos{i}.CVp(j)=std(diff(pkTrueTimes))/mean(diff(pkTrueTimes));
                                pos{i}.pkTrueTimesN{j}=pkTrueTimes; 
                                counter4(k)=counter4(k)+1;  
                            else
                                pos{i}.val(j)=NaN;
                                pos{i}.val2(j)=NaN; 
                                pos{i}.WCV_for(j)=NaN;
                                pos{i}.WCV_back(j)=NaN;
                                Amps_for=[Amps_for NaN ];
                                Amps_back=[Amps_back NaN ];
                                pos{i}.Ampsfor(j, 1:2)=NaN; 
                                pos{i}.Ampsmid(j)=NaN;
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
                        pos{i}.Ampsfor(j, 1:2)=NaN; 
                        pos{i}.Ampsback(j, 1:2)=NaN; 
                        pos{i}.Ampsmid(j)=NaN;
                        pos{i}.ppmeans(j)=NaN; %mean period per cell
                        pos{i}.CVp(j)=NaN;
                        pos{i}.pkTrueTimesN{j}=NaN; 
                    end   
                else
                    pos{i}.val(j)=NaN;
                    pos{i}.val2(j)=NaN; 
                    pos{i}.WCV_for(j)=NaN;
                    pos{i}.WCV_back(j)=NaN;
                    Amps_for=[Amps_for NaN];
                    Amps_back=[Amps_back NaN];       
                    pos{i}.Ampsfor(j, 1:2)=NaN; 
                    pos{i}.Ampsmid(j)=NaN;
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
                     % MinPeakDist and MinPeakDist2 replaced by 18 
                    findpeaks(dtdata, pos{i}.time(1:xs),  'MinPeakDistance',16, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkProm)
                    [a2,b2]=findpeaks(-dtdata, pos{i}.time(1:xs),   'MinPeakDistance',16, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkPromTr);
                    title(['Cell=' num2str(j) '  Pos=' num2str(i)])
                    hold on
                    plot( b2, -a2, 'ro')
                    xlabel('Time (h)')
                    ylabel('Levels of detrended time series')  
                    if isnan(pos{i}.val(j))==1 && isnan(pos{i}.val2(j))==1
                        clf
                    % MinPeakDist and MinPeakDist2 replaced by 18 
                        [a,b]= findpeaks(dtdata, pos{i}.time(1:xs),  'MinPeakDistance',16, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkProm);
                        [a2,b2]=findpeaks(-dtdata, pos{i}.time(1:xs),   'MinPeakDistance',16, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkPromTr);
                        plot(pos{i}.time(1:xs),pos{i}.YFP(1:xs,j),'-', 'Linewidth',3)
                        hold on 
                        clear pkxs trxs
                        pkxs=NaN(1,length(b));trxs=NaN(1,length(b));
                        for s=1: length(b)
                            pkxs(s)=find(pos{i}.time==b(s));
                        end
                        clear s
                        for s=1: length(b2)
                            trxs(s)=find(pos{i}.time==b2(s));
                        end
                        plot(pos{i}.time(pkxs),pos{i}.YFP(pkxs,j),'ko', 'Linewidth',3)
                        plot(pos{i}.time(trxs),pos{i}.YFP(trxs,j),'ro', 'Linewidth',3)
                        ylabel('Levels of time series') 
                        print(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/notenough/Cell_' num2str(j) '.png'],'-dpng','-r0')    
                    else
                        clf
                        [~,b]= findpeaks(dtdata, pos{i}.time(1:xs),  'MinPeakDistance',MinPeakDist, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkProm);
                        [a2,b2]=findpeaks(-dtdata, pos{i}.time(1:xs),   'MinPeakDistance',MinPeakDist2, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkPromTr);
                        plot(pos{i}.time(1:xs),pos{i}.YFP(1:xs,j),'-', 'Linewidth',3)
                        hold on
                        clear pkxs trxs
                        for s=1: length(b)
                            pkxs(s)=find(pos{i}.time==b(s));
                        end
                        clear s
                        for s=1: length(b2)
                            trxs(s)=find(pos{i}.time==b2(s));
                        end
                        plot(pos{i}.time(pkxs),pos{i}.YFP(pkxs,j),'ko', 'Linewidth',3)
                        plot(pos{i}.time(trxs),pos{i}.YFP(trxs,j),'ro', 'Linewidth',3)
                        ylabel('Levels of time series')     
                        print(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/good/Cell_' num2str(j) '.png'],'-dpng','-r0')
                        clf
                        clear s pkss BF1 AF1 
                        [BF1,AF1]=butter(3,0.5); 
                        dtdata2= filtfilt(BF1,AF1,pos{i}.YFP(:,j));
                        [~,b]= findpeaks(dtdata2, pos{i}.time(1:xs),  'MinPeakDistance',16, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkProm); 
                        set(gcf,'DefaultLineLineWidth',0.5)
                        pkss=Nan(1,6); 
                        for s=1: length(b)
                            pkss(s)=find(pos{i}.time==b(s));
                        end
                        plot(pos{i}.time(pkss),pos{i}.YFP(pkss,j),'ko', 'Linewidth',1,'MarkerFaceColor', 'k')
                        hold on
                        findpeaks(pos{i}.YFP(1:xs,j), pos{i}.time(1:xs),  'MinPeakDistance',16, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkProm); 
                        clear a b pkss
                        [a,b]= findpeaks(dtdata, pos{i}.time(1:xs),  'MinPeakDistance',16, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkProm); 
                        set(gcf,'DefaultLineLineWidth',0.5)
                        for s=1: length(b)
                            pkss(s)=find(pos{i}.time==b(s));
                        end
                        plot(pos{i}.time(pkss),pos{i}.YFP(pkss,j),'rd', 'Linewidth',1)
                        hold on
                        if exist(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/good/check/'],'file')~=7
                            mkdir(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/good/check/']);
                        end
                        print(['Doneanalysis/' name  '/figures/'  'individual_traces/dist' num2str(mpkdis2) '/pos' num2str(i) '/good/check/Acell_' num2str(j) '.png'],'-dpng','-r0')
                    end
                end
                clf
            else
                pos{i}.val(j)=NaN;
                pos{i}.val2(j)=NaN; 
                Amps_for=[Amps_for NaN ];
                Amps_back=[Amps_back NaN ];
                pos{i}.Ampsfor(j, 1:2)=NaN; 
                pos{i}.Ampsmid(j)=NaN;
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
            end
            clf
        end
        Amps_forward{i}=Amps_for;
        Amps_backward{i}= Amps_back; 
    end
    save(['Doneanalysis/' name  '/summary/dist' num2str(mpkdis2) '.mat'], 'pos', 'counter', 'counter2', 'counter3','counter4', 'Amps_forward', 'Amps_backward');
else
    a=load([ 'Doneanalysis/' name  '/summary/dist' num2str(mpkdis2) '.mat']);
    pos=a.pos;
    counter=a.counter;
    counter2=a.counter2;
    counter3=a.counter3;
    counter4=a.counter4;
    Amps_forward= a.Amps_forward;
    Amps_backward = a.Amps_backward;
    clear a;  
end
%% Mean trace statistics
for k=start_pos:numpos 
    i=posi(k);
    clear r c
    [r,c]=find(~isnan(pos{i}.val));
    pos{i}.meanlevs=mean(pos{i}.YFP(:,c),2);
end

val=NaN(1,numpos);val2=NaN(1,numpos);

for j=1:numpos
    i=posi(j);
    data_F1(:,i)= pos{i}.meanlevs; 
    mpkdis= 10;
    [a,b]=findpeaks((data_F1(:,i)),  'MinPeakDistance',mpkdis, 'MinPeakProminence',15);  
    [a2,b2]=findpeaks(-(data_F1(:,i)),  'MinPeakDistance',mpkdis, 'MinPeakProminence',15); 
    figure(1)
    findpeaks((data_F1(:,i)),  'MinPeakDistance',mpkdis);
    hold on
    findpeaks(pos{i}.meanlevs,  'MinPeakDistance',mpkdis);
    clf
    figure(2)
    findpeaks(data_F1(:,i), 'MinPeakDistance',10, 'MinPeakProminence',15);
    [ac2,bc2]=findpeaks(-data_F1(:,i),   'MinPeakDistance',10, 'MinPeakProminence',15); %troughsbar
    title([' Mean of Pos=' num2str(i)])
    hold on
    plot( bc2, -ac2, 'ro')
    xlabel('Time (h)')
    ylabel('Levels of detrended mean time series')  
    print(['Doneanalysis/' name  '/figures/pos' num2str(i)  '__forampidentification.png'],'-dpng','-r0')    
    close all
    % need to measure amplitude by two ways trough to peak and peak to trough. 
    if length(b)==length(b2) % b is peak and b2 is trough
        if b(1)<b2(1) 
            diffs_forward=abs(pos{i}.meanlevs(b2(1:end-1),:)-pos{i}.meanlevs(b(2:end),:));
            diffs_back=abs(pos{i}.meanlevs(b2,:)-pos{i}.meanlevs(b,:));
            indTtake=min(find(pos{i}.time(troughTime(1:end-1))>72));          
            pos{i}.Ampsmid(j)=diffs_forward(indTtake);
            pos{i}.AmpsforMean( 1:2)=diffs_forward([ 1 end]); 
            diffs_forward=[diffs_forward; diffs_back]; 
            val(i)=mean(diffs_forward);
            val2(i)=mean(diffs_back);
            pos{i}.AmpsbackMean( 1:2)=diffs_back([ 1 end]); 
        elseif b(1)>b2(1)
            diffs_forward=abs(pos{i}.meanlevs(b2,:)-pos{i}.meanlevs(b,:));
            diffs_back=abs(pos{i}.meanlevs(b2(2:end),:)-pos{i}.meanlevs(b(1:end-1),:));
            pos{i}.AmpsforMean( 1:2)=diffs_forward([ 1 end]); 
            diffs_forward=[diffs_forward; diffs_back]; 
            val(i)=mean(diffs_forward);
            val2(i)=mean(diffs_back);   
            pos{i}.AmpsbackMean( 1:2)=diffs_back([ 1 end]); 
        end   
    elseif length(b)>length(b2)
      if b(1)<b2(1)
          diffs_forward=abs(pos{i}.meanlevs(b2,:)-pos{i}.meanlevs(b(2:end),:));
          diffs_back=abs(pos{i}.meanlevs(b2,:)-pos{i}.meanlevs(b(1:end-1),:));
          pos{i}.AmpsforMean( 1:2)=diffs_forward([ 1 end]); 
          diffs_forward=[diffs_forward; diffs_back]; 
          val(i)=mean(diffs_forward);
          val2(i)=mean(diffs_back); 
          pos{i}.AmpsbackMean( 1:2)=diffs_back([ 1 end]); 
      end
    elseif length(b)<length(b2)
      if b(1)>b2(1)
          diffs_forward=abs(pos{i}.meanlevs(b2(1:end-1),:)-pos{i}.meanlevs(b,:));
          diffs_back=abs(pos{i}.meanlevs(b2(2:end),:)-pos{i}.meanlevs(b,:));
          pos{i}.AmpsforMean( 1:2)=diffs_forward([ 1 end]); 
          diffs_forward=[diffs_forward; diffs_back];   
          val(i)=mean(diffs_forward);
          val2(i)=mean(diffs_back);
          pos{i}.AmpsbackMean( 1:2)=diffs_back([ 1 end]); 
      end
    end
    clear a b a2 b2 data_F1
end
%% Figures to be plotted:
% Fig 2B
% Fig 1 Supplement 5
% Fig 2A
% Fig 2B
%% Figure 2B
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1800 3600]/300)
maxylim=4000;
for j=1:4
    set(gca, 'Fontsize', 7)
    i=posi(j);
    subplot(4,1,5-j)
    hold on 
    y = [0 0 maxylim maxylim];
    x = [12 24 24 12];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
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
    ylim([0 4000])  
    set(gca, 'LineWidth', 1.2)
    set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2); 
    xlim([29 168])  
    clear x y 
    hold on
    if j==1
        cmap=colormap(gray(size(pos{j}.YFP,2)));
        for k= 1:size(pos{j}.YFP,2)
            if k~=6 && mod(k,2)==0
                plot(pos{j}.time, pos{j}.YFP(:,k), 'Color', cmap(randi([1 size(pos{j}.YFP,2)]), :)*diag([1 1 1]), 'Linewidth', 0.05);
                hold on 
                set(gca, 'LineWidth', 1.2)
                set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
            end
        end
        ylim([0 4000]);
    elseif j==2
        concat=[pos{2}.YFP pos{3}.YFP];           
        cmap= colormap(gray(size(concat,2)));
        for k= 1:size(concat,2)
                if mod(k,2)==0
                    if k<=size(pos{2}.YFP,2)
                        plot(pos{2}.time, concat(:,k),  'Color', cmap(randi([1 size(concat,2)]), :)*diag([1 1 1]), 'Linewidth', 0.05)
                        hold on 
                        set(gca, 'LineWidth', 1.2)
                        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
                    else
                            plot( pos{3}.time, concat(:,k),  'Color', cmap(randi([1 size(concat,2)]), :)*diag([1 1 1]), 'Linewidth', 0.05)
                  hold on 
                  set(gca, 'LineWidth', 1.2)
    set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
                    end
                end
        end                
        ylim([0 2000]);
    elseif j==3
        clear concat
        concat=[pos{4}.YFP pos{5}.YFP];
        lev=size(concat,2)*2;
        cc=redblue(lev);
        set(gca, 'LineWidth', 1.2)
        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
        cmap= cc(1:lev/2,:);
        for k= 1:size(concat,2)
            if mod(k,2)==0
                if k<=size(pos{4}.YFP,2)
                    plot(pos{4}.time, concat(:,k),  'Color', cmap(randi([1 lev/2]), :)*diag([1 1 1]), 'Linewidth', 0.05)
                    hold on 
                else
                    plot( pos{5}.time, concat(:,k),  'Color', cmap(randi([1 lev/2]), :)*diag([1 1 1]), 'Linewidth', 0.05)
                    hold on 
                end
             end
         end
        set(gca, 'LineWidth', 1.2)
        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
        ylim([0 2000]);
     elseif j==4
         lev=size(pos{6}.YFP,2)*2;
         clear cc
         cc=redblue(lev);
         cmap= cc(end:-1:lev/2+1,:);
         for k= 1:size(pos{6}.YFP,2)
            if mod(k,2)==0
                plot(pos{6}.time, pos{6}.YFP(:,k), 'Color', cmap(randi([1 lev/2]), :)*diag([1 1 1]), 'Linewidth', 0.05);
                hold on 
            end
            set(gca, 'LineWidth', 1.2)
            set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
         end  
         ylim([0 2000]);
    end   
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
end
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure2B.pdf','file')~=2
    print('-Painters',kkq, 'Figure2B','-dpdf','-r300')
    movefile('Figure2B.pdf', 'Doneanalysis/WT/figures')
end
close all
%% Background we will subtract: 
backg=xlsread([pwd '/Data_singlecell/WT_final_coordinates/WT_pos1-6_RFP_YFP_BF_median_tracked_background_Plot.xlsx']);
meanBG=  mean(backg(1:120, 2:end));
%% Figure 1 Supplement 5 
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1800*3.6 3600*2.3]/300)
maxylim=4000;
for j=1:10
    subplot(22,2,j*2-1)
    hold on 
    y = [0 0 maxylim maxylim];
    x = [12 24 24 12];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
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
    if  j==5
        ylabel('Scaled Intensity', 'Fontsize', 7) 
    end
    xlim([29 168])  
    clear x y 
    hold on
    clear concat  
    concat=[pos{1}.YFP]-meanBG(1);   
    concat(concat<0) = 0;
    cmap=colormap(gray(size(concat,2)));
    lev=size(concat,2)*2;
    clear numcells
    numcells=randi(size(concat,2) , 1, 5); 
    for k= 1:size(numcells,2)
        clear mix
        if numcells(k)<=size(concat,2)
            mix=max(concat(:,numcells(k)));
            plot(pos{1}.time, concat(:,numcells(k))/mix,  'Color', cmap(randi([1 lev/2-40]), :)*diag([1 1 1]), 'Linewidth', 0.05)
            hold on 
            set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
        end
    end
    ylim([0 1]);
    xlim([29 168])  
    if j==1
        title({'\fontsize{7}ROOT TIP',['\fontsize{5} Cells: ' num2str(numcells)]})
    else
        title(['Cells:' num2str(numcells)],'Fontsize', 5)
    end      
    ax = gca;  
    if j==10
        ax.XTick = 0:24:220;
    else
       set(gca,'xtick',[])
    end
    if j==10
        xlabel('Time (h)', 'Fontsize', 7) 
    end
    hold on
end
set(gca, 'LineWidth', 1.2)
rng(1)
for j=1:10
    subplot(22,2,j*2-1+20+2+2 )
    hold on 
    y = [0 0 maxylim maxylim];
    x = [12 24 24 12];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
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
    if j==5
       ylabel('Scaled Intensity', 'Fontsize', 7) 
    end  
    xlim([29 168])  
    clear x y 
    hold on
    clear concat
    concat=[pos{2}.YFP-meanBG(2) pos{3}.YFP-meanBG(3)];
    concat(concat<0) = 0;
    cmap=colormap(gray(size(concat,2)));
    lev=size(concat,2)*2;
    clear numcells     
    numcells=randi(size(concat,2) , 1, 5); 
    for k= 1:size(numcells,2)
        clear mix
        if numcells(k)<=size(concat,2)
            mix=max(concat(:,numcells(k)));
            if numcells(k)<=size(pos{2}.YFP,2)
                plot(pos{2}.time, concat(:,numcells(k))/mix,  'Color', cmap(randi([1 lev/2-40]), :)*diag([1 1 1]), 'Linewidth', 0.05)
                hold on 
                set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
            elseif numcells(k)>size(pos{2}.YFP,2) && numcells(k)<=size(concat,2)
                plot( pos{3}.time, concat(:,numcells(k))/mix,  'Color', cmap(randi([1 lev/2-40]), :)*diag([1 1 1]), 'Linewidth', 0.05)
                hold on 
                set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
            end
        end
    end
    ylim([0 1]);
    xlim([29 168])  
    if j==1
        title({'\fontsize{7}ROOT',['\fontsize{5} Cells:' num2str(numcells)]})
    else
        title(['Cells:' num2str(numcells)],'Fontsize', 5)
    end      
    ax = gca;  
    if j==10
        ax.XTick = 0:24:220;
    else
       set(gca,'xtick',[])
    end
    if j==10
        xlabel('Time (h)', 'Fontsize', 7) 
    end
    hold on
end
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
rng(1)
for j=1:10
    subplot(22,2,j*2)
    hold on 
    y = [0 0 maxylim maxylim];
    x = [12 24 24 12];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
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
    if  j==5
        ylabel('Scaled Intensity', 'Fontsize', 7) 
    end
    xlim([29 168])
    clear x y 
    hold on   
    clear concat
    concat=[pos{4}.YFP-meanBG(4) pos{5}.YFP-meanBG(5)];
    concat(concat<0) = 0;
    cmap=colormap(gray(size(concat,2)));
    lev=size(concat,2)*2;
    clear numcells
    numcells=randi(size(concat,2) , 1, 5);              
    for k= 1:size(numcells,2)
        clear mix
        if numcells(k)<=size(concat,2)
            mix=max(concat(:,numcells(k)));
            if numcells(k)<=size(pos{4}.YFP,2)
                plot(pos{4}.time, concat(:,numcells(k))/mix,  'Color', cmap(randi([1 lev/2-40]), :)*diag([1 1 1]), 'Linewidth', 0.05)
                hold on 
            elseif numcells(k)>size(pos{4}.YFP,2) && numcells(k)<=size(concat,2)
                plot( pos{5}.time, concat(:,numcells(k))/mix,  'Color', cmap(randi([1 lev/2-40]), :)*diag([1 1 1]), 'Linewidth', 0.05)
                hold on 
                set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);    
            end
        end
    end
    ylim([0 1]);
    xlim([29 168])  
    if j==1
        title({'\fontsize{7}HYPOCOTYL',['\fontsize{5} Cells: ' num2str(numcells)]})
    else                   
        title(['Cells:' num2str(numcells)],'Fontsize', 5)
    end      
    ax = gca;  
    if j==10
        ax.XTick = 0:24:220;
    else
       set(gca,'xtick',[])
    end
    if j==10
        xlabel('Time (h)', 'Fontsize', 7) 
    end
    hold on
    set(gca, 'Fontsize', 7)
end
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
for j=1:10
    subplot(22,2,j*2-1+20+2+1+2)
    hold on 
    y = [0 0 maxylim maxylim];
    x = [12 24 24 12];
    patch(x,y,[0.4,0.4,0.4], 'Facealpha',0.1, 'EdgeColor','none')
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
    if  j==5
        ylabel('Scaled Intensity', 'Fontsize', 7) 
    end
    xlim([29 168])  
    clear x y 
    hold on
    clear concat
    concat=[pos{6}.YFP]-meanBG(6);
    concat(concat<0) = 0;
    cmap=colormap(gray(size(concat,2)));
    lev=size(concat,2)*2;
    clear numcells
    numcells=randi(size(concat,2) , 1, 5); 
    for k= 1:size(numcells,2)
        clear mix
        if numcells(k)<=size(concat,2)
            mix=max(concat(:,numcells(k)));
            plot(pos{6}.time, concat(:,numcells(k))/mix,  'Color', cmap(randi([1 lev/2-40]), :)*diag([1 1 1]), 'Linewidth', 0.05)
            hold on 
            set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
        end
    end
    ylim([0 1]);
    xlim([29 168])  
    if j==1
        title({'\fontsize{7}COTYLEDON',['\fontsize{5} Cells:' num2str(numcells)]})
    else
        title(['Cells: ' num2str(numcells)],'Fontsize', 5)
    end       
    ax = gca;  
    if j==10
        ax.XTick = 0:24:220;
    else
       set(gca,'xtick',[])
    end
    if j==10
        xlabel('Time (h)', 'Fontsize', 7) 
    end
    hold on
end
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure1Supp5.pdf','file')~=2
    print('-Painters',kkq, 'Figure1Supp5','-dpdf','-r300')
    movefile('Figure1Supp5.pdf', 'Doneanalysis/WT/figures')
end
close all
%% We want to lump data from sections togetehr if sections are from the same tissue. 
% first since data is sampled at different times we need to interpolate it
% to uniform timepoints.
alltimes=[];
for j=1:numpos
    i=posi(j);
    alltimes=[alltimes pos{i}.time];
end
alltimes=sort(alltimes, 'ascend');
%now interpolate only over the needed region
for k=1:numpos
    i=posi(k);
    for j=1: size(pos{i}.YFP,2)
        pos{i}.YFPinterp(:,j)=interp1(pos{i}.time, pos{i}.YFP(:,j), alltimes);%when outside the range, gives nan
    end
end
% Now lump pos2 and pos3 (root) together nad pos4 and pos5 (hypocotyl) together 
YFPint{1}=pos{1}.YFPinterp;
YFPint{2}=[pos{2}.YFPinterp pos{3}.YFPinterp];
YFPint{3}=[pos{4}.YFPinterp pos{5}.YFPinterp];
YFPint{4}=[pos{6}.YFPinterp];
YFPintot=[YFPint{1} YFPint{2} YFPint{3} YFPint{4}];
%% Figure 1G
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1600 3600]/300)
maxylim=1200;
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
    ylabel('Mean Fluor.', 'Fontsize', 7) ;
    ylim([0 maxylim]) ; 
    xlim([29 168]);  
    clear x y 
    x=alltimes;
    if j==1
        y=nanmean(YFPint{j}, 2);
        plot(x,y, 'Color',[0.6,0.6,0.6]);
        set(gca, 'LineWidth', 1.2);
        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
    elseif j==2
        y=nanmean(YFPint{j}, 2);
        y([1:2 717:720])=nan; 
        plot(x,y, 'k');
        set(gca, 'LineWidth', 1.2);
        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
    elseif j==3
        y=nanmean(YFPint{j}, 2);
        y([1:4 719:720])=nan; 
        plot(x,y, 'b')
        set(gca, 'LineWidth', 1.2)
        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
    elseif j==4
        y=nanmean(YFPint{j}, 2); 
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
end
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure1G.pdf', 'file')~=2
    print('-Painters',kkq, 'Figure1G','-dpdf','-r300')
    movefile('Figure1G.pdf', 'Doneanalysis/WT/figures')
end
close all
%% Percentage rhythmic
xxp=NaN(1,6);
for k=1:6
   i=posi(k); 
   xxp(k)=counter(k)/size(pos{i}.periods,1)*100; %this is percentage of cells that have been used for further analysis
end
%% Figure 2A
Bigcombo=NaN([  6 size((pos{1}.val(pos{1}.index<100))') ]); 
for k=1:6
  Bigcombo(k,1:size((pos{k}.val(pos{k}.index<100))'))=(pos{k}.val(pos{k}.index<100));
  Bigcombo(k,size((pos{k}.val(pos{k}.index<100))')+1:end)=nan;   
end
Bigcombo=Bigcombo';
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  2000 1600]/300)
axes('FontSize',7);
Bigcombox=NaN( length((pos{1}.val(pos{1}.index<100))), 6 ); 
for i=1:numpos
    Bigcombox(:,i)=Bigcombo(:,i)*1./val(i);
end
bh=boxplot(Bigcombox, 'symbol','');
p = prctile(Bigcombox,[9 91]);
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
secto={ 'Root tip';'Lower root';'Upper root';'Lower hypocotyl';'Upper hypocotyl';'Cotyledon'};
ylabel('Normalised amplitude levels')
set(gca, 'XTick', 1:6, 'XTickLabel', secto,'FontSize',7);
hold on
for nums=1:numpos
    line(nums-0.35:0.7:nums+0.35, [1 1], 'Linewidth', 2, 'Color', 'g')
end
ylim([0 5])
box off
xlim([0.5 6+0.5])
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
ax=gca;
ax.XTickLabelRotation=25;
xlim([0.5 6+0.5])
ylim([0 6])
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);

if exist('Figure2A.pdf', 'file')~=2
    print('-Painters',kkq, 'Figure2A','-dpdf','-r300')
    movefile('Figure2A.pdf', 'Doneanalysis/WT/figures')
end
%% Figure 2B
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
%calculating of between and within period variability
btwn=NaN(1,6);
witn=NaN(1,6);
for k=start_pos:numpos
    i=posi(k);   
    btwn(k)=nanstd(pos{i}.ppmeans)/nanmean(pos{i}.ppmeans);
    witn(k)=nanmean(pos{i}.CVp);
    s=s+1;
end
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  2000 1200]/300)
h=bar([btwn' witn'],'BarWidth', 1);
% Figure labels, fonts, limits
set(gca, 'FontSize', 7)
set(gca, 'XTick', 1:numpos, 'XTickLabel', cell(6,1));
set(h(1),'facecolor','k') 
set(h(2),'facecolor','w') 
h(2).LineWidth=1;
h=legend({' Between-cell', ' Within-cell'}, 'Position', [0.75,0.9,0.01,0.01], 'box', 'off');
box off
ylabel({' Period variability (CV)'})
xlim([0.5 6+0.5])
ylim([0 0.225])
set(gca, 'YTick', 0:0.05:2, 'YTickLabel', 0:0.05:2);
secto={ 'Root tip';'Lower root';'Upper root';'Lower hypocotyl';'Upper hypocotyl';'Cotyledon'};
set(gca, 'XTick', 1:6, 'XTickLabel', secto);
set(gca, 'Fontsize', 7)
box off
ax=gca; 
ax.XTickLabelRotation=25;
xlim([0.5 6+0.5])
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure2B.pdf', 'file')~=2
    print('-Painters',kkq, 'Figure2B','-dpdf','-r300')
    movefile('Figure2B.pdf', 'Doneanalysis/WT/figures')
end
close all
%% Figure 1 supp 6 A
posnew{1}.Ampsfor=pos{1}.Ampsfor;
posnew{2}.Ampsfor=[pos{2}.Ampsfor;pos{3}.Ampsfor];
posnew{3}.Ampsfor=[pos{4}.Ampsfor;pos{5}.Ampsfor];
posnew{4}.Ampsfor=pos{6}.Ampsfor;
posnew{1}.Ampsmid=pos{1}.Ampsmid';
posnew{2}.Ampsmid=[pos{2}.Ampsmid';pos{3}.Ampsmid'];
posnew{3}.Ampsmid=[pos{4}.Ampsmid';pos{5}.Ampsmid'];
posnew{4}.Ampsmid=pos{6}.Ampsmid';
sectnew={ {'Root tip'},{'Root (rest)'}, {'Hypocotyl'},{'Cotyledon'}};
kkq=figure;
axes('FontSize',7);
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  5400 3200]/300)
for i=1:4
    subplot(1,4,i)  
    posnew{i}.AmpNew=[ posnew{i}.Ampsmid posnew{i}.Ampsfor(:,2)];
    bh=boxplot(posnew{i}.AmpNew, 'symbol','');
    p = prctile(posnew{i}.AmpNew,[9 91]);
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
    ylim([0 1450])
    set(gca, 'LineWidth', 1.2)
    ylabel('Amplitude levels')
elseif i==2
    ylim([0 900])
    set(gca, 'LineWidth', 1.2)
elseif i==3
    ylim([0 900])
    set(gca, 'LineWidth', 1.2)
elseif i==4
    ylim([0 900])
    set(gca, 'LineWidth', 1.2)
    end
    set(gca, 'Fontsize',7);
    title(sectnew{i})
end
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
if exist('Figure1Supp6A.pdf', 'file')~=2
    print('-Painters',kkq, 'Figure1Supp6A','-dpdf','-r300')
    movefile('Figure1Supp6A.pdf', 'Doneanalysis/WT/figures')
end
NumberNotUsed=zeros(1,4); 
All=zeros(1,4);
for i=1:4 
    NumberNotUsed(i)=sum(isnan(posnew{i}.AmpNew(:,1)));
    All(i)=size(posnew{i}.AmpNew,1);
end
Percentages=(All-NumberNotUsed)./All*100; %percentage of all cells used in the analysis in the last boxplot. 
close all