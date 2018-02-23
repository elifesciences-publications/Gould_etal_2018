
%%  File for analysis of Long period data (CCA1-long) for Gould, Domijan et al.
%  "Coordination of robust single cell rhythms in the Arabidopsis circadian 
%  clock via spatial waves of gene expression" (BioRxiv)
%  doi https://doi.org/10.1101/208900s
% 
%  by Mirela Domijan (U. of Liverpool)
clear all 
%%  
name='CCA1-long'; %user chosen
%% first add in useful paths
addpath([pwd '/Data_singlecell/'])
addpath([pwd '/Data_singlecell/' name '_final_coordinates'])
% to plot Figure 2Supplemenatries2 G you can use the additional toolbox which is
% available from https://uk.mathworks.com/matlabcentral/fileexchange/30884-controllable-tight-subplot
%in the case, uncomment the line below.
addpath([pwd '/subplot_tight/']) 
%% is there a folder where to save all the data?
if exist(['Doneanalysis/' name  '/summary'],'file')~=7
mkdir(['Doneanalysis/' name  '/summary']);
end
%% is there a folder where to save all the figures?
if exist(['Doneanalysis/' name  '/figures'],'file')~=7
mkdir(['Doneanalysis/' name  '/figures']);
end
%% pull out all the period estimates and the data
[aa, texts]=xlsread('Data_singlecell/CCA1-long_final_coordinates/interpolated/Individual results_1h.xlsx');
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

posi=12:18;
for j=1:length(posi)
xs=csvread(['Data_singlecell/CCA1-long_final_coordinates/interpolated/interpolateddata_CCA1-long_pos' num2str(posi(j)) '_spacing1h.csv']);

if j==6
    separators=[separators size(xs,2)-2]; %for pos17 Pete files don'es havt he first trace of pos17 in there for period analysis.
else
separators=[separators size(xs,2)-1];
end
if j==1
    savepos12=xs;
elseif j==2
     savepos13=xs;   
     elseif j==3
     savepos14=xs;   
     elseif j==4
     savepos15=xs;   
     elseif j==5
     savepos16=xs;   
elseif j==6
    savepos17=xs;
elseif j==7
    savepos18=xs;
end
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


%this is to sort out which traces are wrong in pos12: 

[aa2, texts2]=xlsread('Data_singlecell/CCA1-long_final_coordinates/interpolated/metamorph_CCA1-long_1hr_interp.xlsx');

aa2(:, 1:6)=[];
aa2(1,:)=[];

clear texts2

pos12=aa2(:, 1:274); 

rows12torem= find(sum(abs(pos12))==0);


pos12(:, find(sum(abs(pos12)) == 0)) = [];

check12=isequal(pos12, savepos12(:, 2:end));


%rows that we have to remove from the pp etc are: 

 zerorows=find(sum(abs(aa2)) == 0); %*3;
% 
% zerosrow3=zerorows-2;
% 
%  rowsremoving=[];
% for i=1:length(zerorows)
%     rowsremoving=[rowsremoving zerosrow3(i): zerorows(i)];
% end

pp(zerorows,:)=[];
ppgof(zerorows,:)=[];
ppignore(zerorows,:)=[];


for j=1:length(posi)
    k=posi(j);
   
    num1=separates(j)+1;
    num2=separates(j+1);
   pos{k}.periods=pp(num1:num2,:);
   pos{k}.periodsgofs=ppgof(num1:num2,:);
   pos{k}.periodsignore=ppignore(num1:num2,:);
    
end

clearvars -except pos name savepos17 savepos18
%% determine which positions you are looking at
if strcmp(name,'CCA1-long')==1
    posi=[12:18];
    sect={ {'root tip', ''},{'root', '(sect 1)'}, {'root', '(sect 2)'},{'Root/Hyp' , ''}, {'hypocotyl', '(sect 1)'},{'hypocotyl', '(sect 2)'},{'cotyledon', ''}};
end
%% now which positions you will look at 
start_pos=1;
numpos=7;
%% for each file you pull out, write down the positions; the YFP levels and the periods. 
for k=  start_pos:numpos
%     bb=load(a(k).name);
%    %two types of periods: (i) periods(:,1)== FFT unscaled; (ii) periods(:,2)== FFT unscaled RAE; (iii) periods(:,9)== LS acceptance 10%
     i=posi(k); 
%    % i=14; 
%     pos{i}.periods=bb. periods;
%     
     if strcmp(name, 'CCA1-long')==1
      name2=['pos',num2str(i)];
%     name3=['FFT' name2];
%     xa=xlsread('Summary of Period analysis.xls', name3);
%     % Note: xa is the excelt sheet of the FFT periods from Pete
%     % third column are  FFT periods, fifth column are RAEs (called GOFs in the sheet)
%     pos{i}.periods(:,1)=nan;
%     pos{i}.periods(:,2)=nan;
%     xa(isnan(xa(:,1)),:)=[];
%       
%     pos{i}.periods(xa(:,1),1)=xa(:,3); %FFT period
%     pos{i}.periods(xa(:,1),2)=xa(:,5);
%     
%      
%     clear xa

        posfile=dir([pwd , '/Data_singlecell/' name '_final_coordinates/',['pos_' name2 '*tracked.csv']]);

    fid= fopen( posfile.name,'r');
    C=textscan(fid, repmat('%s',1,10), 'delimiter',',', 'CollectOutput',true);
       
    
    if i==12
    rems=[];
    for ink=1:size(C{1})
        if strcmp(C{1}(ink,8), '1000000119')==1 || strcmp(C{1}(ink,8), '1000000128')==1 || strcmp(C{1}(ink,8), '1000000233')==1 || strcmp(C{1}(ink,8), '1000000036')==1
            rems=[rems ink];
        end
    end
   
    C{1}(rems, :)=[];
    end
    
    
    
     if i==17
    rems=[];
    for ink=1:size(C{1})
        if strcmp(C{1}(ink,8), '1000000000')==1
            rems=[rems ink];
        end
    end
   
    C{1}(rems, :)=[];
     end
    
     
    if i==18
    rems=[];
    for ink=1:size(C{1})
        if strcmp(C{1}(ink,8), '1000000271')==1
            rems=[rems ink];
        end
    end
   
    C{1}(rems, :)=[];
    end
    
  triple=str2double(C{1}(3:end,1:3));
  LL=logical(diff(str2double(C{1}(3:end,7))));

  fa=find(LL==1);
 % Note on Labels:
 % triple is the x,y,z triple of the coordinates. 
 % LL keeps clear of when the timing changes 
 % fa identifies where the new time point begins
 % res reshapes the positions so every three columns (i.e. x,y,z) you have
 % a new time point i.e, dimensions are no. of coordinates* time points. 
 % num_cells is the number of cells
 % cell_labels is the labels of the cells 
 
%  if i==16
%      LL(1800:end)=[];
%      triple(1801:end,:)=[];
%      fa(end)=[];
%  end

    if numel(unique(diff([0; fa])))==1 & (size(LL,1)-fa(end))== (fa(1)-1) 
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
    YFPfile=dir([pwd , '/Data_singlecell/'  name '_final_coordinates/',['mean*' name2 '*YFP_*.csv']]);  
    fid= fopen( YFPfile.name,'r');
    clear C    
    if i==12
        num_cells2=num_cells+4;
    elseif i==18
        num_cells2=num_cells+1; 
    elseif i==17
        num_cells2=num_cells+2;
    else 
        num_cells2=num_cells;
    end
    C=textscan(fid, repmat('%s',1,num_cells2+1), 'delimiter',',', 'CollectOutput',true);
    if i==17
        pos{i}.YFP=str2double(C{1}(4:end,3:end-1));
    else
        pos{i}.YFP=str2double(C{1}(4:end,2:end));
    end
    if i==12
        pos{i}.YFP(1:2:end, :)=[];
        pos{i}.YFP( :, [ 120 129 234 38]) =[];
        pos{i}.YFP( :, [37]) =[]; %in the interpolation file this trace is not there. 
        pos{i}.x( :, [37]) =[];
        pos{i}.y( :, [37]) =[];
        pos{i}.z( :, [37]) =[];   
    end

    
  if i==18
      pos{i}.YFP(:, 272)=[];   %cell 1000000271 is missing
  end
    
    cell_labels2=str2double(C{1}(3,2:end));
    
    % verify that everything is pulled out correctly and make sure the cells are in increasing order
    % repeat of above just wiht the new YFP file. 
    diffy2=diff(cell_labels2)-1;
    diffywo2 = diffy2(diffy2~=0);
    allok2(i)=isempty(diffywo2);
    end
 clear cell_labels2  fid YFPfile diffy2 diffywo2 num_cells
% time points 
if strcmp(name, 'CCA1-long')==1
[Times,txt2]=xlsread(['timestamp_CCA1-long.xlsx'], 'B8:AI26');
pos{i}.time= Times(i,:);
end
  clear data FFTrecs periods
end
clear i k  txt2 posfile bb a Times C name2 name3 
clear i k  txt2 posfile bb a Times C name2 name3 

gof_cutoff=0.9;
gof_cutoff2=1;
%% Add in the pad positions  (i.e. connect x and y positions across different square pads)
addpath('pictures');

E26=imread('CCA1-long.png');
imshow(E26)
E26(:, 1:1650, :)=[];
E26(:, 2300:end, :)=[];
E26(1:300, :, :)=[];
E26(2300:end, :, :)=[];
imshow(E26)

bw=im2bw(E26, 0.001);
imshow(bw)
bw2= imfill(bw,'holes');
imshow(bw2)


%make a square: 
imshow(bw2(1500:end, 1700:end))
bw3=bw2(1500:end, 1700:end);
st = regionprops(bw3, 'BoundingBox' );
Col=colormap('jet');

for k=1:size(st,1)
 rectangle('Position',[st(k).BoundingBox(1),st(k).BoundingBox(2),st(k).BoundingBox(3),st(k).BoundingBox(4)],...
'EdgeColor',Col(k*5, :), 'LineWidth',5)
hold all
end

dd=503; 
bb=double(bw2); 
sqs=ones(dd,dd); 

Ba=conv2(bb,sqs,'same');
[ind1,ind2]=find(imregionalmax(Ba)>0);
plot(ind2, ind1, 'o')
set(gca,'YDir','Reverse')
%resizing the whole deal 
vals=354;
factor= vals/ dd; 
J = imresize(bw2,factor);
stJ = regionprops(J, 'BoundingBox' );
Col=colormap('jet');
imshow(J);
for k=1:size(stJ,1)
rectangle('Position',[stJ(k).BoundingBox(1),stJ(k).BoundingBox(2),stJ(k).BoundingBox(3),stJ(k).BoundingBox(4)],...
'EdgeColor',Col(k*2, :), 'LineWidth',5)
hold all
end
bbJ=double(J); 
sqsJ=ones(vals,vals); 
BaJ=conv2(bbJ,sqsJ,'same');
[indJ1,indJ2]=find(imregionalmax(BaJ)>0);
plot(indJ2, indJ1, 'o')
set(gca,'YDir','Reverse') 
ULC=[indJ2-vals/2 indJ1-vals/2];
[sa, sb]=(sort(indJ1, 'ascend'));
indJ1=indJ1(sb);
indJ2=indJ2(sb);
indJ1(6)=[];
indJ2(6)=[];
xx=indJ1([5 4]);
yy=indJ2([5 4]);
indJ1(4:5)=xx;
indJ2(4:5)=yy;
for k=start_pos:numpos
    i=posi(k);
    pos{i}.xmod=indJ2(k)+ (pos{i}.x);
    pos{i}.ymod=indJ1(k)+ (vals-pos{i}.y);  
    plot3(pos{i}.ymod(1, :), pos{i}.xmod(1, :),pos{i}.z(1, :), '.','Markersize',10)
end
% pos{i}.xmod is the correct x position for all cells (taking into account
% square pad info. 
% pos{i}.ymod is the correct x position for all cells (taking into account
% square pad info. 
for k=start_pos:  numpos
    i=posi(k);
    pos{i}.xmod=indJ2(k)+ (pos{i}.x);
    pos{i}.ymod=indJ1(k)+ (vals-pos{i}.y);
    plot(pos{i}.xmod(end, :), pos{i}.ymod(end, :),'.','Markersize',10)
    hold on
     rectangle('position', [indJ2(k) indJ1(k) vals vals])
end
daspect([1 1 1])
%% Now get all data together so you can plot all sections together
xxs=[];
yys=[];
zzs=[];
traces=[];
times=[];
xxrs=[];
yyrs=[];
tracers=[];
timers=[];
yyrxs=[];
tracerxs=[];
yyxs=[];
tracexs=[];
coloring=[];
pers=[];
pergofs=[];
perignore=[];
persmF=[];
pergofsmF=[];
perignoremF=[];
persSR=[];
pergofsSR=[];
perignoreSR=[];
cmap = colormap(jet(128));
cmap=fliplr(cmap);
siz=size(cmap,1);
posix=[  12    13    14    16    17    18];
for k=1: numpos %for posi not posix
    i=posi(k); %posix(k);
    if k==1
    xxs=[ xxs pos{i}.xmod];
    yys=[ yys pos{i}.ymod];
     zzs=[ zzs pos{i}.z];
    traces=[ traces pos{i}.YFP];
    pers=[pers; pos{i}.periods(:,1)];
    pergofs=[pergofs; pos{i}.periodsgofs(:,1)];
    perignore=[perignore; pos{i}.periodsignore(:,1)];
     persmF=[persmF; pos{i}.periods(:,2)];
    pergofsmF=[pergofsmF; pos{i}.periodsgofs(:,2)];
    perignoremF=[perignoremF; pos{i}.periodsignore(:,2)];
     persSR=[persSR; pos{i}.periods(:,3)];
    pergofsSR=[pergofsSR; pos{i}.periodsgofs(:,3)];
    perignoreSR=[perignoreSR; pos{i}.periodsignore(:,3)];
    count(k)=size(pos{i}.YFP,2);
    countkeep(k)=size(pos{i}.YFP,2);
    times=[times repmat(pos{i}.time',1,(size(pos{i}.YFP,2)))];
    elseif k==2
    [ik, ja]=find(pos{i}.ymod<583.65);   
    dontkeep=unique(ja);
    keep=(1:size(pos{i}.xmod,2));
    keep(dontkeep)=[];
    clear ik ja bb  
    xxs=[ xxs pos{i}.xmod(:, keep)];
    yys=[ yys pos{i}.ymod(:, keep)];  
    zzs=[ zzs pos{i}.z(:, keep)];
    traces=[ traces pos{i}.YFP(:, keep)];
    times=[times repmat(pos{i}.time',1,(size(pos{i}.YFP(:, keep),2)))];
    pers=[pers; pos{i}.periods(keep,1)];
    pergofs=[pergofs; pos{i}.periodsgofs(keep,1)];
    perignore=[perignore; pos{i}.periodsignore(keep,1)];
    persmF=[persmF; pos{i}.periods(keep,2)];
    pergofsmF=[pergofsmF; pos{i}.periodsgofs(keep,2)];
    perignoremF=[perignoremF; pos{i}.periodsignore(keep,2)];
    persSR=[persSR; pos{i}.periods(keep,3)];
    pergofsSR=[pergofsSR; pos{i}.periodsgofs(keep,3)];
    perignoreSR=[perignoreSR; pos{i}.periodsignore(keep,3)];
    count(k)=size(pos{i}.YFP,2);
    countkeep(k)=size(pos{i}.YFP(:,keep),2);
    elseif k==3  
    [ik, ja]=find(pos{i}.ymod <923.63);   
     dontkeep=unique(ja);
    keep=(1:size(pos{i}.xmod,2));
    keep(dontkeep)=[];
    clear ik ja bb  
    xxs=[ xxs pos{i}.xmod(:, keep)];
    yys=[ yys pos{i}.ymod(:, keep)];  
    zzs=[ zzs pos{i}.z(:, keep)];
    traces=[ traces pos{i}.YFP(:, keep)];
    times=[times repmat(pos{i}.time',1,(size(pos{i}.YFP(:, keep),2)))];
    pers=[pers; pos{i}.periods(keep,1)];
    pergofs=[pergofs; pos{i}.periodsgofs(keep,1)];
    perignore=[perignore; pos{i}.periodsignore(keep,1)];
    persmF=[persmF; pos{i}.periods(keep,2)];
    pergofsmF=[pergofsmF; pos{i}.periodsgofs(keep,2)];
    perignoremF=[perignoremF; pos{i}.periodsignore(keep,2)];
    persSR=[persSR; pos{i}.periods(keep,3)];
    pergofsSR=[pergofsSR; pos{i}.periodsgofs(keep,3)];
    perignoreSR=[perignoreSR; pos{i}.periodsignore(keep,3)];
    count(k)=size(pos{i}.YFP,2);
    countkeep(k)=size(pos{i}.YFP(:,keep),2);
    elseif k==4
    [ik, ja]=find(pos{i}.ymod <max(max(pos{14}.ymod)) & pos{i}.xmod <max(max(pos{14}.xmod)));  %how originally done  
    dontkeep=unique(ja);
    clear ik ja bb 
    keep=(1:size(pos{15}.xmod,2));
    keep(dontkeep)=[];
    clear ik ja bb 
    xxs=[ xxs pos{i}.xmod(:, keep)];
    yys=[ yys pos{i}.ymod(:, keep)];  
    zzs=[ zzs pos{i}.z(:, keep)];
    traces=[ traces pos{i}.YFP(:, keep)];
    times=[times repmat(pos{i}.time',1,(size(pos{i}.YFP(:, keep),2)))];
    pers=[pers; pos{i}.periods(keep,1)];
    pergofs=[pergofs; pos{i}.periodsgofs(keep,1)];
    perignore=[perignore; pos{i}.periodsignore(keep,1)];
    persmF=[persmF; pos{i}.periods(keep,2)];
    pergofsmF=[pergofsmF; pos{i}.periodsgofs(keep,2)];
    perignoremF=[perignoremF; pos{i}.periodsignore(keep,2)];
    persSR=[persSR; pos{i}.periods(keep,3)];
    pergofsSR=[pergofsSR; pos{i}.periodsgofs(keep,3)];
    perignoreSR=[perignoreSR; pos{i}.periodsignore(keep,3)];
    count(k)=size(pos{i}.YFP,2);
    countkeep(k)=size(pos{i}.YFP(:,keep),2);
    elseif k==6
    [ik, ja]=find(pos{i}.ymod <max(max(pos{16}.ymod)) & pos{i}.xmod <max(max(pos{16}.xmod)));  %how originally done  
    dontkeep=unique(ja);
    clear ik ja bb 
    keep3=(1:size(pos{i}.xmod,2));
    keep3(dontkeep)=[];
    clear ik ja bb 
    [ik, ja]=find(pos{i}.ymod >min(min(pos{i+1}.ymod)) & pos{i}.xmod >min(min(pos{i+1}.xmod)));  %how originally done  
    dontkeep2=unique(ja);
    clear ik ja bb 
    keep2=(1:size(pos{i}.xmod,2));
    keep2(dontkeep2)=[];
    clear ik ja bb    
    keep=intersect(keep3, keep2);
    xxs=[ xxs pos{i}.xmod(:, keep)];
    yys=[ yys pos{i}.ymod(:, keep)];  
    zzs=[ zzs pos{i}.z(:, keep)];
    traces=[ traces pos{i}.YFP(:, keep)];
    times=[times repmat(pos{i}.time',1,(size(pos{i}.YFP(:, keep),2)))];
    pers=[pers; pos{i}.periods(keep,1)];
    pergofs=[pergofs; pos{i}.periodsgofs(keep,1)];
    perignore=[perignore; pos{i}.periodsignore(keep,1)];
    persmF=[persmF; pos{i}.periods(keep,2)];
    pergofsmF=[pergofsmF; pos{i}.periodsgofs(keep,2)];
    perignoremF=[perignoremF; pos{i}.periodsignore(keep,2)];
    persSR=[persSR; pos{i}.periods(keep,3)];
    pergofsSR=[pergofsSR; pos{i}.periodsgofs(keep,3)];
    perignoreSR=[perignoreSR; pos{i}.periodsignore(keep,3)];
    count(k)=size(pos{i}.YFP,2);
    countkeep(k)=size(pos{i}.YFP(:,keep),2);
    else
    xxs=[ xxs pos{i}.xmod];
    yys=[ yys pos{i}.ymod];
    zzs=[ zzs pos{i}.z];
    traces=[ traces pos{i}.YFP];
    pers=[pers; pos{i}.periods(:,1)];
    pergofs=[pergofs; pos{i}.periodsgofs(:,1)];
    perignore=[perignore; pos{i}.periodsignore(:,1)];
    persmF=[persmF; pos{i}.periods(:,2)];
    pergofsmF=[pergofsmF; pos{i}.periodsgofs(:,2)];
    perignoremF=[perignoremF; pos{i}.periodsignore(:,2)];
    persSR=[persSR; pos{i}.periods(:,3)];
    pergofsSR=[pergofsSR; pos{i}.periodsgofs(:,3)];
    perignoreSR=[perignoreSR; pos{i}.periodsignore(:,3)];
    count(k)=size(pos{i}.YFP,2);
    times=[times repmat(pos{i}.time',1,(size(pos{i}.YFP,2)))];
    countkeep(k)=size(pos{i}.YFP,2);
    end
    if i==15
        sizing=length(pers);     
    end 
    if i==17
        sizing2=length(pers);
     end 
   cc=0;
end 
xxs_adjusted=xxs;yys_adjusted=yys;zzs_adjusted=zzs*2/0.692;
%% Which tissues are we looking at:
str={{'Root tip', ''},{' Root', '(sect 1)'}, {' Root', '(sect 2)'},{'Root/', 'Hypocotyl'}, {'Hypocotyl', '(sect 1)'},{'Hypocotyl', '(sect 2)'},{'Cotyledon', ''}}; 

%% Making of the Figures:
% below you can plot for Figure 2 Supplement 2. 
% Fig 2 SI2 D
% Fig 2 SI2 C
% Fig 2 SI2 E
% Fig 2 SI2 F
% Fig 2 SI2 H
% Fig 2 SI2 J
%% Figure2 SI2 D
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1700 2400]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
   
        col2=interp1([min(min(traces)) max(max(traces))],[1 siz],max(traces(:,k)));
        colr= interp1(1:siz,flipud(cmap),col2);
        hold on
        plot( max(traces(:,k)), yys_adjusted(end,k),'.','Markersize',6, 'Color', colr)
        clear col2 colr
      end
end
ylim([indJ1(1) 1800]) %all three lines removed MD2018
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) %all three lines removed MD2018
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5)) %MD 2018. 
xlim([min(min(traces)) max(max(traces))])
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('Maximum expression level ','Fontsize', 7)
print('-Painters',kkq, 'Figure2Supp2D','-dpdf','-r300')
movefile('Figure2Supp2D.pdf', 'Doneanalysis/CCA1-long/figures')
%% Figure2 SI2 C: y positive vs. max value. 
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  3000 2400]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
   
        col2=interp1([min(min(traces)) max(max(traces))],[1 siz],max(traces(:,k)));
        colr= interp1(1:siz,flipud(cmap),col2);
        hold on
        plot( xxs_adjusted(end,k), yys_adjusted(end,k),'.','Markersize',6, 'Color', colr)
        clear col2 colr
      end
end
ylim([indJ1(1) 1800]) %all three lines removed MD2018
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) %all three lines removed MD2018
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5)) %MD 2018.
set(gca,'Xtick',-60:500:6500,'XTickLabel',0:0.5:6.5) %all three lines removed MD2018
 xlim([-60 1590])

set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('x (mm) ','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);


for k=start_pos:numpos
 hold on
     rectangle('position', [indJ2(k) indJ1(k) vals vals], 'Linewidth',0.15)
     set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if k==3 
text((indJ2(k)-200),(indJ1(k)+150),str{k}, 'Fontsize', 7)
elseif k==4
text((indJ2(k)+80),(indJ1(k)+430),str{k}, 'Fontsize', 7)
elseif k==6
text((indJ2(k)+400),(indJ1(k)+120),str{k}, 'Fontsize', 7)
elseif k==7
text((indJ2(k)+355),(indJ1(k)+150),str{k}, 'Fontsize', 7)
else
    text((indJ2(k)+400),(indJ1(k)+150),str{k}, 'Fontsize', 7)
end
end
colormap(jet(128));
caxis([(floor(min(min(traces)))) (ceil(max(max(traces))))]);
h = colorbar( 'Location', 'north','XTick',[(floor(min(min(traces)))) (ceil(max(max(traces))))], 'TickLabels', {'low', 'high'}, 'FontSize',7, 'Position',[[0.65 0.2 0.2 0.01]]);

print('-Painters',kkq, 'Figure2Supp2C','-dpdf','-r300')
movefile('Figure2Supp2C.pdf', 'Doneanalysis/CCA1-long/figures')
%% Figure2 SI2 E:
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  3000 2400]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(xxs_adjusted(end,k), yys_adjusted(end,k),'.','Markersize',6, 'Color', colr)
        clear col2 colr
      end
end
ylim([indJ1(1) 1800]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5)) 
set(gca,'Xtick',0:500:6500,'XTickLabel',0:0.5:6.5) 
 xlim([-60 1590])

colormap(jet(128));
caxis([18 35]);
h = colorbar( 'Location', 'north', 'Direction', 'reverse', 'XTick',[18 35],'TickLabels',{'35h' '18h'}, 'FontSize',7, 'Position', [[0.4 0.05 0.2 0.01]]);


for k=start_pos:numpos
 hold on
     rectangle('position', [indJ2(k) indJ1(k) vals vals], 'Linewidth',0.15)
     set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if k==3 
text((indJ2(k)-200),(indJ1(k)+150),str{k}, 'Fontsize', 7)
elseif k==4
text((indJ2(k)+80),(indJ1(k)+430),str{k}, 'Fontsize', 7)
elseif k==6
text((indJ2(k)+400),(indJ1(k)+120),str{k}, 'Fontsize', 7)
elseif k==7
text((indJ2(k)+355),(indJ1(k)+150),str{k}, 'Fontsize', 7)
else
text((indJ2(k)+400),(indJ1(k)+150),str{k}, 'Fontsize', 7)
end
end
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('x (mm)','Fontsize', 7)

print('-Painters',kkq, 'Figure2Supp2E','-dpdf','-r300')
movefile('Figure2Supp2E.pdf', 'Doneanalysis/CCA1-long/figures')

close all
%% Figure2 SI2 F: y positive vs. max value. 
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1700 2400]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
         col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(pers(k), yys_adjusted(end,k),'.','Markersize',6, 'Color', colr)
        clear col2 colr    
      end
end
ylim([indJ1(1) 1800]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5))
xlim([18 35])
set(gca,'Xtick',16:4:36,'XTickLabel',(16:4:36)) 
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('Period (h) ','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
print('-Painters',kkq, 'Figure2Supp2F','-dpdf','-r300')
movefile('Figure2Supp2F.pdf', 'Doneanalysis/CCA1-long/figures')

%% Figure2 SI3 C: y positive vs. z periods. 
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1400 3600]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
         col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(zzs_adjusted(end,k), yys_adjusted(end,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr     
      end
end
ylim([indJ1(1) 1800]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5)) 
set(gca,'Xtick',0:50:10500,'XTickLabel',(0:0.05:10.5)) 
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('z (mm) ','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);

print('-Painters',kkq, 'Figure2Supp3C','-dpdf','-r300')
movefile('Figure2Supp3C.pdf', 'Doneanalysis/CCA1-long/figures')
close all

%% Figure2 SI3 D: y positive vs. z periods (only roots)
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1400 3600]/300)
cc=[];
allrootcells=sum(countkeep(1:3));
for k=1:allrootcells
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(zzs_adjusted(end,k), yys_adjusted(end,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr     
      end
end
ylim([indJ1(1) 1300]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5)) 
set(gca,'Xtick',0:50:10500,'XTickLabel',(0:0.05:10.5)) 
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('z (mm) ','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
print('-Painters',kkq, 'Figure2Supp3D','-dpdf','-r300')
movefile('Figure2Supp3D.pdf', 'Doneanalysis/CCA1-long/figures')
close all
%% now sort by position in:
[thwo, reloc]=sort(yys(1, :), 'ascend'); 
xxs=xxs(:,reloc); 
yys=yys(:,reloc); 
traces=traces(:,reloc); 
pers=pers(reloc);
pergofs=pergofs(reloc);
perignore=perignore(reloc);
persmF=persmF(reloc);
pergofsmF=pergofsmF(reloc);
perignoremF=perignoremF(reloc);
persSR=persSR(reloc);
pergofsSR=pergofsSR(reloc);
perignoreSR=perignoreSR(reloc);
times=times(:, reloc); 

%% Info for the filmstrip. 
cellused=sum(countkeep(1:4));
tracesmin=traces(:, 1:cellused)- repmat(min(traces(:, 1:cellused)),size(traces(:, 1:cellused),1),1);
tracesmin=tracesmin*diag(1./(max(tracesmin)));
% what do we want to calculate: normalized or raw? 
traces2=traces(:, 1:cellused); 
cmap = colormap(hot); % colour map to use.
siz=size(cmap,1);
s=[ min(min(traces2)) max(max(traces2))]; 

dd6=zeros(size(traces2 ,1), size(traces2,2),3);
full6=zeros(size(traces2,1), size(traces2,2));
for i=1:size(traces2,1)
    for j=1:size(traces2,2)
        full6(i,j)=interp1(s,[1 siz],traces2(i,j));
         dd6(i,j,:)=interp1(1:siz,cmap,full6(i,j));
    end
    
end
traces2x=tracesmin; 
cmap = colormap(hot); % colour map to use.
siz=size(cmap,1);
s=[ min(min(traces2x)) max(max(traces2x))]; 

dd6x=zeros(size(traces2x ,1), size(traces2x,2),3);
full6x=zeros(size(traces2x,1), size(traces2x,2));
for i=1:size(traces2x,1)
    for j=1:size(traces2x,2)
        full6x(i,j)=interp1(s,[1 siz],traces2x(i,j));
         dd6x(i,j,:)=interp1(1:siz,cmap,full6x(i,j));
    end
end
%% Filmstrip; Figure SI2 Fig 2, G.
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  5800 700]/300)
func_hndl={@subplot;@subplot_tight};
for  ti=1:1:17%
       set(gcf, 'Color', 'black')

set(gcf, 'Color', 'black')
 func_hndl{2}(1,17,ti, 0.001); 
 cellused=sum(countkeep(1:4));
   for i=1: sum(countkeep(1:4)) 
              set(gcf, 'Color', 'k')
   plot(xxs_adjusted(ti,i)',yys_adjusted(ti,i)', '.','Markersize',full6x(ti,i)/20+2,'Color',(dd6x(ti,i,1:3)))
  hold on;  
    grid off;
      set(gca, 'color', 'k')
       if ti==1
       line( [max(max(xxs_adjusted(:,1:sum(countkeep(1:4)))))-50 max(max(xxs_adjusted(:,1:sum(countkeep(1:4)))))-50], [min(min(yys_adjusted(:,1:sum(countkeep(1:4)))))+50 min(min(yys_adjusted(:,1:sum(countkeep(1:4)))))+50+250],  'LineStyle','-','Linewidth', 2, 'Color', 'white')
    end
 
    end

    daspect([1.5 1.5 1.5])
    camlight;lighting gouraud;
    set(gca, 'Fontsize', 12)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlim([min(min(xxs_adjusted(:, 1:cellused))) max(max(xxs_adjusted(:, 1:cellused)))])
    ylim([min(min(yys_adjusted(:, 1:cellused))) max(max(yys_adjusted(:, 1:cellused)))])
    hold on
    set(gcf, 'Color', 'k')
end
set(gcf,'color','w');
set(gcf, 'InvertHardCopy', 'off');
print('-Painters',kkq, 'Figure2Supp2G_toppanel','-dpdf','-r300')
movefile('Figure2Supp2G_toppanel.pdf', 'Doneanalysis/CCA1-long/figures')
close all

kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  5800 700]/300)
func_hndl={@subplot;@subplot_tight};
for  ti=18:1:34
    set(gcf, 'Color', 'black')
    set(gcf, 'Color', 'black')
    func_hndl{2}(1,17,ti-17, 0.003); 
   for i=1: sum(countkeep(1:4))
    set(gcf, 'Color', 'k')
    plot(xxs_adjusted(ti,i)',yys_adjusted(ti,i)', '.','Markersize',full6x(ti,i)/20+2,'Color',(dd6x(ti,i,1:3)))
    hold on;  
    grid off;
    set(gca, 'color', 'k')
   end
    daspect([1.5 1.5 1.5])
    camlight;lighting gouraud;
    set(gca, 'Fontsize', 12)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlim([min(min(xxs_adjusted(:, 1:cellused))) max(max(xxs_adjusted(:, 1:cellused)))])
    ylim([min(min(yys_adjusted(:, 1:cellused))) max(max(yys_adjusted(:, 1:cellused)))])
    hold on
    set(gcf, 'Color', 'k')
end
set(gcf,'color','w');
set(gcf, 'InvertHardCopy', 'off');
print('-Painters',kkq, 'Figure2Supp2G_lowerpanel','-dpdf','-r300')
movefile('Figure2Supp2G_lowerpanel.pdf', 'Doneanalysis/CCA1-long/figures')
close all
