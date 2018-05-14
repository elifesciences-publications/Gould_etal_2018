%  File for analysis of WT data for Gould, Domijan et al.
%  "Coordination of robust single cell rhythms in the Arabidopsis circadian 
%  clock via spatial waves of gene expression" (BioRxiv)
%  doi https://doi.org/10.1101/208900s
% 
%  by Mirela Domijan (U. of Liverpool)

clear all 
%%  
name='WT'; %Experiment to run. 
%% add in relevant paths to the data
addpath([pwd '/Data_singlecell/'])
addpath([pwd '/Data_singlecell/' name '_final_coordinates'])
%% is there a folder where to save all the analysis?
if exist(['Doneanalysis/' name  '/summary'],'file')~=7
mkdir(['Doneanalysis/' name  '/summary']);
end
%% is there a folder where to save all the figures?
if exist(['Doneanalysis/' name  '/figures'],'file')~=7
mkdir(['Doneanalysis/' name  '/figures']);
end
%% pull out all the period estimates and the data
[aa, texts]=xlsread([ 'Data_singlecell/' name '_final_coordinates/interpolated/Individual results.xlsx']);
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

posi=1:6;
for j=1:length(posi)
xs=csvread( ['Data_singlecell/' name '_final_coordinates/interpolated/interpolateddata_WT_pos' num2str(posi(j)) '.csv']);
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
    posi=[1:6];
    sect={ {'root tip', ''},{'root', '(sect 1)'}, {'root', '(sect 2)'},{'hypocotyl', '(sect 1)'},{'hypocotyl', '(sect 2)'},{'cotyledon', ''}};
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
        allok=isequal(abs(diffywo),ones(size(res,2)-1,1)*num_cells);
        % allok tells us that there are right number of cells within each time point 
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
    % repeat of above just with the new YFP file. 
    diffy2=diff(cell_labels2)-1;
    diffywo2 = diffy2(diffy2~=0);
    allok2=isempty(diffywo2);
    end
    
   % clear everything not needed: 
   clear cell_labels2 allok allok2 fid YFPfile diffy2 diffywo2 num_cells
% time points 
if strcmp(name, 'WT')==1
[Times,txt2]=xlsread(['timestamp_WT.xlsx'], 'C8:DR15');
pos{i}.time= Times(i,:);
elseif strcmp(name, 'WTrepeat')==1
    
[Times,txt2]=xlsread(['timestamp_WTrepeat.xlsx'], 'B8:AG27');
pos{i}.time= Times(i,2:end);
end
clear data FFTrecs periods
end

clear i k  txt2 posfile bb a Times C name2 name3 

% goodness of Fit cut-offs for accepting data as periodic (used in period
% analysis)
gof_cutoff=0.9;
gof_cutoff2=1;
%% Add in the pad positions 

addpath('pictures');

E26=imread('WT.png');
imshow(E26)
E26(:, 1:1200, :)=[];
E26(:, 800:end, :)=[];
imshow(E26)

bw=im2bw(E26, 0.1);
imshow(bw)
bw2= imfill(bw,'holes');
imshow(bw2)

[B,L] = bwboundaries(bw2);
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
end

bw2=bw2(200:end, :);
bw2(1350:end, :)=[];
imshow(bw2)
st = regionprops(bw2, 'BoundingBox' );
Col=colormap('jet');
for k= 1:size(st,1)
 rectangle('Position',[st(k).BoundingBox(1),st(k).BoundingBox(2),st(k).BoundingBox(3),st(k).BoundingBox(4)],...
'EdgeColor',Col(k*5, :), 'LineWidth',5)
hold all
end

dd=min(st(3).BoundingBox(3:4));

% squares are of dim dd: 
bb=double(bw2); 
sqs=ones(dd,dd); 

Ba=conv2(bb,sqs,'same');
C = imfuse(bw2,imregionalmax(Ba),'blend','Scaling','joint');
imshow(C); 

bx=bw2(143:152, 232:283);

C=corner(bx); 
imshow(bx);
hold on
plot(C(:,1), ones(1,length(C(:,1))), '*', 'Color', 'r')

nums=unique(sort(C(:,1)))-min(unique(sort(C(:,1))))+1;% first pixel must value 1; 

pixs=(nums(end)-nums(3));


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
'EdgeColor',Col(k*5, :), 'LineWidth',5)
hold all
end
for k=3
 rectangle('Position',[stJ(k).BoundingBox(1),stJ(k).BoundingBox(2),stJ(k).BoundingBox(3),stJ(k).BoundingBox(4)],...
'EdgeColor','r', 'LineWidth',5)
hold all
end

bbJ=double(J); 
sqsJ=ones(floor(vals),floor(vals)); %values are they natural numbers?

BaJ=conv2(bbJ,sqsJ,'same');
CJ = imfuse(J,imregionalmax(BaJ),'blend','Scaling','joint');
imshow(CJ); 

[indJ1,indJ2]=find(imregionalmax(BaJ)>0);

hold on
plot(indJ2, indJ1, 'o')
set(gca,'YDir','Reverse')
   
ULC=[indJ2-vals/2 indJ1-vals/2];

[sa, sb]=(sort(indJ1, 'ascend'));
plot(ULC(:,1), ULC(:,2), 'go')

indJ1=indJ1(sb);
indJ2=indJ2(sb);

clf
% % periods on the figure
for k= start_pos:numpos
    i=posi(k);
    pos{i}.xmod=indJ2(k)-vals/2+ (pos{i}.x);
    pos{i}.ymod=indJ1(k)-vals/2+ (vals-pos{i}.y);
    
    plot( pos{i}.xmod(1, :),pos{i}.ymod(1, :), '.','Markersize',10)
    hold on
    rectangle('Position',[indJ2(k)-vals/2 indJ1(k)-vals/2  vals vals])
    daspect([1 1 1])
end


% periods on the figure
for k=start_pos:numpos
    i=posi(k);
    keep=[];
    for ss=1: size(pos{i}.ymod,2) 
        if  pos{i}.periodsgofs(ss, 1)<gof_cutoff && pos{i}.periodsgofs(ss, 2)<gof_cutoff && pos{i}.periodsgofs(ss, 3)<gof_cutoff2  && abs(pos{i}.periods(ss, 1)- pos{i}.periods(ss, 3))<2.5 && abs(pos{i}.periods(ss, 2)- pos{i}.periods(ss, 3))<2.5 && abs(pos{i}.periods(ss, 2)-pos{i}.periods(ss, 1))<2.5 && sum(pos{i}.periodsignore(ss, 1)+pos{i}.periodsignore(ss, 2)+pos{i}.periodsignore(ss, 3))==0 
         keep=[ ss keep];
        
        end
     end 
   pos{i}.di1= [pos{i}.xmod(1, keep); pos{i}.ymod(1, keep); pos{i}.z(1, keep)]; 
   pos{i}.Dist= pdist(  pos{i}.di1');
   pos{i}.Distx= pdist(  pos{i}.di1(1, :)');
   pos{i}.Disty= pdist(  pos{i}.di1(2,:)');
   pos{i}.Distz= pdist(  pos{i}.di1(3,:)');
   pos{i}.DistP=pdist(pos{i}.periods(keep,1)); 
   pos{i}.distancing= pos{i}.DistP./pos{i}.Dist;

end
%% Which tissues are we looking at: 
str={'Root tip', 'Lower root', 'Upper root', 'Lower hypocotyl', 'Upper hypocotyl', 'Cotyledon'};
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

counter=1;
for k= 1:  numpos 
     i=posi(k);
     if k==1 || k==6 || k==4 
    traces=[ traces pos{i}.YFP];
    xxs=[ xxs pos{i}.xmod];
    yys=[ yys pos{i}.ymod];
    zzs=[zzs pos{i}.z];
    times=[times repmat(pos{i}.time',1,(size(pos{i}.YFP,2)))];
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
     elseif k==2  
         
             [ik, ja]=find(pos{i}.ymod<max(max(pos{k-1}.ymod)));   
         dontkeep=unique(ja);
               keep=(1:size(pos{i}.ymod,2));
    keep(dontkeep)=[];
         clear ik ja bb 
     xxs=[ xxs pos{i}.xmod(:, keep)];
    yys=[ yys pos{i}.ymod(:, keep)];  
    zzs=[zzs pos{i}.z(:, keep)];
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
       [ik, ja]=find(pos{i}.ymod >min(min(pos{k+1}.ymod)));   
        dontkeep=unique(ja);
        keep=(1:size(pos{i}.ymod,2));
        keep(dontkeep)=[];
         clear ik ja bb   
       xxs=[ xxs pos{i}.xmod(:, keep)];
        yys=[ yys pos{i}.ymod(:, keep)];  
        zzs=[zzs pos{i}.z(:, keep)];
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
    elseif k==5 
    [ik,ja]=find(pos{i}.ymod <max(max(pos{4}.ymod)) & pos{i}.xmod<max(max(pos{4}.xmod)));  
    dontkeep=unique(ja); 
    keep=(1:size(pos{5}.xmod,2));
    keep(dontkeep)=[];
    clear ik ja bb 
    
    xxs=[ xxs pos{i}.xmod(:, keep)];
    yys=[ yys pos{i}.ymod(:, keep)];  
    zzs=[zzs pos{i}.z(:, keep)];
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
     end    
     counter=counter+1;    
end    

xxs_adjusted=xxs;yys_adjusted=yys;zzs_adjusted=zzs*2/0.692;

%% Making of the Figures:
% below you can plot: 
% Fig 2 F and G:
% Fig 3 A
% Fig 3 B
% Fig 3 C
% Fig 3 D
% Fig 3 F,G,H and
% Fig 3E
%% Figure 2F: scatterplots in Fig F and G:
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1400 3600]/300)
cc=[];
for k=1:size(pers,1)
     % only plot data which is treated as periodic (conditions outlined in  the paper)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        %determine the colour scheme for the cell (from min to max vaues of all traces
        %below is original in the paper
        col2=interp1([min(min(traces)) max(max(traces))],[1 siz],max(traces(:,k)));
        colr= interp1(1:siz,flipud(cmap),col2);
        hold on
        plot( max(traces(:,k)), yys_adjusted(end,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
      end
end
%% Figure 2F sizing and saving options: 
ylim([indJ1(1)-vals/2 2700]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1)-vals/2:500:10500,'YTickLabel',(0:0.5:10.5))  
xlim([min(min(traces)) max(max(traces))])
% make sure lines and axes  are thick enough when plotting.
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
% Fonts
set(gca,'Fontsize', 7)
% Axes labels
ylabel('y (mm)','Fontsize', 10)
xlabel('Maximum expression level ','Fontsize', 10)
%Do you want to save the figure?
if exist('Figure2G.pdf')~=2
print('-Painters', kkq, 'Figure2G','-dpdf','-r300')
movefile('Figure2G.pdf', 'Doneanalysis/WT/figures')
end
%% Figure 2G: y positive vs. max value. 
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1400 3600]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
       % only plot data which is treated as periodic (conditions outlined in  the paper)
        col2=interp1([min(min(traces)) max(max(traces))],[1 siz],max(traces(:,k))); 
        colr= interp1(1:siz,flipud(cmap),col2);
        hold on
        plot( xxs_adjusted(end,k), yys_adjusted(end,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
      end
end
%% Figure 2F sizing and saving options
ylim([indJ1(1)-vals/2 2700]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1)-vals/2:500:10500,'YTickLabel',(0:0.5:10.5)) 
% query limits of x:
xlimits=xlim;
set(gca,'Xtick',xlimits(1):500:14500,'XTickLabel',(0:0.5:10.5)) 
% make sure lines and axes  are thick enough when plotting.
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
% fonts
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 10)
xlabel('x (mm) ','Fontsize', 10)

% plot the imaging squares on top 354 microns in x and y 
% indJ2(k) and indJ1(k) indicate the positions of the tissue sections
% (indexed k) labeled by str{k}. 
% also plot the tissue labels as text. 
for k=start_pos:numpos
 hold on
 if k~=5
 rectangle('position', [indJ2(k)-vals/2 indJ1(k)-vals/2 vals vals])
text((indJ2(k)-vals/2+400),(indJ1(k)-vals/2+150),str{k}, 'Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
 else
 rectangle('position', [indJ2(k)-vals/2 indJ1(k)-vals/2 vals vals])
text((indJ2(k)-vals/2),(indJ1(k)-vals/2+400),str{k}, 'Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
end
end
% set ticker linewidth on the figure
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);

% legend plot 
colormap(jet(128));
caxis([(floor(min(min(traces)))) ceil(max(max(traces)))]);
h = colorbar( 'Location', 'north','XTick',[(floor(min(min(traces)))) (ceil(max(max(traces))))], 'TickLabels', {'low', 'high'}, 'FontSize',7, 'Position',[[0.5 0.02 0.2 0.01]]);
if exist('Figure2F.pdf')~=2
print('-Painters', kkq, 'Figure2F','-dpdf','-r300')
movefile('Figure2F.pdf', 'Doneanalysis/WT/figures')
end
%% Figure 3 A:
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1400 3600]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
     
      for tim=120 
         col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(xxs_adjusted(tim,k), yys_adjusted(tim,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
      end
      end
end
ylim([indJ1(1)-vals/2 2700]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1)-vals/2:500:10500,'YTickLabel',(0:0.5:10.5)) 
set(gca,'Xtick',400:500:14500,'XTickLabel',(0:0.5:10.5)) 
 

colormap(jet(128));
caxis([18 35]);
h = colorbar( 'Location', 'north', 'Direction', 'reverse', 'XTick',[18 35],'TickLabels',{'35h' '18h'}, 'FontSize',7, 'Position', [[0.7 0.9 0.2 0.01]]);


for k=start_pos:numpos
 hold on
 if k~=5
 rectangle('position', [indJ2(k)-vals/2 indJ1(k)-vals/2 vals vals])
text((indJ2(k)-vals/2+400),(indJ1(k)-vals/2+150),str{k}, 'Fontsize', 7)
 else
 rectangle('position', [indJ2(k)-vals/2 indJ1(k)-vals/2 vals vals])
text((indJ2(k)-vals/2),(indJ1(k)-vals/2+400),str{k}, 'Fontsize', 7)
 end
end

set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('x (mm)','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure3A.pdf')~=2
print('-Painters', kkq, 'Figure3A','-dpdf','-r300')
movefile('Figure3A.pdf', 'Doneanalysis/WT/figures')
end

%% Figure 3B: y positive vs. max value. 
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1400 3600]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        for tim=120 %120
         col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(pers(k), yys_adjusted(tim,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
        end
      end
end
ylim([indJ1(1)-vals/2 2700]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1)-vals/2:500:10500,'YTickLabel',(0:0.5:10.5)) 
xlim([18 35])
 set(gca,'Xtick',16:4:36,'XTickLabel',(16:4:36)) 
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('Period (h) ','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure3B.pdf')~=2
print('-Painters', kkq, 'Figure3B','-dpdf','-r300')
movefile('Figure3B.pdf', 'Doneanalysis/WT/figures')
end
%% Figure 3C: y positive vs. z periods. 
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1400 3600]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        for tim=120 %120
         col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(zzs_adjusted(tim,k), yys_adjusted(tim,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
        end
      end
end
ylim([indJ1(1)-vals/2 2700]) 
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1)-vals/2:500:10500,'YTickLabel',(0:0.5:10.5)) 
set(gca,'Xtick',0:50:10500,'XTickLabel',(0:0.05:10.5)) 
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('z (mm) ','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure3C.pdf')~=2
print('-Painters', kkq, 'Figure3C','-dpdf','-r300')
movefile('Figure3C.pdf', 'Doneanalysis/WT/figures')
end
%% Figure 3D: y positive vs. z periods (only roots)
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  700 3600]/300)
cc=[];
allrootcells=sum(countkeep(1:3));
for k=1:allrootcells
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        for tim=120 
         col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(zzs_adjusted(tim,k), yys_adjusted(tim,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
        end
      end
end
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
ylim([indJ1(1)-vals/2 1300]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1)-vals/2:500:10500,'YTickLabel',(0:0.5:10.5)) 
set(gca,'Xtick',0:50:10500,'XTickLabel',(0:0.05:10.5)) 
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('z (mm) ','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
if exist('Figure3D.pdf')~=2
print('-Painters', kkq, 'Figure3D','-dpdf','-r300')
movefile('Figure3D.pdf', 'Doneanalysis/WT/figures')
end
close all
%% now sort by position in: (Commented out Mirela Jan2018). 
% [thwo, reloc]=sort(yys(1, :), 'ascend'); 
% xxs=xxs(:,reloc); 
% yys=yys(:,reloc); 
% traces=traces(:,reloc); 
% pers=pers(reloc);
% pergofs=pergofs(reloc);
% perignore=perignore(reloc);
% persmF=persmF(reloc);
% pergofsmF=pergofsmF(reloc);
% perignoremF=perignoremF(reloc);
% persSR=persSR(reloc);
% pergofsSR=pergofsSR(reloc);
% perignoreSR=perignoreSR(reloc);
% times=times(:, reloc); 
%% Figure out positions where cotyledon, upper hypocotyl and lower hypocotyl to root are in the concantenated data. 
StartNum(1)=1;
StartNum(2)=sum(countkeep(1:4))+1;
StartNum(3)=sum(countkeep(1:4))+ sum(countkeep(5))+1;
StartNum(4)=sum(countkeep)+1;
%% Find projection of lower hypocotyl values  onto the diagonal across the 
% this will be used instead of y values when plotting 
% if you don't want to 

  for j=1: size(xxs_adjusted, 2)
      for st=1:size(xxs_adjusted,1)
      A=[xxs_adjusted(st, j) yys_adjusted(st, j)]; 
      B=[1 1]; 
      yysproj(st, j)=(dot(A,B)/norm(B)^2); 
      end
  end 
  
counting=0; % counting will give us an idea of how many values are there of this type.
    for k=sum(countkeep(1:4))+1:sum(countkeep(1:4))+sum(countkeep(5))
        yys(:,k)=yysproj(:,k);
        counting=counting+1;
    end
% 
%% test for the projection
% y_vec=(vals-pos{i}.y);
% 
% rema=[3     6    40    41    53];
% x_vec(:,rema)=[];
% y_vec(:,rema)=[];
% for j=1:54
%             for st=1:size(xxs_adjusted,1)
%       vec=[x_vec(st, j) y_vec(st, j)]; 
%       B=[1 1]; 
%       yysproja(st, j)=(dot(vec,B)/norm(B)^2); 
%             end
% end 
%   A1=(yys(:,sum(countkeep(1:4))+1:sum(countkeep(1:4))+sum(countkeep(5)))-min(min(yys(:,sum(countkeep(1:4))+1:sum(countkeep(1:4))+sum(countkeep(5))))));
% B1=(yysproja-min(min(yysproja)));
% DD=A1-B1;
% DD=abs(DD);
% max(max(DD))
% %
%%
mina=3000; maxa=0;
for j=  1:size(traces,2)  
      mpkdis2=15; 
      MinPeakDist=mpkdis2*(times(2,j)-times(1,j));
      MinPeakDist2=mpkdis2*(times(2,j)-times(1,j));
             MinPkHeight=-200; 
             MinPkProm=20;
            MinPkPromTr=1;
            % Find the peaks:
            clear dtdata p trend
            p = polyfit(times(:, j),traces(:,j),1);  % find trend as polynomial of degree 'd'
            trend = polyval(p,times(:, j));
            dtdata = traces(:,j) - trend;  
            dtdata= smooth(times(:, j),dtdata,0.05,'rloess');

     if  pergofs(j)<gof_cutoff && pergofsmF(j)<gof_cutoff && pergofsSR(j)<gof_cutoff2  && abs(pers(j)- persSR(j))<2.5 && abs(persmF(j)- persSR(j))<2.5 && abs(persmF(j)-pers(j))<2.5 && sum(perignoremF(j)+perignoreSR(j)+perignore(j))==0%<0.9  & pos{i}.periods(j,5)<=35 & pos{i}.periods(j,5)>=18 & pos{i}.periods(j,9)>0 & abs(pos{i}.periods(j,5)- pos{i}.periods(j,1))<=5% pos{i}.periods(j,2)<0.5
      [a,b]= findpeaks(dtdata, times(:, j),  'MinPeakDistance',MinPeakDist, 'MinPeakHeight', -50, 'MinPeakProminence',MinPkProm); 
        peakscell{j}=b; 
        
        
        for k=1:length(b)
            
            c(k)=find(times(:, j)==b(k));
        end
        
        peakscellval{j}=a; 
        peakscellyy{j}=yys(c,j); 
        peakscelllen{j}=length(b); 
        
        mina=min([mina a']);
        maxa=max([maxa a']);
        clear a b
     else
        peakscell{j}=nan; 
           peakscellyy{j}=nan; 
        peakscelllen{j}=nan; 
     end
         
end 

 for i=1:size(peakscell,2)
     if peakscell{i}(1)>47.5
          peakscellA{i}=[ 0; peakscell{i}];
        peakscellyyA{i}=[ 0 ;peakscellyy{i}];
     else
        peakscellA{i}=peakscell{i};
        peakscellyyA{i}=peakscellyy{i};
     end
         
 end
%% Colour maps 
i=1; 
cmap = colormap(hot); % colour map to use.
siz=size(cmap,1);
colo=1/100*[5.1, 23.9, 33.7; 94.9, 42.7, 6.7;7.8, 58.8, 73.3 ; 	100, 83.5, 0 ;	 80 20 100;	10.2, 67.8, 16.9;87.8, 27.8, 20 ; 10.2, 67.8, 16.9; 45.5, 10.2, 67.8];
%% Figures 3F, 3G  and 3H
for figs=1:3
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  2500 2000]/300)
for j=StartNum(figs):StartNum(figs+1)-1
hold on; 
    if j==StartNum(figs)
        % for every new figure plot the background subjective night and day
        y = [0 0 3000 3000]; 
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
    end
    if isnan(peakscellA{j})==0 
       s=[mina maxa]; 
      for i=1:length( peakscellA{j})
        if peakscellyyA{j}(i)>1000
        plot( peakscellA{j}(i),peakscellyyA{j}(i),'s', 'MarkerFaceColor', colo(i, :), 'MarkerEdgeColor', colo(i, :), 'MarkerSize',2);
        else
        plot( peakscellA{j}(i),peakscellyyA{j}(i),'s', 'MarkerFaceColor', colo(i, :), 'MarkerEdgeColor', colo(i, :), 'MarkerSize',2);
        end
      end
    end;
% Figure details : fonts, axes labels and rest. 
set(gca, 'Fontsize', 7)
%x axis
ax = gca;
ax.XTick = [0:24:144];
xlabel('Time (h)');
xlim([min(min(times)) max(max(times))])
% y axis and labels:
allyys=yys(:, StartNum(figs):StartNum(figs+1)-1);
ylim([floor(min(min(allyys)))  max(max(allyys))])          
set(gca,'Ytick',floor(min(min(allyys))):200:10500)
 % y label for middle figure is y rescaled along the diagonal. 
    if figs==2
    ylabel('rescaled y (mm)');
    else
    ylabel('y (mm)');
    end
clear allyys
end
ax = gca;
set(gca, 'FontSize', 7)
ax.XTick = [0:24:144];

if figs==1
   xx=floor(min(min(yys(:, StartNum(figs):StartNum(figs+1)-1)))):200:10500;% 
   set(gca,'YTickLabel', cellstr(num2str(xx(:)/1000, '%.2f')))% 
   set(gca, 'LineWidth', 1.2)
   set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
   if exist('Figure3H.pdf')~=2
print('-Painters', kkq, 'Figure3H','-dpdf','-r300')
movefile('Figure3H.pdf', 'Doneanalysis/WT/figures')
end
elseif figs==2
    xx=0:200:10500;
   set(gca,'YTickLabel', cellstr(num2str(xx(:)/1000, '%.2f')))
   set(gca, 'LineWidth', 1.2)
    set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
    set(gca, 'LineWidth', 1.2)
    set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
    if exist('Figure3G.pdf')~=2
    print('-Painters', kkq, 'Figure3G','-dpdf','-r300')
    movefile('Figure3G.pdf', 'Doneanalysis/WT/figures')
    end
elseif figs==3
    xx=floor(min(min(yys(:, StartNum(figs):StartNum(figs+1)-1)))):200:10500;
    set(gca,'YTickLabel', cellstr(num2str(xx(:)/1000, '%.2f')))
    set(gca, 'LineWidth', 1.2)
    set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
    if exist('Figure3F.pdf')~=2
print('-Painters', kkq, 'Figure3F','-dpdf','-r300')
movefile('Figure3F.pdf', 'Doneanalysis/WT/figures')
end
end   
end
close all
%%  colour scheme for Fig 3E. 
tracesmin=traces- repmat(min(traces),size(traces,1),1);
tracesmin=tracesmin*diag(1./(max(tracesmin)));

% what do we want to calculate: normalized or raw? 
traces2=traces; 
cmap = colormap(hot); % colour map
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
%% Figure 3E bottom panel space-time plots
addpath([ pwd '/subplot_tight/'])
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  6000 3000]/300)
func_hndl={@subplot;@subplot_tight};
count=0;
for  ti= 16:1:34
       set(gcf, 'Color', 'black')
   func_hndl{2}(1,18,ti-15, 0.005); 
   for i=1: 440 
       set(gcf, 'Color', 'black')
    if  pergofs(i)<gof_cutoff && pergofsmF(i)<gof_cutoff && pergofsSR(i)<gof_cutoff2  && abs(pers(i)- persSR(i))<2.5 && abs(persmF(i)- persSR(i))<2.5 && abs(persmF(i)-pers(i))<2.5 && sum(perignoremF(i)+perignoreSR(i)+perignore(i))==0%<0.9  & pos{i}.periods(j,5)<=35 & pos{i}.periods(j,5)>=18 & pos{i}.periods(j,9)>0 & abs(pos{i}.periods(j,5)- pos{i}.periods(j,1))<=5% pos{i}.periods(j,2)<0.5
     plot(xxs_adjusted(ti,i)',yys_adjusted(ti,i)', '.','Markersize',full6x(ti,i)/20+2,'Color',(dd6x(ti,i,1:3))) %1.2
    if ti==16
       line( [min(min(xxs_adjusted(ti,1:440)))+30 min(min(xxs_adjusted(ti,1:440)))+30], [min(min(yys_adjusted(ti,1:440)))+50 min(min(yys_adjusted(ti,1:440)))+50+250],  'LineStyle','-','Linewidth', 2, 'Color', 'white')
    end
    if ti==16
    count=count+1;
    end
        end
    hold on;  
    grid off;
    set(gca, 'color', 'k')
   end
x = [0 1 1 0];
y = [0 0 1 1];
patch(x,y,'red')

daspect([1.5 1.5 1.5]); camlight;lighting gouraud;
    
set(gca, 'Fontsize', 20)
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([min(min(xxs_adjusted(:,1: 440))) max(max(xxs_adjusted(:,1: 440)))])
ylim([min(min(yys_adjusted(:,1: 440) )) max(max(yys_adjusted(:,1: 440)))])
hold on
  
end
set(gcf,'color','w');
set(gcf, 'InvertHardCopy', 'off');
if exist('Figure3E_lastpanel.pdf')~=2
print('-Painters',kkq, 'Figure3E_lastpanel','-dpdf','-r300')
movefile('Figure3E_lastpanel.pdf', 'Doneanalysis/WT/figures')
end
close all
%% Figure 3E second bottom panel space-time plots
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  6000 3000]/300)
func_hndl={@subplot;@subplot_tight};
count=0;
for  ti= 16:1:34
    set(gcf, 'Color', 'black')
    func_hndl{2}(1,18,ti-15, 0.005);
   for i=441:558  %
       set(gcf, 'Color', 'black')
   if  pergofs(i)<gof_cutoff && pergofsmF(i)<gof_cutoff && pergofsSR(i)<gof_cutoff2  && abs(pers(i)- persSR(i))<2.5 && abs(persmF(i)- persSR(i))<2.5 && abs(persmF(i)-pers(i))<2.5 && sum(perignoremF(i)+perignoreSR(i)+perignore(i))==0%<0.9  & pos{i}.periods(j,5)<=35 & pos{i}.periods(j,5)>=18 & pos{i}.periods(j,9)>0 & abs(pos{i}.periods(j,5)- pos{i}.periods(j,1))<=5% pos{i}.periods(j,2)<0.5
     plot(xxs_adjusted(ti,i)',yys_adjusted(ti,i)', '.','Markersize',full6x(ti,i)/20+2,'Color',(dd6x(ti,i,1:3))) %1.2
    if ti==16
       line( [min(min(xxs_adjusted(ti,441:558)))-23 min(min(xxs_adjusted(ti,441:558)))-23], [min(min(yys_adjusted(ti,441:558)))+10 min(min(yys_adjusted(ti,441:558)))+10+100],  'LineStyle','-','Linewidth', 1, 'Color', 'white')
    end
    if ti==16
    count=count+1;
    end
   end
hold on;  
grid off;
set(gca, 'color', 'k')
end
    
   x = [0 1 1 0];
y = [0 0 1 1];
patch(x,y,'red')

  
    daspect([1.5 1.5 1.5])
    camlight;lighting gouraud;
    
    set(gca, 'Fontsize', 20)

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ytick',[])
     xlim([min(min(xxs_adjusted(:,441:558)))-25 max(max(xxs_adjusted(:,441:558)))])
     ylim([min(min(yys_adjusted(:,441:558) )) max(max(yys_adjusted(:,441:558)))])
hold on


   
end
set(gcf,'color','w');
set(gcf, 'InvertHardCopy', 'off');
if exist('Figure3E_midpanel.pdf')~=2
print('-Painters',kkq, 'Figure3E_midpanel','-dpdf','-r300')
 movefile('Figure3E_midpanel.pdf', 'Doneanalysis/WT/figures')
end
%% Figure 3E next to top panel space-time plots
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  6000 3000]/300)
func_hndl={@subplot;@subplot_tight};
count=0;
for  ti= 16:1:34%
 set(gcf, 'Color', 'black')
 func_hndl{2}(1,18,ti-15, 0.005); %-60);
   for i=559:612 
    set(gcf, 'Color', 'black')
    if  pergofs(i)<gof_cutoff && pergofsmF(i)<gof_cutoff && pergofsSR(i)<gof_cutoff2  && abs(pers(i)- persSR(i))<2.5 && abs(persmF(i)- persSR(i))<2.5 && abs(persmF(i)-pers(i))<2.5 && sum(perignoremF(i)+perignoreSR(i)+perignore(i))==0%<0.9  & pos{i}.periods(j,5)<=35 & pos{i}.periods(j,5)>=18 & pos{i}.periods(j,9)>0 & abs(pos{i}.periods(j,5)- pos{i}.periods(j,1))<=5% pos{i}.periods(j,2)<0.5
    plot(xxs_adjusted(ti,i)',yys_adjusted(ti,i)', '.','Markersize',full6x(ti,i)/20+2,'Color',(dd6x(ti,i,1:3))) %1.2
    if ti==16
       line( [min(min(xxs_adjusted(ti,559:612)))-10 min(min(xxs_adjusted(ti,559:612)))-10], [min(min(yys_adjusted(ti,559:612)))+10 min(min(yys_adjusted(ti,559:612)))+10+100],  'LineStyle','-','Linewidth', 1, 'Color', 'white')
    end
    if ti==16
    count=count+1;
    end
 end
 
hold on;  
grid off;
set(gca, 'color', 'k')
end
    
x = [0 1 1 0];
y = [0 0 1 1];
patch(x,y,'red')

daspect([1.5 1.5 1.5]); camlight;lighting gouraud;
    
set(gca, 'Fontsize', 20)
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ytick',[])
xlim([min(min(xxs_adjusted(:,559:612) ))-25 max(max(xxs_adjusted(:,559:612) ))])
ylim([min(min(yys_adjusted(:,559:612) )) max(max(yys_adjusted(:,559:612)))])
hold on 
end
set(gcf,'color','w');
set(gcf, 'InvertHardCopy', 'off');
if exist('Figure3E_midpanel2.pdf')~=2
print('-Painters',kkq, 'Figure3E_midpanel2','-dpdf','-r300')
movefile('Figure3E_midpanel2.pdf', 'Doneanalysis/WT/figures')
end
%% Figure 3E top panel space-time plots
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  6000 3000]/300)
func_hndl={@subplot;@subplot_tight};
count=0;
for  ti= 16:1:34%
  set(gcf, 'Color', 'black')
 func_hndl{2}(1,18,ti-15, 0.005); 
 for i=613:722  
  set(gcf, 'Color', 'black')
  if  pergofs(i)<gof_cutoff && pergofsmF(i)<gof_cutoff && pergofsSR(i)<gof_cutoff2  && abs(pers(i)- persSR(i))<2.5 && abs(persmF(i)- persSR(i))<2.5 && abs(persmF(i)-pers(i))<2.5 && sum(perignoremF(i)+perignoreSR(i)+perignore(i))==0%<0.9  & pos{i}.periods(j,5)<=35 & pos{i}.periods(j,5)>=18 & pos{i}.periods(j,9)>0 & abs(pos{i}.periods(j,5)- pos{i}.periods(j,1))<=5% pos{i}.periods(j,2)<0.5
    plot(xxs_adjusted(ti,i)',yys_adjusted(ti,i)', '.','Markersize',full6x(ti,i)/20+2,'Color',(dd6x(ti,i,1:3))) %1.2
    if ti==16
    line( [min(min(xxs_adjusted(ti,613:722)))-15 min(min(xxs_adjusted(ti,613:722)))-15], [min(min(yys_adjusted(:,613:722)))+10 min(min(yys_adjusted(:,613:722)))+10+100],  'LineStyle','-','Linewidth', 1, 'Color', 'white')
    end
    if ti==16
    count=count+1;
  end
  end
hold on;  
grid off;
set(gca, 'color', 'k')
end
x = [0 1 1 0];
y = [0 0 1 1];
patch(x,y,'red')
daspect([1.5 1.5 1.5]); camlight;lighting gouraud;

set(gca, 'Fontsize', 20);
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ytick',[])
xlim([min(min(xxs_adjusted(:,613:722)))-25 max(max(xxs_adjusted(:,613:722)))])
ylim([min(min(yys_adjusted(:,613:722) )) max(max(yys_adjusted(:,613:722)))])
hold on
end
set(gcf,'color','w');
set(gcf, 'InvertHardCopy', 'off');
if exist('Figure3E_toppanel.pdf')~=2
print('-Painters',kkq, 'Figure3E_toppanel','-dpdf','-r300')
movefile('Figure3E_toppanel.pdf', 'Doneanalysis/WT/figures')
end
close all