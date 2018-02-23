%  File for analysis of WT data (WTrepeat) for Gould, Domijan et al.
%  "Coordination of robust single cell rhythms in the Arabidopsis circadian 
%  clock via spatial waves of gene expression" (BioRxiv)
%  doi https://doi.org/10.1101/208900s
% 
%  by Mirela Domijan (U. of Liverpool)
clear all 
%%  
name='WTrepeat'; %Experiment to run. 
%% add in relevant paths to the data
addpath([pwd '/Data_singlecell/'])
addpath([pwd '/Data_singlecell/' name '_final_coordinates'])
addpath('/Users/mirela/Documents/MATLAB/Archive/CCA1starline_Sept2015/TESTperiod_results/')
addpath([pwd '/subplot_tight'])
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
    posi=[1:6];
    sect={ {'root tip', ''},{'root', '(sect 1)'}, {'root', '(sect 2)'},{'hypocotyl', '(sect 1)'},{'hypocotyl', '(sect 2)'},{'cotyledon', ''}};
elseif strcmp(name,'WTrepeat')==1
    posi= [13:20];
    sect={{'root tip', ''},{'root', '(sect 1)'}, {'root', '(sect 2)'},{'root', '(sect 3)'}, {'root', '(sect 4)'} ,{'hypocotyl', '(sect 1)'},{'hypocotyl', '(sect 2)'},{'cotyledon', ''}}; 
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

    if numel(unique(diff([0; fa])))==1 & (size(LL,1)-fa(end))== (fa(1)-1) 
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
if strcmp(name, 'Exp26')==1
[Times,txt2]=xlsread(['timestamp_WT.xlsx'], 'C8:DR15');
pos{i}.time= Times(i,:);
elseif strcmp(name, 'WTrepeat')==1
[Times,txt2]=xlsread(['timestamp_WTrepeat.xlsx'], 'B8:AG27');
pos{i}.time= Times(i,2:end);
end
clear data FFTrecs periods
end
clear i k  txt2 posfile bb a Times C name2 name3 

% Goodness of fit cut-off criteria.
gof_cutoff=0.9;
gof_cutoff2=1;
%% Add in the pad positions 
addpath('pictures');
E26=imread('WTrepeat.png');
imshow(E26)
E26(:, 1:1800, :)=[];
E26(:, 1000:end, :)=[];
E26(1:300, :, :)=[];
E26(1980:end, :, :)=[];
imshow(E26);
bw=im2bw(E26, 0.1);
imshow(bw);
bw2= imfill(bw,'holes');
imshow(bw2);
hold on
[B,L] = bwboundaries(bw2);
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
   hold on
end
st = regionprops(bw2, 'BoundingBox' );
Col=colormap('jet');
for k=1:size(st,1)
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
[ind1,ind2]=find(imregionalmax(Ba)>0);
plot(ind2, ind1, 'o')
set(gca,'YDir','Reverse')
%resizing
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
bbJ=double(J); 
sqsJ=ones(vals,vals); 
BaJ=conv2(bbJ,sqsJ,'same');
CJ = imfuse(J,imregionalmax(BaJ),'blend','Scaling','joint');
imshow(CJ); 
[indJ1,indJ2]=find(imregionalmax(BaJ)>0);
plot(indJ2, indJ1, 'o')
set(gca,'YDir','Reverse')  
ULC=[indJ2-vals/2 indJ1-vals/2];
[sa, sb]=(sort(indJ1, 'ascend'));
[ a, xx]=(unique(indJ1));
indJ1u=indJ1(xx);
indJ2u=indJ2(xx);
indJ1=indJ1u(2:9);
indJ2=indJ2u(2:9);
%% quick plots of values
for k=start_pos:numpos
    i=posi(k);
    pos{i}.xmod=indJ2(k)+ (pos{i}.x);
    pos{i}.ymod=indJ1(k)+ (vals-pos{i}.y);   
    plot3(pos{i}.ymod(1, :), pos{i}.xmod(1, :),pos{i}.z(1, :), '.','Markersize',10)
end

for k=start_pos:  numpos
    i=posi(k);
    pos{i}.xmod=indJ2(k)+ (pos{i}.x);
    pos{i}.ymod=indJ1(k)+ (vals-pos{i}.y);
     if i==2
    pos{i}.xmod=pos{i}.xmod+3;
    pos{i}.ymod=pos{i}.ymod-1;  
     end 
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
%% which tissues are we looking at:
str={{'Root tip', ''},{'Root (sect 1)'}, {'Root (sect 2)'},{'Root (sect 3)'}, {'Root (sect 4)'} ,{'Hypocotyl', '(sect 1)'},{'Hypocotyl', '(sect 2)'},{'Cotyledon', ''}}; 
%% if mixing positions together do this: 
for k= 1 : numpos 
     i=posi(k);
 if k==1
    xxs=[ xxs pos{i}.xmod(1:30, :)];
    yys=[ yys pos{i}.ymod(1:30, :)];
    zzs=[ zzs pos{i}.z(1:30, :)];
    traces=[ traces pos{i}.YFP(1:30, :)];
    times=[times repmat(pos{i}.time(1:30)',1,(size(pos{i}.YFP,2)))];
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
    [ik, ja]=find(pos{i}.ymod<724);   
    dontkeep=unique(ja);
    keep=(1:size(pos{i}.ymod,2));
    keep(dontkeep)=[];
    clear ik ja bb  
    xxs=[ xxs pos{i}.xmod(1:30, keep)];
    yys=[ yys pos{i}.ymod(1:30, keep)]; 
    zzs=[ zzs pos{i}.z(1:30, keep)];
    traces=[ traces pos{i}.YFP(1:30, keep)];
    times=[times repmat(pos{i}.time(1:30)',1,(size(pos{i}.YFP(1:30, keep),2)))];
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
    clear keep
    [ik, ja]=find(pos{i}.ymod<994);   
    dontkeep=unique(ja);
    keep=(1:size(pos{i}.ymod,2));
    keep(dontkeep)=[];
    clear ik ja bb  
    xxs=[ xxs pos{i}.xmod(1:30, keep)];
    yys=[ yys pos{i}.ymod(1:30, keep)];  
    zzs=[ zzs pos{i}.z(1:30, keep)];
    traces=[ traces pos{i}.YFP(1:30, keep)];
    times=[times repmat(pos{i}.time(1:30)',1,(size(pos{i}.YFP(1:30, keep),2)))];
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
    clear keep
    [ik, ja]=find(pos{i}.ymod<1303);   
    dontkeep=unique(ja);
    keep=(1:size(pos{i}.ymod,2));
    keep(dontkeep)=[];
    clear ik ja bb 
    xxs=[ xxs pos{i}.xmod(1:30, keep)];
    yys=[ yys pos{i}.ymod(1:30, keep)];  
    zzs=[ zzs pos{i}.z(1:30, keep)];
    traces=[ traces pos{i}.YFP(1:30, keep)];
    times=[times repmat(pos{i}.time(1:30)',1,(size(pos{i}.YFP(1:30, keep),2)))];
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
    clear keep
    [ik, ja]=find(pos{i}.ymod<1529);   
    dontkeep=unique(ja);
    keep=(1:size(pos{i}.ymod,2));
    keep(dontkeep)=[];
    clear ik ja bb
    xxs=[ xxs pos{i}.xmod(1:30, keep)];
    yys=[ yys pos{i}.ymod(1:30, keep)];  
    zzs=[ zzs pos{i}.z(1:30, keep)];
    traces=[ traces pos{i}.YFP(1:30, keep)];
    times=[times repmat(pos{i}.time(1:30)',1,(size(pos{i}.YFP(1:30, keep),2)))];
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
    xxs=[ xxs pos{i}.xmod(1:30, :)];
    yys=[ yys pos{i}.ymod(1:30, :)];
    zzs=[ zzs pos{i}.z(1:30, :)]; 
    traces=[ traces pos{i}.YFP(1:30, :)];
    times=[times repmat(pos{i}.time(1:30)',1,(size(pos{i}.YFP,2)))];
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
 end
end
xxs_adjusted=xxs;yys_adjusted=yys;zzs_adjusted=zzs*2/0.692;

%% Making of the Figures:
% below you can plot: 
% Fig2 Supp1 D
% Fig2 Supp1 C
% Fig2 Supp1 E
% Fig2 Supp1 F
% Fig2 Supp1 G
% Fig2 Supp3 
%  
%% Figure SI2 D
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  2000 3600]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        col2=interp1([min(min(traces)) max(max(traces))],[1 siz],max(traces(:,k)));
        colr= interp1(1:siz,flipud(cmap),col2);
        hold on
        plot( max(traces(:,k)), yys_adjusted(end,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
      end
end
ylim([indJ1(1) 4200]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5))  
xlim([min(min(traces)) max(max(traces))])

set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('Maximum expression level ','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
print('-Painters', kkq, 'Fig2Supp1D','-dpdf','-r300')
movefile('Fig2Supp1D.pdf', 'Doneanalysis/WTrepeat/figures')
%% Figure SI 2C: y positive vs. max value. 
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  2000 3600]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        col2=interp1([min(min(traces)) max(max(traces))],[1 siz],max(traces(:,k)));
        colr= interp1(1:siz,flipud(cmap),col2);
        hold on
        plot( xxs_adjusted(end,k), yys_adjusted(end,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
      end
end
ylim([indJ1(1) 4200]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5)) 
set(gca,'Xtick',0:500:6500,'XTickLabel',0:0.5:6.5) 
xlim([700 max(max(xxs_adjusted))+100])

set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlab=xlabel('x (mm) ','Fontsize', 7);
set(xlab, 'position', get(xlab, 'position')+[0,20,0])

set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
for k=start_pos:numpos
hold on
rectangle('position', [indJ2(k) indJ1(k) vals vals], 'Linewidth', 0.15)
text((indJ2(k)+400),(indJ1(k)+150),str{k}, 'Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
end

colormap(jet(128));
caxis([(floor(min(min(traces)))) (ceil(max(max(traces))))]);
h = colorbar( 'Location', 'north','XTick',[(floor(min(min(traces)))) (ceil(max(max(traces))))], 'TickLabels', {'low', 'high'}, 'FontSize',7, 'Position',[[0.4 0.032 0.2 0.01]]);
daspect([1 1 1]);

print('-Painters', kkq, 'Fig2Supp1C','-dpdf','-r300')
movefile('Fig2Supp1C.pdf', 'Doneanalysis/WTrepeat/figures')
%% Figure SI 2E:
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  2000 3600]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(xxs_adjusted(end,k), yys_adjusted(end,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
      end
end
ylim([indJ1(1) 4200]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5)) 
set(gca,'Xtick',0:500:6500,'XTickLabel',0:0.5:6.5) 
 xlim([700 max(max(xxs_adjusted))+100])
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlab=xlabel('x (mm) ','Fontsize', 7);
set(xlab, 'position', get(xlab, 'position')+[0,20,0])
colormap(jet(128));
caxis([18 35]);
for k=start_pos:numpos
hold on
rectangle('position', [indJ2(k) indJ1(k) vals vals], 'Linewidth', 0.15)
text((indJ2(k)+400),(indJ1(k)+150),str{k}, 'Fontsize', 7)
end
colormap(jet(128));
caxis([18 35]);
h = colorbar( 'Location', 'north','XTick',[18 35], 'TickLabels', {'18h', '35h'}, 'FontSize',7, 'Position',[[0.4 0.032 0.2 0.01]]);daspect([1 1 1])
print('-Painters', kkq, 'Fig2Supp1E','-dpdf','-r300')
movefile('Fig2Supp1E.pdf', 'Doneanalysis/WTrepeat/figures')

%% Figure SI 2F: y positive vs. max value. 
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  2000 3600]/300)
cc=[];
for k=1:size(pers,1)
      if   pergofs(k)<gof_cutoff && pergofsmF(k)<gof_cutoff && pergofsSR(k)<gof_cutoff2  && abs(pers(k)- persSR(k))<2.5 && abs(persmF(k)- persSR(k))<2.5 && abs(persmF(k)-pers(k))<2.5 && sum(perignoremF(k)+perignoreSR(k)+perignore(k))==0 && pers(k)<=35 && pers(k)>=18
        col2=interp1([18 35],[1 siz],pers(k));
        colr= interp1(1:siz,cmap,col2);
        hold on
        plot(pers(k), yys_adjusted(end,k),'.','Markersize',10, 'Color', colr)
        clear col2 colr
      end
end
ylim([indJ1(1) 4200])
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5)) 
colormap(jet(128));
caxis([18 35]);
xlim([18 35])
set(gca,'Xtick',16:4:36,'XTickLabel',(16:4:36)) 

set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('Period (h) ','Fontsize', 7)
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
print('-Painters', kkq, 'Fig2Supp1F','-dpdf','-r300')
movefile('Fig2Supp1F.pdf', 'Doneanalysis/WTrepeat/figures')
%% Figure SI 2H: y positive vs. z periods. 
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  2000 3600]/300)
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
ylim([indJ1(1) 4200]) 
set(gca,'Ytick',1000:500:10500,'YTickLabel',(1:0.5:10.5)) 
set(gca,'Ytick',indJ1(1):500:10500,'YTickLabel',(0:0.5:10.5)) 
set(gca,'Xtick',0:50:10500,'XTickLabel',(0:0.05:10.5)) 
set(gca,'Fontsize', 7)
ylabel('y (mm)','Fontsize', 7)
xlabel('z (mm) ','Fontsize', 7)
colormap(jet(128));
caxis([18 35]);
print('-Painters', kkq, 'Fig2Supp3A','-dpdf','-r300')
movefile('Fig2Supp3A.pdf', 'Doneanalysis/WTrepeat/figures')
%% Figure SI 2J: y positive vs. z periods (only roots)
clf
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1700 2400]/300)
cc=[];
allrootcells=sum(countkeep(1:5));
for k=1:allrootcells
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
colormap(jet(128));
caxis([18 35]);
set(gca, 'LineWidth', 1.2)
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
print('-Painters', kkq, 'Fig2Supp3B','-dpdf','-r300')
movefile('Fig2Supp3B.pdf', 'Doneanalysis/WTrepeat/figures')
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
%% In case you want to plot phases need this info:
mina=3000; maxa=0;
for j=  1:size(traces,2)  
     mpkdis2=8; 
     MinPeakDist=mpkdis2*(times(2,j)-times(1,j));
     MinPeakDist2=mpkdis2*(times(2,j)-times(1,j));
     MinPkHeight=100; 
     MinPkProm=20; 
     if  pergofs(j)<gof_cutoff && pergofsmF(j)<gof_cutoff && pergofsSR(j)<gof_cutoff2  && abs(pers(j)- persSR(j))<2.5 && abs(persmF(j)- persSR(j))<2.5 && abs(persmF(j)-pers(j))<2.5 && sum(perignoremF(j)+perignoreSR(j)+perignore(j))==0%<0.9  & pos{i}.periods(j,5)<=35 & pos{i}.periods(j,5)>=18 & pos{i}.periods(j,9)>0 & abs(pos{i}.periods(j,5)- pos{i}.periods(j,1))<=5% pos{i}.periods(j,2)<0.5
        [a,b]= findpeaks(traces(:,j),times(:, j),  'MinPeakDistance',MinPeakDist, 'MinPeakHeight', MinPkHeight, 'MinPeakProminence',MinPkProm); 
        [d,c]= findpeaks(traces(:,j),  'MinPeakDistance',mpkdis2, 'MinPeakHeight', MinPkHeight, 'MinPeakProminence',MinPkProm); 
        peakscell{j}=b; 
        tracesrhythm(:,j)=traces(:,j);
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
        peakscelllen{j}=nan; 
     end      
 end 
%% Info needed to make filmstrip (below). Some may be redundant. 
[a,b]=find(tracesrhythm==0);tracesrhythm(:,b)=[];
tracesmin=traces- repmat(min(traces),size(traces,1),1);
tracesmin=tracesmin*diag(1./(max(tracesmin)));
tracesminrhythm=tracesrhythm- repmat(min(tracesrhythm),size(tracesrhythm,1),1);
tracesminrhythm=tracesminrhythm*diag(1./(max(tracesminrhythm)));
tracesrhythma=tracesrhythm(2:end,:); 
tracesminrhythma=tracesrhythma- repmat(min(tracesrhythma),size(tracesrhythma,1),1);
tracesminrhythma=tracesminrhythma*diag(1./(max(tracesminrhythma)));
% what do we want to calculate: normalized or raw? 
traces2=traces; 
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
close all
%% Figure2 SI2 G
func_hndl={@subplot;@subplot_tight};
kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  4800 700]/300)
for  ti=1:1:15
set(gcf, 'Color', 'black')
whitebg(1, 'k')
set(gcf, 'Color', 'black')
func_hndl{2}(1,15,ti, 0.001); 
   for i=1: sum(countkeep(1:5))
           set(gcf, 'Color', 'k')
           if  pergofs(i)<gof_cutoff && pergofsmF(i)<gof_cutoff && pergofsSR(i)<gof_cutoff2  && abs(pers(i)- persSR(i))<2.5 && abs(persmF(i)- persSR(i))<2.5 && abs(persmF(i)-pers(i))<2.5 && sum(perignoremF(i)+perignoreSR(i)+perignore(i))==0%
            plot(xxs_adjusted(ti,i)',yys_adjusted(ti,i)', '.','Markersize',full6x(ti,i)/25+1,'Color',(dd6x(ti,i,1:3)))
           end
    hold on;  
    grid off;
    set(gca, 'color', 'k')
   end
    daspect([1.5 1.5 1.5])
    camlight;lighting gouraud;   
    set(gca, 'Fontsize', 12)   
    if ti==1
       line( [max(max(xxs_adjusted(:,1:sum(countkeep(1:5)))))-50 max(max(xxs_adjusted(:,1:sum(countkeep(1:5)))))-50], [min(min(yys_adjusted(:,1:sum(countkeep(1:5)))))+50 min(min(yys_adjusted(:,1:sum(countkeep(1:5)))))+50+250],  'LineStyle','-','Linewidth', 2, 'Color', 'white')
    end
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([min(min(xxs_adjusted(:,1:sum(countkeep(1:5))))) max(max(xxs_adjusted(:,1:sum(countkeep(1:5)))))])
ylim([min(min(yys_adjusted(:,1:sum(countkeep(1:5))))) max(max(yys_adjusted(:,1:sum(countkeep(1:5)))))])
hold on
set(gcf, 'Color', 'k') 
end
set(gcf,'color','w');
set(gcf, 'InvertHardCopy', 'off');
print('-Painters', kkq, 'Fig2Supp1G_upperpanel','-dpdf','-r300')
movefile('Fig2Supp1G_upperpanel.pdf', 'Doneanalysis/WTrepeat/figures')
close all

kkq=figure;
set(kkq,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  4800 700]/300)
for  ti=16:1:30%
       set(gcf, 'Color', 'black')
 whitebg(1, 'k')
set(gcf, 'Color', 'black')
 func_hndl{2}(1,15,ti-15, 0.001); 
   for i=1: sum(countkeep(1:5))
           set(gcf, 'Color', 'k')
           if  pergofs(i)<gof_cutoff && pergofsmF(i)<gof_cutoff && pergofsSR(i)<gof_cutoff2  && abs(pers(i)- persSR(i))<2.5 && abs(persmF(i)- persSR(i))<2.5 && abs(persmF(i)-pers(i))<2.5 && sum(perignoremF(i)+perignoreSR(i)+perignore(i))==0%
           plot(xxs_adjusted(ti,i)',yys_adjusted(ti,i)', '.','Markersize',full6x(ti,i)/25+1,'Color',(dd6x(ti,i,1:3)))
           end
    hold on;  
    grid off;
    set(gca, 'color', 'k')
   end
    daspect([1.5 1.5 1.5])
    camlight;lighting gouraud;   
    set(gca, 'Fontsize', 12)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlim([min(min(xxs_adjusted(:,1:sum(countkeep(1:5))))) max(max(xxs_adjusted(:,1:sum(countkeep(1:5)))))])
    ylim([min(min(yys_adjusted(:,1:sum(countkeep(1:5))))) max(max(yys_adjusted(:,1:sum(countkeep(1:5)))))])
    hold on
    set(gcf, 'Color', 'k') 
end
set(gcf,'color','w');
set(gcf, 'InvertHardCopy', 'off');
print('-Painters', kkq, 'Fig2Supp1G_lowerpanel','-dpdf','-r300')
movefile('Fig2Supp1G_lowerpanel.pdf', 'Doneanalysis/WTrepeat/figures')
close all
