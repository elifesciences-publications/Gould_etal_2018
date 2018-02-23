%% Model from  "Coordination of robust single cell rhythms in the Arabidopsis circadian clock via spatial waves of gene expression" 
%  This is a script for the Kuramoto phase oscillator model described in Gould, Domijan et al.
%  "Coordination of robust single cell rhythms in the Arabidopsis circadian 
%  clock via spatial waves of gene expression" (BioRxiv)
%  doi https://doi.org/10.1101/208900s
% 
%  by Mirela Domijan (U. of Liverpool)
%%
clear all 
%% initialise  
% T =end time (144h)
% dt= time step (dt): adjusted to 1/24.  In 24h will go through 2*pi angle.  
% kappa = couping constant. 
T= 144; 
dt=0.0417; 
kappa= 1;  
%% initalise model plant template: 
% template from ModelTemplate.mat. Matrix of 0's and 1's
% plot it and record the number of cells (n)
% Nx= number of pixels in x dimension
% Ny= number of pixels in y dimension
load('ModelTemplate')  
Fig1=figure;
set(Fig1,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  1200 3000]/300)
spy(template);
xlabel('')
title({'Arabidopsis Template',['(n=' num2str(nnz(template)) ')']})
Nx=size(template,1);
Ny=size(template,2);
print('-Painters', Fig1, 'PlantTemplate','-dpdf','-r300')
%% initalise data arrays:
% X= data array. Initial phase (X(:, :, 1) is zero (all cells are synched
% by LD cycles)
% omega = intrinsic period (2*pi) with differing values it root tip & shoot. 
% S = coupling information (initialised as matrix of 0's). 

X=zeros(Nx,Ny,T); 
omega=2*pi*ones(Nx,Ny); 
omega(end-59:end, :)=2*pi*0.9385; 
omega(end-15:end, :)=2*pi*1.0588; 
S= zeros(Nx,Ny); 
%% Simulations (using Euler method)
% k is time in hours. 
for k=2:T
    %nearest neighbour (coupling) information: 
    for i=1:Nx
        for j=1:Ny
            if template(i,j)~=0 
            A= sin(X(:,:,k-1)- X(i,j, k-1));
            AA=A.*template;
            K=[1 1 1; 1 0 1; 1 1 1]; 
            C=conv2(AA,K,'same');
            S(i,j)=C(i,j); 
            else 
               S(i,j)=0; 
            end
        end
    end
    % Euler method
    X(:,:,k)=X(:,:,k-1).*template+dt*(omega.*template+ kappa*S).*template; 
    X(:,:,k)=mod(X(:,:,k), 2*pi); % modulus 2*pi. 
    clear S; 
end
%% Plot: check of phases for  3 cell types over 144h. 
Fig2=figure;
set(Fig2,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  3000 1200]/200)
plot(squeeze(X(25, 22, :)), 'o-r')
hold on
plot(squeeze(X(60, 22, :)), 'o-b')
hold on
plot(squeeze(X(94, 22, :)), 'o-c')
xlabel('time (h)'); ylabel('angle (rad)'); 
ylim([0 10])
legend({'hypocotyl cell', 'root cell', 'root tip cell'}, 'Location','best')
print('-Painters', Fig2, 'ThreeCellPhases','-dpdf','-r300')
clf
%% Kymograph (raster plot) of the data
for k=1:T
    cc=(cos(X(:,:,k))+1);
    cc=cc.*template;
    Expression(:,k)=sum(cc,2); %sum expression along the plant longitudinally
    clear cc   
end

Expression_longtudinal=diag(1./max(Expression(:,1:144)'))*Expression(:,1:144); %scale by maximum expression in each longitudinal sections
%remove images of the plot that don't contain the plant: 
Expression_longtudinal(1:5, :)=[];
Expression_longtudinal(end-5:end, :)=[];

Fig3=figure;
set(Fig3,'PaperUnits', 'centimeters',  'PaperPosition', [0 0  3000 1200]/300)
imagesc(Expression_longtudinal(end:-1:1,:))
ylabel('Distance from root tip')
xlabel('Time (h)')
set(gca, 'Fontsize', 7); set(gca, 'FontName', 'Helvetica')
ax = gca;
ax.XTick = [0:24:144];
ax.YTick = [0:20:100];
set(gca,'YDir','normal')
print('-Painters', Fig3, 'ModelKymograph','-dpdf','-r300')