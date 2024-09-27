close all
clc
clear all
%% 00 - Pre-Definitions
% Total Number of DoFs
  NDoFs=10;
% The DoF to Draw the Plots for
  DoFtoDraw=10;
% External Load Increments Vector
  ExternalLoadIncrements=[0; 0; 0; 0; 0; 0;0;0; 3; 2.5];
%% 01 - Calculations
SCAN='%f'; for i=1:NDoFs-1; SCAN=[SCAN '%f']; end
system('OpenSees.exe Example1OpenSees_Example_File.tcl');
fid=fopen('OSNodalDisplacements.txt','r');
a=textscan(fid,SCAN,'CollectOutput',1); Displacements=a{1}; clear a;
fclose(fid);
delete('OSNodalDisplacements.txt');
OpenSeesForces=zeros(NDoFs,1);
OpenSeesDisplacements=zeros(NDoFs,1);
for i=1:size(Displacements,1)
OpenSeesDisplacements(:,i+1)=transpose(Displacements(i,:));
OpenSeesForces(:,i+1)=(i)*ExternalLoadIncrements;
end
OpenSeesDisplacements=abs(OpenSeesDisplacements);
OpenSeesForces=abs(OpenSeesForces);
clear Displacements; clear fid; clear i;
figure('units','normalized','outerposition',[0 0 1 1]); 
ylim([0, 1.05*max(OpenSeesForces(DoFtoDraw,size(OpenSeesForces,2)))]); 
xlim([0, 1.05*max(OpenSeesDisplacements(DoFtoDraw,size(OpenSeesDisplacements,2)))]);
grid on; grid minor; ax=gca; ax.GridLineStyle='--'; ax.GridAlpha=0.6; ax.GridColor=['k']; ax.FontSize=12; ax.LineWidth=0.8; ax.TickLength=[0.01 0.01];    
Title=['Load-Dispacement Plot for DoF ' num2str(DoFtoDraw) '']; title(Title,'fontsize',15,'fontweight','bold'); hold on;
xlabel('Displacement [in]','fontsize',13,'fontweight','bold'); ylabel('Load [lb]','fontsize',13,'fontweight','bold');  
plot(OpenSeesDisplacements(DoFtoDraw,:),OpenSeesForces(DoFtoDraw,:),'b','LineWidth',3);
Title=['DoF ' num2str(DoFtoDraw) ' Load-Dispacement Plot']; hold on;
Plot_Name=['OpenSees Result'];
h1=suptitle(Plot_Name); set(h1,'FontSize',17,'FontWeight','bold'); hold off;
clear DoFtoDraw;
%% 02 - Save the Results
save('Example1Results_OpenSees')