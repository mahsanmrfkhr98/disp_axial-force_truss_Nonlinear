close all
clc
clear all
%% 00 - Pre-Definitions
% Total Number of DoFs
  NDoFs=9;
% Total Number of Elements
  NElements=4;
% The DoF to Draw the Plots for
  DoFtoDraw=9;
% Axial Force to Nodal Force Transformation Vector
  QtoP=[-0.529998940003180,0,0;-0.847998304005088,0,0;0,0.242535625036333,0;0,-0.970142500145332,0;0,0,0.868243142124459;0,0,-0.496138938356834;0.529998940003180,-0.242535625036333,-0.868243142124459;0.847998304005088,0.970142500145332,0.496138938356834];
% The Model to draw the Plot for
  Multilinear=0; Parallel=0;
%% 01 - Calculations
SCAN1='%f'; for i=1:NDoFs-1;     SCAN1=[SCAN1 '%f']; end
SCAN2='%f'; for i=1:NElements-1; SCAN2=[SCAN2 '%f']; end
if     Multilinear==1
   system('OpenSees.exe Example2OpenSees_A_Multilinear_Material.tcl');
elseif Parallel==1
   system('OpenSees.exe Example2OpenSees_B_Bilinear_Prallel_Materials.tcl');
end
fid=fopen('OSNodalDisplacements.txt','r');
a=textscan(fid,SCAN1,'CollectOutput',1); Displacements=a{1}; clear a;
fclose(fid);
fid=fopen('OSElementsForce.txt','r');
a=textscan(fid,SCAN2,'CollectOutput',1); Forces=a{1}; clear a;
fclose(fid);
delete('OSNodalDisplacements.txt');
delete('OSElementsForce.txt');
OpenSeesForces=zeros(NDoFs,1);
OpenSeesDisplacements=zeros(NDoFs,1);
for i=1:size(Displacements,1)
OpenSeesDisplacements(:,i+1)=transpose(Displacements(i,:));
OpenSeesForces(:,i+1)=QtoP*transpose(Forces(i,:));
end
clear b; clear Displacements; clear fid; clear i;
figure('units','normalized','outerposition',[0 0 1 1]); 
ylim([1.05*min((OpenSeesForces(DoFtoDraw,:))), 1.05*max(abs(OpenSeesForces(DoFtoDraw,:)))]); 
xlim([1.05*min((OpenSeesDisplacements(DoFtoDraw,:))), 1.05*max(abs(OpenSeesDisplacements(DoFtoDraw,:)))]);
grid on; grid minor; ax=gca; ax.GridLineStyle='--'; ax.GridAlpha=0.6; ax.GridColor=['k']; ax.FontSize=12; ax.LineWidth=0.8; ax.TickLength=[0.01 0.01];    
Title=['DoF ' num2str(DoFtoDraw) ' Load-Dispacement Plot']; title(Title,'fontsize',15,'fontweight','bold'); hold on;
xlabel('Displacement [in]','fontsize',13,'fontweight','bold'); ylabel('Load [lb]','fontsize',13,'fontweight','bold');  
plot(OpenSeesDisplacements(DoFtoDraw,:),OpenSeesForces(DoFtoDraw,:),'b','LineWidth',3);
if     Multilinear==1
       Plot_Name=['Multilinear Material - OpenSees']; hold on;
elseif Parallel==1
       Plot_Name=['Bilinear Parallel Materials - OpenSees']; hold on;
end
h1=suptitle(Plot_Name); set(h1,'FontSize',17,'FontWeight','bold'); hold off;