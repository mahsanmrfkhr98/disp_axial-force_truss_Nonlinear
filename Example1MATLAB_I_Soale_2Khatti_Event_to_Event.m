clc;clear all;close all
%% units: lb, psi, in
%% 00 - Pre-Definitions
% Pause Durations when Plotting
  Pause1=0.15; Pause2=1; 
% Disable Plotting Animations
  NoAnimation=0;
% Plot Axis Scale
  XScale=0.4; YScale=0.7;
% The DoF to Draw the Plots for
  DoFtoDraw=9;                             % Shekastegie nemoodar lozooman too plote hameye darajate azadi moshakhas nist. masalan vaghi nemoodara ro baraye daraje azadie 7 bekeshid, shekastegie noghteye akahar asan be cheshm moshakhas nist.         
% Draw the Plots for All Free DoFs or Not
  DrawAllFreeDoFsPlots=0;
%% 01 - Input
% Nodes Data Matrix:                               Node#       x[in]      y[in]
                                     Nodes=[
                                                    1,           0,         12
                                                    2,           0,         0
                                                    3,          17,         0
                                                    4,          24,         12
                                                    5,           9,         17];
                           
% Elements Data Matrix:                          Element#   Start_Node   End_Node  A[in^2]   E[psi]     Fy[psi]    alpha                                 
                                  Elements=[
                                                    1,           1          5,       5,     36000000     25000     0.03;
                                                    2,           2          5,       4,     36000000     25000     0.03;
                                                    3,           3          5,       5,     36000000     25000     0.03;
                                                    4,           4          5,       6,     36000000     25000     0.03];
                           
% Restrains Vector:                                DoFs
                                 Restrains=[
                                                    1;
                                                    2;
                                                    3;
                                                    4;
                                                    5;
                                                    6;
                                                    7;
                                                    8];
                           
% External Loads Vector:                      Load_Value[lb]   DoF#  
                                     Loads=[
                                                    0           %1
                                                    0           %2
                                                    0           %3
                                                    0           %4
                                                    0           %5
                                                    0           %6
                                                    0           %7
                                                    0           %8
                                                    3           %9
                                                  2.5];         %10 
                                              
% Target DoFs to Draw Load Displacement Plot:     DoFs
                                    Target=[
                                                    9
                                                    10];
%% 02 - Complie Inputs
% Nodes Data                                                                                               
  % Available Data                                                                        
    x=Nodes(:,2);                                                                              
    y=Nodes(:,3);                                                                              
  % New Data                                                                               
    % Number of Nodes                                                                
      NNodes=size(Nodes,1);                                                                   
      NDoFs=2*NNodes;                                                                   
    % Free DoFs                                                                           
      FreeDoFs=(1:NDoFs)';                                                                      
      FreeDoFs(Restrains)=[];                                                                   
% Elements Data
  % Available Data                                                                     
    Start=Elements(:,2);                                                                        
    End=Elements(:,3);                                                                          
    A=Elements(:,4);                                                                            
    E=Elements(:,5);
    Fy=Elements(:,6);                                                                                    
    alpha=Elements(:,7);                                                                                 
  % New Data                                                                           
    % Number of Elements                                                                   
      NElements=size(Elements,1);     
    % Geometric and Mechanical Info                                                       
      DtoV=zeros(NElements,NDoFs);                                                                                                                    
      QtoP=zeros(NDoFs,NElements);
      for i=1:NElements                                                                                                                                                              
          DoFs(i,:)=[2*Start(i)-1  2*Start(i)  2*End(i)-1  2*End(i)];                                              
          invL(i)=1/((x(End(i))-x(Start(i)))^2+(y(End(i))-y(Start(i)))^2)^0.5;                                                               
          C(i)=(x(End(i))-x(Start(i)))*invL(i);                                                                                              
          S(i)=(y(End(i))-y(Start(i)))*invL(i);                                                                                              
          a(i,:)=[-C(i) -S(i) C(i) S(i)];                                            
          DtoV(i,DoFs(i,:))=a(i,:);                                                                                        
          QtoP(DoFs(i,:),i)=transpose(a(i,:));                                                                                        
          Kel(i,i)=A(i)*E(i)*invL(i);                                                                                                                                                      
          Qy(i,1)=Fy(i)*A(i);                                                                
          Vy(i,1)=Qy(i)/Kel(i,i);                                                                                                                   
      end 
      Kelt=Kel;                                                                                         
%% 03 - Calculations
 % 00 - Pre-Definitions
        D=zeros(NDoFs,1);                                                                           % Avale kar bordaraii ke gharare ye size moshakhas date bashan va faghat ye seri derayehaye khaseshoon por she ro ta'rif mikonim.
        ResultNodalDisplacements=zeros(NDoFs,1);                                                    % In bordar baraye ja'm avari result estefade mishe (derayehaye radifaye 2*j-1 va 2j too setoone i+1, be tartib meghdare taghiire makane node j, be tartib dar jahate x va y too lahzeye taslime elemane i hastan) va be tore ma'mool bayad hamoon entehaye har step behesh adad dade beshe va niazi naboode inja ta'rif she. Amma mikhaim az in bordar too toole hal too hamoon stepe avval estefade konim, pas bayad inja taarif she chie vagarna error migirim.
        ResultMemberDeformations=zeros(NElements,1);                                                % In bordar baraye ja'm avari result estefade mishe (derayeye radife j va setoone i+1, meghdare taghiire toole mehvarie too elemane j too lahzeye taslime elemane i hast) va be tore ma'mool bayad hamoon entehaye har step behesh adad dade beshe va niazi naboode inja ta'rif she. Amma mikhaim az in bordar too toole hal too hamoon stepe avval estefade konim, pas bayad inja taarif she chie vagarna error migirim.
        ResultMemberAxialForces=zeros(NElements,1);                                                 % In bordar baraye ja'm avari result estefade mishe (derayeye radife j va setoone i+1, meghdare nirooye mehvarie too elemane j too lahzeye taslime elemane i hast) va be tore ma'mool bayad hamoon entehaye har step behesh adad dade beshe va niazi naboode inja ta'rif she. Amma mikhaim az in bordar too toole hal too hamoon stepe avval estefade konim, pas bayad inja taarif she chie vagarna error migirim.
        Yielded=[];                                                                                 % In ye bordare ke derayeye i behemoon mige too stepe i kodoom eleman taslim shode. Man ye doone "if" neveshtam ke agar too stepe aval hastim, az in bordar va bordare khate bala estefade nashe. Amma pish-ta'rife in do ta bordaro neveshtam ke betoonim bedoone oon "if" ham code bezanim.
 % 01 - Solve
        tic;
        for i=1:NElements                                                                           % Too Event-to-Event baraye eleman ba materiale 2khatti, be te'dade elemanaye saze step darim.      
            % Assemble the Global Stifness Matrix
              Kt=zeros(NDoFs,NDoFs);                                                                % Matrice sakhti inja update mishe.                                       
              for j=1:NElements
                  Kt(DoFs(j,:),DoFs(j,:))=Kt(DoFs(j,:),DoFs(j,:))+(transpose(a(j,:))*Kelt(j,j)*a(j,:));
              end
            % Calculate the Nodal Displacements Vectors
              D(FreeDoFs)=Kt(FreeDoFs,FreeDoFs)\Loads(FreeDoFs); 
            % Calculate the Element Deformations Vectors
              V=DtoV*D;
            % Calculate the Elements Internal Forces Vector
              Q=Kelt*V;
            % Calculate the deltaQ Vector                                                            % Inja ba tavajoh be in ke too stepe aval hastim ya na, baraye mohasebeye deltaQ ye formool entekhab mishe.                     
              if i==1
                 deltaQ=Qy; 
              else
                 deltaQ=abs(Qy-abs(ResultMemberAxialForces(:,i)));                   %%1             % DeltaQ = Fasele ta noghteye taslim.                 
              end
%                deltaQ=abs(Qy-abs(ResultMemberAxialForces(:,i)));                   %%1             % Choon ghabl aza shoroo' matrice "ResultMemberAxialForces" ro ta'rif kardim va choon khoroojiaye stepe i ro too setoone i+1 matricaye "Result" mirizim, mitoonim kolle "if" ro pak konim va faghat hamin khatto bezarim.     
            % Calculate the Lambda Vector                                                            % Inja baraye har eleman, [Lambda=deltaQ/Q] hesab mishe. | Lambda = ye matrice ke dereyeye i setoone 1, shmoareye elemane i hast, va derayeye i setoone 2, meghdare Lambdaye mohasebe shode baraye oon eleman hast.   
              Lambda(:,1)=(1:NElements)';
              for j=1:NElements
                  Lambda(j,2)=abs(deltaQ(j)/Q(j));                                   %%2 
              end
%                 Lambda(:,2)=abs(deltaQ./Q);                                        %%2             
            % Remove Previously Yielded Elements from Lambda Vector                                  % Inja satraye marboot be elemanaii ke ghablan taslim shodan, az matrice Lambda hazf mishan ta dige Lambdahashoon too peyda kardane Lambdaye minimum dekhalat nakone. 
              if i>1 
                 Lambda(Yielded(:),:)=[];                                            %%3 
              end
%                Lambda(ResultYielded(:),:)=[];                                      %%3             % Choon ghabl az shroo' bordare "ResultYielded" ro taarif kardim (va avale kar khalie) mitoonim kollem "if" ro pak konim va faghat ham ye satro bezarim.             
            % Identify the Most-Probable-to-Yield Element
              Lambda=sortrows(Lambda,2);                                                             % Aval ba in dastoor redifaye matrice "Lambda" ro ba tavajoh be setoone 2 (minimum be maximum) sort mikonim. In dastoor mese dokmeye "Filter" too "Excel" kar mikone.
              Yielded(i)=Lambda(1,1);                                                                % Derayeye avale setoone avale matrice "Lambda" mishe shomare elemani ke too in step taslim shode.
              minLambda(i)=Lambda(1,2);                                                              % Derayeye avale setoone dovome matrice "Lambda" mishe minimum Lambdaye mohasebe shode (meghdare Lambdaye mohasebe shode baraye elemani ke taslim shode).
              Kelt(Lambda(1,1),Lambda(1,1))=alpha(Lambda(1,1))*Kel(Lambda(1,1),Lambda(1,1));         % Sakhtie locale ozvi ke too in marhale taslim shode update mishe. 
              clear Lambda;                                                                          % Lambda pak mishe ke too stepe ba'd dobare sakhte she. Te'dade setoonaye lambda barabare "NElements" Boode va too khatte 119 ye seri az radifash pak shode. Agar inja pakesh nakonim, too stepe ba'd az hamin matrie naghes estefade mishe va too khatte 112 error migirim ke te'dade radifaye "Lambda" kamtar az "NElements" hast.
            % Collect the Results                                                                    % Inja moshabehe formulaye to jozve khroojiaro mohasebe mikonim. 
              ResultYieldings(i)=Yielded(i);
              ResultLambda(i)=minLambda(i);
              ResultNodalForces(:,i+1)=sum(ResultLambda(:))*Loads;
              if i==1
                 ResultNodalDisplacements(:,i+1)=ResultLambda(i)*D; 
                 ResultMemberDeformations(:,i+1)=ResultLambda(i)*V; 
                 ResultMemberAxialForces(:,i+1)=ResultLambda(i)*Q;
              else 
                 ResultNodalDisplacements(:,i+1)=ResultNodalDisplacements(:,i)+ResultLambda(i)*D; 
                 ResultMemberDeformations(:,i+1)=ResultMemberDeformations(:,i)+ResultLambda(i)*V;
                 ResultMemberAxialForces(:,i+1)=ResultMemberAxialForces(:,i)+ResultLambda(i)*Q; 
              end                
%                ResultNodalDisplacements(:,i+1)=ResultNodalDisplacements(:,i)+ResultLambda(i)*D;    % Choon khoroojiaye stepe i ro too setoone i+1 matricaye "Result" mirizim, va choon ghabl az shoroo' be setoone avae in matrica 0 dadim, mitoonim kolle "if" ro pak konim va faghat hamin khatto bezarim.     
%                ResultMemberDeformations(:,i+1)=ResultMemberDeformations(:,i)+ResultLambda(i)*V;
%                ResultMemberAxialForces(:,i+1)=ResultMemberAxialForces(:,i)+ResultLambda(i)*Q; 
        end
        Duration=toc;
%% 04 - Visualization
close all;
ResultNodalForces=abs(ResultNodalForces);
ResultNodalDisplacements=abs(ResultNodalDisplacements);
ResultMemberDeformations=abs(ResultMemberDeformations);
ResultMemberAxialForces=abs(ResultMemberAxialForces);
load('Example1Results_OpenSees')
figure('units','normalized','outerposition',[0 0 1 1]); 
if DrawAllFreeDoFsPlots==1
    ii=1:size(Target,1);
else
    ii=1;
end
   for i=ii
       % Predefinition
         if DrawAllFreeDoFsPlots==1
            subplot(1,2,i);
            DoFtoDraw=Target(i);
         end
         ylim([0, OpenSeesForces(DoFtoDraw,size(OpenSeesDisplacements,2))]); 
         xlim([0, OpenSeesDisplacements(DoFtoDraw,size(OpenSeesDisplacements,2))]);
         grid on; grid minor; ax=gca; ax.GridLineStyle='--'; ax.GridAlpha=0.6; ax.GridColor=['k']; ax.FontSize=12; ax.LineWidth=0.8; ax.TickLength=[0.01 0.01];    
         xlabel('Displacement [in]','fontsize',13,'fontweight','bold'); ylabel('Load [lb]','fontsize',13,'fontweight','bold');      
         Title=['DoF ' num2str(DoFtoDraw) ' Load-Dispacement Plot']; title(Title,'fontsize',15,'fontweight','bold'); hold on;
       % Plotting - Real Behaviour
         P1=plot(OpenSeesDisplacements(DoFtoDraw,:),OpenSeesForces(DoFtoDraw,:),'b','LineWidth',2);
         if DrawAllFreeDoFsPlots==0 && NoAnimation==0
            pause(Pause2);
         end
       % Plotting - Matlab Results
         scatter(ResultNodalDisplacements(DoFtoDraw,1),ResultNodalForces(DoFtoDraw,1),'m','filled'); hold on;
         if DrawAllFreeDoFsPlots==0 && NoAnimation==0
            pause(0.5*Pause2);
         end
         for j=1:NElements
             if j==NElements
                plot(ResultNodalDisplacements(DoFtoDraw,j:j+1),ResultNodalForces(DoFtoDraw,j:j+1),'g','LineWidth',1.5); hold on;
                P2=scatter(ResultNodalDisplacements(DoFtoDraw,j+1),ResultNodalForces(DoFtoDraw,j+1),'m','filled'); hold on;
                TXT=['Element ' num2str(ResultYieldings(j)) ' yields in this point.'];
                text(ResultNodalDisplacements(DoFtoDraw,j+1)+0.0010,ResultNodalForces(DoFtoDraw,j+1)-1500,TXT,'HorizontalAlignment','left','VerticalAlignment','top','Color','Magenta','FontSize',13,'FontWeight','bold');
                clear TXT;
             else
                plot(ResultNodalDisplacements(DoFtoDraw,j:j+1),ResultNodalForces(DoFtoDraw,j:j+1),'g','LineWidth',1.5); hold on;
                scatter(ResultNodalDisplacements(DoFtoDraw,j+1),ResultNodalForces(DoFtoDraw,j+1),'m','filled'); hold on;
                TXT=['Element ' num2str(ResultYieldings(j)) ' yields in this point.'];
                text(ResultNodalDisplacements(DoFtoDraw,j+1)+0.0010,ResultNodalForces(DoFtoDraw,j+1)-1500,TXT,'HorizontalAlignment','left','VerticalAlignment','top','Color','Magenta','FontSize',13,'FontWeight','bold');
                clear TXT; 
                if DrawAllFreeDoFsPlots==0 && NoAnimation==0
                   pause(0.5*Pause2);
                end
             end
         end 
                if NoAnimation==1
                scatter(ResultNodalDisplacements(DoFtoDraw,:),ResultNodalForces(DoFtoDraw,:),'m','filled'); hold on;
                end
                pause(0.5*Pause1);
                leg=legend([P1, P2],{'Real Behaviour from OpenSees', 'MATLAB Events'},'Location','southeast','fontsize',15);  

   end
Plot_Name=['Bilinear Material - Event-to-Event Method - ' num2str(NElements) ' Steps - Duration ' num2str(Duration) ' Seconds']; hold on;
h1=suptitle(Plot_Name); set(h1,'FontSize',17,'FontWeight','bold'); hold off;