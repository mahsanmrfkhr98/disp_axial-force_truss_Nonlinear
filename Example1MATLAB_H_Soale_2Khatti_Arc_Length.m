clear all;clc
%% units: lb, psi, in
%% 00 - Pre-Definitions
% Pause Durations when Plotting
  Pause1=0.15; Pause2=1; 
% Disable Plotting Animations
  NoAnimation=0;
% Plot Axis Scale
  XScale=0.45; YScale=0.75; YScalePr=1.1;
% Total Number of Steps
  Radius=0.3655; Betta=0.005;                 % Radius = 0~1. Harchi "Radius" bozorgtar, archaye bozorg tar. "Radius" agar az ye haddi bishtar she momkene "Pext" hadafo rad konim.
% Radius=0.17; Betta=0.16;                 % Meghdare "Betta" ba sa'y va khata peida mishe. Harchi "Betta" bozorg tar, enhenaye arca kamtar mishe.     
% Tolerance for the Error of Calculations
  Tolerance=0.00001;
% The DoF to Draw the Plots for
  DoFtoDraw=10;
% Conect the Last Iteration Points of the Steps
  ConnectSteps=0;
% Draw Residual Nodal Forces Plot or Not
  DrawPrPlots=0;
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
                                                    0
                                                    0
                                                    0           
                                                    0          
                                                    0           
                                                    0           
                                                    0           
                                                    0        
                                                250000           
                                                208300];        
                                              
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
        Pext=Loads;                                                                                         % Bordare "Pext" ta'rif shode ta tebghe notatione kelas pish berim.     
        Arc=0;                                                                                              % Meghdare tooye shomarandeye "Arc" te'dade arc'haye zade shode too har lahze ro neshoon mide. Aval ke hal shoroo' nashode meghdaresh ba 0 barabare.             
        Lambda(1)=0;                                                                                        % Ye bordar ta'rif mikonim ke derayeye k+1, meghdare Lambdaye mohasebe shode too Iteratione "k" az arce "Arc" hast.         
        deltaD1=zeros(NDoFs,1);                                                                             % "deltaD1" ye bordare ke derayeye i, meghdare deltaD1 mohasebe shode baraye daraje azadie i, too Iteratione "k" az arce "Arc" hast. Choon "deltaD1" size moshakhas dare bayad inja ta'rif beshe. Too in mesal moshkeli baraye ma pish nemiad, amma farz konid too ye mesale dige darajate 5 va 6 azadan. Agar "deltaD1" az ghabl ta'rif nashode bashe, vaghti code be khatte 121 mirese, faghat ye bordare 1x6 misaze.  
        deltaD2=zeros(NDoFs,1);                                                                             % "deltaD2" ye bordare ke derayeye i, meghdare deltaD2 mohasebe shode baraye daraje azadie i, too Iteratione "k" az arce "Arc" hast. Choon "deltaD2" size moshakhas dare bayad inja ta'rif beshe. Too in mesal moshkeli baraye ma pish nemiad, amma farz konid too ye mesale dige darajate 5 va 6 azadan. Agar "deltaD2" az ghabl ta'rif nashode bashe, vaghti code be khatte 122 mirese, faghat ye bordare 1x6 misaze.              
        Pint=zeros(NDoFs,1);                                                                                % Agar inja bordare "Pint" ta'rif nashe, too khatte 105 va 122 error migirim.          
        D0=zeros(NDoFs,1);                                                                                  % D0 ye bordare ke derayehaye 2i-1 va 2i, meghdare kolle jabejaiie node i dar jahate x va y ta akhare arce "Arc-1" hastan. Aagar inja ta'rif nashe too khatte 124 error migirim.   
        D=zeros(NDoFs,1);                                                                                   % D ye bordare ke derayehaye 2i-1 va 2i, meghdare kolle jabejaiie node i dar jahate x va y ta akhare Iteratione "k" az arce "Arc" hastan. Aagar inja ta'rif nashe too khatte 124 error migirim.   
        ResultNodalResidualForces(:,1)=Pext;                                                                % In matrice baraye rasme nemoodare Presidual-Displacement taarif shode. Choon avale kar saze saken hast, bayad setoone aval ba niroohaye khrejie avale kar barabar bashe.
        ResultNodalResidualForcesDisplacement(:,1)=zeros(NDoFs,1);                                          % In matrice baraye rasme nemoodare Presidual-Displacement taarif shode. Choon avale kar saze saken hast, bayad setoone aval sefr bashe.
        Lambda0=0;                                                                                          % Lambda0 ye moteghayere ke meghdare kolle Lambdaye mohasebe shode ta akhare arce "Arc-1" ro neshoon mide. Agar inja ta'rif nashe too khatte 155 error migirim. Meghdare Lambda0 bayad payane kar be 1 berese, agar bishtar she neshoon mide az "Pext" rad shodim.
        k0=0;                                                                                               % k0 ye shomarandas ke meghdare kolle Iterationaye anjam shode ta akhare arce "Arc-1" ro neshoon mide. Agar inja ta'rif nashe too khatte 157 error migirim.         
        Reached(1)=1;                                                                                       % Reached = ye bordar ke derayeye a+1, beramoon moshakhas mikone kodoom setoon az matricaye "Result..." motenazer ba noghte payanie arce a hast. Felan ye derayas, too ravande hal bordar mishe.
 % 01 - Solve                                                                                            %% Pext az aval ta akhare kar bayad ba kolle "Loads" barabar bashe va be hich vajh nabayad update she.        
   tic;
   Pr=Pext-Pint;                                                                                            % Avale kar baraye in ke betoonim varede loop beshim ye "Pr" mohasebe mikonim. Choon avalesh "Pint" 0 hast, meghdare "Pr" ba "Pext" barabar mishe.
   while norm(Pr(FreeDoFs))>=Tolerance && Lambda0<1                                                         % Loope aval ba check kardane "Pr" baraye ma moshakhas mikone ke be "Pext" residim ya na. Sharte dovome [Lambda0<1] tazmin mikone az "Pext" rad nashim. 
         Arc=Arc+1;                                                                                         % Har bar ke varede loope aval beshim code ye arc mizane, pas be shomarandeye "Arc" yki ezaf mikonim.
         k=0;                                                                                               % Avale har arc, shomarandeye "k" ro sefr mikonim ke te'dade Iterationaye anjam shode too arce jadido beshmorim.
         Lambda(1)=0;                                                                                       % Avale har arc, "Lambda" ro sefr mikonim, choon too formoole Arc-Length az Lambdaye kolli estefade nemishe. 
         deltaLambda=10;                                                                                    % Baraye vorood be arc, be "deltaLambda" alaki ye adade bozorg midim ke sharte vorood be arc baravorde she.
         while abs(deltaLambda)>=Tolerance                                                                  % Hadaf too har arc, koochick kardane "deltaLambda" hast. Ta vaghti in shart baravorde nashe, be Iteration too hamoon arci ke hastim edame midim. 
               k=k+1;                                                                                       % Har bar ke "deltaLambda" be andazeye kafi koochik nabashe, be tedade iterationaye too oon arc (k) ye vahed ezafe mishe.
             % Assemble the Global Stifness Matrix
               if k0==0 && k==1                                                                             % Raveshe jadid baraye estefade az Initial Stiffness Method, faghat vaghti too avalin Iteration az kolle Iterationa bashim (hanooz naraftim arce ba'di ke "k0" chizi gheir az 0 beshe va taze "k" barabare 1 shode) matrice sakhti update mishe.
                  Kt=zeros(NDoFs,NDoFs);
                  for i=1:NElements
                      Kt(DoFs(i,:),DoFs(i,:))=Kt(DoFs(i,:),DoFs(i,:))+(transpose(a(i,:))*Kelt(i,i)*a(i,:));
                  end                                                                                                                                                                                                                                           
               end                                                                                                                                                                                                                                              
             % Calculate deltaD1 and deltaD2                                                                                                                                                                                                                    
               deltaD1(FreeDoFs)=Kt(FreeDoFs,FreeDoFs)\Pext(FreeDoFs);                                                                                                                                                                                          
               deltaD2(FreeDoFs)=Kt(FreeDoFs,FreeDoFs)\(((Lambda0+Lambda(k))*Pext(FreeDoFs))-Pint(FreeDoFs));                                                                                                                                                   
             % Calculate deltaLambda                                                                       
               fun = @(dL)((Radius^2)-((Lambda(k)+dL)^2)-(((D(FreeDoFs)-D0(FreeDoFs)+deltaD2(FreeDoFs)+(dL*deltaD1(FreeDoFs)))'*(D(FreeDoFs)-D0(FreeDoFs)+deltaD2(FreeDoFs)+(dL*deltaD1(FreeDoFs))))/Betta));
               options=optimset('MaxFunEvals',100,'MaxIter',100,'Display','off');                          
               [deltaLambda,~]=fsolve(fun,1000,options);                                                    % Inja choon midoonim "deltaLambda" risheye mosbate mo'adeleye beyzi has ke beyn 0 va 1 ham hast, be onvane hadse avalie bayad hatman ye adade mosbate bozorg tar az 1 bedim. Man 1000 gozashtam, ba 1 ham hal mishe. Amma ba gozashtane a'dade manfi, a'dade mosbate kamtar az 1, va hata 0 be onvane hadse avalie, javabe ghalat mide.
             % Calculate the New Lambda
               Lambda(k+1)=Lambda(k)+deltaLambda; 
             % Calculate Total Nodal Displacements Until the Curent Iteration
               D=D+deltaD2+(deltaLambda*deltaD1); 
             % Calculate Total Element Deformations Until the Curent Iteration
               V=DtoV*D; 
             % Calculate Elements Total Internal Forces Until the Curent Iteration
               for i=1:NElements
                   if abs(V(i))<Vy(i)
                      Q(i,1)=Kel(i,i)*V(i);
                   else
                      Q(i,1)=sign(V(i))*(Qy(i)+(alpha(i)*Kel(i,i)*(abs(V(i))-Vy(i)))); 
                      Kelt(i,i)=alpha(i)*Kel(i,i);                                                                          
                   end                      
               end
             % Calculate the Total Internal Nodal Loads Until the Curent Iteration
               Pint=QtoP*Q;
             % Update the Residual Nodal Loads Vector
               Pr=Pext-Pint;
             % Collect the Results
               ResultNodalForces(:,(2*(k0+k)))=((Lambda0+Lambda(k+1))*Pext);   ResultNodalForces(:,(2*(k0+k))+1)=Pint;
               ResultNodalDisplacements(:,(2*(k0+k)))=D;                       ResultNodalDisplacements(:,(2*(k0+k))+1)=D; 
               ResultNodalResidualForces(:,k0+k+1)=Pr; 
               ResultNodalResidualForcesDisplacement(:,k0+k+1)=D;
               ResultMemberDeformations(:,k0+k+1)=V; 
               ResultMemberAxialForces(:,k0+k+1)=Q;
         end
         D0=D;                                                                                              % Ba'd az etmame har arc, "D0" update mishe ke nedoonim ta akhare oon arc noda cheghad jabeja shodan. 
         Lambda0=Lambda0+Lambda(k);                                                                         % Ba'd az etmame har arc, Lambdaye naha'iie arc (Lambda(k)) ba Lambdaye nahaii ta arce ghabl (Lambda0) jam' mishe va "Lambda0" update mishe.
         clear Lambda;                                                                                      % Ba'd az etmame har arc, Lambdahaye hesab shode tooye arce tamoom shode pak mishan, choon estefade'ii nadaran.        
         k0=k0+k;                                                                                           % Ba'd az etmame har arc, teedade iterationaye anjam shode too oon arc (k) ba teedade kolle iterationaye anjam shode ta akhare arce ghabl (k0) jam' mishe va "k0" update mishe.      
         Reached(Arc+1)=k0;                                                                                       
   end
   Duration=toc; 
%% 04 - Visualization
close all;
ResultNodalForces=abs(ResultNodalForces);
ResultNodalDisplacements=abs(ResultNodalDisplacements);
ResultNodalResidualForcesDisplacement=abs(ResultNodalResidualForcesDisplacement);
ResultMemberDeformations=abs(ResultMemberDeformations);
ResultMemberAxialForces=abs(ResultMemberAxialForces);
if DrawAllFreeDoFsPlots==1 || DrawPrPlots==1; Font=12; else; Font=15; end
load('Example1Results_OpenSees')
figure('units','normalized','outerposition',[0 0 1 1]); 
if     DrawAllFreeDoFsPlots==1 && DrawPrPlots==1
       ii=1:2*size(Target,1);
elseif DrawAllFreeDoFsPlots==1 && DrawPrPlots==0
       ii=1:size(Target,1);
elseif DrawAllFreeDoFsPlots==0 && DrawPrPlots==1
       ii=1:2;  
elseif DrawAllFreeDoFsPlots==0 && DrawPrPlots==0
       ii=1;
end
   for i=ii
       % Predefinition
         if     DrawAllFreeDoFsPlots==1 && DrawPrPlots==1 && i-(2*floor((i/2)))==1
                subplot(2,2,floor(i/2)+1); DoFtoDraw=Target(floor(i/2)+1);
                DrawPr=0;
         elseif DrawAllFreeDoFsPlots==1 && DrawPrPlots==1 && i-(2*floor((i/2)))==0
                if   i==2
                     subplot(2,2,i+1); 
                else
                     subplot(2,2,i);   
                end
                     DrawPr=1;         DoFtoDraw=Target(i/2);
         elseif DrawAllFreeDoFsPlots==1 && DrawPrPlots==0
                     subplot(1,2,i);   DoFtoDraw=Target(i);
                     DrawPr=0;
         elseif DrawAllFreeDoFsPlots==0 && DrawPrPlots==1
                     subplot(2,1,i);
                if   mode(i,2)==1
                     DrawPr=0;
                else
                     DrawPr=1;
                end
         elseif DrawAllFreeDoFsPlots==0 && DrawPrPlots==0
                     DrawPr=0;
         end


         grid on; grid minor; ax=gca; ax.GridLineStyle='--'; ax.GridAlpha=0.6; ax.GridColor=['k']; ax.FontSize=12; ax.LineWidth=0.8; ax.TickLength=[0.01 0.01];    
         xlabel('Displacement [in]','fontsize',13,'fontweight','bold');
         if     DrawPr==0
                Title=['DoF ' num2str(DoFtoDraw) ' Load-Dispacement Plot'];
                ylabel('Load [lb]','fontsize',13,'fontweight','bold');      
         elseif DrawPr==1
                Title=['DoF ' num2str(DoFtoDraw) ' Residual Load-Dispacement Plot'];
                ylabel('Residual Load [lb]','fontsize',13,'fontweight','bold');    
         end
         xlim([0, XScale*OpenSeesDisplacements(DoFtoDraw,size(OpenSeesDisplacements,2))]);
         if     DrawPr==0
                ylim([0, YScale*OpenSeesForces(DoFtoDraw,size(OpenSeesDisplacements,2))]); 
         elseif DrawPr==1 
                ylim([min(0,YScalePr*min(ResultNodalResidualForces(DoFtoDraw,:))), YScalePr*max(ResultNodalResidualForces(DoFtoDraw,:))]);     
         end
         title(Title,'fontsize',15,'fontweight','bold'); hold on;
       % Plotting - Real Behaviour
         if     DrawPr==0
                P2=plot(OpenSeesDisplacements(DoFtoDraw,:),OpenSeesForces(DoFtoDraw,:),'b','LineWidth',2);
                if DrawPrPlots==1 || DrawAllFreeDoFsPlots==1 || ConnectSteps==1 || NoAnimation==1
                else; pause(Pause2);
                end
                
         end
       % Plotting - Target Load
         if     DrawPr==0
                P1=plot([0,XScale*OpenSeesDisplacements(DoFtoDraw,size(OpenSeesDisplacements,2))],[Pext(DoFtoDraw),Pext(DoFtoDraw)],'r','LineWidth',1.5);
                if DrawPrPlots==1 || DrawAllFreeDoFsPlots==1 || ConnectSteps==1 || NoAnimation==1
                else; pause(Pause1);
                end
         end
       % Plotting - Matlab Results 
         if     DrawPr==0
                for a=1:Arc
                    if a==1
                       scatter(ResultNodalDisplacements(DoFtoDraw,2*Reached(a)-1),ResultNodalForces(DoFtoDraw,2*Reached(a)-1),'m','filled'); hold on;
                       if DrawPrPlots==1 || DrawAllFreeDoFsPlots==1 || ConnectSteps==1 || NoAnimation==1
                       else; pause(Pause1);
                       end
                    end
                    if a==Arc
                       P5=plot(ResultNodalDisplacements(DoFtoDraw,2*Reached(a)+2:2:2*Reached(a+1)),ResultNodalForces(DoFtoDraw,2*Reached(a)+2:2:2*Reached(a+1)),'c','LineWidth',1.5); hold on;
                       if DrawPrPlots==1 || DrawAllFreeDoFsPlots==1 || ConnectSteps==1 || NoAnimation==1
                       else; pause(0.5*Pause1);
                       end
                       P3=plot(ResultNodalDisplacements(DoFtoDraw,2*Reached(a)-1:2*Reached(a+1)+1),ResultNodalForces(DoFtoDraw,2*Reached(a)-1:2*Reached(a+1)+1),'g','LineWidth',1.5); hold on;
                       P4=scatter(ResultNodalDisplacements(DoFtoDraw,2*Reached(a)+3:2:2*Reached(a+1)+1),ResultNodalForces(DoFtoDraw,2*Reached(a)+3:2:2*Reached(a+1)+1),'m','filled'); hold on;
                       P6=scatter(ResultNodalDisplacements(DoFtoDraw,2*Reached(a+1)+1),ResultNodalForces(DoFtoDraw,2*Reached(a+1)+1),[],[0.4940 0.1840 0.5560],'filled'); hold on;
                    else
                       plot(ResultNodalDisplacements(DoFtoDraw,2*Reached(a)+2:2:2*Reached(a+1)),ResultNodalForces(DoFtoDraw,2*Reached(a)+2:2:2*Reached(a+1)),'c','LineWidth',1.5); hold on;
                       if DrawPrPlots==1 || DrawAllFreeDoFsPlots==1 || ConnectSteps==1 || NoAnimation==1
                       else; pause(0.5*Pause1);
                       end
                       plot(ResultNodalDisplacements(DoFtoDraw,2*Reached(a)-1:2*Reached(a+1)+1),ResultNodalForces(DoFtoDraw,2*Reached(a)-1:2*Reached(a+1)+1),'g','LineWidth',1.5); hold on;
                       scatter(ResultNodalDisplacements(DoFtoDraw,2*Reached(a)+3:2:2*Reached(a+1)+1),ResultNodalForces(DoFtoDraw,2*Reached(a)+3:2:2*Reached(a+1)+1),'m','filled'); hold on;
                       scatter(ResultNodalDisplacements(DoFtoDraw,2*Reached(a+1)+1),ResultNodalForces(DoFtoDraw,2*Reached(a+1)+1),[],[0.4940 0.1840 0.5560],'filled'); hold on;
                       if DrawPrPlots==1 || DrawAllFreeDoFsPlots==1 || ConnectSteps==1 || NoAnimation==1
                       else; pause(Pause1);
                       end
                    end

                end
            PLOT(1)=1;
            for a=1:Arc
                PLOT(2*a:2*a+1)=[2*Reached(a), 2*Reached(a+1)+1];
            end
            
         if ConnectSteps==1
            scatter(ResultNodalDisplacements(DoFtoDraw,1:2:2*Reached(a+1)+1),ResultNodalForces(DoFtoDraw,1:2:2*Reached(a+1)+1),'m','filled'); hold on;
            scatter(ResultNodalDisplacements(DoFtoDraw,PLOT),ResultNodalForces(DoFtoDraw,PLOT),[],[0.4940 0.1840 0.5560],'filled'); hold on;
            if DrawPrPlots==1 || DrawAllFreeDoFsPlots==1
            else; pause(Pause2);
            end
            P7=plot(ResultNodalDisplacements(DoFtoDraw,PLOT),ResultNodalForces(DoFtoDraw,PLOT),'color',[0.4660 0.6740 0.1880],'LineWidth',3); hold on;
            scatter(ResultNodalDisplacements(DoFtoDraw,PLOT),ResultNodalForces(DoFtoDraw,PLOT),[],[0.4940 0.1840 0.5560],'filled'); hold on;
            if DrawPrPlots==1 || DrawAllFreeDoFsPlots==1 || NoAnimation==1
            scatter(ResultNodalDisplacements(DoFtoDraw,1:2:2*Reached(a+1)+1),ResultNodalForces(DoFtoDraw,1:2:2*Reached(a+1)+1),'m','filled'); hold on;
            scatter(ResultNodalDisplacements(DoFtoDraw,PLOT),ResultNodalForces(DoFtoDraw,PLOT),[],[0.4940 0.1840 0.5560],'filled'); hold on;
            pause(0.5*Pause1);
               leg=legend([P2, P1, P5, P3, P4, P6, P7],{'Real Behaviour from OpenSees', 'Target Load', 'MATLAB Arcs', 'MATLAB Iterations', 'MATLAB Trial Points', 'MATLAB Final Points', 'MATLAB Result Curve'},'Location','southeast','fontsize',Font); 
            else
                pause(0.5*Pause1);
               leg=legend([P2, P1, P5, P3, P4, P6, P7],{'Real Behaviour from OpenSees', 'Target Load', 'MATLAB Arcs', 'MATLAB Iterations', 'MATLAB Trial Points', 'MATLAB Final Points', 'MATLAB Result Curve'},'Location','southeast','fontsize',Font); 
            end
         else
            if DrawPrPlots==1 || DrawAllFreeDoFsPlots==1 || NoAnimation==1
               scatter(ResultNodalDisplacements(DoFtoDraw,1:2:2*Reached(a+1)+1),ResultNodalForces(DoFtoDraw,1:2:2*Reached(a+1)+1),'m','filled'); hold on;
               scatter(ResultNodalDisplacements(DoFtoDraw,PLOT),ResultNodalForces(DoFtoDraw,PLOT),[],[0.4940 0.1840 0.5560],'filled'); hold on;
               pause(0.5*Pause1);
               leg=legend([P2, P1, P5, P3, P4, P6],{'Real Behaviour from OpenSees', 'Target Load', 'MATLAB Arcs', 'MATLAB Iterations', 'MATLAB Trial Points', 'MATLAB Final Points'},'Location','southeast','fontsize',Font);
            else
                pause(0.5*Pause1);
               leg=legend([P2, P1, P5, P3, P4, P6],{'Real Behaviour from OpenSees', 'Target Load', 'MATLAB Arcs', 'MATLAB Iterations', 'MATLAB Trial Points', 'MATLAB Final Points'},'Location','southeast','fontsize',Font);
            end
         end
         elseif DrawPr==1
                plot(ResultNodalResidualForcesDisplacement(v,:),ResultNodalResidualForces(DoFtoDraw,:),'g','LineWidth',1.5); hold on; 
         end

   end
Plot_Name=['Bilinear Material - Arc-Length Method - ' num2str(Arc) ' Arcs - Duration ' num2str(Duration) ' Seconds - ' num2str(k0) ' Iterations']; hold on;
h1=suptitle(Plot_Name); set(h1,'FontSize',17,'FontWeight','bold'); hold off;