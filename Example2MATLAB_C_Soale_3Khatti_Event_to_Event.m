clear all
clc
%% units: lb, psi, in
%% 00 - Pre-Definitions
% Pause Durations when Plotting
  Pause1=0.15; Pause2=1; 
% Disable Plotting Animations
  NoAnimation=0;
% Plot Axis Scale
  XScale=1.1; YScale=1.1; 
% Pause Durations when Plotting
  NSparedElements=2;
% The DoF to Draw the Plots for
  DoFtoDraw=10;
% Plot Arc-Length Results
  ArcLength=0;
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
                           
                           
% Elements Data Matrix:                          Element#   Start_Node   End_Node  A[in^2]   E[psi]     Fy[psi]   alpha1    ec/ey    ec/ey                               
                                  Elements=[
                                                    1,           1          5,       5,     36000000     25000     0.035       6       10;
                                                    2,           2          5,       4,     36000000     25000     0.035       6       10;
                                                    3,           3          5,       5,     36000000     25000     0.035       6       10;
                                                    4,           4          5,       6,     36000000     25000     0.035       6       10];
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
    alpha1=Elements(:,7);                                                                                   % Zarib Shibaye ba'd az taslime elemana, az setoone 7 matrice "Elements" rikhte mishe too bordare "alpha1". | alpha1 = Ye brodar ke derayeye i, shibe nemoodare tanesh-korneshe elemane i baad az noghteye taslim hast.
    VytoVc=Elements(:,8);                                                                                   % Zarib Tabdile "Vy" be "Vc" elemana, az setoone 8 matrice "Elements" rikhte mishe too bordare "VytoVc".    | VytoVc = Ye brodar ke derayeye i, Zarib Tabdile "Vy" be "Vc" baraye elemane i hast.
    VytoVf=Elements(:,9);                                                                                   % Zarib Tabdile "Vy" be "Vf" elemana, az setoone 9 matrice "Elements" rikhte mishe too bordare "VytoVf".    | VytoVf = Ye brodar ke derayeye i, Zarib Tabdile "Vy" be "Vf" baraye elemane i hast.
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
          Qc(i,1)=((VytoVc(i)*alpha1(i))+(1-alpha1(i)))*Qy(i);                               %%1            % Baraye har eleman, nirooye mehvarie motenazer ba noghteye capping mohasebe mishe.                           | Qy     = ye bordar ke derayeye i, nirooye cappe elemane i hast. Formulation, safheye 13 az file "Session 02.pdf".
          Vy(i,1)=Qy(i)/Kel(i,i);                                        
          Vc(i,1)=VytoVc(i)*Vy(i);                                                           %%2            % Ba estefade az zariba, taghiire toole mehvarie motenazer ba noghteye cappe har eleman mohasebe mishe.       | Vc     = ye bordar ke derayeye i, taghiire tooli has ke ba'ese capp kardane elemane i mishe. Too hal estefade nemishe va faghat baraye yadgiri hesab shode. Formulation, safheye 15 az file "Session 02.pdf".
          Vf(i,1)=VytoVc(i)*Vy(i);                                                           %%3            % Ba estefade az zariba, taghiire toole mehvarie motenazer ba noghteye faile har eleman mohasebe mishe.       | Vf     = ye bordar ke derayeye i, taghiire tooli has ke ba'ese fail kardane elemane i mishe. Too hal estefade nemishe va faghat baraye yadgiri hesab shode. Formulation, safheye 15 az file "Session 02.pdf".
          alpha2(i,1)=-((VytoVc(i)*alpha1(i))+(1-alpha1(i)))/(VytoVf(i)-VytoVc(i));          %%4            % Shibe khatte sevome nemoodare niroo-taghiire shekle har eleman (khatte ba'd az noghteye capp). hesab mishe. | alpha2 = ye bordar ke derayeye i, shibe khatte ba'd az capp too nemoodare niroo-taghiire shekle elemane i hast. Formulation, safheye 15 az file "Session 02.pdf".
      end 
      Kelt=Kel;
%         Qc=((VytoVc.*alpha1)+(1-alpha1)).*Qy;                                              %%1         
%         Vc=VytoVc.*Vy;                                                                     %%2            
%         Vf=VytoVc.*Vy;                                                                     %%3               
%         alpha2=((VytoVc.*alpha1)+(1-alpha1))./(VytoVf-VytoVc);                             %%4      
%% 03 - Calculations
 % 00 - Pre-Definitions
        f=0;                                                                                                % Ye shomarandas ke neshoon mide chand ta eleman fail kardan.
        i=0;                                                                                                % Ye shomarandas ke stepi ke toosh hastimo neshoon mide.
        D=zeros(NDoFs,1); 
        Yielded=zeros(NElements,1);                                                                         % Yielded = ye bordar ke be andazeye te'dade elemana deraye dare va hame sefran. Agar elemane i taslim she, derayeye i az in bordar 1 mishe. Agar inja ta'rif nashe too khatte 127 error migirim.                
        Capped=zeros(NElements,1);                                                                          % Capped  = ye bordar ke be andazeye te'dade elemana deraye dare va hame sefran. Agar elemane i varede khatte nozooli she, derayeye i az in bordar 1 mishe. Agar inja ta'rif nashe too khatte 127 error migirim. 
        Failed=zeros(NElements,1);                                                                          % Failed  = ye bordar ke be andazeye te'dade elemana deraye dare va hame sefran. Agar elemane i fail kone, derayeye i az in bordar 1 mishe. Agar inja ta'rif nashe too khatte 127 error migirim.   
        ResultNodalDisplacements=zeros(NDoFs,1);
        ResultMemberDeformations=zeros(NElements,1);
        ResultMemberAxialForces=zeros(NElements,1);                                              
        Events=zeros(1,4);                                                                                  % Events  = ye matrice ke too evente i (stepe i), rdife i por mishe, derayeye setoone 1 shomareye elemane oon evente, agar eleman taslim she derayeye setoone 2 barabare 1 mishe, agar eleman capp kone derayeye setoone 3 barabare 1 mishe, agar eleman fail she derayeye setoone 4 barabare 1 mishe. Agar inja ta'rif nashe too khatte 168 error migirim.
 % 01 - Solve
        tic;
        while f<NElements-NSparedElements                                                                   % Event-to-Event baraye eleman ba materiale 3khatti, ta vaghti too nemoodare Force-Displacement bargardan be nirooye 0 edame peyda mikone. Bayad be soorate dasti barresi she ba'd az fail kardane chan ta eleman, farayand motevaghghef she.
              i=i+1;                                                                                        % Varede stepe jadid ke mishim ye doone be i ezafe mikonim.                          
              % Assemble the Global Stifness Matrix
                Kt=zeros(NDoFs,NDoFs);                                                                                                      
                for j=1:NElements
                    Kt(DoFs(j,:),DoFs(j,:))=Kt(DoFs(j,:),DoFs(j,:))+(transpose(a(j,:))*Kelt(j,j)*a(j,:));
                end
              % Calculate the Nodal Displacements Vectors
                D(FreeDoFs)=Kt(FreeDoFs,FreeDoFs)\Loads(FreeDoFs); 
              % Calculate the Element Deformations Vectors
                V=DtoV*D;
              % Calculate the Elements Internal Forces Vector
                Q=Kelt*V;
              % Calculate deltaQ Vectors                                                                 %% 3 no' deltaQ mohasebe mikonim.         
                deltaQ1=zeros(NElements,1);                                                                 % DeltaQ1 = Ye bordare ke derayeye i, meghdare deltaQ baraye elemane i ke hanooz taslim nashode hast.         
                deltaQ2=zeros(NElements,1);                                                                 % DeltaQ2 = Ye bordare ke derayeye i, meghdare deltaQ baraye elemane i ke taslim shode amma hanooz capp nakarde hast.         
                deltaQ3=zeros(NElements,1);                                                                 % DeltaQ3 = Ye bordare ke derayeye i, meghdare deltaQ baraye elemane i ke capp karde amma hanooz fail nakarde hast.      
                for j=1:NElements  
                    if Yielded(j)==0 && Capped(j)==0 && Failed(j)==0                                        % Barresi mishe, faghat dar soorati ke too stepe aval bashim va ya agar elemane j taslim nashode bashe "DeltaQ1" barash mohasebe mishe, vagarna derayeye j az "deltaQ1" sefr mimoone.
                       if i==1
                           deltaQ1(j)=Qy(j);                                                                % Too stepe 1 "deltaQ1" baraye haeye elemana be ye soorat hsab mishe.
                       else
                           deltaQ1(j)=Qy(j)-abs(ResultMemberAxialForces(j,i));                              % Fasele ta "Qy" va taslim kardan mohasebe mishe.
                       end
                    elseif Yielded(j)==1 && Capped(j)==0 && Failed(j)==0                                    % Barresi mishe, faghat dar soorati ke elemane j tslim shode bashe amma hanooz capp nakarde bashe "DeltaQ2" barash mohasebe mishe, vagarna derayeye j az "deltaQ2" sefr mimoone.           
                           deltaQ2(j)=Qc(j)-abs(ResultMemberAxialForces(j,i));                              % Fasele ta "Qc" va capp kardan mohasebe mishe.
                    elseif Yielded(j)==1 && Capped(j)==1 && Failed(j)==0                                    % Barresi mishe, faghat dar soorati ke elemane j capp karde bashe amma hanooz fail nakarde bashe "DeltaQ3" barash mohasebe mishe, vagarna derayeye j az "deltaQ2" sefr mimoone.                 
                           deltaQ3(j)=abs(ResultMemberAxialForces(j,i));                                    % Fasele ta 0 shodane nirooye mehvari va fail kardan mohasebe mishe. 
                    end
                end
              % Calculate the Landa Vectors                                                              %% 3 no' deltaQ mohasebe mikonim.
                Lambda1(:,1)=(1:NElements)';  Lambda1(:,2)=10^24;                                           % Lambda1 = Ye matrice 2 setoone ke too radife i, derayeye setoone 1. shomareye elemane i ke hanooz taslim nashode hast va derayeye setoone 2, meghdare Lambdaye mohasebe shode baraye oon eleman hast. Fe'lam be hameye derayehaye setoone 2 ye adade kheili bozorg dadim. 
                Lambda2(:,1)=(1:NElements)';  Lambda2(:,2)=10^24;                                           % Lambda2 = Ye matrice 2 setoone ke too radife i, derayeye setoone 1. shomareye elemane i ke taslim shode amma hanooz capp nakarde hast va derayeye setoone 2, meghdare Lambdaye mohasebe shode baraye oon eleman hast. Fe'lam be hameye derayehaye setoone 2 ye adade kheili bozorg dadim.   
                Lambda3(:,1)=(1:NElements)';  Lambda3(:,2)=10^24;                                           % Lambda3 = Ye matrice 2 setoone ke too radife i, derayeye setoone 1. shomareye elemane i ke capp karde amma hanooz fail nakarde hast va derayeye setoone 2, meghdare Lambdaye mohasebe shode baraye oon eleman hast. Fe'lam be hameye derayehaye setoone 2 ye adade kheili bozorg dadim.  
                for j=1:NElements
                    if deltaQ1(j,1)~=0                                                                      % Braye har eleman, agar "deltaQ1" adad dashte bashe va sefr nabashe (eleman hanooz taslim nashode bashe), "Lambda1" barash mohasebe mishe. Vagarna oon adade bozorg sare jash mimoone.
                       Lambda1(j,2)=abs(deltaQ1(j,1)/Q(j)); 
                    end
                    if deltaQ2(j,1)~=0                                                                      % Braye har eleman, agar "deltaQ2" adad dashte bashe va sefr nabashe (eleman tslim shode bashe amma hanooz capp nakarde bashe), "Lambda2" barash mohasebe mishe. Vagarna oon adade bozorg sare jash mimoone.      
                       Lambda2(j,2)=abs(deltaQ2(j,1)/Q(j)); 
                    end
                    if deltaQ3(j,1)~=0                                                                      % Braye har eleman, agar "deltaQ3" adad dashte bashe va sefr nabashe (capp karde bashe amma hanooz fail nakarde bashe), "Lambda3" barash mohasebe mishe. Vagarna oon adade bozorg sare jash mimoone.
                       Lambda3(j,2)=abs(deltaQ3(j,1)/Q(j));  
                    end  
                end
              % check if Any Element Was Capped During Previous Event and is Failing in Current Step      
                fs=0;                                                                                       % Ye shomarande ta'rif mikonim ke be tore default meghdaresh 0 hast. Amma agar too evente ghabl elemani capp karde va shoroo karde be fail kardan, in shomarande 1 beshe.
                if i>1 && Events(i-1,3)==1                                                               %% Evente ghabl barresi mishe, agar elemani too evente ghabl capp karde bashe va dar hale fail shodan bashe, too evente jadidemoon dige be Lambdaye elemanaye dige kari nadarim. Lambdaye hamoon elemano estefade mikonim ke kamel fail kone.
                   fs=1;                                                                                    % Shomarande ro 1 mikonim.
                   Failed(Events(i-1,1))=1;                                                                 % Too bordare "Failed" too derayeye marboot be oon eleman ye 1 mizarim ke moshakhas she fail karde.
                   Kelt(Lambda3(1,1),Lambda3(1,1))=0;                                                       % Sakhtie locale elemano sefr mikonim.
                   minLambda(i)=(-1)*Lambda3(Events(i-1,1),2);                                              % Lambdaye marboot be oon elemano az matrice "Lambda3" ke marboot be elemanaye dar hale fail hast zakhire mikonim.
                   Events(i,1)=Events(i-1,1);                                                               % Too matrice "Events" shomare elemano too redife stepemoon va setoone 1 minevisim.
                   Events(i,4)=1;                                                                           % Too matrice "Events" too redife stepemoon va setoone 4 ye 1 mizarim ke bedoonim in eleman fail karde.
                   f=f+1;                                                                                   % Be shomarandeye elemanaye fail shode ye doone ezaf mikonim.
                end
              % Identify the Most-Probable-to-Yield/Capp Element
                Lambda1=sortrows(Lambda1,2);                                                                % Redifaye matrice "Lambda1" ro ba tavajoh be setoone 2 (minimum be maximum) sort mikonim. Adadye bozorgi ke ghablan dade boodim, inja ba'es mishe elemanaye namarboot beran kenar.
                Lambda2=sortrows(Lambda2,2);                                                                % Redifaye matrice "Lambda2" ro ba tavajoh be setoone 2 (minimum be maximum) sort mikonim. Adadye bozorgi ke ghablan dade boodim, inja ba'es mishe elemanaye namarboot beran kenar.
                Lambda3=sortrows(Lambda3,2);                                                                % Redifaye matrice "Lambda3" ro ba tavajoh be setoone 2 (minimum be maximum) sort mikonim. Adadye bozorgi ke ghablan dade boodim, inja ba'es mishe elemanaye namarboot beran kenar.
                if Lambda1(1,2)<Lambda2(1,2) && Lambda1(1,2)<Lambda3(1,2) && fs==0                       %% Agar derayeye 2 az radife 1 matrice "Lambda1" az hameye lambdahaye dige kamtar bashe va too evente ghabl elemani capp nakarde bashe (dar hale fail nabashe), too evente jadidemoon elemane derayeye 1 radife 1 az "Lambda1" taslim mishe.
                   Yielded(Lambda1(1,1))=1;                                                                 % Too bordare "Yielded" too derayeye marboot be oon eleman ye 1 mizarim ke moshakhas she taslim shode.           
                   Kelt(Lambda1(1,1),Lambda1(1,1))=alpha1(Lambda2(1,1))*Kel(Lambda1(1,1),Lambda1(1,1));     % Sakhtie locale elemano too "alpha1" zarb mikonim.
                   minLambda(i)=Lambda1(1,2);                                                               % Lambdaye marboot be oon elemano az matrice "Lambda1" ke marboot be elemanaye taslim nashodeye dar hale taslim hast, zakhire mikonim.          
                   Events(i,1)=Lambda1(1,1);                                                                % Too matrice "Events" shomare elemano too redife stepemoon va setoone 1 minevisim.                  
                   Events(i,2)=1;                                                                           % Too matrice "Events" too redife stepemoon va setoone 2 ye 1 mizarim ke bedoonim in eleman taslim shode.
                end
                if Lambda2(1,2)<Lambda1(1,2) && Lambda2(1,2)<Lambda3(1,2) && fs==0                       %% Agar derayeye 2 az radife 1 matrice "Lambda2" az hameye lambdahaye dige kamtar bashe va too evente ghabl elemani capp nakarde bashe (dar hale fail nabashe), too evente jadidemoon elemane derayeye 1 radife 1 az "Lambda2" capp mikone.
                   Capped(Lambda2(1,1))=1;                                                                  % Too bordare "Capped" too derayeye marboot be oon eleman ye 1 mizarim ke moshakhas she capp karde. 
                   Kelt(Lambda2(1,1),Lambda2(1,1))=alpha2(Lambda2(1,1))*Kel(Lambda2(1,1),Lambda2(1,1));     % Sakhtie locale elemano too "alpha2" zarb mikonim.  
                   minLambda(i)=Lambda2(1,2);                                                               % Lambdaye marboot be oon elemano az matrice "Lambda1" ke marboot be elemanaye taslim shodeye capp nakardeye dar hale capp hast, zakhire mikonim.                                                                         
                   Events(i,1)=Lambda2(1,1);                                                                % Too matrice "Events" shomare elemano too redife stepemoon va setoone 1 minevisim.  
                   Events(i,3)=1;                                                                           % Too matrice "Events" too redife stepemoon va setoone 3 ye 1 mizarim ke bedoonim in eleman capp karde.
                end 
                clear Lambda1; clear Lambda2; clear Lambda3;                                                % Niazi nist hatman Lambdaha pak shan. Baraye tamizie kar pak shodan.          
              % Collect the Results   
                ResultEvents(i,:)=Events(i,:);
                ResultLambda(i)=minLambda(i);
                ResultNodalForces(:,i+1)=sum(ResultLambda(:))*Loads;
                if i==1
                   ResultNodalDisplacements(:,i+1)=abs(ResultLambda(i)*D); 
                   ResultMemberDeformations(:,i+1)=ResultLambda(i)*V; 
                   ResultMemberAxialForces(:,i+1)=(ResultLambda(i)*Q);
                else
                   ResultNodalDisplacements(:,i+1)=ResultNodalDisplacements(:,i)+abs(ResultLambda(i)*D); 
                   ResultMemberDeformations(:,i+1)=ResultMemberDeformations(:,i)+ResultLambda(i)*V;
                   ResultMemberAxialForces(:,i+1)=ResultMemberAxialForces(:,i)+(ResultLambda(i)*Q);
                end
        end
        Duration=toc;                     
%% 04 - Visualization
close all;
ResultNodalDisplacements=abs(ResultNodalDisplacements);
ResultNodalForces=abs(ResultNodalForces);
ResultMemberDeformations=abs(ResultMemberDeformations);
ResultMemberAxialForces=abs(ResultMemberAxialForces);
% load('Example2Results_B_3Khati_Arc_Length')
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
         ylim([0, YScale*max(abs(ResultNodalForces(DoFtoDraw,:)))]); 
         xlim([0, XScale*ResultNodalDisplacements(DoFtoDraw,size(ResultNodalDisplacements,2))]);
         grid on; grid minor; ax=gca; ax.GridLineStyle='--'; ax.GridAlpha=0.6; ax.GridColor=['k']; ax.FontSize=12; ax.LineWidth=0.8; ax.TickLength=[0.01 0.01];    
         xlabel('Displacement [in]','fontsize',13,'fontweight','bold'); ylabel('Load [lb]','fontsize',13,'fontweight','bold');      
         Title=['DoF ' num2str(DoFtoDraw) ' Load-Dispacement Plot']; title(Title,'fontsize',15,'fontweight','bold'); hold on;
       % Plotting - Arc-Length Results
         if ArcLength==1
            P1=plot(ArcLengthDisplacements(DoFtoDraw,:),ArcLengthForces(DoFtoDraw,:),'b','LineWidth',3);
            if DrawAllFreeDoFsPlots==0 && NoAnimation==0
               pause(Pause2);
            end
         end
       % Plotting - Event-to-Event Results
         if DrawAllFreeDoFsPlots==0 && NoAnimation==0 && ArcLength==0
            pause(Pause2);
         end
         P2=scatter(ResultNodalDisplacements(DoFtoDraw,1),ResultNodalForces(DoFtoDraw,1),'m','filled'); hold on;
         if DrawAllFreeDoFsPlots==0 && NoAnimation==0 && ArcLength==0
            pause(0.5*Pause2);
         end
         for j=2:size(ResultNodalForces,2)
             plot(ResultNodalDisplacements(DoFtoDraw,j-1:j),ResultNodalForces(DoFtoDraw,j-1:j),'g','LineWidth',2); hold on;
             scatter(ResultNodalDisplacements(DoFtoDraw,j),ResultNodalForces(DoFtoDraw,j),'m','filled'); hold on;
             if Events(j-1,2)==1 && Events(j-1,3)==0 && Events(j-1,4)==0
                TXT=['Element ' num2str(Events(j-1,1)) ' yields here.'];
             end
             if Events(j-1,2)==0 && Events(j-1,3)==1 && Events(j-1,4)==0
                TXT=['Element ' num2str(Events(j-1,1)) ' Capps here.'];
             end                
             if Events(j-1,2)==0 && Events(j-1,3)==0 && Events(j-1,4)==1
                TXT=['Element ' num2str(Events(j-1,1)) ' Fails here.'];
             end 
             if j==size(ResultNodalForces,2)
                text(ResultNodalDisplacements(DoFtoDraw,j)-0.0055,ResultNodalForces(DoFtoDraw,j)+2000,TXT,'HorizontalAlignment','left','VerticalAlignment','bottom','Color','Magenta','FontSize',13,'FontWeight','bold');
             elseif  j~=size(ResultNodalForces,2) && Events(j-1,4)==1
                text(ResultNodalDisplacements(DoFtoDraw,j)-0.0035,ResultNodalForces(DoFtoDraw,j)-1500,TXT,'HorizontalAlignment','left','VerticalAlignment','top','Color','Magenta','FontSize',13,'FontWeight','bold');
             else 
                text(ResultNodalDisplacements(DoFtoDraw,j)+0.0006,ResultNodalForces(DoFtoDraw,j)-1500,TXT,'HorizontalAlignment','left','VerticalAlignment','top','Color','Magenta','FontSize',13,'FontWeight','bold');
             end
             clear TXT;
             if j~=size(ResultNodalForces,2) && DrawAllFreeDoFsPlots==0 && NoAnimation==0 && ArcLength==0
                pause(0.5*Pause2);
             end
             if ArcLength==1 && j==size(ResultNodalForces,2)
                scatter(ResultNodalDisplacements(DoFtoDraw,:),ResultNodalForces(DoFtoDraw,:),'m','filled'); hold on;
                pause(0.5*Pause1);
                leg=legend([P1, P2],{'Arc-Length Results', 'Event-to-Event Results'},'Location','northeast','fontsize',15);
             elseif ArcLength~=1 && j==size(ResultNodalForces,2)
                 if NoAnimation==1
                scatter(ResultNodalDisplacements(DoFtoDraw,:),ResultNodalForces(DoFtoDraw,:),'m','filled'); hold on;
                 end
                 pause(0.5*Pause1);
                leg=legend([P2],{'Events'},'Location','northeast','fontsize',15); 
             end

         end

   end
Plot_Name=['Trilinear Material - Event-to-Event Method - ' num2str(size(ResultNodalForces,2)-1) ' Steps - Duration ' num2str(Duration) ' Seconds']; hold on;
h1=suptitle(Plot_Name); set(h1,'FontSize',17,'FontWeight','bold'); hold off;