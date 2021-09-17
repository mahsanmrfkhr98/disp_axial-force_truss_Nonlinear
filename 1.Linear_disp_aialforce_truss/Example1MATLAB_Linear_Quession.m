clc;clear;close all
%% units: kN, GPa, mm
%% 01 - Input
% Nodes Data Matrix:          Node#       x[mm]         y[mm]
                Nodes=[
                               1,           0,             0
                               2,        5000,             0
                               3,        5000,         12000    
                               4,           0,         24000];
% Elements Data Matrix:     Element#   Start_Node     End_Node  A[mm^2]  E[GPa]
             Elements=[
                               1,           1             2,     4000,     58;
                               2,           2             3,     4000,     58;
                               3,           3             4,     4000,     58;
                               4,           1             3,     4000,     58;
                               5,           2             4,     4000,     58];
% Restrains Vector:           DoFs
            Restrains=[
                               1;
                               2;
                               4;
                               7;
                               8];
% External Loads Vector:  Load_Value[kN)   DoF#  
                Loads=[
                               0           %1
                               0           %2
                             680           %3
                               0           %4
                            -560           %5
                             980           %6
                               0           %7
                               0];         %8
%% 02 - Complie Inputs
% Nodes Data                                                                                               
  % Available Data                                                                                   %%% Etela'ate vared shode dar matrice "Nodes" inja be bordarhaye koochik tar shekaste mishe.
    x=Nodes(:,2);                                                                                           % Mokhtasate x noghat, az setoone 2 matrice "Nodes" rikhte mishe too bordare "x". | x = Ye brodar ke derayeye i, mokhtasate x node i hast.
    y=Nodes(:,3);                                                                                           % Mokhtasate y noghat, az setoone 3 matrice "Nodes" rikhte mishe too bordare "y". | y = Ye brodar ke derayeye i, mokhtasate y node i hast.
  % New Data                                                                                         %%% Ye seri etela'at mostaghiman be matricaye voroodi dade nashode.           
    % Number of Nodes                                                                                    %% Az size matrice "Nodes" ye seri etelaat estekhraj mishe.     
      NNodes=size(Nodes,1);                                                                                 % Ba tavajoh be te'dade radifaye matrice "Nodes", te'dade kolle nodha (NNodes) moshakhas mishe.
      NDoFs=2*NNodes;                                                                                       % Te'dade kolle darajate azadi (NDoFs) = 2 barabare Te'dade kolle noda (NNodes).
    % Free DoFs                                                                                          %% Ma ye bordar darim ke toosh dade shode kodoom darajate azadi baste shodan, mikhaim ye bordar dashte bashim ke toosh faghat shomare darajate azadi azad neveshte shode.
      FreeDoFs=(1:NDoFs)';                                                                                  % Ye bordar sakhtam ke derayehash (har radif) az 1 shoroo mishe va to teedade kolle drajate azadi (NDoFs) edame dare.
      FreeDoFs(Restrains)=[];                                                                               % Az bordari ke shakhtam, derayehaye marboot be darajate azadi baste (Restrains) ro hazf mikonam.
% Elements Data
  % Available Data                                                                                   %%% Etela'ate vared shode dar matrice "Elements" inja be bordarhaye koochik tar shekaste mishe.    
    Start=Elements(:,2);                                                                                    % Node ebtedaye elemana, az setoone 2 matrice "Elements" rikhte mishe too bordare "Start". | Start = Ye brodar ke derayeye i, node ebtedaye elemane i hast.
    End=Elements(:,3);                                                                                      % Node entehaye elemana, az setoone 3 matrice "Elements" rikhte mishe too bordare "End".   | End   = Ye brodar ke derayeye i, node entehaye elemane i hast.  
    A=Elements(:,4);                                                                                        % Sathe maghta'e elemana, az setoone 4 matrice "Elements" rikhte mishe too bordare "A".    | A     = Ye brodar ke derayeye i, sathe maghta'e elemane i hast.  
    E=Elements(:,5);                                                                                        % Module elasticite elemana, az setoone 5 matrice "Elements" rikhte mishe too bordare "E". | E     = Ye brodar ke derayeye i, module elasticiteye elemane i hast.  
  % New Data                                                                                         %%% Ye seri etela'at mostaghiman be matricaye voroodi dade nashode. 
    % Number of Elements                                                                                 %% Az size matrice "Elements" ye seri etelaat estekhraj mishe.     
      NElements=size(Elements,1);                                                                           % Ba tavajoh be te'dade radifaye matrice "Elements", te'dade kolle elemana (NElements) moshakhas mishe.
      % Geometric and Mechanical Info                                                                       % Inja az moshakhasate noda va elemana ba ham estefade mishe                        
      K=zeros(NDoFs,NDoFs);                                                                                 % Midoonim matrice sakhtie koll (K) ye matrice moraba'ii hast ke te'dade satra va setoonash barabare teedade kolle darajate azadie.
%     DoFs=[2*Start-1,  2*Start,  2*End-ones(NElements,1),  2*End];                          %%1  
%     invL=1./((x(End)-x(Start)).^2+(y(End)-y(Start)).^2).^0.5;                              %%2                        
%     C=(x(End)-x(Start)).*invL;                                                             %%3
%     S=(y(End)-y(Start)).*invL;                                                             %%4
%     a=[-C, -S, C, S];                                                                      %%5
%     Kel=A.*E.*invL;                                                                        %%6         
      for i=1:NElements                                                                                  %% Tamame Bordara va matricaye morede niaz baraye sakhte matrice sakhti, too in loop baraye har eleman hesab mishan. 
          DoFs(i,:)=[2*Start(i)-1  2*Start(i)  2*End(i)-1  2*End(i)];                        %%1            % Ba tavajoh be shomareye nodaye sar o tahe har eleman, darajate azadi sar va tahe eleman peida mishe.         | DoFs = Ye matrice ke radife i, bordare darajate azadie sar o tahe elemane i hast.
          invL(i)=1/((x(End(i))-x(Start(i)))^2+(y(End(i))-y(Start(i)))^2)^0.5;               %%2            % Az mokhtasate nodhaye sar o tahe har eleman, toole eleman va 1 rooye toole eleman (invL) hesab mishe.        | invL = Ye brodar ke derayeye i, 1 rooye toole elemane i hast.
          C(i)=(x(End(i))-x(Start(i)))*invL(i);                                              %%3            % Az mokhtasate nodhaye sar o tahe har eleman va invL, cosinuse zavieye eleman ba khate ofogh (C) hesab mishe. | C    = Ye brodar ke derayeye i, cosinuse zavieye elemane i ba khate ofogh hast.
          S(i)=(y(End(i))-y(Start(i)))*invL(i);                                              %%4            % Az mokhtasate nodhaye sar o tahe har eleman va invL, sinuse zavieye eleman ba khate ofogh (S) hesab mishe.   | S    = Ye brodar ke derayeye i, sinuse zavieye elemane i ba khate ofogh hast.
          a(i,:)=[-C(i) -S(i) C(i) S(i)];                                                    %%5            % Ba sinusa va cosinusaye hesab shode bara har eleman, bordare "a" tashkil mishe.                              | a    = Ye matrice ke radife i, bordare "a" marboot be elemane i hast. Formulation, safheye 2 va 4 az file "Session 01.pdf".
          Kel(i)=A(i)*E(i)*invL(i);                                                          %%6            % Az rabeteye [AE/L ya A*E*invL], sakhtie locale har eleman (Kel) be dast miad.                                | Kel  = Ye bordar ke derayeye i, sakhtie locale elemane i hast. Formulation, safheye 3 az file "Session 01.pdf".
          K(DoFs(i,:),DoFs(i,:))=K(DoFs(i,:),DoFs(i,:))+(transpose(a(i,:))*Kel(i)*a(i,:));   %%7            % Ba tavajoh be darajate azadie sar o tahe har eleman, matrice sakhtie kolle saze (K) assemble mishe.          | K    = Matrice sakhtie kolle saze. Formulation, safheye 5 va 6 az file "Session 01.pdf".                                                
%         K(DoFs(i,:),DoFs(i,:))=K(DoFs(i,:),DoFs(i,:))+((a(i,:)')*Kel(i)*a(i,:));           %%7
      end
%% 03 - Calculate the Nodal Displacements 
  Pext=Loads;                                                                                               % In bordaro ta'rif kardim ke tebghe notatione kelas pish berim.                                                         | Pext = ye bordar ke derayehaye 2i-1 va 2i, niroohaye kharejie rooye node i dar jahate x va y hastan.       
  D=zeros(NDoFs,1);                                                                                         % Midoonim bordare kolle jabeja'iia ye bordare ke sizesh (te'dade radifash) barabare ba teedade kolle darajate azadi.
  D(FreeDoFs)=K(FreeDoFs,FreeDoFs)\Pext(FreeDoFs);                                           %%8            % Jabeja'ii faghat baraye dajate azadie baz mohasebe mishe.                                                              | D    = ye bordar ke derayehaye 2i-1 va 2i, jabeja'iihaye node i dar jahate x va y hastan. Formulation, safheye 6 az file "Session 01.pdf".
% D(FreeDoFs)=inv(K(FreeDoFs,FreeDoFs))*Pext(FreeDoFs);                                      %%8
%% 04 - Calculate Members Deformation
for i=1:NElements
    V(i)=a(i,:)*D(DoFs(i,:));                                                                               % Ba zarbe "a" har eleman dar jabeja'iie darajate azadie sar o tahe eleman, taghiire toole mehvarie eleman be dast miad. | V    = ye bordar ke derayeye i, taghiire toole mehvarie elemane i hast. Formulation, safheye 2 az file "Session 01.pdf".
end
%% 05 - Calculate Members Axial Forces
for i=1:NElements
    Q(i)=Kel(i)*V(i);                                                                        %%9            % Az zarbe sakhtie locale har eleman dar taghiire toole mehvari, nirooye mehvarie ijad shode dar eleman be dast miad.    | Q    = ye bordar ke derayeye i, nirooye mehvarie ijad shode dar elemane i hast. Formulation, safheye 3 az file "Session 01.pdf".
end
%   Q=Kel.*V;                                                                                %%9
%% 06 - Calculate the Internal Nodal Loads Vector
Pint=K*D;                                                                                                   % Az zarbe matrice sakhti dar bordare jabejaiia, bordare niroo shamele tamame niroohaye dakheli va khareji be dast miad. | Pint = ye bordar ke derayehaye 2i-1 va 2i, niroohaye rooye node i dar jahate x va y hastan (niroohaye khareji va niroohaye tekyegah). Formulation, safheye 6 az file "Session 01.pdf".
%% 07 - Collect the Results
ResultNodalDisplacements=D;
ResultMemberDeformations=V;
ResultMemberAxialForces=Q;
ResultNodalForces=Pint;
%% 08 - Visualization
i=1:NNodes;
T1=table(i',ResultNodalDisplacements(1:2:2*NNodes),ResultNodalDisplacements(2:2:2*NNodes),ResultNodalForces(1:2:2*NNodes),ResultNodalForces(2:2:2*NNodes),'VariableNames',{'Node_Number' 'X_Displacement_mm' 'Y_Displacement_mm' 'X_Force_kN' 'Y_Force_kN'})
j=1:NElements;
T2=table(j',ResultMemberDeformations',ResultMemberAxialForces','VariableNames',{'Member_Number' 'Axial_Deformation_mm' 'Axial_Force_kN'})