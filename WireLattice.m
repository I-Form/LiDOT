function [N,S,ElVol,M] = WireLattice(varargin)
% [NODES,STRUTS,ElVol, MAXWELL NUMBER] = WireLattice(TYPE,DIM,SIZE,SHAPE, name-value pairs)
% Creates a lattice of desired type, number of unit cells, and size.
% -TYPE -  one of  ["BCC","BCCZ","FCC","FCCZ","OCTET","Dual","Edge",
%       "tetDual","tetEdge"]
% -DIM - unit cell count in [x y z] directions - 
%       Cylinder [radial axial Inner_diameter] - must be 1x3 array of integers
% -SIZE - size of a unit cell in [x y z] or 
%       [radial hoop(approx) axial] direction - 1x3 array real values
% -SHAPE - Base mesh - Box or cylinder - SquareSymetrical for symmetry FEA models -
%       ["Square","SquareSymetrical","Cylinder","InputMesh"]
% Name - value pairs:
% -"Plot" - 1 visualise lattice, 0 no visualisation
% -"hex2tetType" - int 1-6 how hexs are split to tets for dual and edge lattice in
%       hex2tet function- 
% -"ZDir" - ["Radial","Hoop","Axial"] which direction has Z struts for BCCZ and FCCZ in cyclinders for BCCZ and FCCZ
% 
% 

%% Parse Inputs
% Defaults
lattice_types = ["BCC","BCCZ","FCC","FCCZ","Octet","Dual","Edge","tetDual","tetEdge"];
shape_types = ["Square","Cylinder","SquareSymetrical"];
Z_dirs = ["Radial","Hoop","Axial"];

p=inputParser;
addOptional(p,'TYPE',"BCC", @(x) isstring(validatestring(x,lattice_types))...
                                || ischar(validatestring(x,lattice_types)))
addOptional(p,'DIM',[1 1 1], @(x) nnz(size(x) == [1,3])==2)    %&& (nnz(x == round(x)) == 3)
addOptional(p,'SIZE',[1 1 1], @(x) nnz(size(x) == [1,3])==2)
addOptional(p,'SHAPE',"Square", @(x) isstring(validatestring(x,shape_types))...
                                || ischar(validatestring(x,shape_types)))
addParameter(p,'ZDir',"Axial", @(x) isstring(validatestring(x,Z_dirs))...
                                || ischar(validatestring(x,Z_dirs)))
addParameter(p,'Plot',0);
addParameter(p,'hex2tetType',1);
addParameter(p,'InputMesh',[]);

parse(p,varargin{:})

TYPE = validatestring(p.Results.TYPE,lattice_types);
DIM = p.Results.DIM;
SIZE = p.Results.SIZE;
SHAPE = validatestring(p.Results.SHAPE,shape_types);
PlotOn =p.Results.Plot;
Z_dir = p.Results.ZDir;
hex2tetType = p.Results.hex2tetType;
meshStruct = p.Results.InputMesh;

if isempty(meshStruct) == 0
    SHAPE = "Input";
end

%% Create Base Mesh
switch SHAPE
    case {'Square',"SquareSymetrical"}
        % Create Cube with HEX Mesh using hexMeshBox
        cubeDimensions=DIM.*SIZE; %Dimensions
        [meshStruct]=hexMeshBox(cubeDimensions,DIM,2);
        V1=meshStruct.nodes; %The nodes (vertices)
        E1=meshStruct.elements;  %the elements
        F1=meshStruct.faces;
    
    case {'Cylinder'} 
        % Settings
        numRings=DIM(1);     %number of cells radially
        CellsH=DIM(2)+1;       %Number of Cells in Height +1
        di=DIM(3);          %Inner Diameter
        do=di + DIM(1)*SIZE(1)*2;   %Outer Diameter
        Height=DIM(2)*SIZE(3);       %Height
        
        Ds=linspace(di,do,numRings+1);  %Diameters at each ring
        CellsCirc=round((pi()*mean([di do]))/SIZE(2));        %Number of cells along circumfrence
       
        %Initilaise Cell Arrays
        F=cell(numRings,1); V=cell(numRings,1);
        %Create first cylinder - inside surF1ce
        [F{1},V{1}]=patchcylinder(Ds(1)/2,CellsCirc,Height,CellsH,'quad');
        E1=zeros(numRings*size(F{1},1),8);  %Initialise Element Array
    
        for i=1:numRings   
            %Creates Next cylinder - one step to the outside
            [F{i+1},V{i+1}]=patchcylinder(Ds(i+1)/2,CellsCirc,Height,CellsH,'quad');
            %Correct indices
            F{i+1}=F{i+1} + i*size(V{1},1);
            %Create Element array based on Faces
            switch Z_dir
                case {'Radial'}
                %Radial Z struts
                E1( ((i-1)*size(F{1},1)+1) : (i)*size(F{1},1) ,:)=[F{i} F{i+1}];
                case {'Hoop'}
                %Hoop Z struts
                E1( ((i-1)*size(F{1},1)+1) : (i)*size(F{1},1) ,:)=...
                    [F{i}(:,1) F{i+1}(:,1) F{i+1}(:,4) F{i}(:,4) ...
                    F{i}(:,2) F{i+1}(:,2) F{i+1}(:,3) F{i}(:,3)];
                case {'Axial'}
                %Axial Z struts
                E1( ((i-1)*size(F{1},1)+1) : (i)*size(F{1},1) ,:)=...
                    [F{i}(:,1:2) F{i+1}(:,2) F{i+1}(:,1) ...
                    F{i}(:,4) F{i}(:,3) F{i+1}(:,3:4)];
            end
            
        end
        %Combine all nodes
        V1 = vertcat(V{:});
        %Create final patch data
        [F1]=element2patch(E1,[],'hex8');  
    case {'Input'}
        E1 = meshStruct.elements;
        V1  =meshStruct.nodes;
        [F1]=element2patch(E1);       

end

if TYPE == "tetDual" || TYPE =="tetEdge"
    [E1,V1,~]=hex2tet(E1,V1,[],hex2tetType);     %convert hex to tet elements
    [F1]=element2patch(E1);            %Face data
end

%% Create Lattice

switch TYPE
    case {'BCC','BCCZ'}
        %Get Midpoints of all unit cells
%         Vb=mean(pagetranspose(reshape(transpose(V1(E1',:)),3,8,[])));
%         Vb=reshape(Vb,3,size(E1,1))';
        Vb=zeros(size(E1,1),3);
        for k = 1:size(E1,1)    
            Vb(k,:)=mean(V1(E1(k,:),:)) ;   
        end
        
        N=[V1;Vb]; %Combine arrays to get =lattice nodes
        V2Ind = (1:length(Vb)) + length(V1);
        
        S=zeros(8*size(E1,1),2);
     
        for k = 1:size(E1,1)   %create BCC struts
            RowIndex=((1+8*(k-1)):8*k);
            S( RowIndex, [1 2])=[ V2Ind(k) E1(k,1);...
                                   V2Ind(k) E1(k,2);...
                                   V2Ind(k) E1(k,3);...
                                   V2Ind(k) E1(k,4);...
                                   V2Ind(k) E1(k,5);...
                                   V2Ind(k) E1(k,6);...
                                   V2Ind(k) E1(k,7);...
                                   V2Ind(k) E1(k,8)];
        end
        
        
        if TYPE=="BCCZ"
            SV=zeros(4*size(E1,1),2);
            for k = 1:size(E1,1)   %create BCC struts
                RowIndex=(1)+4*(k-1):(4*k);
                 SV( RowIndex , [1 2])=[ E1(k,1) E1(k,5);...
                                        E1(k,2) E1(k,6);...
                                        E1(k,3) E1(k,7);...
                                        E1(k,4) E1(k,8)];
               
            end
            
            S=[S;SV];
        end
    
        
    case {'FCC','FCCZ'}
       
        S=zeros(8*size(E1,1),2);
        FMid=zeros(size(F1,1),3);
        % Get Diagonals for each element/unit cell
        % get rid of Top bottom F1ces
        %F1=F1((2*size(E1,1)+1):end,:);  %Removes XY F1ces
        %F1=F1([1:(2*size(E1,1)) ((4*size(E1,1))+1):6*size(E1,1)],:);  %Removes ZX F1ces
        F1=F1([((2*size(E1,1))+1):(6*size(E1,1))],:);  %Removes XY F1ces
        
        for k = 1:size(F1,1)    
            FMid(k,:)=mean(V1(F1(k,:),:)) ;   
        end
        
       [FMidu,~,ic]=unique(FMid,'stable','rows');
        MidVs=length(V1)+ic;
         
        for k = 1:size(F1,1)   %get struts on all F1ces
            RowIndex=(1+4*(k-1)):4*k;
            S(RowIndex,[1 2])= [ MidVs(k) F1(k,1) ;...   
                                 MidVs(k) F1(k,2) ;...
                                 MidVs(k) F1(k,3) ;...
                                 MidVs(k) F1(k,4) ];
        end
        
        if TYPE=="FCCZ"
            SV=zeros(4*size(E1,1),2);
            for k = 1:size(E1,1)   %Get vertical struts
                 RowIndex=1+4*(k-1):(4*k);
                 SV(RowIndex, [1 2])=[  E1(k,1) E1(k,5);...
                                        E1(k,2) E1(k,6);...
                                        E1(k,3) E1(k,7);...
                                        E1(k,4) E1(k,8)];
                
            end
        end
    
        % Ammend Strut and Node Arrays
        N = [V1; FMidu];                     %Join node sets
            
        if TYPE=="FCCZ"   
            S=[S ; SV];  %Strut Array

        end
    case 'Octet' 
        %Creates all struts for each Octet Truss unit cell
        S=zeros(12*size(E1,1),2);
        FMid=zeros(size(F1,1),3);
        
        for k = 1:size(F1,1)    
            FMid(k,:)=mean(V1(F1(k,:),:)) ;   
        end
        
        [FMidu,~,ic]=unique(FMid,'stable','rows');
        MidVs=length(V1)+ic;
         
        for k = 1:size(F1,1)    %get struts on all F1ces
            RowIndex=(1+4*(k-1)):4*k;
            S(RowIndex,[1 2])= [ MidVs(k) F1(k,1) ;...   
                                 MidVs(k) F1(k,2) ;...
                                 MidVs(k) F1(k,3) ;...
                                 MidVs(k) F1(k,4) ];
        end
     
         %get struts in octahedron core using the F1ce midpoints
         SOctet=zeros(size(E1,1)*12,2);   %initialize array
         for k=1:size(E1,1)
             RowIndex=(1+12*(k-1)):12*k;
             SOctet( RowIndex, [1 2])=[...
                  MidVs(k+(1-1)*size(E1,1)) MidVs(k+(3-1)*size(E1,1));...
                  MidVs(k+(1-1)*size(E1,1)) MidVs(k+(4-1)*size(E1,1));...
                  MidVs(k+(2-1)*size(E1,1)) MidVs(k+(3-1)*size(E1,1));...
                  MidVs(k+(2-1)*size(E1,1)) MidVs(k+(4-1)*size(E1,1));...

                  MidVs(k+(3-1)*size(E1,1)) MidVs(k+(5-1)*size(E1,1));...
                  MidVs(k+(3-1)*size(E1,1)) MidVs(k+(6-1)*size(E1,1));...
                  MidVs(k+(4-1)*size(E1,1)) MidVs(k+(5-1)*size(E1,1));...
                  MidVs(k+(4-1)*size(E1,1)) MidVs(k+(6-1)*size(E1,1));...

                  MidVs(k+(5-1)*size(E1,1)) MidVs(k+(1-1)*size(E1,1));...
                  MidVs(k+(5-1)*size(E1,1)) MidVs(k+(2-1)*size(E1,1));...
                  MidVs(k+(6-1)*size(E1,1)) MidVs(k+(1-1)*size(E1,1));...
                  MidVs(k+(6-1)*size(E1,1)) MidVs(k+(2-1)*size(E1,1))];
         end

        % Ammend Strut and Node Arrays
        N = [V1; FMidu];                       %Join node sets
        S = [S;SOctet]  ;

    case {'Dual','tetDual'}   %based on Gibbon Duallattice function - connects element centroids
       FsPerE=size(F1,1)/size(E1,1);      %Number of Faces per element

       X=V1(:,1); Y=V1(:,2); Z=V1(:,3);

       Vm=[mean(X(E1),2) mean(Y(E1),2) mean(Z(E1),2)]; %mid element points
       
       Vf=[mean(X(F1),2) mean(Y(F1),2) mean(Z(F1),2)];  %get Face mid-points (from hex2tet)

       [VfUnique, ~ ,ic]=unique(Vf,'rows');  %Gets unique Face mid-points

       N=[Vm ; VfUnique];               %Centroids and F1ce mid-points
       VInd=[(1:length(Vm))' ; ic];     %Indices of all centroids F1ce points mapped to new ones 
        
       %Row is an elements, columns Face-mid point for each
       E2FInd=size(Vm,1) + reshape(ic,size(E1,1),FsPerE);
        
       S0n=zeros(size(ic,1),2);
        for k = 1:size(E1,1)     %Creates struts from centroid and F1cepoints
            RowIndex=(1+FsPerE*(k-1)):FsPerE*k;
            if FsPerE == 4
                S0n(RowIndex,1:2) = [k E2FInd(k,1)
                                    k E2FInd(k,2)
                                    k E2FInd(k,3)
                                    k E2FInd(k,4)];
            elseif FsPerE==6
                S0n(RowIndex,1:2) = [k E2FInd(k,1)
                                    k E2FInd(k,2)
                                    k E2FInd(k,3)
                                    k E2FInd(k,4)
                                    k E2FInd(k,5)
                                    k E2FInd(k,6)];
            end
        end
    
        S0Sort=sort(S0n,2,'descend');  %Sort along rows left to right(brings IND for Face min-points to left)
        S0Sort=sortrows(S0Sort);       %Sort rows top to bottom
        
        % Find first and 2nd occurance of the same Face mid-point- 
        % intersections where 2 struts meet at a Face mid-point
        [~, i1 ,~] = unique(S0Sort(:,1),'first'); 
        [~, i2 ,~] = unique(S0Sort(:,1),'last');
        
         combinedind=[i1 i2];
         %if Face midpoint is unique - Keep strut from centroid->F1ce midpoint
         %if Face midpoint is not unique, connects element centroids to skip intersection at F1ce midpoint
         S_F1ce= combinedind(:,1) == combinedind(:,2);  
         S=[S0Sort(combinedind(S_F1ce,1),:)   ;...
         S0Sort(combinedind(~S_F1ce,1),2) S0Sort(combinedind(~S_F1ce,2),2)];
      
     
%         %Remove overhanging struts on outer regions of the lattice 
%         OuterNodes=  (N(:,3)==min(N(:,3)))  + (N(:,3)==max(N(:,3)))... 
%                 + (N(:,2)==min(N(:,2)))  + (N(:,2)==max(N(:,2)))...
%                 + (N(:,1)==min(N(:,1)))  + (N(:,1)==max(N(:,1))) ;
% 
%         N=N(~OuterNodes,:);
%         InnerNodes=find(~OuterNodes);   %Indices of outer nodes
%         OuterNodes=find(OuterNodes);   %Indices of outer nodes
% 
% 
%         OuterStruts= ismember(S(:,1),OuterNodes);  %Find outer stuts
%         InnerStruts= S(~OuterStruts,:);
%         InnerStruts=changem(InnerStruts,1:length(InnerNodes),InnerNodes);
%         S=InnerStruts;
%           
    case {'Edge','tetEdge'}
        [S]=patchEdges(F1);
        N = V1;

end
if SHAPE == "SquareSymetrical"  % Remove struts inner cells to make a symmetrical FE Model May only work for some unit cells
    % Get index of Nodes at minimum X and Y
    N_rm = find((N(:,1) == min(N(:,1))) + (N(:,2) == min(N(:,2)))) ;  
    % Get Struts contatining these nodes
    S_rm = logical(sum(ismember(S,N_rm),2));   
    % Remove these struts from strut array
    S =  S(~S_rm,:) ;
  

end

% Remove Unused nodes
[S,N]=patchCleanUnused(S,N);

%Merge duplicate nodes~
[S,N]=mergeVertices(S,N);

%Remove duplicate struts
S = sort(S,2);
S = unique(S,'stable','rows');    


% Apply displacement to make all co-ordinates positive
N(:,1) = N(:,1) + (0 - min(N(:,1)));
N(:,2) = N(:,2) + (0 - min(N(:,2)));
N(:,3) = N(:,3) + (0 - min(N(:,3)));

%Calculate Maxwell Number
M = size(S,1) - 3*size(N,1) +6;

%Calculate Fill Volume
if size(E1,2) == 4 
    [ElVol,~]=tetVol(E1,V1);
    ElVol=sum(ElVol);
elseif size(E1,2) == 8 
    [ElVol,~]=hexVol(E1,V1);
    ElVol=sum(ElVol);
end

if PlotOn==1
    %% Plot Lattice
    Xstruts= [N(S(:,1),1) N(S(:,2),1)]';
    Ystruts= [N(S(:,1),2) N(S(:,2),2)]';
    Zstruts= [N(S(:,1),3) N(S(:,2),3)]';
    
    cFigure;
    subplot(1,2,1)
    title('Base Mesh','FontSize',12);
    gpatch(F1,V1,'kw','black',0.5);
    axisGeom(gca,15);
    
    subplot(1,2,2)   
    plot3(Xstruts,Ystruts,Zstruts,'-k','LineWidth',1,'Marker','o','MarkerSize',2)
    title('Lattice','FontSize',12);
    axisGeom(gca,15);
    camlight headlight;
    drawnow;
end
end

%%
% Lattice Inverse Design & Optimisation Tool 
% 06/12/2023 - Brian McDonnell - University of Galway
% GNU AFFERO GENERAL PUBLIC LICENSE - See LICENSE file details