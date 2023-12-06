% Demo of exporting lattice to .3mf file


%% Define custom hex mesh to fill with lattice
% (lattice variables controlling cell size/aspect ratio are overwritten)
% optionStruct.sphereRadius=5;
% optionStruct.coreRadius=3;
% optionStruct.numElementsMantel=2;
% optionStruct.numElementsCore=3;
% optionStruct.makeHollow=1;
% optionStruct.outputStructType=2;
% 
% %Creating sphere
% [meshStruct]=hexMeshSphere(optionStruct);

%%
name = 'Cubic1cell' ;
IN = [5,1,1 1];
export=0;   % Export = 1 to write to file, 0 to view.
%IN = x;
%lat_opts.meshStruct = meshStruct;
%% Lattice Geometry Settings
lat_opts.Shape='Cylinder'; % 'Square' or 'Cylinder'; 
lat_opts.ZDir='Axial';        %["Radial","Hoop","Axial"] orientation of vertical struts for cylinders

% Number of unit cells in [x,y,z] [radial,axial,inner-diameter]
lat_opts.Dim=[2 2 10];

% base unit cell size [x,y,z] or [radial,hoop,axial]
%Size=[5 5 5];
%Size=[1000 1000 1000*IN(1)];

% Change aspect ratio depending on number of cells:
lat_opts.Size=[4 4 4];
%lat_opts.Size=[4 4 4*(lat_opts.Dim(3)/IN(1))];
%lat_opts.Dim(3)=IN(1);

%lat_opts.nominal_diam = IN(3);
lat_opts.nominal_diam = 0.7;

taper_on =0; gradient_on =0; AR_grad_on=0;  

%lat_opts.taper=[1 IN(2) 1];          % Global Strut taper
lat_opts.taper=[1 0.8 0.8 0.8 1];          % Global Strut taper

%lat_opts.gradient =[lat_opts.nominal_diam IN(4:end)];    % Z-Graded strut thicknes
lat_opts.gradient =[1 1 1 1 0.6 0.6 0.6];    % Z-Graded strut thicknes
lat_opts.AR_gradient=[1 0.5 1];   % Z-Graded Angle Ratio


lat_opts.el_per_strut=1;
lat_opts.el_lengths=[];
lat_opts.element_type=1;
%-TYPE -  one of  ["BCC","BCCZ","FCC","FCCZ","OCTET","Dual","Edge", "tetDual","tetEdge"]
lat_opts.lattice_type="Edge"; 

if taper_on==0 
    lat_opts.taper=[]; 
end

if gradient_on==0 
    lat_opts.gradient =[];    % Z-Graded strut thicknes 
end

if AR_grad_on==0 
    lat_opts.AR_gradient=[];   % Z-Graded Angle Ratio 
end

Geometry = GenerateLatticeExport(lat_opts);


%%
%write3mf('C:\Users\16438722\Documents\GitHub\Lattice\Matlab\3mf\test2\FCC_2mm.3mf', N, [], S,[],DiamU,[]);

if export == 1
    file_name = ['C:\Users\16438722\Documents\GitHub\All-PhD-Codes\Lattice\3mf\'...
                  name '.3mf' ];

    write3mf(file_name, Geometry.V, [], Geometry.E(:,[1 2]),[],Geometry.Diameters);
else
    ViewLattice(Geometry,'Generate',0,'n_sides',20)
end