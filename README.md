# LiDOT
Lattice Inverse Design &amp; Optimisation Tool

## Description

Scripts developed for optimising lattice structures based on ABAQUS FEA. Lattices are generated in a beam-elemnt format which can be written to an ABAQUS .INP file or exported to .3MF to be printed (beam-element extension supported mainly in software/slicers with a focus on metal AM, tested in 3DXpert, Netfabb, nTopology)

## Getting Started

### Dependencies
* [GIBBON](https://www.gibboncode.org)
* Matlab Optimisation Toolbox / Global Optimisation Toolbox 
* ABAQUS for FEA

### Installing
* Download and install the Matlab toolbox [GIBBON](https://www.gibboncode.org)
* Clone repository / download and extract files, ensure GIBBON is added to path as per the install instructions.
* Add subfolders to path
* Check file path for ABAQUS (set in ExplicitBeamCompression/StandardBeamCompression)

### Executing program
Scripts are set up as applied for publication....

Key functions:
Optimisation
* Script to run optimisation algorithm and visualise results
* Set varialbes, bounds, constraints, optimisation algorithm parameters etc...

OptRun (Objective function)
* Generate lattice based on input variables, run FE with ABAQUS and return objective function.
* The variable func_mode allows you to visualise lattices or check the target curve instead of running ABAQUS
* For sample inputs, un-comment line 13:
```
%IN = [3 70 100 90 80 70];    % Variables Example
```

GenerateLattice / GenerateLatticeExport
*Function generates beam-element lattices based on an input structure of lattice parameters.
*Input an empty structure for defaults/demo. See function OptRun for details on all settings.
*e.g. 
```
%% Generate example lattice to view output structure
GenerateLattice((struct))

%% 3D visualisation of lattice
ViewLattice(GenerateLattice((struct)))
```

WireLattice
* Generates lattices based on unit cells of a set size and grid (e.g BCC, FCC, or OCTET), or on the edges or dual lattice of an input hex or tet mesh.
* Lattices of unit cells can also be generated based on an input hex mesh.
* Output is: N - nodes of lattices [x y z] co-ordinates, and S - strut array with nodes in each strut.

LatticeExport_DEMO
* Demo of exporting lattice to .3mf with write3mf

## Version History

* 1.0
    * Initial Release

## License

This project is licensed under the GNU Affero General Public License version 3 License - see the LICENSE.md file for details

## Acknowledgments

* [GIBBON](https://www.gibboncode.org)
* [write3mf](https://github.com/cvergari/write3mf)

