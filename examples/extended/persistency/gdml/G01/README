-------------------------------------------------------------------

     =========================================================
     Geant4 - an Object-Oriented Toolkit for Simulation in HEP
     =========================================================

                          GDML read/write
                          ---------------

This example demonstrates the usage of the GDML reader and writer. It allows
to export geometry descriptions in an application independent format (GDML,
Geometry Description Markup Language).
The GDML files can be then used to interchange geometries between different
applications and users.

The detector construction consists of a call to GDMLProcessor which parses a
GDML file and returns the pointer to the world volume. The user can also write
her/his own GDML file and use it as the primary input format for her/his Geant4
application.

Several simple GDML files are provided:
- axes.gdml,   showing loading and orientation of Cartesian axes;
- solids.gdml, list of all supported solids with placement; 
- scale.gdml,  a simple diamond structure made of extruded solids;
- divisionvol.gdml, a divided box;
- parameterized.gdml, a parameterised box;
- pTube.gdml, a parameterised tube;
- auxiliary.gdml, showing association of volume with auxiliary information;
- etc...

HOW TO BUILD THE EXAMPLE ?

- You need to have built the persistency/gdml module by having
  set the -DGEANT4_USE_GDML=ON flag during the CMAKE configuration step, 
  as well as the -DXERCESC_ROOT_DIR=<path_to_xercesc> flag pointing to 
  the path where the XercesC XML parser package is installed in your system.

- Compile and link to generate the executable (in your CMAKE build directory):
 	      % make

- Execute the application.
  o For reading and visualize interactively a GDML file:
 	      % load_gdml [GDML-file-in].gdml

  o For reading, writing and visualize interactively a GDML file:
 	      % load_gdml [GDML-file-in].gdml [GDML-file-out].gdml

  o For reading, writing a GDML file and running in batch a macro:
 	      % load_gdml [GDML-file-in].gdml [GDML-file-out].gdml [macro].in
