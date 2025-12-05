\page ExampleP03  Example P03

This example illustrates the use of the text geometry. 

## GEOMETRY EXAMPLES

Several examples of text geometries are provided:

g4geom_simple.txt :
- Simple construction of materials and single placements 

g4geom_matemixt.txt :
- Isotopes, elements and materials

g4geom_boolean.txt :
- Boolean solids
 
g4geom_reflections.txt :
- Reflections

g4geom_replicas.txt :
- Replicas

g4geom_divisions.txt :
- Divisions

g4geom_paramLinear.txt :
- Linear parameterisations

g4geom_paramSquare.txt :
- Square parameterisations

g4geom_assembly.txt : 
- Assembly placements

     	
## HOW TO START ?
 
Execute textGeom in 'batch' mode from macro file 
```
% ./textGeom run.mac
```
 		
It will read the geometry from a file named 'g4geom.txt', and it will create a
VRML2 file to visualise the geometry.
Therefore if you want to try any of the above-mentioned files, copy it to a
file with this name,

 	
## DEFINING A SENSITIVE DETECTOR

The detector construction class ExTGDetectorConstructionWithSD shows how to
access a volume of the text geometry and assign to it a sensitive detector.
To use it, replace at exampleTextGeom.cc the line

```cpp
runManager->SetUserInitialization(new ExTGDetectorConstruction);
```

by

```cpp
runManager->SetUserInitialization(new ExTGDetectorConstructionWithSD);
```

and recompile.

It will use the geometry from the file 'g4geom_SD.txt'


## MIXING TEXT AND C++ GEOMETRIES

The detector construction class ExTGDetectorConstructionWithCpp shows how to
create a volume with C++ and place it in the world or inside a volume defined
in the text geometry.
To use it, replace at exampleTextGeom.cc the line

```cpp
runManager->SetUserInitialization(new ExTGDetectorConstruction);
```

by

```cpp
runManager->SetUserInitialization(new ExTGDetectorConstructionWithCpp);
```

and recompile.

It will use the geometry from the file 'g4geom_simple.txt'


## CREATING NEW TAGS IN THE TEXT GEOMETRY: DEFINING CUTS PER REGION

The detector construction class ExTGDetectorConstructionWithCuts, together with
ExTGRCLineProcessor, ExTGRCDetectorBuilder, ExTGRCRegionCutsMgr and
ExTGRCRegionData show how to add a couple of tags, ':REGION' and ':CUT', that
allow to define cuts per region in your input geometry text file. 
To use it, replace at exampleTextGeom.cc the line

```cpp
runManager->SetUserInitialization(new ExTGDetectorConstruction);
```

by

```cpp
runManager->SetUserInitialization(new ExTGDetectorConstructionWithCuts);
```

and recompile.

It will use the geometry from the file 'g4geom_cutsPerRegion.txt'


## DUMP THE IN-MEMORY GEOMETRY TO TEXT FILE

The run action, ExTGRunAction, triggers the writing of the in-memory Geant4
geometry to a text file. 
To use it you just have to uncomment in exampleTextGeom.cc the line

```cpp
runManager->SetUserAction(new ExTGRunAction);
```

and it will read the geometry from the file 'g4geom.txt' and will write the
geometry in a file named 'geom.txt'
