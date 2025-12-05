\page ExampleuserPrimaryGenerator Example userPrimaryGenerator

 The example shows how to create a primary event including several vertices and 
 several primary particles per vertex

## Geometry construction

 It is a simple box which represents an 'infinite' homogeneous medium.
  
## Physics list

 PhysicsList defines only geantino and transportation process.

## Primary generator : several vertices and particles per event

 - vertex A and particle 1 : a geantino uniformly randomized on a cylinder surface.
 - vertex B and particles 2 and 3 : symetric to vertex A.

## Visualisation

 Visualization Manager is set in the main().
 Initialisation of the drawing is done via the commands
 `/vis/..`` in the macro vis.mac. This macro is automatically read from the main 
 in case of interactive running mode.

## How to start ?

  - execute userPrimaryGenerator in 'batch' mode from macro files
```
	% ./userPrimaryGenerator run1.mac
```

  - execute userPrimaryGenerator in 'interactive mode' with visualization
```
	% ./userPrimaryGenerator
	....
	Idle>  ---> type your commands. For instance:
	Idle> /run/beamOn 1
	....
	Idle> exit
```
