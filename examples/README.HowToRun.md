\page README_HowToRun How to build and run an example

All basic and most of extended examples have a similar structure.
The `main()` function is included in `exampleXYZ.cc` in the example top directory
and the example source code is structered in `include` and `src` subdirectories.
When the example is built, the executable takes the same name as the file with
`main()` function without the `.cc` extension, `exampleXYZ`.

Then several macros are provided to run the example with various start-up
conditions. These macros have usually the `.mac` extension. Besides these macros,
there is often a macro `exampleXYZ.in` (note its different extension)
which is used in Geant4 testing and which output, `exampleXYZ.out`, can be also included
in the distribution.

You can find all details about building the examples in the <a href= "https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/index.html"> Geant4 Installation Guide </a>, in the section <a href="https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/buildtools.html">  How to Use the Geant4 Toolkit Libraries </a>.
Here we recall only the basics. 

## Compile and link to generate an executable

```
% cd path_to_exampleXYZ     # go to directory which contains your example
% mkdir exampleXYZ_build
% cd exampleXYZ_build
% cmake -DGeant4_DIR=path_to_Geant4_installation/lib[64]/cmake/Geant4/ ../exampleXYZ
% make -j N exampleXYZ      # "N" is the number of processes 
```

## Execute exampleXYZ in 'batch' mode from macro files
 		
```
% cd exampleXYZ_build
% ./exampleXYZ  xyz.mac
```
 		
## Execute exampleXYZ in the 'interactive mode' with visualization
```
% cd exampleXYZ_build
% ./exampleXYZ
....
Idle> type your commands
....
Idle> exit
```

See also the explicit instructions for building TestEm1 example at 
\ref README_HowToRunTestEm1 page.

\page README_HowToRunTestEm1 How to build and run TestEm1 example

Below we give the explicit instructions (described in general at \ref README_HowToRun page) for the example extended/electromagnetic/TestEm1.

Let's suppose that the `TestEm1` directory is available in `$HOME` and Geant4 installation in `/usr/local` and we work within the `bash` shell on a `64-bit machine.

## Compile and link TestEm1 to generate an executable

```
% cd $HOME
% mkdir TestEm1_build
% cd TestEm1_build
% cmake -DGeant4_DIR=/usr/local/lib64/Geant4-11.3.0/ ../TestEm1
% make -j 2 TestEm1
```

## Execute TestEm1 in 'batch' mode from macro files
```
% cd $HOME/TestEm1_build
% ./TestEm1 annihil.mac
% ./TestEm1 brem.mac
% ./TestEm1 TestEm1.in >& myTestEm1.out   # redirecting output in a file
```
 		
## Execute TestEm1 in 'interactive mode' with visualization
```
% cd $HOME/TestEm1_build
% ./TestEm1
PreInit> /run/initialize 
Idle>    /run/beamOn 1
...
Idle>    exit
```
