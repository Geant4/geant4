     =========================================================
     Geant4 - an Object-Oriented Toolkit for Simulation in HEP
     =========================================================

                            fastAerosol
                            -----------

 This example demonstrates the fastAerosol and fastAerosolSolid
 classes for accurately and efficiently simulating aerosols with
 many droplets.  

 Using the fastAerosol geometry classes, it is possible to efficiently
 simulate clouds containing billions of randomly-positioned droplets
 in GEANT4.  

 This example demonstrates the use of fastAerosol.  It is based
 off exampleB1.

 0- SAMPLE EVALUATION

   It is recommended to run the test.mac script for a test/sample
   evaluation. This will simulate transport of 50 MeV through an
   aerosol (modelled as a fastAerosol object) and save the deposited
   energy in a detector as a csv file.

   By default, droplets are dynamically populated in the aerosol. One
   can populate all droplets at the beginning by setting prePopulate
   to true (this saves the generated distribution of droplet centers).
   One can also run the same experiment, except modelling the aerosol
   as a parameterised solid by setting parameterisedCloud to true and
   fastAerosolCloud to false.

 1- GEOMETRY DEFINITION

   FADetectorConstruction builds a user-defined-shaped aerosol of
   arbitarily-shaped droplets (placed like spheres), centered at the
   origin, with a box-shaped aluminum detector at one end along the
   z-axis.

   The background material was set to be air at 14km and the cloud
   is allowed to consist of either water or ice droplets.

   The example allows three methods for building the cloud: 'granular',
   'parameterised', and 'smooth'.  'granular' demonstrates building
   the cloud using the new fastAerosol class.  'parameterized' and 
   'smooth' are included for comparision and benchmarking purposes.
   The 'smooth' method treats the aerosol as containing water at its
   average density; the 'parameterised' method models the aerosol
   as spherical droplets using G4PVParameterised.

   The 'smooth' cloud is included to demonstrate the need to model
   some aerosols at a droplet level for accurate physics results.
   Using the example, one can demonstrate that the fastAerosol class
   returns nearly the same results as 'paramaterised,' but using 
   significantly less memory and computation time.
   	
 2- PHYSICS LIST
 
   This simulation uses the QGSP_BIC physics list.  A global step length
   limiter physics process is also included, because this can signifiantly
   speed up calculations using fastAerosol when particle trajectories 
   curve due to the application of fields.
   
 3- ACTION INITALIZATION

   Nothing special.
  	 
 4- PRIMARY GENERATOR
  
   The primary generator is a disk behind the cloud shooting towards
   the detector (+Z) through the cloud. The example uses 50 MeV protons.
     
 5- DETECTOR RESPONSE

   Scoring is done primarily with a scoring grid, to allow histograms
   of the energy deposited to allow comparison of the results between
   geometry modeling methods. 