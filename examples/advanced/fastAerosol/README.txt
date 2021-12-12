 =========================================================
 		      Geant4 - FastAerosol advanced example
 =========================================================

                             README
                      ---------------------
Authors:

Ara Knaian    : ara@nklabs.com
Nate MacFadden: natemacfadden@gmail.com
NK Labs, LLC (http://www.nklabs.com)

Related Publication:

MacFadden, N., Knaian, A., 2020, "Efficient Modeling of Particle
Transport through Aerosols in GEANT4", Manuscript in preparation.

----------------------

 Using the FastAerosol geometry classes, it is possible to efficiently
 and accurately simulate particle transport through aerosols containing
 billions of randomly-positioned droplets, using an ordinary workstation.
 This example demonstrates the use of these classes. It is based off
 exampleB1.

 1- GEOMETRY DEFINITION

   FADetectorConstruction builds a user-defined-shape aerosol of
   non-intersecting, arbitrarily-shaped droplets. The default
   bulk shape is a 500 x 500 x 5000 mm box centerred at the origin;
   the default droplet shape is a r = 1 mm sphere. Intersection checking
   for non-spherical droplets is performed with the droplet's bounding
   raidus. These droplets are randomly placed in the bulk.

   A 500 x 500 x 50 mm box-shaped aluminum detector is placed at
   (0,0,2625) mm so that particles shot along +Z from (0,0,-2612.5) mm
   travel through the aerosol before being captured in the detector.

   The background material is atmospheric air at an altitude of 14 km.
   The droplet material is liquid water. 

 2- AEROSOL MODELLING

   The example can be configured to allow modeling of the aerosol using 
   one of three methods, to allow for benchmarking and testing.
   The methods are (a) using the demonstrated FastAerosol class,
   (b) using a single low-density volume, and (c) using parameterized
   volumes. Alternating between these build methods can be done by
   setting one of "/geometry/fastAerosolCloud", "/geometry/smoothCloud",
   and "/geometry/parameterisedCloud" true for a, b, and c respectively.

   By default, FastAerosol dynamically populates droplets in the bulk
   as particles transport through it. This reduces memory consumption,
   especially when the number of primaries is small compared to the
   number of droplets. Pre-population is slower and uses more memory,
   but allows saving the droplet distribution over the full volume.
   Pre-population may be enabled using the "/geometery/prePopulate true"
   UI command. The example saves the droplet positions in files called
   "distribution_r???mm_n???mm-3.csv" where ??? denotes the droplet
   radius and number density respectively).

   User-specified distibution functions for droplet position and rotation
   may be specified for FastAerosol simulations. Examples of this can be
   found in lines 311-317 and 326-330 of FADetectorConstruction.cc.

   The example can be configured to model the aerosol using a
   G4PVParameterized object, instead of using FastAerosol object.
   This can be used to demonstrate that FastAerosol gives nearly
   identical transport results as parameterized geometery, but
   achieves order-of-magnitude performance gains.
   
   To build the aerosol as a G4PVParameterised object, a droplet
   distribution is required. The example is set up to use the droplet
   distribution automatically generated and saved by a previously-run
   pre-populated FastAerosol simulation, as this allows better comparison
   between the two classes. By default, the example looks for this
   distribution under the file name "distribution_r???mm_n???mm-3.csv",
   so the generated pre-populated FastAerosol distribution can be used
   by simply running a pre-populated FastAerosol simulation with the 
   same geometry before running the parameterised simulation.
   	
   For comparison, the example can also be run using a single volume
   containg air and water mixed together at the correct average density.
   For small-enough droplets, this works well. However, it becomes
   increasingly innacurate as the droplets appoach a critical size range,
   which is about r = 1 mm for the materials and energies used in this
   example. This demonstraties the need to model some aerosols at a
   droplet level for accurate physics results.

 2- PHYSICS LIST
 
   This example uses the QGSP_BIC physics list. A global step length
   limiter physics process is also included, because this can significantly
   speed up calculations using FastAerosol when particle trajectories 
   curve due to the application of fields.
   
 3- ACTION INITALIZATION

    Nothing special.

 4- PRIMARY GENERATOR
  
   Particles are shot from (0,0,-2612.5) mm with a spread in X and Y of
   +/-55 mm with momentum in the +z direction. The example shoots 50 MeV
   protons by default.

 5- DETECTOR RESPONSE

   Scoring is done with a scoring grid, to allow histograms of the energy
   deposited to allow comparison of the results between geometry modeling
   methods. 

 6- SAMPLE EVALUATION

   It is recommended to run the test.mac script for a test/sample
   evaluation. This will simulate the transport of 100 protons (50 MeV)
   through an aerosol. The default bulk shape is "box"; other allowed
   shapes are "ellipsoid", "cylinder", and "pipe". All shapes are aligned
   along the z axis and maximally sized to fit in the 500 x 500 x 5000 mm
   bounding box. 

   By default, droplets are dynamically populated in the aerosol. One
   can populate all droplets at the beginning of the simulation by
   setting prePopulate to true. This saves the generated distribution
   of droplet centers. One can also run the same experiment, except
   modelling the aerosol as a parameterised solid, by setting
   parameterisedCloud to true and FastAerosolCloud to false. This
   requires a distribution of droplet centers; the one generated by a
   pre-opulated FastAerosol simulation automatically works. One can
   also simulate the cloud as a single average-density object by setting
   FastAerosolCloud and parameterisedCloud to false and smoothCloud to
   true.

 7 - OUTPUT

   The energy deposited into the aluminum detector by any particles is
   measured with command line scoring by a 20 x 20 x 1 scoring grid
   and saved via "/score/dumpQuantityToFile" as a csv file named
   "eDep_FastAerosol_r1p0mm_n1E-3p7mm-3_bulkbox_dropletsphere.csv"
   for the default aerosol build parameters.

   The GNU program "/usr/bin/time" may be used measure the simulation
   time and memory load. If simulating a pre-populated FastAerosol cloud,
   the time to populate the cloud is printed at population time in the
   program. This population time is also saved as
   "popTime_r???mm_n???mm-3.csv".

   The distribution of droplet centers is saved as
   "distribution_r???mm_n???mm-3.csv" where each row of this csv file
   corresponds to a droplet with the 1st, 2nd, and 3rd columns
   corresponding the the x, y, and z position of it's center respectively.
   
8 - HOW TO RUN THE EXAMPLE
    Batch mode: fastAerosol test.mac
    Interactive mode: fastAerosol (init_vis.mac is executed to set-up visualisation)
     
   