Note that this is a *generic* guide, allowing the use of the `FLUKA` interface in *any* G4 application.   
If you are interested in studying the XS and final states at the interaction level,   
or if you have never used the interface to `FLUKA` before,    
please directly follow the `README.md` files included in     
`ProcessLevel/CrossSection` and `ProcessLevel/FinalState`.   

If you want to integrate a physics list including FLUKA hadron-nucleus inelastic models into your G4 application,   
please follow `FlukaCern/FlukaInterface/README.md`.    
You will find there a generic guide, allowing you to use the `FLUKA` interface in *any* G4 application.    
It notably contains a list of dependencies, how to setup your environment,   
and how to adapt your build system.       

In order to include the `G4_HP_CernFLUKAHadronInelastic_PhysicsList` into the source of your application,   
you can then follow what has been done in `FlukaCern/ProcessLevel/CrossSection/HadronNucleusXS.cc`.   
Note that the integration of that `PhysicsList` is similar to the integration of any other custom PhysicsList,   
with, in addition, the `G4_USE_FLUKA` treatment, and the necessary call to `fluka_particle_table::initialize()`.   