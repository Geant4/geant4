#ifndef G4ParallelTransport_hh
#define G4ParallelTransport_hh G4ParallelTransport_hh

#include "G4VProcess.hh"
#include "globals.hh"
#include "g4std/strstream"

class G4VPGeoDriver;
class G4VParallelStepper;

class G4ParallelTransport : public G4VProcess {
public:
  G4ParallelTransport(G4VPGeoDriver &pgeodriver, 
		      G4VParallelStepper &aStepper,
		      const G4String &aName = "ParallelTransport");
  virtual ~G4ParallelTransport();

  // the post step pair of functions may be overwritten
  // in certain derived classes like the importance process
  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);
  
  virtual G4VParticleChange * 
  PostStepDoIt(const G4Track& ,
	       const G4Step&);
  
  //  no operation in  AtRestDoIt and  AlongStepDoIt
  G4double 
  AlongStepGetPhysicalInteractionLength(const G4Track&,
					G4double  ,
					G4double  ,
					G4double& ,
					G4GPILSelection*
					){ return -1.0; };
  
  G4double 
  AtRestGetPhysicalInteractionLength(const G4Track& ,
				     G4ForceCondition* 
				     ) { return -1.0; };
  
  //  no operation in  AtRestDoIt and  AlongStepDoIt
  G4VParticleChange* 
  AtRestDoIt(const G4Track& ,
	     const G4Step&
	     ) {return 0;};
  G4VParticleChange* AlongStepDoIt(const G4Track& ,
				   const G4Step& 
				   ) {return 0;};

  G4VParallelStepper &GetPStepper(){return fPStepper;}

protected:
  virtual void Error(const G4String &m){
    G4cout << "ERROE: in G4ParallelTransport: " << m << G4endl;
    exit(1);
  }
  virtual void Warning(const G4String &m);

  G4ParticleChange *fParticleChange;


private:
  G4ParallelTransport(const G4ParallelTransport &);
  G4ParallelTransport &operator=(const G4ParallelTransport &);
  G4VPGeoDriver &fPgeodriver;
  G4VParallelStepper &fPStepper;
  G4bool fCrossBoundary;

};


#endif
