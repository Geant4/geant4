#ifndef Beam_H
#define Beam_H 1

#include "G4LorentzVector.hh"
#include "G4String.hh"

class Beam
{

   public:
   
      Beam() : fBeamPartName(""), fBeamPartMass(0.), fBeamEnergy(0.),
                fLabV(G4LorentzVector()), fLabP(G4LorentzVector()) {}
      ~Beam() {}
      
      void SetBeam( G4String name, G4double mass, G4double energy ) { fBeamPartName=name;
                                                                      fBeamPartMass=mass;
						                      fBeamEnergy=energy;
						                      return; }
						 
      void SetLabV( G4LorentzVector lv ) { fLabV=lv; return; }
      void SetLabP( G4LorentzVector lp ) { fLabP=lp; return; }
      
      G4String         GetBeamPartName()   const { return fBeamPartName; }
      G4double         GetBeamPartMass()   const { return fBeamPartMass; }
      G4double         GetBeamEnergy()     const { return fBeamEnergy; }
      const G4LorentzVector&  GetLabV()    const { return fLabV; }
      const G4LorentzVector&  GetLabP()    const { return fLabP; }

   private:
   
      G4String        fBeamPartName;
      G4double        fBeamPartMass;
      G4double        fBeamEnergy;
      G4LorentzVector fLabV;
      G4LorentzVector fLabP;

};

#endif
