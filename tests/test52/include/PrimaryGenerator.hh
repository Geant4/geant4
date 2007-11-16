#ifndef PRIMARYGENERATOR_HH
#define PRIMARYGENERATOR_HH 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4Event;
class G4ParticleGun;
class PrimaryGeneratorMessenger;


class PrimaryGenerator : public G4VUserPrimaryGeneratorAction {

 public:
   PrimaryGenerator();
   ~PrimaryGenerator();

   void GeneratePrimaries(G4Event*);

   void SetPrimaryKineticEnergy(G4double);
   void SetSigmaKineticEnergy(G4double);
   void SetSigmaSpatialPlacement(G4double);
   void SetIncidentAngle(G4double);

 private:
   G4ParticleGun* particleGun;
   PrimaryGeneratorMessenger* messenger;

   G4double primaryKineticEnergy;
   G4double sigmaKineticEnergy;
   G4double sigmaSpatialPlacement;
   G4double incidentAngle;
};

#endif // PRIMARYGENERATOR_HH
