#ifndef PRIMARYGENERATORMESSENGER_HH
#define PRIMARYGENERATORMESSENGER_HH 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGenerator;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;


class PrimaryGeneratorMessenger : public G4UImessenger {

 public:
   PrimaryGeneratorMessenger(PrimaryGenerator*);
   ~PrimaryGeneratorMessenger();

   void SetNewValue(G4UIcommand*, G4String);

 private:
   PrimaryGenerator* primaryGenerator;

   G4UIdirectory* sourceDirectory;
   G4UIcmdWithADoubleAndUnit* primEnergyCmd;
   G4UIcmdWithADoubleAndUnit* sigmaEnergyCmd;
   G4UIcmdWithADoubleAndUnit* sigmaSpatialCmd;
   G4UIcmdWithADoubleAndUnit* incidAngleCmd;
};

#endif // PRIMARYGENERATORMESSENGER_HH
