
#ifndef exrdmHistoMessenger_h
#define exrdmHistoMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class exrdmHisto;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class exrdmHistoMessenger: public G4UImessenger
{
  public:

   exrdmHistoMessenger(exrdmHisto* );
  ~exrdmHistoMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   exrdmHisto*                  histo;
   
   G4UIdirectory*          histoDir;   
   G4UIcmdWithAString*     factoryCmd;
   G4UIcmdWithAString*     fileCmd;
   G4UIcommand*            histoCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
