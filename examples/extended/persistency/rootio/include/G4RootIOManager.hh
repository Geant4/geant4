// $Id: G4RootIOManager.hh,v 1.1 2002-12-04 02:44:28 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4RootIOManager.hh
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#ifndef ROOT_IO_MANAGER_HH
#define ROOT_IO_MANAGER_HH 1

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4VPhysicalVolume.hh"
#include <string>
#include <map>

class G4VPHitIO;
class G4VPDigitIO;
class VPHepMCIO;
class VPMCTruthIO;
class G4RootTransManager;
class G4HitRootIO;
class G4DigitRootIO;
class G4EventRootIO;
class G4HepMCRootIO;
class G4MCTruthRootIO;

// Class inherited:
#include "G4PersistencyManager.hh"

// Class Description:
//   Implementation of G4PersistencyManager with ROOT I/O package

class G4RootIOManager
 : public G4PersistencyManager
{
    friend class G4RootTransManager;

    public: // With description
      G4RootIOManager(G4PersistencyCenter* pc, std::string n);
      // Constructor

      virtual ~G4RootIOManager();
      // Destructor

    public: // With description
      G4VPHitIO* HitIO() { return (G4VPHitIO*) f_G4HitRootIO; };
      // returns the hit I/O handling manager.

      G4VPDigitIO* DigitIO() { return (G4VPDigitIO*) f_G4DigitRootIO; };
      // returns the digit I/O handling manager.

      G4VPEventIO* EventIO() { return (G4VPEventIO*) f_G4EventRootIO; };
      // returns the event I/O handling manager.

      G4VHepMCIO* HepMCIO() { return (G4VHepMCIO*) f_G4HepMCRootIO; };
      // returns the HepMC event I/O handling manager.

      G4VMCTruthIO* MCTruthIO() { return (G4VMCTruthIO*) f_G4MCTruthRootIO; };
      // returns the HepMC event I/O handling manager.

      G4VTransactionManager* TransactionManager() { return (G4VTransactionManager*) f_G4RootTransManager; };
      // returns the event I/O handling manager.

      void Initialize();
      // Initialize the ROOT I/O package.

    private:
      G4PersistencyManager* Create() {return 0;};
      // Dummy in this class (used by PersistentManagerT).

    private:
      G4RootTransManager* f_G4RootTransManager;
      G4HitRootIO*        f_G4HitRootIO;
      G4DigitRootIO*      f_G4DigitRootIO;
      G4EventRootIO*      f_G4EventRootIO;
      G4HepMCRootIO*      f_G4HepMCRootIO;
      G4MCTruthRootIO*    f_G4MCTruthRootIO;

}; // End of class G4RootIOManager

#endif

