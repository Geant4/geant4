#include "G4ios.hh"
#include "globals.hh"
#include "g4templates.hh"

#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"
#include "G4Cerenkov.hh"
#include "G4OpRayleigh.hh"
#include "G4OpAbsorption.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

void LoopUntilPressEnter(void);

void main()
{

/////////////////////////
// Material Definition
/////////////////////////

  G4double a, z, density;
  G4String name, symbol;
  G4int nel;

// Water
// -----

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., a);

  density = 1.0*g/cm3;
  G4Material* Water = new G4Material(name="Water", density, nel=2);

  Water->AddElement(elH, 2);
  Water->AddElement(elO, 1);

/////////////////////////////////////////
// Material Properties Table Definition
/////////////////////////////////////////

  const G4int NUMENTRIES = 32;

  G4double PPCKOV[NUMENTRIES] =
            { 2.038E-9*GeV, 2.072E-9*GeV, 2.107E-9*GeV, 2.143E-9*GeV,
              2.181E-9*GeV, 2.220E-9*GeV, 2.260E-9*GeV, 2.302E-9*GeV,
              2.346E-9*GeV, 2.391E-9*GeV, 2.438E-9*GeV, 2.486E-9*GeV,
              2.537E-9*GeV, 2.590E-9*GeV, 2.645E-9*GeV, 2.702E-9*GeV,
              2.763E-9*GeV, 2.825E-9*GeV, 2.891E-9*GeV, 2.960E-9*GeV,
              3.032E-9*GeV, 3.108E-9*GeV, 3.188E-9*GeV, 3.271E-9*GeV,
              3.360E-9*GeV, 3.453E-9*GeV, 3.552E-9*GeV, 3.656E-9*GeV,
              3.767E-9*GeV, 3.884E-9*GeV, 4.010E-9*GeV, 4.144E-9*GeV };

  G4double RINDEX[NUMENTRIES] =
            { 1.33, 1.33, 1.33, 1.33, 1.33, 1.33, 1.33,
              1.33, 1.33, 1.34, 1.34, 1.34, 1.34, 1.34,
              1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34,
              1.34, 1.34, 1.35, 1.35, 1.35, 1.35, 1.35,
              1.35, 1.35, 1.35, 1.35 };

  G4MaterialPropertiesTable myMPT;
  myMPT.AddProperty("RINDEX", PPCKOV, RINDEX, NUMENTRIES);

  Water->SetMaterialPropertiesTable(&myMPT);

////////////////////////
// Particle Definition
////////////////////////

  G4OpticalPhoton* theOpticalPhoton =
			G4OpticalPhoton::OpticalPhotonDefinition();

  G4MuonPlus* theMuonPlus = G4MuonPlus::MuonPlusDefinition();

///////////////////////
// Process Definition
///////////////////////

  G4Cerenkov theCerenkovProcess;
  G4OpRayleigh theRayleighProcess;
  G4OpAbsorption theAbsorptionProcess;

//////////////////////////////////
// Test G4Cerenkov Physics Table
//////////////////////////////////

  G4cout << "G4Cerenkov's Physics Table" << G4endl;
  G4cout << "--------------------------" << G4endl;

  theCerenkovProcess.DumpPhysicsTable();
  LoopUntilPressEnter();

  G4cout << "G4OpAbsorption's Physics Table" << G4endl;
  G4cout << "----------------------------" << G4endl;

  theAbsorptionProcess.DumpPhysicsTable();
  LoopUntilPressEnter();

  G4cout << "G4OpRayleigh's Physics Table" << G4endl;
  G4cout << "--------------------------" << G4endl;

  theRayleighProcess.DumpPhysicsTable();
  LoopUntilPressEnter();

}

void LoopUntilPressEnter()
{
        char ch;

        G4cout << "Press <Enter> to continue ... ";
        while ( G4cin.get(ch) )
        {
                if (ch == '\n') break;
        }
        G4cout << G4endl;
}

