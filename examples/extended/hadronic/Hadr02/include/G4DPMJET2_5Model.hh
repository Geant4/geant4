//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 19770/06/NL/JD (Technology Research Programme).         *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/include/G4DPMJET2_5Model.hh
/// \brief Definition of the G4DPMJET2_5Model class
//
// $Id: G4DPMJET2_5Model.hh 81932 2014-06-06 15:39:45Z gcosmo $
//

#ifndef G4DPMJET2_5Model_h
#define G4DPMJET2_5Model_h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4DPMJET2_5Model.hh
//
// Version:             0.B
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
//
//
// Class Description - End
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"

#include "G4ParticleTable.hh"

#include "G4HadronicInteraction.hh"
#include "G4WilsonAblationModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4VPreCompoundModel.hh"
#include "G4HadFinalState.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"
#include "G4Fragment.hh"
#include "G4HadProjectile.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include "G4GlaubAADataSetHandler.hh"
#include <sstream>

class G4GlaubAADataSetHandler;
enum G4DPMJET2_5InitialisationType { DEFAULT=1, DPM2_5=2, DPM3=3, CORSIKA=4 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DPMJET2_5Model : public G4HadronicInteraction
{
public:
//
//
// Standard constructor, destructor, copy etc declarations.
//
  G4DPMJET2_5Model ();
  G4DPMJET2_5Model (const G4DPMJET2_5InitialisationType);
  G4DPMJET2_5Model (G4ExcitationHandler *,
                    const G4DPMJET2_5InitialisationType initType = DEFAULT);
  G4DPMJET2_5Model (G4VPreCompoundModel *,
                    const G4DPMJET2_5InitialisationType initType = DEFAULT);
  virtual ~G4DPMJET2_5Model ();

  G4DPMJET2_5Model(const G4DPMJET2_5Model &right);

  const G4DPMJET2_5Model& operator=(G4DPMJET2_5Model &right);
    
  virtual G4bool IsApplicable(const G4HadProjectile &theTrack, 
                              G4Nucleus &theTarget);
    
  virtual G4HadFinalState *ApplyYourself(const G4HadProjectile &, 
                                         G4Nucleus &);

  virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;
    
  inline G4GlaubAADataSetHandler *GetGlauberDataSetHandler ();

  inline void SetVerboseLevel (G4int);
    
  void SetNoDeexcitation ();
  void SetDefaultDeexcitation ();
  inline void SetExcitationHandler (G4ExcitationHandler *);
  inline G4ExcitationHandler *GetExcitationHandler () const;
    
  void SetNoPreCompoundModel ();
  void SetDefaultPreCompoundModel ();    
  inline void SetPreCompoundModel(G4VPreCompoundModel* value);
  inline G4VPreCompoundModel* GetPreCompoundModel() const;
    
  inline void SetDPMInitialRandomSeeds (const G4int seed1, const G4int seed2);
  inline G4int GetDPMInitialRandomSeeds (const G4int i) const;

  inline G4double GetMinEnergy( const G4Material *aMaterial,
                                const G4Element *anElement ) const;
  inline G4double GetMaxEnergy( const G4Material *aMaterial,
                                const G4Element *anElement ) const;

  inline G4bool SetVerboseFortranOutput (const G4String filename);
  G4String GetVerboseFortranOutput () const;
    
  inline void SetDPMVariablesTAUFOR   (const G4double TAUFOR_P,
                                       const G4int    KTAUGE_P,
                                       const G4int    ITAUVE_P);
                                  
private:

  void DumpVerboseInformation1 (const G4int n) const;
  void DumpVerboseInformation2 (const G4String particleName,
                                const G4ThreeVector p, 
                                const G4double E, 
                                const G4double T,
                                const G4ThreeVector pinit) const;
  void DumpVerboseInformation3 (const G4int i, const G4int A, 
                                const G4int Z,
                                const G4ThreeVector p, 
                                const G4double E, 
                                const G4double T,
                                const G4ThreeVector pinit) const;
  void DumpVerboseInformation4 (const G4int i, 
                                const G4String particleName,
                                const G4ThreeVector p, 
                                const G4double E, 
                                const G4double T,
                                const G4ThreeVector pinit) const;
  void PrintWelcomeMessage () const;
  void Initialise ();
    
  G4DPMJET2_5InitialisationType theInitType;
                             
  G4GlaubAADataSetHandler *theGlauberDataSetHandler;
    
  G4ParticleTable         *theParticleTable;
  G4IonTable              *theIonTable;

  G4bool                   debug;
  G4int                    debug_level;
  G4int                    lunber;
  G4double                 dpmver;
  G4String                 defaultDirName;
    
  G4ExcitationHandler     *theExcitationHandler;
  G4VPreCompoundModel     *thePreComp;
    
  G4double                 TAUFOR;
  G4int                    KTAUGE;
  G4int                    ITAUVE;
  G4double                 UNON;
  G4double                 UNOM;
  G4double                 UNOSEA;
  G4double                 CVQ;
  G4double                 CDQ;
  G4double                 CSEA;
  G4double                 SSMIMA;
  G4double                 VVMTHR;
  G4double                 SEASQ;
  G4int                    MKCRON;
  G4double                 CRONCO;
  G4int                    ISINGD;
  G4int                    ISINGX;
  G4int                    IDUBLD;
  G4double                 SDFRAC;  
    
  ftnlogical LTRUE;
  ftnlogical LFALSE;

  G4String                 verboseFortranFile;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4GlaubAADataSetHandler *G4DPMJET2_5Model::GetGlauberDataSetHandler ()
  {return theGlauberDataSetHandler;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4DPMJET2_5Model::SetExcitationHandler
  (G4ExcitationHandler *aExcitationHandler)
  {theExcitationHandler = aExcitationHandler;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4ExcitationHandler *G4DPMJET2_5Model::GetExcitationHandler () const
  {return theExcitationHandler;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4DPMJET2_5Model::SetPreCompoundModel
  (G4VPreCompoundModel *aPreCompoundModel)
  {thePreComp = aPreCompoundModel;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4VPreCompoundModel *G4DPMJET2_5Model::GetPreCompoundModel () const
  {return thePreComp;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4DPMJET2_5Model::SetVerboseLevel (G4int verboseLevel1)
  {verboseLevel = verboseLevel1;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4DPMJET2_5Model::GetMinEnergy( const G4Material *,
  const G4Element * ) const
  {return theMinEnergy;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4DPMJET2_5Model::GetMaxEnergy( const G4Material *,
  const G4Element * ) const
  {return theMaxEnergy;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// SetVerboseFortranOutput
//
inline G4String G4DPMJET2_5Model::GetVerboseFortranOutput () const
  {return verboseFortranFile;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4DPMJET2_5Model::SetDPMInitialRandomSeeds (const G4int seed1,
  const G4int seed2)
{
  G4int s1 = seed1;
  G4int s2 = seed2;
  rd2in_ (&s1, &s2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4DPMJET2_5Model::GetDPMInitialRandomSeeds (const G4int i)
  const
{
  G4int seed1, seed2;
  rd2out_ (&seed1, &seed2);
  if (i == 0)      return seed1;
  else if (i == 2) return seed2;
  else             return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4DPMJET2_5Model::SetDPMVariablesTAUFOR (const G4double TAUFOR_P,
  const G4int KTAUGE_P, const G4int ITAUVE_P)
{
  TAUFOR = TAUFOR_P;
  KTAUGE = KTAUGE_P;
  ITAUVE = ITAUVE_P;
  taufo_.taufor = TAUFOR;
  taufo_.ktauge = KTAUGE;
  taufo_.itauve = ITAUVE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
