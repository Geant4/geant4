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

////////////////////////////////////////////////////////////////////////////////
//
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
    ~G4DPMJET2_5Model ();

    G4DPMJET2_5Model(const G4DPMJET2_5Model &right);

    const G4DPMJET2_5Model& operator=(G4DPMJET2_5Model &right);

    G4GlaubAADataSetHandler *GetGlauberDataSetHandler ();
    
    G4bool IsApplicable (const G4HadProjectile &theTrack, G4Nucleus &theTarget);
    
    virtual G4HadFinalState *ApplyYourself
      (const G4HadProjectile &, G4Nucleus &);
    
    void SetVerboseLevel (G4int);
    
    void SetExcitationHandler (G4ExcitationHandler *);
    void SetNoDeexcitation ();
    void SetDefaultDeexcitation ();
    G4ExcitationHandler *GetExcitationHandler () const;
    
    void SetPreCompoundModel(G4VPreCompoundModel* value);
    void SetNoPreCompoundModel ();
    void SetDefaultPreCompoundModel ();    
    G4VPreCompoundModel* GetPreCompoundModel() const;
    
    void SetDPMInitialRandomSeeds (const G4int seed1, const G4int seed2);
    G4int GetDPMInitialRandomSeeds (const G4int i) const;

    G4double GetMinEnergy( const G4Material *aMaterial,
                                  const G4Element *anElement ) const;
    G4double GetMaxEnergy( const G4Material *aMaterial,
                                  const G4Element *anElement ) const;

    G4bool SetVerboseFortranOutput (const G4String filename);
    G4String GetVerboseFortranOutput () const;
    
    void SetDPMVariablesTAUFOR   (const G4double TAUFOR_P,
                                  const G4int    KTAUGE_P,
                                  const G4int    ITAUVE_P);
/*    void SetDPMVariablesXCUTS    (const G4double CVQ_P,
                                  const G4double CDQ_P,
                                  const G4double CSEA_P,
                                  const G4double SSMIMA_P);
    void SetDPMVariablesCRONINPT (const G4int    MKCRON_P,
                                  const G4double CRONCO_P);
    void SetDPMVariablesSEADISTR (const G4double 
                                  const G4double UNON_P,
                                  const G4double UNOM_P,
                                  const G4double UNOSEA_P);*/
                                  
  private:
    void DumpVerboseInformation1 (const G4int n) const;
    void DumpVerboseInformation2 (const G4String particleName,
      const G4ThreeVector p, const G4double E, const G4double T,
      const G4ThreeVector pinit) const;
    void DumpVerboseInformation3 (const G4int i, const G4int A, const G4int Z,
      const G4ThreeVector p, const G4double E, const G4double T,
      const G4ThreeVector pinit) const;
    void DumpVerboseInformation4 (const G4int i, const G4String particleName,
      const G4ThreeVector p, const G4double E, const G4double T,
      const G4ThreeVector pinit) const;
    void PrintWelcomeMessage () const;
    void Initialise ();
    
  private:
    G4DPMJET2_5InitialisationType
                             theInitType;
                             
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
////////////////////////////////////////////////////////////////////////////////
//
inline G4GlaubAADataSetHandler *G4DPMJET2_5Model::GetGlauberDataSetHandler ()
  {return theGlauberDataSetHandler;}
////////////////////////////////////////////////////////////////////////////////
//
inline void G4DPMJET2_5Model::SetExcitationHandler
  (G4ExcitationHandler *aExcitationHandler)
  {theExcitationHandler = aExcitationHandler;}
////////////////////////////////////////////////////////////////////////////////
//
inline G4ExcitationHandler *G4DPMJET2_5Model::GetExcitationHandler () const
  {return theExcitationHandler;}
////////////////////////////////////////////////////////////////////////////////
//
inline void G4DPMJET2_5Model::SetPreCompoundModel
  (G4VPreCompoundModel *aPreCompoundModel)
  {thePreComp = aPreCompoundModel;}
////////////////////////////////////////////////////////////////////////////////
//
inline G4VPreCompoundModel *G4DPMJET2_5Model::GetPreCompoundModel () const
  {return thePreComp;}
////////////////////////////////////////////////////////////////////////////////
//
inline void G4DPMJET2_5Model::SetVerboseLevel (G4int verboseLevel1)
  {verboseLevel = verboseLevel1;}
////////////////////////////////////////////////////////////////////////////////
//
inline G4double G4DPMJET2_5Model::GetMinEnergy( const G4Material *,
  const G4Element * ) const
  {return theMinEnergy;}
////////////////////////////////////////////////////////////////////////////////
//
inline G4double G4DPMJET2_5Model::GetMaxEnergy( const G4Material *,
  const G4Element * ) const
  {return theMaxEnergy;}
////////////////////////////////////////////////////////////////////////////////
//
// SetVerboseFortranOutput
//
inline G4String G4DPMJET2_5Model::GetVerboseFortranOutput () const
  {return verboseFortranFile;}
////////////////////////////////////////////////////////////////////////////////
//
inline void G4DPMJET2_5Model::DumpVerboseInformation2
  (const G4String particleName, const G4ThreeVector p,
  const G4double E, const G4double T, const G4ThreeVector pinit) const
{
  G4cout <<"Name = " <<particleName <<G4endl;
  G4cout <<"            Momentum          = " <<p/MeV <<" MeV/c" <<G4endl;
  G4cout <<"            T. Energy         = " <<E/MeV <<" MeV"   <<G4endl;
  G4cout <<"            K. Energy         = " <<T/MeV <<" MeV"   <<G4endl;
  if (verboseLevel >= 3)
  {
    G4ThreeVector axis = pinit.unit();
    G4double pz = p.dot(axis);
    G4cout <<"            Transverse mass   = " <<std::sqrt(E*E-pz*pz)/MeV <<" MeV"
           <<G4endl;
    G4cout <<"            Rapidity          = "
           <<0.5*std::log((E+pz)/(E-pz)) <<G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
inline void G4DPMJET2_5Model::DumpVerboseInformation3 (const G4int i,
  const G4int A, const G4int Z, const G4ThreeVector p,
  const G4double E, const G4double T, const G4ThreeVector pinit) const
{
  G4cout <<"----------------------------------------" 
         <<"----------------------------------------" <<G4endl;
  G4cout <<"The nuclear fragment #" <<i <<" before" <<G4endl;
  G4cout <<"----------------------------------------"
         <<"----------------------------------------" <<G4endl;

  std::ostringstream tmpStream;
  tmpStream <<"(A = " <<A <<", Z = " <<Z <<")";
  G4String AZ = tmpStream.str();
  
  DumpVerboseInformation2(AZ, p, E, T, pinit);
}
////////////////////////////////////////////////////////////////////////////////
//
inline void G4DPMJET2_5Model::DumpVerboseInformation4 (const G4int i,
  const G4String particleName, const G4ThreeVector p,
  const G4double E, const G4double T, const G4ThreeVector pinit) const
{
  G4cout <<"----------------------------------------" 
         <<"----------------------------------------" <<G4endl;
  G4cout <<"The nuclear fragment #" <<i <<" after" <<G4endl;
  G4cout <<"----------------------------------------" 
         <<"----------------------------------------" <<G4endl;

  DumpVerboseInformation2(particleName, p, E, T, pinit);
}
////////////////////////////////////////////////////////////////////////////////
//
inline void G4DPMJET2_5Model::SetDPMInitialRandomSeeds (const G4int seed1,
  const G4int seed2)
{
  G4int s1 = seed1;
  G4int s2 = seed2;
  rd2in_ (&s1, &s2);
}
////////////////////////////////////////////////////////////////////////////////
//
inline G4int G4DPMJET2_5Model::GetDPMInitialRandomSeeds (const G4int i)
  const
{
  G4int seed1, seed2;
  rd2out_ (&seed1, &seed2);
  if (i == 0)      return seed1;
  else if (i == 2) return seed2;
  else             return 0;
}
////////////////////////////////////////////////////////////////////////////////
//
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
////////////////////////////////////////////////////////////////////////////////
//
#endif
