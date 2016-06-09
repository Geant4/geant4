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
/// \file hadronic/Hadr02/src/G4DPMJET2_5Model.cc
/// \brief Implementation of the G4DPMJET2_5Model class
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4DPMJET2_5Model.cc
//
// Version:             0.B
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#ifdef G4_USE_DPMJET


#include "G4DPMJET2_5Model.hh"
#include "G4GlaubAADataSetHandler.hh"

#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4FermiBreakUp.hh"
#include "G4PhotonEvaporation.hh"
#include "G4PreCompoundModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include "G4Fragment.hh"
#include "G4VNuclearDensity.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4FermiMomentum.hh"
#include "G4ReactionProductVector.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Poisson.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4LorentzVector.hh"
#include "G4HadTmpUtil.hh"
#include "globals.hh"

#include <fstream>

#include "G4DPMJET2_5Interface.hh"

/////////////////////////////////////////////////////////////////////////////////
//
// Constructor without arguments
//
// This constructor uses a default pre-compound.  It initialises the
// variables (including the de-excitation), but note that much of the work is
// done in the member function Initialise(), which is dedicated to
// initialising variables in DPMJET-II.5 to class-default values.
//
G4DPMJET2_5Model::G4DPMJET2_5Model () : G4HadronicInteraction("DPMJET2_5")
{
//
// Set the default verbose level to 0 - no output.
//
  SetVerboseLevel(1);
//
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout <<"G4DPMJET2_5Model default constructor" <<G4endl;
  }
#endif
//
//
// Send message to stdout to advise that the G4DPMJET2_5Model model is being used.
//
  theInitType = DEFAULT;
  PrintWelcomeMessage();
//
// Use the precompound model for nuclear de-excitation.
//
  theExcitationHandler = 0;
  SetDefaultPreCompoundModel();
//
//
// Set the minimum and maximum range for the model (despite nomanclature, this
// is in energy per nucleon number).
//
  SetMinEnergy(5.0*GeV);
  SetMaxEnergy(1000.0*TeV);
//
//
// Initialise the DPMJET model - this effectively does what the DPMJET
// subroutine DMINIT does without reading an input file.
//
  debug       = false;
  debug_level = 0;
  lunber      = 14;
  dpmver      = 2.5;
  
  LFALSE = 0;
  LTRUE  = 1;

  Initialise ();
//
//
// Next bit directs how and how many Glauber data sets are loaded
// or created.
//
  theGlauberDataSetHandler = G4GlaubAADataSetHandler::getInstance();
  theGlauberDataSetHandler->SetMaxGlauberDataSets(-1);
  
  theParticleTable = G4ParticleTable::GetParticleTable();
  theIonTable      = const_cast <G4IonTable *> (theParticleTable->GetIonTable());
//
//
}
/////////////////////////////////////////////////////////////////////////////////
//
// Constructor with DPMJET-II.5 initialisation type
//
// This constructor uses a default pre-compound.  It initialises the
// variables (including the de-excitation), but note that much of the work is
// done in the member function Initialise(), which is dedicated to
// initialising variables in DPMJET-II.5.  The user is able to define whether to
// use the default values, DPMJET-II.5 settings or DPMJET-III settings.
//
G4DPMJET2_5Model::G4DPMJET2_5Model (const G4DPMJET2_5InitialisationType initType)
{
//
// Set the default verbose level to 0 - no output.
//
  SetVerboseLevel(1);
//
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout <<"G4DPMJET2_5Model constructor #1 " <<G4endl;
  }
#endif
//
//
// Send message to stdout to advise that the G4DPMJET2_5Model model is being used.
//
  theInitType = initType;
  PrintWelcomeMessage();
//
// Use the precompound model for nuclear de-excitation.
//
  theExcitationHandler = 0;
  SetDefaultPreCompoundModel();
//
//
// Set the minimum and maximum range for the model (despite nomanclature, this
// is in energy per nucleon number).
//
  SetMinEnergy(5.0*GeV);
  SetMaxEnergy(1000.0*TeV);
//
//
// Initialise the DPMJET model - this effectively does what the DPMJET
// subroutine DMINIT does without reading an input file.
//
  debug       = false;
  debug_level = 0;
  lunber      = 14;
  dpmver      = 2.5;
  
  LFALSE = 0;
  LTRUE  = 1;

  Initialise ();
//
//
// Next bit directs how and how many Glauber data sets are loaded
// or created.
//
  theGlauberDataSetHandler = G4GlaubAADataSetHandler::getInstance();
  theGlauberDataSetHandler->SetMaxGlauberDataSets (-1);
  
  theParticleTable = G4ParticleTable::GetParticleTable();
  theIonTable      = theParticleTable->GetIonTable();
//
//
}
////////////////////////////////////////////////////////////////////////////////
//
// Constructor with de-excitation handler.
//
// This constructor uses the user-provided de-excitation handler.  However, it
// is possible for the use to provide a NULL pointer, in which case, it is
// assumed that the user doesn't want to simulate de-excitation - USER, BEWARE!
//
// The member function initialises the variables (including the de-excitation),
// but note that much of the work is done in the member function Initialise(),
// which is dedicated to initialising variables in DPMJET-II.5.
//
G4DPMJET2_5Model::G4DPMJET2_5Model (G4ExcitationHandler *aExcitationHandler,
  const G4DPMJET2_5InitialisationType initType) 
{
//
// Set the default verbose level to 0 - no output.
//
  SetVerboseLevel(1);
//
// Send message to stdout to advise that the G4DPMJET2_5Model model is being used.
//
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout <<"G4DPMJET2_5Model constructor #2 " <<G4endl;
  }
#endif
  theInitType = initType;
  PrintWelcomeMessage();
//
// The user is able to provide the excitation handler.
//
  theExcitationHandler = aExcitationHandler;
  thePreComp           = 0;
//
//
// Set the minimum and maximum range for the model (despite nomanclature, this
// is in energy per nucleon number).  
//
  SetMinEnergy(5.0*GeV);
  SetMaxEnergy(1000.0*TeV);
//
//
// Initialise the DPMJET model - this effectively does what the DPMJET
// subroutine DMINIT does without reading an input file.
//
  debug       = false;
  debug_level = 0;
  lunber      = 14;
  dpmver      = 2.5;
  
  LFALSE = 0;
  LTRUE  = 1;
  
  Initialise ();
//
//
// Next bit directs how and how many Glauber data sets are loaded
// or created.
//
  theGlauberDataSetHandler = G4GlaubAADataSetHandler::getInstance();
  theGlauberDataSetHandler->SetMaxGlauberDataSets (-1);

  theParticleTable = G4ParticleTable::GetParticleTable();
  theIonTable      = theParticleTable->GetIonTable();
//
//
}
////////////////////////////////////////////////////////////////////////////////
//
// Constructor with pre-compound model.
//
// This constructor uses the user-provided pre-equilibrium model.  However, it
// is possible for the use to provide a NULL pointer, in which case, it is
// assumed that the user doesn't want to simulate pre-equilibrium. - USER, BEWARE!
//
// The member function initialises the variables (including the de-excitation),
// but note that much of the work is done in the member function Initialise(),
// which is dedicated to initialising variables in DPMJET-II.5.
//
G4DPMJET2_5Model::G4DPMJET2_5Model (G4VPreCompoundModel *aPreComp,
  const G4DPMJET2_5InitialisationType initType) 
{
//
// Set the default verbose level to 0 - no output.
//
  SetVerboseLevel(1);
//
// Send message to stdout to advise that the G4DPMJET2_5Model model is being used.
//
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout <<"G4DPMJET2_5Model constructor #3 " <<G4endl;
  }
#endif
  theInitType = initType;
  PrintWelcomeMessage();
//
// The user is able to provide the pre-compound model.
//
  theExcitationHandler = 0;
  thePreComp           = aPreComp;
//
//
// Set the minimum and maximum range for the model (despite nomanclature, this
// is in energy per nucleon number).  
//
  SetMinEnergy(5.0*GeV);
  SetMaxEnergy(1000.0*TeV);
//
//
// Initialise the DPMJET model - this effectively does what the DPMJET
// subroutine DMINIT does without reading an input file.
//
  debug       = false;
  debug_level = 0;
  lunber      = 14;
  dpmver      = 2.5;
  
  LFALSE = 0;
  LTRUE  = 1;
  
  Initialise ();
//
//
// Next bit directs how and how many Glauber data sets are loaded
// or created.
//
  theGlauberDataSetHandler = G4GlaubAADataSetHandler::getInstance();
  theGlauberDataSetHandler->SetMaxGlauberDataSets (-1);

  theParticleTable = G4ParticleTable::GetParticleTable();
  theIonTable      = theParticleTable->GetIonTable();
//
//
}
////////////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4DPMJET2_5Model::~G4DPMJET2_5Model ()
{
//
//
// The destructor doesn't have to do a great deal!
//
  if (theExcitationHandler) delete theExcitationHandler;
  if (thePreComp)           delete thePreComp;
//  delete theGlauberDataSetHandler;
  theGlauberDataSetHandler->UnloadAllGlauberData();
}
////////////////////////////////////////////////////////////////////////////////
//
// IsApplicable
//
// This member function simply determines whether there is relevant information
// in Glauber data for this projectile and target, and if the nucleus is
// sensible.
//
G4bool G4DPMJET2_5Model::IsApplicable (
  const G4HadProjectile &theTrack, G4Nucleus &theTarget)
{
//
//
// Get relevant information about the projectile and target (A, Z)
//
  const G4ParticleDefinition *definitionP = theTrack.GetDefinition();
  G4int AP   = definitionP->GetBaryonNumber();
  G4int ZP   = G4int(definitionP->GetPDGCharge()/eplus + 0.5);
  G4int AT   = theTarget.GetA_asInt();
  G4int ZT   = theTarget.GetZ_asInt();
  
  if (AP >= 2 && ZP >= 1 && AT >= 2 && ZT >=1) {
    return theGlauberDataSetHandler->IsGlauberDataSetAvailable(AP,AT);
  }
  else {
    return false;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
// ApplyYourself
//
// Member function to process an event, and get information about the products.
// The phases are:
//
// (1) Determine the information about the projectile and target.
//
// (2) Identify to the Glauber data set handler which data need to be used.  If 
// the GDSH finds there are no Glauber profile data for the collision,the 
// product is the unchanged projectile.
//
// (3) Initialise further common-block variables in DPMJET-II.5, and perform
// FORTRAN calls (note this is taken largely from an interpretation of CORSIKA
// and DPMJET-II.5 FORTRAN ... I think there's duplication in some of the
// initialisation on top of what Initialise() does, but CORSIKA has similar
// duplication.  I'm not confident at removing this out at the moment.
//
// (4) Call the DPMJET-II.5 FORTRAN subroutine DPMEVT.
//
// (5) Pick out the final state particles and nuclei.  In the case of nuclei
// use the de-excitation handler if one has been defined.  Transfer all these
// particles to the final-state vector.
//
// (6) if very-verbose output is demanded by the user, there is a print-out
// of the total energy and total momentum before and after collision.
//
G4HadFinalState *G4DPMJET2_5Model::ApplyYourself (
  const G4HadProjectile &theTrack, G4Nucleus &theTarget)
{
//
//
// The secondaries will be returned in G4HadFinalState &theParticleChange -
// initialise this.  The original track will always be discontinued and
// secondaries followed.
//
  theParticleChange.Clear();
  theParticleChange.SetStatusChange(stopAndKill);
//
//
// Get relevant information about the projectile and target (A, Z, energy/nuc,
// momentum, etc).
//
  const G4ParticleDefinition *definitionP = theTrack.GetDefinition();
  G4int AP   = definitionP->GetBaryonNumber();
  G4int ZP   = G4int(definitionP->GetPDGCharge()/eplus+0.5);
  G4double M          = definitionP->GetPDGMass();
  G4ThreeVector pP    = theTrack.Get4Momentum().vect();
  G4double T          = theTrack.GetKineticEnergy()/G4double(AP);   // Units are MeV/nuc 
  G4double E          = theTrack.GetTotalEnergy()/G4double(AP);            // Units are MeV/nuc
  G4int AT         = theTarget.GetA_asInt();
  G4int ZT         = theTarget.GetZ_asInt();
  G4double mpnt  = theTarget.AtomicMass(AT, ZT);
  G4double TotalEPre  = theTrack.GetTotalEnergy() + mpnt;
    //    theTarget.AtomicMass(AT, ZT) + theTarget.GetEnergyDeposit();
//  G4LorentzRotation transformToLab =
//          (const_cast <G4HadProjectile*> (&theTrack))->GetTrafoToLab();
//
//
// Output relevant information on initial conditions if verbose.  Note that
// most of the verbsse output is dealt with through private member function
// calls.
//
  if (verboseLevel >= 2)
  {
    G4cout <<"########################################"
           <<"########################################"
           <<G4endl;
    G4cout.precision(6);
    G4cout <<"IN G4DPMJET2_5Model::ApplyYourself" <<G4endl;
    G4cout <<"START OF EVENT" <<G4endl;
    G4cout <<"Initial projectile (A,Z) = (" <<AP <<", " <<ZP <<")" <<G4endl;
    G4cout <<"Initial target     (A,Z) = (" <<AT <<", " <<ZT <<")" <<G4endl;
    G4cout <<"Projectile momentum      = " <<pP/MeV <<" MeV/c" <<G4endl;
    G4cout <<"Total energy             = " <<TotalEPre/MeV <<" MeV" <<G4endl;
    G4cout <<"Kinetic energy/nuc       = " <<T/MeV <<" MeV" <<G4endl;
  }
//
//
// Setup variables and call the DPMJET model.  There is a significant amount of
// initialisation which still needs to be done, based on DPMJET-II.5 main
// program and card reader subroutine, and the CORSIKA implementation.
//
  G4int AP1 = AP;
  G4int ZP1 = ZP;
  G4int AT1 = AT;
  G4int ZT1 = ZT;
//
//
// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
// ***** WARNING *****
// The following is a provisional "catch" for ions with A==Z.  The de-excitation
// and precompound model can produce such nuclei, although they should decay
// into free protons.  For the moment, DPMJET-II.5 doesn't treat them ... i.e.
// the FORTRAN code would crash.  Therefore, return such ions without 
// nuclear interactions.
//
  if (AP1 > 1 && AP1 == ZP1) {
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(theTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(theTrack.Get4Momentum().vect().unit());
    if (verboseLevel >= 2) {
      G4cout <<"PROJECTILE WITH AP = " <<AP1 <<"  ==  ZP = " <<ZP1 
             <<" REJECTED" <<G4endl;
      G4cout <<"########################################"
             <<"########################################"
             <<G4endl;
    }
    return &theParticleChange;
  }
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  nucc_.ibproj   = 1;                           // IBPROJ = 1
  nucc_.ijproj   = 1;                           // IJPROJ = 1
  collis_.ijprox = 1;                           // IJPROX = 1
  nucc_.ip       = AP1;                         // IP     = IP_P
  nucc_.ipz      = ZP1;                         // IPZ    = IPZ_P
  nucc_.it       = AT1;                         // IT     = IT_P
  nucc_.itz      = ZT1;                         // ITZ    = ITZ_P
  collis_.ijtar  = 1;                           // IJTAR=1
//
//
// Note that epn is deliberately given units of GeV/nuc, because of the units
// used in DPMJET-II.5.
//
  G4double epn   = E / GeV;
  G4double mpn   = M / (GeV*AP);                // Projectile mass per nucleon in GeV/(c2*nuc)
//  G4double gamma = epn / mpn;
//  G4double elab  = epn * M / GeV;             // ELAB   = EPN * AMPRO_P

  diffra_.isingd = ISINGD;
  user2_.isingx  = ISINGX;
  user2_.idubld  = IDUBLD;
  user2_.sdfrac  = SDFRAC;

//  G4double amu_c2GeV = amu_c2 / GeV;
  G4double ppn   = std::sqrt((epn-mpn)*(epn+mpn));
                                                // Units of GeV/(c*nuc)
                                                // PPN    = SQRT( (EPN-AMPROJ)*(EPN+AMPROJ) )
//
//
// Before setting the remainder of the variables for DPMJET-II.5, check for
// appropriate Glauber data.  If the value returned is false, then set the change
// object to the source particle.
//
  if (!(theGlauberDataSetHandler->SetCurrentGlauberDataSet(AP1,AT1,ppn))) {
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(theTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(theTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }
  dtumat_.ntaxx[0]  = AT1;
  dtumat_.nztaxx[0] = ZT1;
  dtumat_.nprxx[0]  = AP1;
  dtumat_.nzprxx[0] = ZP1;
//
//
// Set the remainder of the variables for DPMJET-II.5 FORTRAN.
//
  nncms_.pproj   = ppn;                         // PPROJ  = PPN
  nncms_.eproj   = epn;                         // EPROJ  = EPN
  mpnt          /= (AT * GeV);                  // Mass per nuclon of target in GeV/(c2*nuc)
  nncms_.umo     = std::sqrt(mpn*mpn + mpnt*mpnt + 2.0*mpnt*epn);
                                                // UMO    = SQRT( AMPROJ**2 + AMTAR**2 +2.D0*AMTAR*EPROJ )
                                                // Note I believe this equation is only correct
                                                // if the subsequent equations (for pTthr) 
                                                // needs the Ecm for the NUCLEON-NUCLEON system 
  user2_.cmener  = nncms_.umo;                  // CMENER = UMO
  collis_.s      = nncms_.umo * nncms_.umo;
                                                // SS     = UMO**2
  collis_.ptthr  = 3.0;
  if (strufu_.istrut == 1)
  {
    collis_.ptthr = 2.1 + 0.15*std::pow(std::log10(user2_.cmener/50.),3.0);
                                                // PTTHR  = 2.1D0+0.15D0*(LOG10(CMENER/50.))**3
  }
  else if (strufu_.istrut == 2)
  {
    collis_.ptthr = 2.5 + 0.12*std::pow(std::log10(user2_.cmener/50.),3.0);
                                                // PTTHR  = 2.5D0+0.12D0*(LOG10(CMENER/50.))**3
  }
  collis_.ptthr2 = collis_.ptthr;               // PTTHR2 = PTTHR
  nncms_.gamcm   = (epn + mpnt) / nncms_.umo;
                                                // GAMCM  = (EPROJ+AMTAR)/UMO
                                                // Note I believe this equation is only correct
                                                // if the subsequent equations (for pTthr)
                                                // need the Ecm for the NUCLEON-NUCLEON system 
  nncms_.bgcm    = ppn / nncms_.umo;                // PPROJ/UMO
  nncms_.pcm     = nncms_.gamcm*ppn - nncms_.bgcm*epn;
                                                // PCM    = GAMCM*PPROJ - BGCM*EPROJ
  sigma_.sigsof  = 37.8 * std::pow(collis_.s,0.076);
                                                // ALFA   = 1.076D0
                                                // A      = 37.8D0
                                                // SIGSOF = A * SS**(ALFA-1.D0)
  seasu3_.seasq  = SEASQ;                       // SEASQ  = 0.50D0
  xseadi_.ssmima = SSMIMA;                      // SSMIMA = 1.201D0
  xseadi_.ssmimq = xseadi_.ssmima * xseadi_.ssmima;
                                                // SSMIMQ = SSMIMA**2
  taufo_.taufor  = TAUFOR;
  taufo_.ktauge  = KTAUGE;
  
  if ( theInitType == DEFAULT ) {
    final_.ifinal  = 0;                         // IFINAL = 1
    evappp_.ievap  = 0;                         // IEVAP  = 0
    parevt_.levprt = LTRUE;                     // LEVPRT = .FALSE.
    parevt_.ilvmod = 1;                         // ILVMOD = 1
    parevt_.ldeexg = LFALSE;                    // LDEEXG = .FALSE.
    parevt_.lheavy = LFALSE;                    // LHEAVY = .FALSE.
    frbkcm_.lfrmbk = LFALSE;                    // LFRMBK = .FALSE.
    inpflg_.ifiss  = 0;                         // IFISS  = 0
  } else if ( theInitType == CORSIKA ) {
    final_.ifinal  = 0;                         // IFINAL = 1
    evappp_.ievap  = 0;                         // IEVAP  = 0
    parevt_.levprt = LTRUE;                     // LEVPRT = .FALSE.
    parevt_.ilvmod = 1;                         // ILVMOD = 1
    parevt_.ldeexg = LFALSE;                    // LDEEXG = .FALSE.
    parevt_.lheavy = LTRUE;                     // LHEAVY = .FALSE.
    frbkcm_.lfrmbk = LFALSE;                    // LFRMBK = .FALSE.
    inpflg_.ifiss  = 0;                         // IFISS  = 0
  }
  else if ( theInitType == DPM2_5 ) {
    final_.ifinal  = 0;                         // IFINAL = 0
    evappp_.ievap  = 0;                         // IEVAP  = 0
    parevt_.levprt = LTRUE;                     // LEVPRT = .TRUE. NOTE: THIS IS AT ODDS WITH WHAT'S IN DPMJET-II.5, BUT IF NOT SET, ALL EVENTS GET REJECTED.
    parevt_.ilvmod = 1;                         // ILVMOD = 1
    parevt_.ldeexg = LFALSE;                    // LDEEXG = .FALSE.
    parevt_.lheavy = LFALSE;                    // LHEAVY = .FALSE.
    frbkcm_.lfrmbk = LFALSE;                    // LFRMBK = .FALSE.
    inpflg_.ifiss  = 0;                         // IFISS  = 0
  }
  else if ( theInitType == DPM3 ) {
    final_.ifinal  = 0;                         // IFINAL = 0
    evappp_.ievap  = 0;                         // IEVAP  = 0
    parevt_.levprt = LTRUE;                     // LEVPRT = .TRUE. NOTE: THIS IS AT ODDS WITH WHAT'S IN DPMJET-II.5, BUT IF NOT SET, ALL EVENTS GET REJECTED.
    parevt_.ilvmod = 1;                         // ILVMOD = 1
    parevt_.ldeexg = LFALSE;                    // LDEEXG = .FALSE.
    parevt_.lheavy = LFALSE;                    // LHEAVY = .FALSE.
    frbkcm_.lfrmbk = LFALSE;                    // LFRMBK = .FALSE.
    inpflg_.ifiss  = 0;                         // IFISS  = 0
  }
  xsecpt_.ptcut  = collis_.ptthr;               // PTCUT = PTTHR
  
  G4double dsig1[maxpro+1];
  csj1mi_ (&xsecpt_.ptcut, &dsig1[0]);          // CALL CSJ1MI(PTCUT,DSIG1)
  xsecpt_.dsigh   = dsig1[0];                   // SIG1  = DSIG1(0)
                                                // DSIGH = SIG1
  G4int i         = 0;
  G4double pt     = 0.0;
  samppt_ (&i,&pt);                             // SAMPPT(0,PT)
  collap_.s3      = collis_.s;                  // S3      = SS
  collap_.ijproj1 = collis_.ijprox;             // IJPROJ1 = IJPROX
  collap_.ijtar1  = collis_.ijtar;              // IJTAR1  = IJTAR
  collap_.ptthr1  = collis_.ptthr;              // PTTHR1  = PTTHR
  collap_.iophrd1 = collis_.iophrd;             // IOPHRD1 = IOPHRD
  collap_.ijprlu1 = collis_.ijprlu;             // IJPRLU1 = IJPRLU
  collap_.ijtalu1 = collis_.ijtalu;             // IJTALU1 = IJTALU
  collap_.ptthr3  = collis_.ptthr2;             // PTTHR3  = PTTHR2

  G4int iiipro          = nucc_.ijproj;         // IIPROJ = IJPROJ & IIIPRO = IIPROJ
  G4int iiitar          = nucc_.ijtarg;         // IITARG = IJTARG
  G4int kkmat           = 1;
  G4int nhkkh1          = 1;                    // NHKKH1 = 1
  
  for (i=0; i<8; i++) user1_.projty[i] = paname_.btype[iiipro][i];
                                                // PROJTY=BTYPE(IPROJ)
  G4int irej   = 1;
  G4int evtcnt = 0;
  
  do {
/*    bufueh_.annvv   = 0.001;                   // ANNVV = 0.001
    bufueh_.annss   = 0.001;                    // ANNSS = 0.001
    bufueh_.annsv   = 0.001;                    // ANNSV = 0.001
    bufueh_.annvs   = 0.001;                    // ANNVS = 0.001
    bufueh_.anncc   = 0.001;                    // ANNCC = 0.001
    bufueh_.anndv   = 0.001;                    // ANNDV = 0.001
    bufueh_.annvd   = 0.001;                    // ANNVD = 0.001
    bufueh_.annds   = 0.001;                    // ANNDS = 0.001
    bufueh_.annsd   = 0.001;                    // ANNSD = 0.001
    bufueh_.annhh   = 0.001;                    // ANNHH = 0.001
    bufueh_.annzz   = 0.001;                    // ANNZZ = 0.001
    bufueh_.anndi   = 0.001;                    // ANNDI = 0.001
    bufueh_.annzd   = 0.001;                    // ANNZD = 0.001
    bufueh_.anndz   = 0.001;                    // ANNDZ = 0.001
    bufueh_.ptvv    = 0.0;                      // PTVV = 0.
    bufueh_.ptss    = 0.0;                      // PTSS = 0.
    bufueh_.ptsv    = 0.0;                      // PTSV = 0.
    bufueh_.ptvs    = 0.0;                      // PTVS = 0.
    bufueh_.ptcc    = 0.0;                      // PTCC = 0.
    bufueh_.ptdv    = 0.0;                      // PTDV = 0.
    bufueh_.ptvd    = 0.0;                      // PTVD = 0.
    bufueh_.ptds    = 0.0;                      // PTDS = 0.
    bufueh_.ptsd    = 0.0;                      // PTSD = 0.
    bufueh_.pthh    = 0.0;                      // PTHH = 0.
    bufueh_.ptzz    = 0.0;                      // PTZZ = 0.
    bufueh_.ptdi    = 0.0;                      // PTDI = 0.
    bufueh_.ptzd    = 0.0;                      // PTZD = 0.
    bufueh_.ptdz    = 0.0;                      // PTDZ = 0.
    bufueh_.eevv    = 0.0;                      // EEVV = 0.
    bufueh_.eess    = 0.0;                      // EESS = 0.
    bufueh_.eesv    = 0.0;                      // EESV = 0.
    bufueh_.eevs    = 0.0;                      // EEVS = 0.
    bufueh_.eecc    = 0.0;                      // EECC = 0.
    bufueh_.eedv    = 0.0;                      // EEDV = 0.
    bufueh_.eevd    = 0.0;                      // EEVD = 0.
    bufueh_.eeds    = 0.0;                      // EEDS = 0.
    bufueh_.eesd    = 0.0;                      // EESD = 0.
    bufueh_.eehh    = 0.0;                      // EEHH = 0.
    bufueh_.eezz    = 0.0;                      // EEZZ = 0.
    bufueh_.eedi    = 0.0;                      // EEDI = 0.
    bufueh_.eezd    = 0.0;                      // EEZD = 0.
    bufueh_.eedz    = 0.0;                      // EEDZ = 0.
    ncouch_.acouvv  = 0.0;                      // ACOUVV = 0.
    ncouch_.acouss  = 0.0;                      // ACOUSS = 0.
    ncouch_.acousv  = 0.0;                      // ACOUSV = 0.
    ncouch_.acouvs  = 0.0;                      // ACOUVS = 0.
    ncouch_.acouzz  = 0.0;                      // ACOUZZ = 0.
    ncouch_.acouhh  = 0.0;                      // ACOUHH = 0.
    ncouch_.acouds  = 0.0;                      // ACOUDS = 0.
    ncouch_.acousd  = 0.0;                      // ACOUSD = 0.
    ncouch_.acoudz  = 0.0;                      // ACOUDZ = 0.
    ncouch_.acouzd  = 0.0;                      // ACOUZD = 0.
    ncouch_.acoudi  = 0.0;                      // ACOUDI = 0.
    ncouch_.acoudv  = 0.0;                      // ACOUDV = 0.
    ncouch_.acouvd  = 0.0;                      // ACOUVD = 0.
    ncouch_.acoucc  = 0.0;                      // ACOUCC = 0.*/
    if (evtcnt > 0)
      G4cout <<"REJECTED KKINC EVENT.  RETRY # = " <<evtcnt <<G4endl;
//
//
// Generate an event using DPMJET II-5.  NOTE this is a call to dpmevt, and
// could possibly be replaced by call to kkinc_.  After the call,
// ResetCurrentGlauberDataSet is used to delete any temporary Glauber data 
// set if one had to be set up.
//
    //G4cout << "Call to kkinc_" << G4endl;
    kkinc_ (&epn, &AT1, &ZT1, &AP1, &ZP1, &iiipro, &kkmat, &iiitar, &nhkkh1,
      &irej);
                        //      CALL KKINC(EPN,IIT,IITZ,IIP,IIPZ,IIPROJ,KKMAT,
                        //     * IITARG,NHKKH1,IREJ)
//    dpmevt_ (&elabt, &iiipro, &AP1, &ZP1, &AT1, &ZT1, &kkmat, &nhkkh1);
  } while (irej == 1 && ++evtcnt <100);
  //G4cout << "Call to reset G-data" << G4endl;

  theGlauberDataSetHandler->ResetCurrentGlauberDataSet();
//
//
// If the event has been rejected more than 100 times, then set the track
// as still active and return to the calling routine.
//
  if (irej == 1) {
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(theTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(theTrack.Get4Momentum().vect().unit());
    if (verboseLevel >= 2) {
      G4cout <<"Event rejected and original track maintained" <<G4endl;
      G4cout <<"########################################"
             <<"########################################"
             <<G4endl;
    }
    return &theParticleChange;
  }
//
//
// Determine number of final state particles (including nuclear fragments,
// and load into the particle change if you can identify the particles.
//  
  G4int n              = hkkevt_.nhkk;
  G4int m              = 0;
  G4Fragment *fragment = 0;
  if (verboseLevel >= 2) DumpVerboseInformation1 (n);
//
//
// Now go through each of the secondaries and add to theParticleChange.
//
  for (G4int i=0; i<n; i++)
  {
    if (hkkevt_.isthkk[i]==1 || hkkevt_.isthkk[i]==-1)
    {
//
// Particle is a final state secondary and not a nucleus.
// Determine what this secondary particle is, and if valid, load dynamic
// parameters.
//
      G4ParticleDefinition* theParticle =
        theParticleTable->FindParticle(hkkevt_.idhkk[i]);
      if (theParticle)
      {
        G4double px        = hkkevt_.phkk[i][0] * GeV;
        G4double py        = hkkevt_.phkk[i][1] * GeV;
        G4double pz        = hkkevt_.phkk[i][2] * GeV;
        G4double et        = hkkevt_.phkk[i][3] * GeV;
//        G4LorentzVector lv = transformToLab * G4LorentzVector(px,py,pz,et);
        G4LorentzVector lv = G4LorentzVector(px,py,pz,et);
        
        G4DynamicParticle *theDynamicParticle = 
          new G4DynamicParticle(theParticle,lv);
        theParticleChange.AddSecondary (theDynamicParticle);
        
        if (verboseLevel >= 2)
          DumpVerboseInformation2 (theParticle->GetParticleName(),
          lv.vect(), et, theDynamicParticle->GetKineticEnergy(),pP);
      }
    }
    else if (hkkevt_.idhkk[i]==80000 && hkkevt_.isthkk[i]==1001)
    {
//
//
// Particle is a secondary nucleus. Determine the details of the nuclear
// fragment prior to de-excitation. (Note that the 1 eV in the total energy
// is a safety factor to avoid any possibility of negative rest mass energy.)
// Note also that we don't full trust the energy provided by the DPMJET-II.5,
// and there it's based on the Geant4-determined ion rest-mass.
//
      G4int nucA = extevt_.idres[i];
      G4int nucZ = extevt_.idxres[i];
      if (nucA>0 && nucZ>0) {
        m++;
        fragment           = 0;
        G4double px        = hkkevt_.phkk[i][0] * GeV;
        G4double py        = hkkevt_.phkk[i][1] * GeV;
        G4double pz        = hkkevt_.phkk[i][2] * GeV;
        G4double ionMass   = theIonTable->GetIonMass(nucZ,nucA);
//                         GetIonMass(nucZ,nucA) + nucex[i];  // check how to get this energy
        G4double dpmMass   = hkkevt_.phkk[i][4] * GeV;
        if (dpmMass > ionMass) ionMass = dpmMass;
        G4double et        = std::sqrt(px*px + py*py + pz*pz + ionMass*ionMass);
//        G4LorentzVector lv = transformToLab * G4LorentzVector(px,py,pz,et+1.0*eV);
        G4LorentzVector lv = G4LorentzVector(px,py,pz,et+1.0*eV);
        fragment           = new G4Fragment(nucA, nucZ, lv);
        if (verboseLevel >= 2)
          DumpVerboseInformation3 (m, nucA, nucZ, lv.vect(), et, et-ionMass, pP);
//
//
// Now we can decay the nuclear fragment if present.  The secondaries are
// collected and boosted as well.  The priority is to use a pre-compound
// de-excitation, otherwise the standard excitation-handler is used.
//
        if (fragment && (thePreComp || theExcitationHandler))
        {
          G4ReactionProductVector *products = 0;
          if (thePreComp && fragment->GetA() > 1.5)
            products = thePreComp->DeExcite(*fragment);
          else
            products = theExcitationHandler->BreakItUp(*fragment);
          delete fragment;
          fragment = 0;
          G4ReactionProductVector::iterator iter;
          for (iter = products->begin(); iter != products->end(); ++iter)
          {
            G4DynamicParticle *secondary =
              new G4DynamicParticle((*iter)->GetDefinition(),
              (*iter)->GetTotalEnergy(), (*iter)->GetMomentum());
            theParticleChange.AddSecondary (secondary);
            G4String particleName = (*iter)->GetDefinition()->GetParticleName();
            delete (*iter);
            if (verboseLevel >= 2) {
              if (particleName.find("[",0) < particleName.size())
                DumpVerboseInformation4 (m, particleName, secondary->GetMomentum(),
                secondary->GetTotalEnergy(), secondary->GetKineticEnergy(), pP);
              else DumpVerboseInformation2 (particleName, secondary->GetMomentum(),
                secondary->GetTotalEnergy(), secondary->GetKineticEnergy(), pP);
            }
          }
          delete products;
        }
        if (fragment != 0)
        {
//
//
// Add the excited fragment to the product vector.  Note that this is temporary
// since we should at the atomic excitation to strip all electrons ... i.e. it's
// actually a bare nucleus not an atom in the ground state.
//
          G4ParticleDefinition *theParticleDefinition = theIonTable->
            GetIon(nucZ,nucA);
          G4DynamicParticle *theDynamicParticle = 
            new G4DynamicParticle(theParticleDefinition,lv);
          theParticleChange.AddSecondary (theDynamicParticle);
          delete fragment;
          fragment = 0;
        }
      }
    }
  }
  
  if (verboseLevel >= 3) {
//
//
// Calculate and display the energy and momenta before and after the collision.
// Everything is calculated for the lab frame.
//
    G4double TotalEPost = 0.0;
    G4ThreeVector TotalPPost;
    G4double charge     = 0.0;
    G4int baryon        = 0;
    G4int lepton        = 0;
//    G4int parity        = 0;
    G4int nSecondaries  = theParticleChange.GetNumberOfSecondaries();
    for (G4int j=0; j<nSecondaries; j++) {
      TotalEPost += theParticleChange.GetSecondary(j)->
        GetParticle()->GetTotalEnergy();
      TotalPPost += theParticleChange.GetSecondary(j)->
        GetParticle()->GetMomentum();
      G4ParticleDefinition *theParticle = theParticleChange.GetSecondary(j)->
        GetParticle()->GetDefinition();
      charge += theParticle->GetPDGCharge();
      baryon += theParticle->GetBaryonNumber();
      lepton += theParticle->GetLeptonNumber();
//      parity += theParticle->GetPDGiParity();
    }
    G4cout <<"----------------------------------------"
           <<"----------------------------------------"
           <<G4endl;
    G4cout <<"Total energy before collision   = " <<TotalEPre/MeV
           <<" MeV" <<G4endl;
    G4cout <<"Total energy after collision    = " <<TotalEPost/MeV
           <<" MeV" <<G4endl;
    G4cout <<"Total momentum before collision = " <<pP/MeV
           <<" MeV/c" <<G4endl;
    G4cout <<"Total momentum after collision  = " <<TotalPPost/MeV
           <<" MeV/c" <<G4endl;
    if (verboseLevel >= 4) {
      G4cout <<"Total charge before collision   = " <<(ZP+ZT)*eplus
             <<G4endl;
      G4cout <<"Total charge after collision    = " <<charge
             <<G4endl;
      G4cout <<"Total baryon number before collision = "<<AP+AT
             <<G4endl;
      G4cout <<"Total baryon number after collision  = "<<baryon
             <<G4endl;
      G4cout <<"Total lepton number before collision = 0"
             <<G4endl;
      G4cout <<"Total lepton number after collision  = "<<lepton
             <<G4endl;
    }
    G4cout <<"----------------------------------------"
           <<"----------------------------------------"
           <<G4endl;
  }
  
  if (verboseLevel >= 2)
     G4cout <<"########################################"
            <<"########################################"
            <<G4endl;
  
  return &theParticleChange;
}
////////////////////////////////////////////////////////////////////////////////
//
// SetNoDeexcitation
//
// Deletes an exiting de-excitation handlers and zeros the pointer.  This
// allows the simulation to run without any de-excitation, or just pre-compound
// model only.
//
// Note that you need to separately SetNoPreCompoundModel, if you REALLY don't
// want any nuclear de-excitation.  But please only do this to understand the
// contribution of de-excitation/pre-equilibrium to the full simulation.
// Running without de-excitation or pre-compound is physically unrealistic!
//
//
void G4DPMJET2_5Model::SetNoDeexcitation ()
{
  if (theExcitationHandler)
  {
    delete theExcitationHandler;
    theExcitationHandler = 0;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
// SetNoPreCompoundModel
//
// Deletes an exiting pre-equilibrium model and zeros the pointer.  This
// allows the simulation to run without any pre-compound, or just -de-excitation
// model only.
//
// Note that you need to separately SetNoDeexcitation, if you REALLY don't want
// any nuclear de-excitation.  But please only do this to understand the
// contribution of de-excitation/pre-equilibrium to the full simulation.
// Running without de-excitation or pre-compound is physically unrealistic!
//
//
void G4DPMJET2_5Model::SetNoPreCompoundModel ()
{
  if (thePreComp)
  {
    delete thePreComp;
    thePreComp = 0;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
// SetDefaultDeexcitation
//
// Note that this is used by the default constructor, as well as can be called
// directly by the user.
//
void G4DPMJET2_5Model::SetDefaultDeexcitation ()
{
//  SetNoDeexcitation();
  
  theExcitationHandler               = new G4ExcitationHandler;
  G4Evaporation * theEvaporation     = new G4Evaporation;
  G4FermiBreakUp * theFermiBreakUp   = new G4FermiBreakUp;
  G4PhotonEvaporation* thePhotonEvap = new G4PhotonEvaporation;
  theExcitationHandler->SetEvaporation(theEvaporation);
  theExcitationHandler->SetFermiModel(theFermiBreakUp);
  theExcitationHandler->SetMaxAandZForFermiBreakUp(17, 9);
  theExcitationHandler->SetFermiModel(theFermiBreakUp);
  theExcitationHandler->SetPhotonEvaporation(thePhotonEvap);
}
////////////////////////////////////////////////////////////////////////////////
//
// SetDefaultPreCompoundModel
//
// Note that this is used by the default constructor, as well as can be called
// directly by the user.
//
void G4DPMJET2_5Model::SetDefaultPreCompoundModel ()
{
//  SetNoPreCompoundModel();
  
  G4ExcitationHandler *anExcitationHandler = new G4ExcitationHandler;
  G4Evaporation * theEvaporation           = new G4Evaporation;
  G4FermiBreakUp * theFermiBreakUp         = new G4FermiBreakUp;
  G4PhotonEvaporation* thePhotonEvap       = new G4PhotonEvaporation;
  anExcitationHandler->SetEvaporation(theEvaporation);
  anExcitationHandler->SetFermiModel(theFermiBreakUp);
  anExcitationHandler->SetMaxAandZForFermiBreakUp(17, 9);
  anExcitationHandler->SetPhotonEvaporation(thePhotonEvap);
  
  thePreComp = new G4PreCompoundModel(anExcitationHandler);
}
////////////////////////////////////////////////////////////////////////////////
//
// SetVerboseFortranOutput
//
G4bool G4DPMJET2_5Model::SetVerboseFortranOutput (const G4String filename)
{
  g4dpmjet_close_fort6_ ();
  if (filename == ""       || filename == "stdo" ||
      filename == "stdout" || filename == "std::out" )
  {
    verboseFortranFile = "std::out";
    return true;
  }
  else
  {
    G4int namelen     = filename.length();
    char *ptr         = new char[namelen+1];
    filename.copy(ptr,namelen,0);
    ptr[namelen]      = '\0';
//    char *ptr         = 0;
//    ptr               = const_cast<char*> (filename.c_str());
    ftnlogical opened = LFALSE; 
    g4dpmjet_open_fort6_ (&namelen, &opened, ptr);
    delete [] ptr;
    if (opened == LTRUE) {
      verboseFortranFile = filename;
      return true;
    } else {
      verboseFortranFile = "std::out";
      return false;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
// PrintWelcomeMessage
//
void G4DPMJET2_5Model::PrintWelcomeMessage () const
{
  G4cout <<G4endl;
  G4cout <<" *****************************************************************"
         <<G4endl;
  G4cout <<" Interface to DPMJET2.5 for nuclear-nuclear interactions activated"
         <<G4endl;
  G4cout <<" Version number : 00.00.0B          File date : 23/05/08" <<G4endl;
  G4cout <<" (Interface written by QinetiQ Ltd for the European Space Agency)"
         <<G4endl;
  G4cout <<G4endl;
  G4cout <<" Initialisation of DPMJET-II.5 variables will be according to "
         <<theInitType <<G4endl;
  G4cout <<" *****************************************************************"
         <<G4endl;
  G4cout << G4endl;

  return;
}
////////////////////////////////////////////////////////////////////////////////
//
// DumpVerboseInformation1
//
// Dumps raw information about the DPMJET-II.5 simulation if verbosity set
// to 4 or more.
//
void G4DPMJET2_5Model::DumpVerboseInformation1 (const G4int n) const
{
  G4cout <<"----------------------------------------"
         <<"----------------------------------------" <<G4endl;
  G4cout <<n <<" INTERMEDIATE AND FINAL-STATE SECONDARIES PRODUCED" <<G4endl;
  if (verboseLevel >= 4)
  {
    G4cout <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
           <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" <<G4endl;
    G4cout <<"ORIGINAL DPMJET-II.5 OUTPUT FOR EVENT:" <<G4endl;
    G4cout <<"Note that (1) the particles are yet to be transformed according"
           <<G4endl;
    G4cout <<"              to incident particle direction" <<G4endl;
    G4cout <<"          (2) the units of energy, momentum and mass are GeV,"
           <<G4endl;
    G4cout <<"              GeV/c and GeV/c^2 respectively" <<G4endl;
    G4cout <<"    I"
           <<"    ISTHKK"
           <<"     IDHKK"
           <<"     IDRES"
           <<"    IDXRES"
           <<"             PX"
           <<"             PY"
           <<"             PZ"
           <<"   TOTAL ENERGY"
           <<"           MASS"
           <<G4endl;
    for (G4int i=0; i<n; i++)
    {
      G4cout.unsetf(std::ios::scientific);
      G4cout.setf(std::ios::fixed|std::ios::right|std::ios::adjustfield);
      G4cout.precision(0);
      G4cout <<std::setw(5)  <<i
             <<std::setw(10) <<hkkevt_.isthkk[i]
             <<std::setw(10) <<hkkevt_.idhkk[i]
             <<std::setw(10) <<extevt_.idres[i]
             <<std::setw(10) <<extevt_.idxres[i];
      G4cout.unsetf(std::ios::fixed);
      G4cout.setf(std::ios::scientific|std::ios::right|std::ios::adjustfield);
      G4cout.precision(7);
      G4cout <<std::setw(15) <<hkkevt_.phkk[i][0]
             <<std::setw(15) <<hkkevt_.phkk[i][1]
             <<std::setw(15) <<hkkevt_.phkk[i][2]
             <<std::setw(15) <<hkkevt_.phkk[i][3]
             <<std::setw(15) <<hkkevt_.phkk[i][4]
             <<G4endl;
    }
    G4cout <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
           <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" <<G4endl;
  }
  G4cout.setf(std::ios::fixed);
  G4cout <<" THE FOLLOWING LISTS ONLY THE FINAL-STATE SECONDARIES" <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
// Initialise
// This is intended to do exactly what the main program and subroutine DMINIT
// attempt to do with variable initialisation.  i've included most of the
// FORTRAN source code for reference.  Note that I'm certain that many of these
// lines relate to tallying of events for output by the standalone version
// of DPMJET-II.5, but for the moment, I'd rather initialise those variables,
// just in case not doing so has an adverse effect.
//
void G4DPMJET2_5Model::Initialise ()
{
//
//
// This first line is intended to make sure the block data statements are
// executed, since we're not running from a FORTRAN main program.
//
  g4dpmjet_initialise_block_data_ ();
//
//
  dpar_.aam[4]   = 0.001;                        // AAM(5)=0.001D0
  dpar_.aam[5]   = 0.001;                        // AAM(6)=0.001D0
  dpar_.aam[132] = 0.001;                        // AAM(133)=0.001D0
  dpar_.aam[133] = 0.001;                        // AAM(134)=0.001D0
  dpar_.aam[134] = 0.001;                        // AAM(135)=0.001D0
  dpar_.aam[135] = 0.001;                        // AAM(136)=0.001D0

  vxsvd_.nxsp  = 0;                                // NXSP=0
  vxsvd_.nxst  = 0;                                // NXST=0
  vxsvd_.nxsap = 0;                                // NXSAP=0
  vxsvd_.nxsat = 0;                                // NXSAT=0
  vxsvd_.nxvp  = 0;                                // NXVP=0
  vxsvd_.nxvt  = 0;                                // NXVT=0
  vxsvd_.nxdp  = 0;                                // NXDP=0
  vxsvd_.nxdt  = 0;                                // NXDT=0

  for (G4int i=0; i<50; i++)
  {
    vxsvd_.vxsp[i]  = 1.0E-08;                        // VXSP(II)=1.D-8
    vxsvd_.vxst[i]  = 1.0E-08;                        // VXST(II)=1.D-8
    vxsvd_.vxsap[i] = 1.0E-08;                        // VXSAP(II)=1.D-8
    vxsvd_.vxsat[i] = 1.0E-08;                        // VXSAT(II)=1.D-8
    vxsvd_.vxvp[i]  = 1.0E-08;                        // VXVP(II)=1.D-8
    vxsvd_.vxvt[i]  = 1.0E-08;                        // VXVT(II)=1.D-8
    vxsvd_.vxdp[i]  = 1.0E-08;                        // VXDP(II)=1.D-8
    vxsvd_.vxst[i]  = 1.0E-08;                        // VXST(II)=1.D-8
  }

  if (debug)
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout <<"AT G4DPMJET2_5Model::Initialise:"  <<G4endl;
    }
#endif
    dprin_.ipri      = debug_level;              // IPRI  = LEVLDB
    dprin_.ipev      = debug_level;              // IPEV  = LEVLDB
    dprin_.ippa      = debug_level;              // IPPA  = LEVLDB
    dprin_.ipco      = debug_level;              // IPCO  = LEVLDB
    dprin_.init      = debug_level;              // INIT  = LEVLDB
    dprin_.iphkk     = debug_level;              // IPHKK = LEVLDB
    pydat1_.mstu[25] = 10;                       // MSTU(26) = 10
  }
  else
  {
    dprin_.ipri      = 0;                        // IPRI  = 0
    dprin_.ipev      = 0;                        // IPEV  = 0
    dprin_.ippa      = 0;                        // IPPA  = 0
    dprin_.ipco      =-2;                        // IPCO  = -2
    dprin_.init      = 0;                        // INIT  = 0
    dprin_.iphkk     = 0;                        // IPHKK = 0
    pydat1_.mstu[25] = 0;                        // MSTU(26) = 10
  }

  for (G4int i=0; i<7; i++)
  {
    diqrej_.idiqre[i] = 0;
    diqrej_.idiqrz[i] = 0;
  }

  for (G4int i=0; i<3; i++)
  {
    diqrej_.idvre[i] = 0;
    diqrej_.ivdre[i] = 0;
    diqrej_.idsre[i] = 0;
    diqrej_.isdre[i] = 0;
    diqrej_.idzre[i] = 0;
    diqrej_.izdre[i] = 0;
  }
  
  diqsum_.ndvuu  = 0;                            // NDVUU     = 0
  diqsum_.ndvus  = 0;                            // NDVUS     = 0
  diqsum_.ndvss  = 0;                            // NDVSS     = 0
  diqsum_.nvduu  = 0;                            // NVDUU     = 0
  diqsum_.nvdus  = 0;                            // NVDUS     = 0
  diqsum_.nvdss  = 0;                            // NVDSS     = 0
  diqsum_.ndsuu  = 0;                            // NDSUU     = 0
  diqsum_.ndsus  = 0;                            // NDSUS     = 0
  diqsum_.ndsss  = 0;                            // NDSSS     = 0
  diqsum_.nsduu  = 0;                            // NSDUU     = 0
  diqsum_.nsdus  = 0;                            // NSDUS     = 0
  diqsum_.nsdss  = 0;                            // NSDSS     = 0
  diqsum_.ndzuu  = 0;                            // NDZUU     = 0
  diqsum_.ndzus  = 0;                            // NDZUS     = 0
  diqsum_.ndzss  = 0;                            // NDZSS     = 0
  diqsum_.nzduu  = 0;                            // NZDUU     = 0
  diqsum_.nzdus  = 0;                            // NZDUS     = 0
  diqsum_.nzdss  = 0;                            // NZDSS     = 0
  diqsum_.nadvuu = 0;                            // NADVUU    = 0
  diqsum_.nadvus = 0;                            // NADVUS    = 0
  diqsum_.nadvss = 0;                            // NADVSS    = 0
  diqsum_.navduu = 0;                            // NAVDUU    = 0
  diqsum_.navdus = 0;                            // NAVDUS    = 0
  diqsum_.navdss = 0;                            // NAVDSS    = 0
  diqsum_.nadsuu = 0;                            // NADSUU    = 0
  diqsum_.nadsus = 0;                            // NADSUS    = 0
  diqsum_.nadsss = 0;                            // NADSSS    = 0
  diqsum_.nasduu = 0;                            // NASDUU    = 0
  diqsum_.nasdus = 0;                            // NASDUS    = 0
  diqsum_.nasdss = 0;                            // NASDSS    = 0
  diqsum_.nadzuu = 0;                            // NADZUU    = 0
  diqsum_.nadzus = 0;                            // NADZUS    = 0
  diqsum_.nadzss = 0;                            // NADZSS    = 0
  diqsum_.nazduu = 0;                            // NAZDUU    = 0
  diqsum_.nazdus = 0;                            // NAZDUS    = 0
  diqsum_.nazdss = 0;                            // NAZDSS    = 0
  hdjase_.nhse1  = 0;                            // NHSE1     = 0
  hdjase_.nhse2  = 0;                            // NHSE2     = 0
  hdjase_.nhse3  = 0;                            // NHSE3     = 0
  hdjase_.nhase1 = 0;                            // NHASE1    = 0
  hdjase_.nhase2 = 0;                            // NHASE2    = 0
  hdjase_.nhase3 = 0;                            // NHASE3    = 0
//
//
// Parton pt distribution.
//
  G4int i      = 1;
  G4double pt1 = 0.0;
  G4double pt2 = 0.0;
  G4int ipt    = 0;
  G4int nevt   = 0;
  parpt_ (&i,&pt1,&pt2,&ipt,&nevt);              // CALL PARPT(1,PT1,PT2,IPT,NEVT)
//
//
// Initialise BAMJET, DECAY and HADRIN.
//
  ddatar_ ();                                    // CALL DDATAR
  dhadde_ ();                                    // CALL DHADDE
  dchant_ ();                                    // CALL DCHANT
  dchanh_ ();                                    // CALL DCHANH

  G4double epn = 0.0;
  G4double ppn = 0.0;
  defaul_ (&epn,&ppn);                           // CALL DEFAUL(EPN,PPN)
  defaux_ (&epn,&ppn);                           // CALL DEFAUX(EPN,PPN)

  coulo_.icoul   = 1;                            // ICOUL  = 1
  nuclea_.icoull = 1;                            // ICOULL = 1
  edens_.ieden   = 0;                            // IEDEN = 0
  dprin_.itopd   = 0;                            // ITOPD = 0

  if ( theInitType == DEFAULT ) {
    TAUFOR  = 5.0E+00;
    KTAUGE  = 25;
  } else if ( theInitType == CORSIKA ) {
    TAUFOR  = 5.0E+00;                           // TAUFOR = 5.D0
    KTAUGE  = 25;                                // KTAUGE = 25
  } else if ( theInitType == DPM2_5 ) {
    TAUFOR  = 105.0E+00;                         // TAUFOR = 105.D0
    KTAUGE  = 10;                                // KTAUGE = 10
  } else if ( theInitType == DPM3 ) {
    TAUFOR  = 3.5E+00;
    KTAUGE  = 10;
  }
  taufo_.taufor = TAUFOR;
  taufo_.ktauge = KTAUGE;

  ITAUVE        = 1;                             // ITAUVE = 1
  taufo_.itauve = ITAUVE;
  taufo_.incmod = 1;                             // INCMOD = 1
//
//
// Definition of soft quark distributions, Fermi, Pauli
//
  G4double xseaco = 1.0;                         // XSEACO = 1.00D0
  xseadi_.xseacu  = 1.05 - xseaco;               // XSEACU = 1.05D0-XSEACO

  if ( theInitType == DEFAULT ) {
    UNON           = 3.50;
    UNOM           = 1.11;
    UNOSEA         = 5.0;
    droppt_.fermp  = LTRUE;
    nucimp_.fermod = 0.6;
  } else if ( theInitType == CORSIKA ) {
    UNON           = 3.50;                       // UNON   = 3.50D0
    UNOM           = 1.11;                       // UNOM   = 1.11D0
    UNOSEA         = 5.0;                        // UNOSEA = 5.0D0
    droppt_.fermp  = LTRUE;                      // FERMP  = .TRUE.
    nucimp_.fermod = 0.6;                        // FERMOD = 0.6D0
  } else if ( theInitType == DPM2_5 ) { 
    UNON           = 3.50;                       // UNON   = 3.50D0
    UNOM           = 1.11;                       // UNOM   = 1.11D0
    UNOSEA         = 5.0;                        // UNOSEA = 5.0D0
    droppt_.fermp  = LTRUE;                      // FERMP  = .TRUE.
    nucimp_.fermod = 0.6;                        // FERMOD = 0.6D0
  } else if ( theInitType == DPM3 ) { 
    UNON           = 2.00;
    UNOM           = 1.5;
    UNOSEA         = 5.0;
    droppt_.fermp  = LTRUE;
    nucimp_.fermod = 0.55;
  }
  xseadi_.unon   = UNON;
  xseadi_.unom   = UNOM;
  xseadi_.unosea = UNOSEA;
  
  nuclea_.fermdd = 0.6;                          // FERMDD = 0.6D0
  ferfor_.iferfo = 1;                            // IFERFO = 1
  dprin_.ipaupr  = 0;                            // IPAUPR = 0
  droppt_.lpauli = LTRUE;                        // LPAULI = .TRUE.
//
//
// Definition of cuts for x-sampling.
//
  if ( theInitType == DEFAULT ) {
    CVQ    = 1.8;
    CDQ    = 2.0;
    CSEA   = 0.5;
    SSMIMA = 0.9;
  } else if ( theInitType == CORSIKA ) {
    CVQ    = 1.8;                                // CVQ  = 1.8D0
    CDQ    = 2.0;                                // CDQ  = 2.0D0
    CSEA   = 0.5;                                // CSEA = 0.5D0
    SSMIMA = 0.901;                              // SSMIMA = 0.901D0
  } else if ( theInitType == DPM2_5 ) {
    CVQ    = 1.8;                                // CVQ  = 1.8D0
    CDQ    = 2.0;                                // CDQ  = 2.0D0
    CSEA   = 0.5;                                // CSEA = 0.5D0
    SSMIMA = 1.201;                              // SSMIMA = 1.201D0
  } else if ( theInitType == DPM3 ) {
    CVQ    = 1.0;
    CDQ    = 2.0;
    CSEA   = 0.1;
    SSMIMA = 0.14;
  }
  xseadi_.cvq    = CVQ;
  xseadi_.cdq    = CDQ;
  xseadi_.csea   = CSEA;
  xseadi_.ssmima = SSMIMA;
  
  xseadi_.ssmimq = xseadi_.ssmima * xseadi_.ssmima;
                                                 // SSMIMQ = SSMIMA**2
  if ( theInitType == DEFAULT ) {
    VVMTHR = 0.0;  
  } else if ( theInitType == CORSIKA ) {
    VVMTHR = 0.0;                                // VVMTHR = 0.D0
  } else if ( theInitType == DPM2_5 ) {
    VVMTHR = 0.0;                                // VVMTHR = 0.D0
  } else if ( theInitType == DPM3 ) {
    VVMTHR = 2.0;
  }
  xseadi_.vvmthr = VVMTHR;
//
//
// There is a final call.  Set recombin, seasu3, coninpt, allpart and interdpm
//
  final_.ifinal  = 0;                            // IFINAL = 0
  recom_.irecom  = 0;                            // IRECOM = 0
  seadiq_.lseadi = LTRUE;                        // LSEADI = .TRUE.

  if ( theInitType == DEFAULT ) {
    SEASQ  = 0.5;
    MKCRON = 1;
    CRONCO = 0.64;
  } else if ( theInitType == CORSIKA ) {
    SEASQ  = 0.5;                                // SEASQ  = 0.50D0
    MKCRON = 0;                                  // MKCRON = 0
    CRONCO = 0.0;                                // CRONCO = 0.00D0
  } else if ( theInitType == DPM2_5 ) {
    SEASQ  = 0.5;                                // SEASQ  = 0.50D0
    MKCRON = 1;                                  // MKCRON = 1
    CRONCO = 0.64;                               // CRONCO = 0.64D0
  } else if ( theInitType == DPM3 ) {
    SEASQ  = 1.0;
    MKCRON = 1;
    CRONCO = 0.64;
  }
  seasu3_.seasq  = SEASQ;
  cronin_.mkcron = MKCRON;
  cronin_.cronco = CRONCO;
  
  droppt_.ihada  = LTRUE;                        // IHADA  = .TRUE.
  
//  inxdpm_.intdpm = 0;                          // INTDPM = 0
// Note FORTRAN initialises the variable IROEH  = 0, but this isn't used
// nor is it contained within a common block.
//
//
// Definition for popcork, casadiqu, popcorse
//
  popcck_.pdbck  = 0.0;                         // PDBCK  = 0.D0
  popcck_.ijpock = 0;                           // IJPOCK = 0
  casadi_.icasad = 1;                           // ICASAD = 1
  casadi_.casaxx = 0.05;                        // CASAXX = 0.05D0             ! corrected Nov. 2001
  popcck_.pdbse  = 0.45;                        // PDBSE  = 0.45D0             ! with baryon stopping
  popcck_.pdbseu = 0.45;                        // PDBSEU = 0.45D0             ! with baryon stopping
//                                                C     PDBSE  = 0.D0               ! without baryon stopping
//                                                C     PDBSEU = 0.D0               ! without baryon stopping
  popcck_.irejck = 0;                           // IREJCK = 0
  popcck_.irejse = 0;                           // IREJSE = 0
  popcck_.irejs3 = 0;                           // IREJS3 = 0
  popcck_.irejs0 = 0;                           // IREJS0 = 0
  popcck_.ick4   = 0;                           // ICK4   = 0
  popcck_.ise4   = 0;                           // ISE4   = 0
  popcck_.ise43  = 0;                           // ISE43  = 0
  popcck_.ihad4  = 0;                           // IHAD4  = 0
  popcck_.ick6   = 0;                           // ICK6   = 0
  popcck_.ise6   = 0;                           // ISE6   = 0
  popcck_.ise63  = 0;                           // ISE63  = 0
  popcck_.ihad6  = 0;                           // IHAD6  = 0
  popcck_.irejsa = 0;                           // IREJSA = 0
  popcck_.ireja3 = 0;                           // IREJA3 = 0
  popcck_.ireja0 = 0;                           // IREJA0 = 0
  popcck_.isea4  = 0;                           // ISEA4  = 0
  popcck_.isea43 = 0;                           // ISEA43 = 0
  popcck_.ihada4 = 0;                           // IHADA4 = 0
  popcck_.isea6  = 0;                           // ISEA6  = 0
  popcck_.isea63 = 0;                           // ISEA63 = 0
  popcck_.ihada6 = 0;                           // IHADA6 = 0

  if ( theInitType == DEFAULT || theInitType == CORSIKA ) {
    popcor_.pdb    = 0.1;                       // PDB    = 0.10D0
  } else if ( theInitType == DPM2_5 ) {
    popcor_.pdb    = 0.1;                       // PDB    = 0.10D0
  } else if ( theInitType == DPM3 ) {
    popcor_.pdb    = 0.15;
  }
  popcor_.ajsdef = 0.0;                         // AJSDEF = 0.D0
//
//
// Definition of fluctuat, intpt, hadroniz, diquarks, singlech, evapor
// (Charmed particles set to decay : IHADRINZ>=2
//
  fluctu_.ifluct = 0;                           // IFLUCT = 0
  droppt_.intpt  = LTRUE;                       // INTPT  = .TRUE.
  colle_.ihadrz  = 2;                           // IHADRZ = 2
  ifragm_.ifrag  = 1;                           // IFRAG  = 1
  promu_.ipromu  = 1;                           // IPROMU = 1
  if (colle_.ihadrz >= 2)
  {
    ifragm_.ifrag  = colle_.ihadrz - 1;
                                                // IFRAG = IHADRZ-1
    lundin_ ();                                 // CALL LUNDIN
  }
  diquax_.idiqua  = 1;                          // IDIQUA = 1
  diquax_.idiquu  = 1;                          // IDIQUU = 1
  diquax_.amedd   = 0.9;                        // AMEDD  = 0.9D0
  sincha_.isicha  = 0;                          // ISICHA = 0
  evappp_.ievap   = 0;                          // IEVAP = 0
  seaqxx_.seaqx   = 0.5;                        // SEAQX = 0.5D0
  seaqxx_.seaqxn  = 0.5;                        // SEAQXN = 0.5D0
  kglaub_.jglaub  = 2;                          // JGLAUB = 2
  hadthr_.ehadth  = 5.0;                        // EHADTH = 5.D0
//
//
// Definitions for hbook, pomtable, cmhisto, central, strucfun
//
  hboo_.ihbook    = 1;                          // IHBOOK = 1
  pomtab_.ipomta  = 1;                          // IPOMTA = 1
  cmhico_. cmhis  = 0.0;                        // CMHIS = 0.0D+00           !   Lab System
  zentra_.icentr  = 0;                          // ICENTR = 0
  user2_.istruf   = 222;                        // ISTRUF = 222
  strufu_.istrum  = 0;                          // ISTRUM = 0
  strufu_.istrut  = user2_.istruf / 100;
                                                // ISTRUT = ISTRUF/100
  user2_.istruf   = user2_.istruf - strufu_.istrut*100;
                                                // ISTRUF = ISTRUF-ISTRUT*100
  strufu_.istrum  = user2_.istruf;              // ISTRUM = ISTRUF

  if ( theInitType == DEFAULT ) {
    ISINGD = 1;
    ISINGX = 1;
    IDUBLD = 0;
    SDFRAC = 1.0;
  } else if ( theInitType == CORSIKA ) {
    ISINGD  = 0;                                // ISINGD = 0
    ISINGX  = 0;                                // ISINGX = 0
    IDUBLD  = 0;                                // IDUBLD = 0
    SDFRAC  = 0.0;                              // SDFRAC = 0.
  } else if ( theInitType == DPM2_5 ) {
    ISINGD  = 1;                                // ISINGD = 1
    ISINGX  = 1;                                // ISINGX = 1
    IDUBLD  = 0;                                // IDUBLD = 0
    SDFRAC  = 1.0;                              // SDFRAC = 1.
  } else if ( theInitType == DPM3 ) {
    ISINGD  = 0;
    ISINGX  = 1;
    IDUBLD  = 0;
    SDFRAC  = 1.0;
  }
  diffra_.isingd = ISINGD;
  user2_.isingx  = ISINGX;
  user2_.idubld  = IDUBLD;
  user2_.sdfrac  = SDFRAC;
//
//
// Definitions for start, inforeje, gluxplit, partev, sampt
//  NEVNTS passed to DTMAI as argument, NEVHAD in COMMON
//
  colle_.nfile    = 0;                          // NFILE  = 0
  nstari_.nstart  = 1;                          // NSTART = 1
  stars_.istar2   = 0;                          // ISTAR2 = 0
  stars_.istar3   = 0;                          // ISTAR3 = 0
  user2_.ptlar    = 2.0;                        // PTLAR  = 2.D0

//  G4int iglaub    = 0;                        // IGLAUB = 0 !! Note that this is just a local variable

//  infore_.ifrej   = 0;                        // IFREJ  = 0  Note rejection diagnostics not required
  gluspl_.nugluu  = 1;                          // NUGLUU = 1
  gluspl_.nsgluu  = 0;                          // NSGLUU = 0
  colle_.nvers    = 1;                          // NVERS  = 1
  ptsamp_.isampt  = 4;                          // ISAMPT = 4
//
//
// Definitions for selhard, sigmapom, pshow, secinter.
//
  dropjj_.dropjt  = 0.0;                        // DROPJT = 0.D0
  collis_.iophrd  = 2;                          // IOPHRD = 2
  collis_.ptthr   = 3.0;                        // PTTHR  = 3.D0
//  collis_.ptthr2  = collis_.ptthr;            // PTTHR2 = PTTHR
  user2_.cmener   = 100.0;                      // CMENER = 100.D0
  if (strufu_.istrut == 1)
  {
    collis_.ptthr = 2.1 + 0.15*std::pow(std::log10(user2_.cmener/50.),3.0);
                                                // PTTHR  = 2.1D0+0.15D0*(LOG10(CMENER/50.))**3
  }
  else if (strufu_.istrut == 2)
  {
    collis_.ptthr = 2.5 + 0.12*std::pow(std::log10(user2_.cmener/50.),3.0);
                                                // PTTHR  = 2.5D0+0.12D0*(LOG10(CMENER/50.))**3
  }
  collis_.ptthr2  = collis_.ptthr;              // PTTHR2 = PTTHR
  pomtyp_.ipim    = 2;                          // IPIM   = 2
  pomtyp_.icon    = 48;                         // ICON   = 48
  pomtyp_.isig    = 10;                         // ISIG   = 10
  pomtyp_.lmax    = 30;                         // LMAX   = 30
  pomtyp_.mmax    = 100;                        // MMAX   = 100
  pomtyp_.nmax    = 2;                          // NMAX   = 2
  pomtyp_.difel   = 0.0;                        // DIFEL  = 0.D0
  pomtyp_.difnu   = 1.0;                        // DIFNU  = 1.D0
  pshow_.ipshow   = 1;                          // IPSHOW = 1
  secint_.isecin  = 0;                          // ISECIN = 0
//
//
// This next bit is associated with evaporation.  I'm not sure if it's needed
// as IEVAP = 0, but will initialise in any case.
//
  if ( theInitType == DEFAULT ) {
    parevt_.levprt  = LTRUE;
    parevt_.ilvmod  = 1;
    parevt_.ldeexg  = LFALSE;
    parevt_.lheavy  = LFALSE;
    frbkcm_.lfrmbk  = LFALSE;
    inpflg_.ifiss   = 0;
  } else if ( theInitType == CORSIKA ) {
    parevt_.levprt  = LTRUE;                    // LEVPRT = .TRUE.
    parevt_.ilvmod  = 1;                        // ILVMOD = 1
    parevt_.ldeexg  = LTRUE;                    // LDEEXG = .TRUE.
    parevt_.lheavy  = LTRUE;                    // LHEAVY = .TRUE.
    frbkcm_.lfrmbk  = LTRUE;                    // LFRMBK = .TRUE.
    inpflg_.ifiss   = 0;                        // IFISS  = 0
  } else if ( theInitType == DPM2_5 ) {
    parevt_.levprt  = LFALSE;                   // LEVPRT = .FALSE.
    parevt_.ilvmod  = 1;                        // ILVMOD = 1
    parevt_.ldeexg  = LFALSE;                   // LDEEXG = .FALSE.
    parevt_.lheavy  = LFALSE;                   // LHEAVY = .FALSE.
    frbkcm_.lfrmbk  = LFALSE;                   // LFRMBK = .FALSE.
    inpflg_.ifiss   = 0;                        // IFISS  = 0
  } else if ( theInitType == DPM3 ) {
    parevt_.levprt  = LFALSE;
    parevt_.ilvmod  = 1;
    parevt_.ldeexg  = LFALSE;
    parevt_.lheavy  = LFALSE;
    frbkcm_.lfrmbk  = LFALSE;
    inpflg_.ifiss   = 0;
  }

  hettp_.nbertp   = lunber;                     // NBERTP = LUNBER

  verboseFortranFile = "fort.6";
  G4int namelen      = verboseFortranFile.length();
  char *ptr1         = new char[namelen+1];
  verboseFortranFile.copy(ptr1,namelen,0);
  ptr1[namelen]      = '\0';
//  ptr1               = const_cast<char*> (verboseFortranFile.c_str());
  ftnlogical opened  = LFALSE; 
  g4dpmjet_open_fort6_ (&namelen, &opened, ptr1);
  delete [] ptr1;
  if (opened == LFALSE)
  {
    G4cout <<"ATTEMPTED TO OPEN fort.6 TO OUTPUT VERBOSE FORTRAN TEXT" <<G4endl;
    G4cout <<"HOWEVER THIS WAS NOT POSSIBLE" <<G4endl;
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout <<"AT G4DPMJET2_5Model::Initialise: before NUCLEAR.BIN"  <<G4endl;
    G4cout <<"OPENING NUCLEAR.BIN ON FILE UNIT " <<lunber <<G4endl;
  }
#endif
  if ( !getenv("G4DPMJET2_5DATA") )
  {
    G4cout <<"ENVIRONMENT VARIABLE G4DPMJET2_5DATA NOT SET " <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, 
      "Please setenv G4DPMJET2_5DATA to point to the dpmjet2.5 data files.");
  }
  defaultDirName = G4String(getenv("G4DPMJET2_5DATA")) + "/NUCLEAR.BIN";
  namelen        = defaultDirName.length();
  ptr1           = new char[namelen+1];
  defaultDirName.copy(ptr1,namelen,0);
  ptr1[namelen]  = '\0';
//  ptr1           = const_cast<char*> (defaultDirName.c_str());
  opened         = LFALSE;
  g4dpmjet_open_nuclear_bin_ (&namelen, &lunber, &opened, ptr1);
  delete [] ptr1;
  if (opened == LFALSE)
  {
//
//
// Problems with locating NUCLEAR.BIN file.
//
    G4cout <<"NUCLEAR.BIN FILE NOT FOUND IN DIRECTORY " <<defaultDirName
           <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__,
    "NUCLEAR.BIN file not present.");
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0)
    G4cout <<"AT G4DPMJET2_5Model::Initialise: after NUCLEAR.BIN"  <<G4endl;
#endif
  G4cout << "CALL BERTTP" << G4endl;
  berttp_ ();                                        // CALL BERTTP
  G4cout << "CALL BERTTP done" << G4endl;
  if (evappp_.ievap == 1)
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
      G4cout <<"AT G4DPMJET2_5Model::Initialise: before INCINI"  <<G4endl;
#endif
    incini_ ();                                      // CALL INCINI
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
      G4cout <<"AT G4DPMJET2_5Model::Initialise: after INCINI"  <<G4endl;
#endif
  }
  g4dpmjet_close_nuclear_bin_ (&lunber);
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0)
    G4cout <<"AT G4DPMJET2_5Model::Initialise: NUCLEAR.BIN closed"  <<G4endl;
#endif

  ptlarg_.xsmax = 0.8;                              // XSMAX  = 0.8D0
  dprin_.itopd  = 0;                                // ITOPD  = 0

  G4int iseed1 = 0;
  G4int iseed2 = 0;
  rd2out_ (&iseed1,&iseed2);
//
//
// This next bit outputs important variables if debug is switched on.
//
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout <<"AT G4DPMJET2_5Model::Initialise:" <<G4endl;
    G4cout <<"Printout of important Parameters before DPMJET2.5" <<G4endl;
    G4cout <<"Please note for DPMJET input all numbers are floating point!" <<G4endl;
    G4cout <<"PROJPAR  " <<nucc_.ip       <<" "
                         <<nucc_.ipz      <<G4endl;
    G4cout <<"TARPAR   " <<nucc_.it       <<" "
                         <<nucc_.itz      <<G4endl;
    G4cout <<"MOMENTUM " <<ppn            <<G4endl;
    G4cout <<"ENERGY   " <<epn            <<G4endl;
    G4cout <<"CMENERGY " <<nncms_.umo     <<G4endl;
    G4cout <<"NOFINALE " <<final_.ifinal  <<G4endl;
    G4cout <<"EVAPORAT " <<evappp_.ievap  <<G4endl;
    G4cout <<"OUTLEVEL " <<dprin_.ipri    <<" "
                         <<dprin_.ipev    <<" "
                         <<dprin_.ippa    <<" "
                         <<dprin_.ipco    <<" "
                         <<dprin_.init    <<" "
                         <<dprin_.iphkk   <<G4endl;
    G4cout <<"RANDOMIZ " <<iseed1         <<" "
                         <<iseed2         <<G4endl;
    G4cout <<"STRUCFUN " <<user2_.istruf+100*strufu_.istrut <<G4endl;
    G4cout <<"SAMPT    " <<ptsamp_.isampt <<G4endl;
    G4cout <<"SELHARD  " <<0              <<" "
                         <<collis_.iophrd <<" "
                         <<0              <<" "
                         <<dropjj_.dropjt <<" "
                         <<collis_.ptthr  <<" "
                         <<collis_.ptthr2 << G4endl;
    G4cout <<"SIGMAPOM " <<0              <<" "
                         <<pomtyp_.isig   <<" "
                         <<pomtyp_.ipim + 10*pomtyp_.icon <<" "
                         <<pomtyp_.lmax   <<" "
                         <<pomtyp_.mmax   <<" "
                         <<pomtyp_.nmax   <<G4endl;
    G4cout <<"PSHOWER  " <<pshow_.ipshow  <<G4endl;
    G4cout <<"CENTRAL  " <<zentra_.icentr <<G4endl;
    G4cout <<"CMHISTO  " <<cmhico_.cmhis  <<G4endl;
    G4cout <<"SEASU3   " <<seasu3_.seasq  <<G4endl;
    G4cout <<"RECOMBIN " <<recom_.irecom  <<G4endl;
    G4cout <<"SINGDIFF " <<diffra_.isingd <<G4endl;
    G4cout <<"TAUFOR   " <<taufo_.taufor  <<" "
                         <<taufo_.ktauge  <<" "
                         <<taufo_.itauve  <<G4endl;
    G4cout <<"POPCORN  " <<popcor_.pdb    <<G4endl;
    G4cout <<"POPCORCK " <<popcck_.ijpock <<" "
                         <<popcck_.pdbck  <<G4endl;
    G4cout <<"POPCORSE " <<popcck_.pdbse  <<" "
                         <<popcck_.pdbseu <<G4endl;
    G4cout <<"CASADIQU " <<casadi_.icasad <<" "
                         <<casadi_.casaxx <<G4endl;
    G4cout <<"DIQUARKS " <<diquax_.idiqua <<" "
                         <<diquax_.idiquu <<" "
                         <<diquax_.amedd  <<G4endl;
    G4cout <<"HADRONIZ " <<colle_.ihadrz  <<G4endl;
    G4cout <<"INTPT    " <<droppt_.intpt  <<G4endl;
    G4cout <<"PAULI    " <<droppt_.lpauli <<G4endl;
    G4cout <<"FERMI    " <<droppt_.fermp  <<" "
                         <<nucimp_.fermod <<G4endl;
    G4cout <<"CRONINPT " <<cronin_.mkcron <<" "
                         <<cronin_.cronco <<G4endl;
    G4cout <<"SEADISTR " <<xseadi_.xseacu+0.95 <<" "
                         <<xseadi_.unon   <<" "
                         <<xseadi_.unom   <<" "
                         <<xseadi_.unosea <<G4endl;
    G4cout <<"SEAQUARK " <<seaqxx_.seaqx  <<" "
                         <<seaqxx_.seaqxn <<G4endl;
    G4cout <<"SECINTER " <<secint_.isecin <<G4endl;
    G4cout <<"XCUTS    " <<xseadi_.cvq    <<" "
                         <<xseadi_.cdq    <<" "
                         <<xseadi_.csea   <<" "
                         <<xseadi_.ssmima <<G4endl;
  }
#endif

  bufues_.bnnvv   = 0.001;                      // BNNVV=0.001
  bufues_.bnnss   = 0.001;                      // BNNSS=0.001
  bufues_.bnnsv   = 0.001;                      // BNNSV=0.001
  bufues_.bnnvs   = 0.001;                      // BNNVS=0.001
  bufues_.bnncc   = 0.001;                      // BNNCC=0.001
  bufues_.bnndv   = 0.001;                      // BNNDV=0.001
  bufues_.bnnvd   = 0.001;                      // BNNVD=0.001
  bufues_.bnnds   = 0.001;                      // BNNDS=0.001
  bufues_.bnnsd   = 0.001;                      // BNNSD=0.001
  bufues_.bnnhh   = 0.001;                      // BNNHH=0.001
  bufues_.bnnzz   = 0.001;                      // BNNZZ=0.001
  bufues_.bnndi   = 0.001;                      // BNNDI=0.001
  bufues_.bnnzd   = 0.001;                      // BNNZD=0.001
  bufues_.bnndz   = 0.001;                      // BNNDZ=0.001
  bufues_.bptvv   = 0.0;                        // BPTVV=0.
  bufues_.bptss   = 0.0;                        // BPTSS=0.
  bufues_.bptsv   = 0.0;                        // BPTSV=0.
  bufues_.bptvs   = 0.0;                        // BPTVS=0.
  bufues_.bptcc   = 0.0;                        // BPTCC=0.
  bufues_.bptdv   = 0.0;                        // BPTDV=0.
  bufues_.bptvd   = 0.0;                        // BPTVD=0.
  bufues_.bptds   = 0.0;                        // BPTDS=0.
  bufues_.bptsd   = 0.0;                        // BPTSD=0.
  bufues_.bpthh   = 0.0;                        // BPTHH=0.
  bufues_.bptzz   = 0.0;                        // BPTZZ=0.
  bufues_.bptdi   = 0.0;                        // BPTDI=0.
  bufues_.bptzd   = 0.0;                        // BPTZD=0.
  bufues_.bptdz   = 0.0;                        // BPTDZ=0.
  bufues_.beevv   = 0.0;                        // BEEVV=0.
  bufues_.beess   = 0.0;                        // BEESS=0.
  bufues_.beesv   = 0.0;                        // BEESV=0.
  bufues_.beevs   = 0.0;                        // BEEVS=0.
  bufues_.beecc   = 0.0;                        // BEECC=0.
  bufues_.beedv   = 0.0;                        // BEEDV=0.
  bufues_.beevd   = 0.0;                        // BEEVD=0.
  bufues_.beeds   = 0.0;                        // BEEDS=0.
  bufues_.beesd   = 0.0;                        // BEESD=0.
  bufues_.beehh   = 0.0;                        // BEEHH=0.
  bufues_.beezz   = 0.0;                        // BEEZZ=0.
  bufues_.beedi   = 0.0;                        // BEEDI=0.
  bufues_.beezd   = 0.0;                        // BEEZD=0.
  bufues_.beedz   = 0.0;                        // BEEDZ=0.
  ncoucs_.bcouvv  = 0.0;                        // BCOUVV=0.
  ncoucs_.bcouss  = 0.0;                        // BCOUSS=0.
  ncoucs_.bcousv  = 0.0;                        // BCOUSV=0.
  ncoucs_.bcouvs  = 0.0;                        // BCOUVS=0.
  ncoucs_.bcouzz  = 0.0;                        // BCOUZZ=0.
  ncoucs_.bcouhh  = 0.0;                        // BCOUHH=0.
  ncoucs_.bcouds  = 0.0;                        // BCOUDS=0.
  ncoucs_.bcousd  = 0.0;                        // BCOUSD=0.
  ncoucs_.bcoudz  = 0.0;                        // BCOUDZ=0.
  ncoucs_.bcouzd  = 0.0;                        // BCOUZD=0.
  ncoucs_.bcoudi  = 0.0;                        // BCOUDI=0.
  ncoucs_.bcoudv  = 0.0;                        // BCOUDV=0.
  ncoucs_.bcouvd  = 0.0;                        // BCOUVD=0.
  ncoucs_.bcoucc  = 0.0;                        // BCOUCC=0.
//
//
// Initialisation of
//   ANNVV, ANNSS ... ANNDZ
//   PTVV, PTSS ... PTDZ
//   EEVV, EESS ... EEDZ
//   ACOUVV, ACOUSS, ... ACOUCC
// now all moved to ApplyYourself member function.
//
//  droppt_.ipadis = LFALSE;                        // IPADIS = .FALSE.
//  droppt_.ihadvv = LFALSE;                        // IHADVV = .FALSE.
//  droppt_.ihadsv = LFALSE;                        // IHADSV = .FALSE.
//  droppt_.ihadvs = LFALSE;                        // IHADVS = .FALSE.

  nucc_.ijtarg   = 1;                               // IJTARG=1

//
//
// The following commented out since it seems to have more to do with histogram
// generation.
//
//  i = 1;
//  G4int idummy;
//  distr_ (&i,&nucc_.ijproj,&ppn,&idummy);        // CALL DISTR( 1,IJPROJ,PPN,IDUMMY )
  if (pomtyp_.ipim == 2) {prblm2_ (&user2_.cmener);}
//
//
// Initialise hard scattering & transverse momentum for soft scattering
//
  i = 0;
  jtdtu_ (&i);                                     // CALL JTDTU( 0 )
  i = 0;
  G4double pt;
  samppt_ (&i,&pt);                                // CALL SAMPPT(0,PT)
}
////////////////////////////////////////////////////////////////////////////////
//
#endif
