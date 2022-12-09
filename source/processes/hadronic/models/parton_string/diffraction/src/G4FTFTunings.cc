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
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//---------------------------------------------------------------------------
//
// ClassName:      G4FTFTunings
//
// Author:         2022 Alberto Ribon
//
// Description:    Singleton to keep sets of parameters, called "tunes",
//                 for the FTF model.
//
//                 Please NOTE that, as of now (Fall 2022) ONLY ONE tune
//                 can be selected/applied; attempt to select multiple tunes
//                 will not results in any error messages, however further
//                 down the workflow only the FIRST of the activated tunes
//                 will be used.
//
//                 To use one of the tunes of this class, there is no need to
//                 change anything in this class, and use instead one of the
//                 following two UI commands, before initialization:
//                     /process/had/models/ftf/selectTuneByIndex integerIndex
//                 or  /process/had/models/ftf/selectTuneByName  stringName
//                 for instance:
//                     /process/had/models/ftf/selectTuneByIndex 1
//                     or
//                     /process/had/models/ftf/selectTuneByIndex 2
//                     or
//                     /process/had/models/ftf/selectTuneByIndex 3
//                 or
//                     /process/had/models/ftf/selectTuneByName baryon-tune2022-v0
//                     or
//                     /process/had/models/ftf/selectTuneByName pion-tune2022-v0
//                     or
//                     /process/had/models/ftf/selectTuneByName combined-tune2022-v0
//
//                 If you want to create a new tune, then you need to modify
//                 this class as follows: look for the first "dummy" tune
//                 available; if you find it, then specify its name in the
//                 std::array fNameOfTunes and the values of the parameters
//                 in the methods: G4FTFParamCollection::SetTuneN()
//                                 G4FTFParamCollBaryonProj::SetTuneN()
//                                 G4FTFParamCollMesonProj::SetTuneN()
//                                 G4FTFParamCollPionProj::SetTuneN
//                 Note that you need to set explicitly only the parameters
//                 with non-default values - all the others inherit the
//                 corresponding default values.
//                 If you don't find available "dummy" tune, then you need
//                 to increase by (at least) 1 the number of tunes, and add
//                 the corresponding "SetTuneN()" methods in the 4 classes
//                   G4FTFParamCollection, G4FTFParamCollBaryonProj,
//                   G4FTFParamCollMesonProj, G4FTFParamCollPionProj
//
//                 In order to explore some variations of FTF parameters
//                 (for instance to find out a new tune), please select
//                 (via UI command, as explained above) the existing tune
//                 from which you want to start with as "baseline", and
//                 then set the values of the parameters you want to change
//                 via the following C++ code (to used before initialization):
//                   G4HadronicDeveloperParameters::GetInstance()->Set(...)
//
//                 Note: in its current, first version, of this class,
//                       any FTF tune is applied "globally", i.e. for all
//                       projectile hadrons and regardless of their kinetic
//                       energy.
//                       In future versions, we might try to have tunes that
//                       are meant for specific projectile type and/or for
//                       intervals of kinetic energy (e.g. low-energy,
//                       medium-energy, high-energy).
//
//                 Note: a few classes (written by Julia Yarba) used only in
//                       G4FTFParameters, related to the set of parameters of
//                       the FTF models, have been moved from the header and
//                       source files of the class G4FTFParameters to this
//                       (G4FTFTunings) class, with minimal modifications.
//
// Modified:
//
//----------------------------------------------------------------------------

#include "G4FTFTunings.hh"
#include <CLHEP/Units/PhysicalConstants.h>
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4HadronicDeveloperParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4FTFTuningsMessenger.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"



G4FTFTunings* G4FTFTunings::sInstance = nullptr;

namespace { G4Mutex paramMutex = G4MUTEX_INITIALIZER; }


G4FTFTunings* G4FTFTunings::Instance() {
  if ( sInstance == nullptr ) {
    G4AutoLock l( &paramMutex );
    if ( sInstance == nullptr ) {
      static G4FTFTunings theFTFTuningsObject;
      sInstance = &theFTFTuningsObject;
    }
    l.unlock();
  }
  return sInstance;
}


G4FTFTunings::~G4FTFTunings() {
  delete fMessenger;
}


G4FTFTunings::G4FTFTunings() {
  fMessenger = new G4FTFTuningsMessenger;
}


G4bool G4FTFTunings::IsLocked() const {
  return ( ! G4Threading::IsMasterThread() ||
           G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit );
}


void G4FTFTunings::SetTuneApplicabilityState( const G4int index, const G4int state ) {
  if ( IsLocked() || index < 0 || index >= sNumberOfTunes ) return;
  fApplicabilityOfTunes[index] = state;
}


G4int G4FTFTunings::GetIndexTune( const G4ParticleDefinition* /* particleDef */, const G4double /* ekin */ ) const {
  // For the time being, select the first alternative (i.e. with index > 0) tune which is switched on.
  // If nothing is found, then returns 0 (which corresponds to the default set of parameters).
  G4int indexTune = 0;
  for ( G4int i = 1; i < sNumberOfTunes; ++i ) {
    if ( GetTuneApplicabilityState(i) != 0 ) {  // tune i-th is switched on
      indexTune = i;
      break;
    }
  }
  
  /* For the future
  G4int pdgCode = std::abs( particleDef->GetPDGEncoding() );
  if ( particleDef != nullptr  &&  ekin >= 0.0  &&  pdgCode != 0 ) {
    G4bool isLowEnergy = ( ekin < fLowEnergyThreshold );
    G4bool isHighEnergy = ( ekin > fHighEnergyThreshold );
    G4bool isMediumEnergy = ( ( ! isLowEnergy ) && ( ! isHighEnergy ) );
    G4bool isMeson = ( pdgCode < 1000 );
    G4bool isPion = ( pdgCode == 211 || pdgCode == 111 );
    G4bool isKaon = ( pdgCode == 321 || pdgCode == 311 || pdgCode == 130 || pdgCode == 310 );
    G4bool isBaryon = ( pdgCode > 1000 );
    G4bool isNucleon = ( pdgCode == 2212 || pdgCode == 2112 );
    G4bool isAntiBaryon = particleDef->GetBaryonNumber() < 0;
    // Based on the projectile type, its kinetic energy, and the "applicability" flag
    // of each tune, find the right tune to be applicable in this case.
    // ...
  }
  */
  
  //G4cout << "G4FTFTunings::GetIndexTune : projectile=" << particleDef->GetParticleName()
  //       << "  ekin[MeV]=" << ekin << " -> indexTune=" << indexTune
  //       << "  " << fNameOfTunes[indexTune] << G4endl;
  
  return indexTune;
}


//============================================================================

G4HadronicDeveloperParameters& HDP = G4HadronicDeveloperParameters::GetInstance();


class G4FTFSettingDefaultHDP {
  public:   
    G4FTFSettingDefaultHDP() {
      // Cross sections for elementary processes
      //
      // these are for Inelastic interactions, i.e. Xinelastic=(Xtotal-Xelastix)>0.
      // for elastic, all the A's & B's, Atop & Ymin are zeros
      // general formula: Pp = A1*exp(B1*Y) + A2*exp(B2*Y) + A3
      // but if Y<Ymin, then Pp=max(0.,Atop)
      // for details, see also G4FTFParameters::GetProcProb( ProcN, y )
      //
      // baryons
      /* JVY, Oct. 31, 2017: Per Alberto R. & Vladimir U., keep this group of parameters FIXED */
      /* JVY, June 11, 2020: try to open up... */
      // Process=0 --> Qexchg w/o excitation
      HDP.SetDefault( "FTF_BARYON_PROC0_A1",  13.71 );
      HDP.SetDefault( "FTF_BARYON_PROC0_B1",   1.75 );
      HDP.SetDefault( "FTF_BARYON_PROC0_A2", -30.69 ); 
      HDP.SetDefault( "FTF_BARYON_PROC0_B2",   3.0  ); 
      HDP.SetDefault( "FTF_BARYON_PROC0_A3",   0.0  );
      HDP.SetDefault( "FTF_BARYON_PROC0_ATOP", 1.0  ); 
      HDP.SetDefault( "FTF_BARYON_PROC0_YMIN", 0.93 ); 
      // Process=1 --> Qexchg w/excitation
      HDP.SetDefault( "FTF_BARYON_PROC1_A1",  25.0  );
      HDP.SetDefault( "FTF_BARYON_PROC1_B1",   1.0  );
      HDP.SetDefault( "FTF_BARYON_PROC1_A2", -50.34 );
      HDP.SetDefault( "FTF_BARYON_PROC1_B2",   1.5  );
      HDP.SetDefault( "FTF_BARYON_PROC1_A3",   0.0  );
      HDP.SetDefault( "FTF_BARYON_PROC1_ATOP", 0.0  );
      HDP.SetDefault( "FTF_BARYON_PROC1_YMIN", 1.4  );
      // NOTE: Process #2 & 3 are projectile & target diffraction
      //       they have more complex definition of A1 & A2 
      //      (see example below)
      // SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0, 0.93 );  // Projectile diffraction
      // SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0, 0.93 );  // Target diffraction
      // Also, for ( AbsProjectileBaryonNumber > 10 ||  NumberOfTargetNucleons > 10 )
      // projectile and/or target diffraction (dissociation) may be switched ON/OFF 
      HDP.SetDefault( "FTF_BARYON_DIFF_DISSO_PROJ", false );
      HDP.SetDefault( "FTF_BARYON_DIFF_DISSO_TGT",  false );
      /* JVY, Oct. 31, 2017: Per Alberto R. & Vladimir U., keep this group of parameters FIXED */
      /* JVY, June 11, 2020: try to open up... */
      // Process=4 --> Qexchg w/additional multiplier in excitation 
      HDP.SetDefault( "FTF_BARYON_PROC4_A1",   0.6 ); 
      HDP.SetDefault( "FTF_BARYON_PROC4_B1",   0.0 );
      HDP.SetDefault( "FTF_BARYON_PROC4_A2",  -1.2 ); 
      HDP.SetDefault( "FTF_BARYON_PROC4_B2",   0.5 );
      HDP.SetDefault( "FTF_BARYON_PROC4_A3",   0.0 );
      HDP.SetDefault( "FTF_BARYON_PROC4_ATOP", 0.0 );
      HDP.SetDefault( "FTF_BARYON_PROC4_YMIN", 1.4 );
      // Parameters of participating hadron (baryon) excitation
      HDP.SetDefault( "FTF_BARYON_DELTA_PROB_QEXCHG", 0.0 );
      HDP.SetDefault( "FTF_BARYON_PROB_SAME_QEXCHG", 0.0 );
      HDP.SetDefault( "FTF_BARYON_DIFF_M_PROJ", 1.16, 1.16, 3.0 );     // it's supposed to be in GeV but do NOT do (*CLHEP::GeV) 
                                                                       // because it'll be done in the G4FTFParameters::SetProjMinDiffMass
      HDP.SetDefault( "FTF_BARYON_NONDIFF_M_PROJ", 1.16, 1.16, 3.0 );  // do NOT (*CLHEP::GeV) - same as above
      HDP.SetDefault( "FTF_BARYON_DIFF_M_TGT", 1.16, 1.16, 3.0 );      // do NOT (*CLHEP::GeV) - same as above
      HDP.SetDefault( "FTF_BARYON_NONDIFF_M_TGT", 1.16, 1.16, 3.0 );   // do NOT (*CLHEP::GeV) - same as above
      HDP.SetDefault( "FTF_BARYON_AVRG_PT2", 0.3, 0.08, 1.0 );         // do NOT (*CLHEP::GeV*CLHEP::GeV) 
      // JVY, Oct. 6, 2017: Per Alberto R., keep these two settings fixed (for now)
      // HDP.SetDefault( "FTF_BARYON_PROB_DISTR_PROJ", 0.3 ); 
      // HDP.SetDefault( "FTF_BARYON_PROB_DISTR_TGT", 0.3 );       
      // pions
      // JVY, Aug.8, 2018 --> Feb.14, 2019 --> June 25, 2019: 
      // Parameters of participating hadron (pions) excitation
      /* JVY, June 25, 2019: For now, keep this group of parameters FIXED */
      // Process=0 --> Qexchg w/o excitation
      HDP.SetDefault( "FTF_PION_PROC0_A1",  150.0  );
      HDP.SetDefault( "FTF_PION_PROC0_B1",    1.8  );
      HDP.SetDefault( "FTF_PION_PROC0_A2", -247.3  ); 
      HDP.SetDefault( "FTF_PION_PROC0_B2",    2.3  ); 
      HDP.SetDefault( "FTF_PION_PROC0_A3",    0.0  );
      HDP.SetDefault( "FTF_PION_PROC0_ATOP",  1.0  ); 
      HDP.SetDefault( "FTF_PION_PROC0_YMIN",  2.3  ); 
      // Process=1 --> Qexchg w/excitation
      HDP.SetDefault( "FTF_PION_PROC1_A1",   5.77 );
      HDP.SetDefault( "FTF_PION_PROC1_B1",   0.6  );
      HDP.SetDefault( "FTF_PION_PROC1_A2",  -5.77 );
      HDP.SetDefault( "FTF_PION_PROC1_B2",   0.8  );
      HDP.SetDefault( "FTF_PION_PROC1_A3",   0.0  );
      HDP.SetDefault( "FTF_PION_PROC1_ATOP", 0.0  );
      HDP.SetDefault( "FTF_PION_PROC1_YMIN", 0.0  );
      /*
      // NOTE: Process #2 & 3 are projectile & target diffraction
      // Process=2 --> Projectile diffraction
      // Q: Would it even make sense to make these configurable ?
      //    The following is hadrcoded:
      //    Projectile Baryon Number > 10 (AbsProjectileBaryonNumber > 10)
      //    ... which is "strange" because projectile is a pion !!!... so it's always OFF    
      //    (see also lines 1007-1016)
      HDP.SetDefault( "FTF_PION_PROC2_A1",      2.27 );
      HDP.SetDefault( "FTF_PION_PROC2_B1",      0.5  );
      HDP.SetDefault( "FTF_PION_PROC2_A2", -98052.0  );
      HDP.SetDefault( "FTF_PION_PROC2_B2",      4.0  );
      HDP.SetDefault( "FTF_PION_PROC2_A3",      0.0  );
      HDP.SetDefault( "FTF_PION_PROC2_ATOP",    0.0  );
      HDP.SetDefault( "FTF_PION_PROC2_YMIN",    3.0  );
      */
      // Process=3 --> Target diffraction
      HDP.SetDefault( "FTF_PION_PROC3_A1",    7.0  );
      HDP.SetDefault( "FTF_PION_PROC3_B1",    0.9  );
      HDP.SetDefault( "FTF_PION_PROC3_A2",  -85.28 );
      HDP.SetDefault( "FTF_PION_PROC3_B2",    1.9  );
      HDP.SetDefault( "FTF_PION_PROC3_A3",    0.08 );
      HDP.SetDefault( "FTF_PION_PROC3_ATOP",  0.0  );
      HDP.SetDefault( "FTF_PION_PROC3_YMIN",  2.2  );
      // projectile and/or target diffraction (dissociation) may be switched ON/OFF 
      // NOTE: Both projectile and target diffraction are turned OFF if
      // a) Number of Target Nucleons > 10 (NumberOfTargetNucleons > 10)
      //    OR
      // b) Projectile Baryon Number > 10 (AbsProjectileBaryonNumber > 10)
      //    ... which is "strange" because projectile is a pion !!!... so it's always OFF
      HDP.SetDefault( "FTF_PION_DIFF_DISSO_PROJ", false );
      HDP.SetDefault( "FTF_PION_DIFF_DISSO_TGT",  false ); 
      /* JVY, June 25, 2019: For now keep this group of parameters FIXED */
      /* JVY, June 11, 2020: try to open up... */
      // Process=4 --> Qexchg w/additional multiplier in excitation 
      HDP.SetDefault( "FTF_PION_PROC4_A1",   1.0  ); 
      HDP.SetDefault( "FTF_PION_PROC4_B1",   0.0  );
      HDP.SetDefault( "FTF_PION_PROC4_A2", -11.02 ); 
      HDP.SetDefault( "FTF_PION_PROC4_B2",   1.0  );
      HDP.SetDefault( "FTF_PION_PROC4_A3",   0.0  );
      HDP.SetDefault( "FTF_PION_PROC4_ATOP", 0.0  );
      HDP.SetDefault( "FTF_PION_PROC4_YMIN", 2.4  );
      //    
      HDP.SetDefault( "FTF_PION_DELTA_PROB_QEXCHG", 0.56 );
      HDP.SetDefault( "FTF_PION_DIFF_M_PROJ", 1.0, 0.5, 3.0 );
      HDP.SetDefault( "FTF_PION_NONDIFF_M_PROJ", 1.0, 0.5, 3.0 );
      HDP.SetDefault( "FTF_PION_DIFF_M_TGT", 1.16, 1.16, 3.0 );   // All (NON)DIFF_M's are supposed to be in GeV but do NOT do (*CLHEP::GeV) 
                                                                  // because it'll be done in the G4FTFParameters::SetProjMinDiffMass
      HDP.SetDefault( "FTF_PION_NONDIFF_M_TGT", 1.16, 1.16, 3.0 );
      HDP.SetDefault( "FTF_PION_AVRG_PT2", 0.3, 0.08, 1.0 );      //  do NOT (*CLHEP::GeV*CLHEP::GeV)      
      // nuclear destruction 
      // NOTE: Settings of most of these parameters are the same
      //       for different types of projectile hadron
      //       However, we decided to introduce separate variables
      //       and configuration cards for each type of projectile
      // baryons
      // projectile destruction
      HDP.SetDefault( "FTF_BARYON_NUCDESTR_P1_PROJ", 1.0, 0.0, 1.0 ); // in principle, it should be 1./NBRN - FIXME later !
      HDP.SetDefault( "FTF_BARYON_NUCDESTR_P1_NBRN_PROJ", false );
      // for now, keep fixed p2 & p3 for the proj destruction
      // they're defined explicitly in G4FTFParamCollection ctor
      // target destruction
      HDP.SetDefault( "FTF_BARYON_NUCDESTR_P1_TGT", 1.0, 0.0, 1.0 );   
      HDP.SetDefault( "FTF_BARYON_NUCDESTR_P1_ADEP_TGT", false );          
      HDP.SetDefault( "FTF_BARYON_NUCDESTR_P2_TGT", 4.0, 2.0, 16.0 );
      HDP.SetDefault( "FTF_BARYON_NUCDESTR_P3_TGT", 2.1, 0.0, 4.0 );
      //
      HDP.SetDefault( "FTF_BARYON_PT2_NUCDESTR_P1", 0.035, 0.0, 0.25 ); 
      HDP.SetDefault( "FTF_BARYON_PT2_NUCDESTR_P2", 0.04, 0.0, 0.25 ); 
      HDP.SetDefault( "FTF_BARYON_PT2_NUCDESTR_P3", 4.0, 2.0, 16.0 ); 
      HDP.SetDefault( "FTF_BARYON_PT2_NUCDESTR_P4", 2.5, 0.0, 4.0 ); 
      //
      HDP.SetDefault( "FTF_BARYON_NUCDESTR_R2", 1.5*CLHEP::fermi*CLHEP::fermi, 0.5*CLHEP::fermi*CLHEP::fermi, 2.0*CLHEP::fermi*CLHEP::fermi  );
      HDP.SetDefault( "FTF_BARYON_EXCI_E_PER_WNDNUCLN", 40.0*CLHEP::MeV, 0.0, 100.0*CLHEP::MeV );
      HDP.SetDefault( "FTF_BARYON_NUCDESTR_DISP", 0.3, 0.1, 0.4 );
      // JVY, Oct. 6, 2017: Per Alberto R., this is just a technical parameter,
      //                    and it should NOT be changed
      // HDP.SetDefault( "FTF_BARYON_NUCDESTR_MAXPT2", 1. * CLHEP::GeV*CLHEP::GeV  ); 	 
      // mesons - these parameters are common for pions, kaons, etc. (per original code)
      // NOTE: *NO* projectile destruction for mesons !!!
      // target destruction
      HDP.SetDefault( "FTF_MESON_NUCDESTR_P1_TGT", 0.00481, 0.0, 1.0 );    
      HDP.SetDefault( "FTF_MESON_NUCDESTR_P1_ADEP_TGT", true );           
      HDP.SetDefault( "FTF_MESON_NUCDESTR_P2_TGT", 4.0, 2.0, 16.0 );
      HDP.SetDefault( "FTF_MESON_NUCDESTR_P3_TGT", 2.1, 0.0, 4.0 );
      //
      HDP.SetDefault( "FTF_MESON_PT2_NUCDESTR_P1", 0.035, 0.0, 0.25 ); 
      HDP.SetDefault( "FTF_MESON_PT2_NUCDESTR_P2", 0.04, 0.0, 0.25 ); 
      HDP.SetDefault( "FTF_MESON_PT2_NUCDESTR_P3", 4.0, 2.0, 16.0 ); 
      HDP.SetDefault( "FTF_MESON_PT2_NUCDESTR_P4", 2.5, 0.0, 4.0 ); 
      //
      HDP.SetDefault( "FTF_MESON_NUCDESTR_R2", 1.5*CLHEP::fermi*CLHEP::fermi, 
                                               0.5*CLHEP::fermi*CLHEP::fermi, 
      					     2.0*CLHEP::fermi*CLHEP::fermi );
      HDP.SetDefault( "FTF_MESON_EXCI_E_PER_WNDNUCLN", 40.0*CLHEP::MeV, 0.0, 100.0*CLHEP::MeV );
      HDP.SetDefault( "FTF_MESON_NUCDESTR_DISP", 0.3, 0.1, 0.4 );
    }
};


G4FTFSettingDefaultHDP FTFDefaultsHDP;  


G4FTFParamCollection::G4FTFParamCollection() {
  // zero out everything
  // parameters of excitation
  // Proc=0 --> Qexchg w/o excitation
  fProc0A1 = 0.0;
  fProc0B1 = 0.0;
  fProc0A2 = 0.0;
  fProc0B2 = 0.0;
  fProc0A3 = 0.0;
  fProc0Atop = 0.0;
  fProc0Ymin = 0.0;
  // Proc=1 --> Qexchg w/excitation
  fProc1A1 = 0.0;
  fProc1B1 = 0.0;
  fProc1A2 = 0.0;
  fProc1B2 = 0.0;
  fProc1A3 = 0.0;
  fProc1Atop = 0.0;
  fProc1Ymin = 0.0;
  //
  fProjDiffDissociation = false;
  fTgtDiffDissociation = false;
  // Proc=2 --> Projectile diffraction
  fProc2A1 = 0.0; 
  fProc2B1 = 0.0; 
  fProc2A2 = 0.0; 
  fProc2B2 = 0.0; 
  fProc2A3 = 0.0; 
  fProc2Atop = 0.0; 
  fProc2Ymin = 0.0;
  // Proc=3 --> Target diffraction
  fProc3A1 = 0.0; 
  fProc3B1 = 0.0; 
  fProc3A2 = 0.0; 
  fProc3B2 = 0.0; 
  fProc3A3 = 0.0; 
  fProc3Atop = 0.0; 
  fProc3Ymin = 0.0; 
  // Proc=4 --> Qexchg w/additional multiplier in excitation  
  fProc4A1 = 0.0;
  fProc4B1 = 0.0;
  fProc4A2 = 0.0;
  fProc4B2 = 0.0; 
  fProc4A3 = 0.0;
  fProc4Atop = 0.0;
  fProc4Ymin = 0.0;
  // parameters of participating baryon excitation
  fDeltaProbAtQuarkExchange = 0.0; 
  fProbOfSameQuarkExchange = 0.0;
  fProjMinDiffMass = 0.0;
  fProjMinNonDiffMass = 0.0;
  fTgtMinDiffMass = 0.0;
  fTgtMinNonDiffMass = 0.0;
  fAveragePt2 = 0.0;
  fProbLogDistrPrD = 0.0;
  fProbLogDistr = 0.0;
  // parameters of nuclear distruction
  // COMMONs
  fNuclearProjDestructP1 = 0.0;
  fNuclearProjDestructP1_NBRNDEP = false;
  fNuclearTgtDestructP1 = 0.0;
  fNuclearTgtDestructP1_ADEP = false;
  fNuclearProjDestructP2 = 0.0;
  fNuclearProjDestructP3 = 0.0;
  fNuclearTgtDestructP2 = 0.0;
  fNuclearTgtDestructP3 = 0.0;
  fPt2NuclearDestructP1 = 0.0;
  fPt2NuclearDestructP2 = 0.0;
  fPt2NuclearDestructP3 = 0.0;
  fPt2NuclearDestructP4 = 0.0;
  // baryons
  fR2ofNuclearDestruct = 0.0;
  fExciEnergyPerWoundedNucleon = 0.0;
  fDofNuclearDestruct = 0.0;
  fMaxPt2ofNuclearDestruct = 0.0;
  // keep the 2 parameters below fixed for now (i.e. do not take them from HDP)
  fNuclearProjDestructP2 = 4.0;
  fNuclearProjDestructP3 = 2.1;  
}


G4FTFParamCollBaryonProj::G4FTFParamCollBaryonProj() : G4FTFParamCollection() {
  // parameters of participating hadron (baryon) excitation
  // baryons projectile
  // Proc=0 --> Qexchg w/o excitation
  HDP.DeveloperGet( "FTF_BARYON_PROC0_A1",   fProc0A1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC0_B1",   fProc0B1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC0_A2",   fProc0A2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC0_B2",   fProc0B2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC0_A3",   fProc0A3 );
  HDP.DeveloperGet( "FTF_BARYON_PROC0_ATOP", fProc0Atop ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC0_YMIN", fProc0Ymin );
  // Proc=1 --> Qexchg w/excitation
  HDP.DeveloperGet( "FTF_BARYON_PROC1_A1",   fProc1A1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC1_B1",   fProc1B1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC1_A2",   fProc1A2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC1_B2",   fProc1B2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC1_A3",   fProc1A3 );
  HDP.DeveloperGet( "FTF_BARYON_PROC1_ATOP", fProc1Atop ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC1_YMIN", fProc1Ymin ); 
  // Proc=2 & Proc=3 for the case ( AbsProjectileBaryonNumber > 10 ||  NumberOfTargetNucleons > 10 )
  // (diffraction dissociation)
  // NOTE-1: used to be ( AbsProjectileBaryonNumber > 1 ||  NumberOfTargetNucleons > 1 )...
  // NOTE-2: As of 10.5, both are set to false (via HDP)
  HDP.DeveloperGet( "FTF_BARYON_DIFF_DISSO_PROJ", fProjDiffDissociation );
  HDP.DeveloperGet( "FTF_BARYON_DIFF_DISSO_TGT",  fTgtDiffDissociation );
  // Proc=4 --> Qexchg "w/additional multiplier" in excitation 
  HDP.DeveloperGet( "FTF_BARYON_PROC4_A1",   fProc4A1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC4_B1",   fProc4B1 );
  HDP.DeveloperGet( "FTF_BARYON_PROC4_A2",   fProc4A2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC4_B2",   fProc4B2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC4_A3",   fProc4A3 );
  HDP.DeveloperGet( "FTF_BARYON_PROC4_ATOP", fProc4Atop ); 
  HDP.DeveloperGet( "FTF_BARYON_PROC4_YMIN", fProc4Ymin );
  //
  HDP.DeveloperGet( "FTF_BARYON_DELTA_PROB_QEXCHG", fDeltaProbAtQuarkExchange );
  HDP.DeveloperGet( "FTF_BARYON_PROB_SAME_QEXCHG", fProbOfSameQuarkExchange );
  HDP.DeveloperGet( "FTF_BARYON_DIFF_M_PROJ", fProjMinDiffMass );
  HDP.DeveloperGet( "FTF_BARYON_NONDIFF_M_PROJ", fProjMinNonDiffMass );
  HDP.DeveloperGet( "FTF_BARYON_DIFF_M_TGT", fTgtMinDiffMass );
  HDP.DeveloperGet( "FTF_BARYON_NONDIFF_M_TGT", fTgtMinNonDiffMass );
  HDP.DeveloperGet( "FTF_BARYON_AVRG_PT2", fAveragePt2 );
  // JVY - Per Alberto R., we're curretly keeping these two settings fixed,
  // thus they're defined here explicitly, rather than via HDP
  // HDP.DeveloperGet( "FTF_BARYON_PROB_DISTR_PROJ", fProbLogDistrPrD );
  // HDP.DeveloperGet( "FTF_BARYON_PROB_DISTR_TGT", fProbLogDistr );
  fProbLogDistrPrD = 0.55; 
  fProbLogDistr    = 0.55;   
  // nuclear destruction
  // ---> LATER !!! ---> fBaryonMaxNumberOfCollisions = 2.;
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P1_PROJ", fNuclearProjDestructP1 );
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P1_NBRN_PROJ",fNuclearProjDestructP1_NBRNDEP );
  // keep the 2 parameters below fixed for now (i.e. do not take them from HDP)
  fNuclearProjDestructP2 = 4.0;
  fNuclearProjDestructP3 = 2.1;
  //
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P1_TGT", fNuclearTgtDestructP1 );
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P1_ADEP_TGT", fNuclearTgtDestructP1_ADEP );
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P2_TGT", fNuclearTgtDestructP2 );
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_P3_TGT", fNuclearTgtDestructP3 );
  //
  HDP.DeveloperGet( "FTF_BARYON_PT2_NUCDESTR_P1", fPt2NuclearDestructP1 ); 
  HDP.DeveloperGet( "FTF_BARYON_PT2_NUCDESTR_P2", fPt2NuclearDestructP2 ); 
  HDP.DeveloperGet( "FTF_BARYON_PT2_NUCDESTR_P3", fPt2NuclearDestructP3 ); 
  HDP.DeveloperGet( "FTF_BARYON_PT2_NUCDESTR_P4", fPt2NuclearDestructP4 ); 
  //
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_R2", fR2ofNuclearDestruct );
  HDP.DeveloperGet( "FTF_BARYON_EXCI_E_PER_WNDNUCLN", fExciEnergyPerWoundedNucleon );
  //
  HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_DISP", fDofNuclearDestruct ); // NOTE: "Dof" means "Dispersion of..."
  // 
  // NOTE-1: this parameter has changed from 1. to 9. between 10.2 and 10.3.ref07 !!!
  //         ... then it went back to 1. for the 10.4-candidate... 
  // NOTE-2: this is a "technical" parameter, it should not be changed; this is why
  //         it is defined explicitly rather than via HDP
  // --> HDP.DeveloperGet( "FTF_BARYON_NUCDESTR_MAXPT2", fMaxPt2ofNuclearDestruct );
  fMaxPt2ofNuclearDestruct = 9.0 * CLHEP::GeV*CLHEP::GeV; 
}


G4FTFParamCollMesonProj::G4FTFParamCollMesonProj() : G4FTFParamCollection() {
  // nuclear destruction
  // These parameters are common for all mesons
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_P1_TGT", fNuclearTgtDestructP1 );
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_P1_ADEP_TGT", fNuclearTgtDestructP1_ADEP );
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_P2_TGT", fNuclearTgtDestructP2 );
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_P3_TGT", fNuclearTgtDestructP3 );
  //
  HDP.DeveloperGet( "FTF_MESON_PT2_NUCDESTR_P1", fPt2NuclearDestructP1 ); 
  HDP.DeveloperGet( "FTF_MESON_PT2_NUCDESTR_P2", fPt2NuclearDestructP2 ); 
  HDP.DeveloperGet( "FTF_MESON_PT2_NUCDESTR_P3", fPt2NuclearDestructP3 ); 
  HDP.DeveloperGet( "FTF_MESON_PT2_NUCDESTR_P4", fPt2NuclearDestructP4 ); 
  //
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_R2", fR2ofNuclearDestruct );
  HDP.DeveloperGet( "FTF_MESON_EXCI_E_PER_WNDNUCLN", fExciEnergyPerWoundedNucleon );
  HDP.DeveloperGet( "FTF_MESON_NUCDESTR_DISP", fDofNuclearDestruct ); // NOTE: "Dof" means "Dispersion of..." 
  // NOTE: it is a "technical" parameter, it should not be changed; 
  //       this is why it is defined explicitly rather than via HDP
  fMaxPt2ofNuclearDestruct = 1.0 * CLHEP::GeV*CLHEP::GeV; 
}


G4FTFParamCollPionProj::G4FTFParamCollPionProj() : G4FTFParamCollMesonProj() {
  // parameters of participating pion excitation (pi+/- or pi0)
  // Proc=0 --> Qexchg w/o excitation
  HDP.DeveloperGet( "FTF_PION_PROC0_A1",   fProc0A1 );
  HDP.DeveloperGet( "FTF_PION_PROC0_B1",   fProc0B1 );
  HDP.DeveloperGet( "FTF_PION_PROC0_A2",   fProc0A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC0_B2",   fProc0B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC0_A3",   fProc0A3 );
  HDP.DeveloperGet( "FTF_PION_PROC0_ATOP", fProc0Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC0_YMIN", fProc0Ymin );
  // Proc=1 --> Qexchg w/excitation
  HDP.DeveloperGet( "FTF_PION_PROC1_A1",   fProc1A1 );
  HDP.DeveloperGet( "FTF_PION_PROC1_B1",   fProc1B1 );
  HDP.DeveloperGet( "FTF_PION_PROC1_A2",   fProc1A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC1_B2",   fProc1B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC1_A3",   fProc1A3 );
  HDP.DeveloperGet( "FTF_PION_PROC1_ATOP", fProc1Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC1_YMIN", fProc1Ymin ); 
  // Proc=2 --> Projectile diffraction
  // Q: Would it even make sense to make these configurable ?
  //    The following is hadrcoded:
  //    Projectile Baryon Number > 10 (AbsProjectileBaryonNumber > 10)
  //    ... which is "strange" because projectile is a pion !!!... so it's always OFF    
  //    (see also lines 1007-1016)
  /* As of Oct. 31, 2017 keep these fixed 
  HDP.DeveloperGet( "FTF_PION_PROC2_A1",   fProc2A1 );
  HDP.DeveloperGet( "FTF_PION_PROC2_B1",   fProc2B1 );
  HDP.DeveloperGet( "FTF_PION_PROC2_A2",   fProc2A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC2_B2",   fProc2B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC2_A3",   fProc2A3 );
  HDP.DeveloperGet( "FTF_PION_PROC2_ATOP", fProc2Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC2_YMIN", fProc2Ymin ); 
  */
  // keep fixed so far; see note above
  fProc2A1 =      2.27;
  fProc2B1 =      0.5;
  fProc2A2 = -98052.0;
  fProc2B2 =      4.0;
  fProc2A3 =      0.0;
  fProc2Atop =    0.0;
  fProc2Ymin =    3.0;
  //
  // Proc=3 --> Target diffraction
  HDP.DeveloperGet( "FTF_PION_PROC3_A1",   fProc3A1 );
  HDP.DeveloperGet( "FTF_PION_PROC3_B1",   fProc3B1 );
  HDP.DeveloperGet( "FTF_PION_PROC3_A2",   fProc3A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC3_B2",   fProc3B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC3_A3",   fProc3A3 );
  HDP.DeveloperGet( "FTF_PION_PROC3_ATOP", fProc3Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC3_YMIN", fProc3Ymin ); 
  // for Proc2 & Proc3, pprojectile or target diffraction can be turned ON/OFF
  // if num.baryons >10 (which is strange for projectile which is pion !!!)
  HDP.DeveloperGet( "FTF_PION_DIFF_DISSO_PROJ", fProjDiffDissociation );
  HDP.DeveloperGet( "FTF_PION_DIFF_DISSO_TGT",  fTgtDiffDissociation );
  // Proc=4 --> Qexchg "w/additional multiplier" in excitation 
  HDP.DeveloperGet( "FTF_PION_PROC4_A1",   fProc4A1 );
  HDP.DeveloperGet( "FTF_PION_PROC4_B1",   fProc4B1 );
  HDP.DeveloperGet( "FTF_PION_PROC4_A2",   fProc4A2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC4_B2",   fProc4B2 ); 
  HDP.DeveloperGet( "FTF_PION_PROC4_A3",   fProc4A3 );
  HDP.DeveloperGet( "FTF_PION_PROC4_ATOP", fProc4Atop ); 
  HDP.DeveloperGet( "FTF_PION_PROC4_YMIN", fProc4Ymin );
  //
  HDP.DeveloperGet( "FTF_PION_DELTA_PROB_QEXCHG", fDeltaProbAtQuarkExchange );
  HDP.DeveloperGet( "FTF_PION_DIFF_M_PROJ", fProjMinDiffMass );
  HDP.DeveloperGet( "FTF_PION_NONDIFF_M_PROJ", fProjMinNonDiffMass );
  HDP.DeveloperGet( "FTF_PION_DIFF_M_TGT", fTgtMinDiffMass );
  HDP.DeveloperGet( "FTF_PION_NONDIFF_M_TGT", fTgtMinNonDiffMass );
  HDP.DeveloperGet( "FTF_PION_AVRG_PT2", fAveragePt2 );
  //
  fProbOfSameQuarkExchange = 0.0; // This does NOT seem to apply to the pion case 
  // currently keep these two parameters fixed
  // thus they're defined here explicitly, rather than via HDP
  fProbLogDistrPrD = 0.55; 
  fProbLogDistr    = 0.55; 
}


void G4FTFParamCollection::SetTune( const G4int tuneIndex ) {
  if ( tuneIndex <= 0 || tuneIndex >= G4FTFTunings::sNumberOfTunes ) return;
  switch ( tuneIndex ) {
    case 1 :
      if ( G4FTFTunings::Instance()->GetTuneApplicabilityState(1) != 0 ) SetTune1();
      break;
    case 2 :
      if ( G4FTFTunings::Instance()->GetTuneApplicabilityState(2) != 0 ) SetTune2();
      break;
    case 3 :
      if ( G4FTFTunings::Instance()->GetTuneApplicabilityState(3) != 0 ) SetTune3();
      break;
    case 4 :
      if ( G4FTFTunings::Instance()->GetTuneApplicabilityState(4) != 0 ) SetTune4();
      break;
    case 5 :
      if ( G4FTFTunings::Instance()->GetTuneApplicabilityState(5) != 0 ) SetTune5();
      break;
    case 6 :
      if ( G4FTFTunings::Instance()->GetTuneApplicabilityState(6) != 0 ) SetTune6();
      break;
    case 7 :
      if ( G4FTFTunings::Instance()->GetTuneApplicabilityState(7) != 0 ) SetTune7();
      break;
    case 8 :
      if ( G4FTFTunings::Instance()->GetTuneApplicabilityState(8) != 0 ) SetTune8();
      break;
    case 9 :
      if ( G4FTFTunings::Instance()->GetTuneApplicabilityState(9) != 0 ) SetTune9();
      break;
    // Add here other cases
    default:
      G4ExceptionDescription ed;
      ed << " tuneIndex= " << tuneIndex << G4endl; 
      G4Exception( "G4FTFParamCollection::SetTune", "FTF_PARAM_COLLECTION_001", FatalException,
                   "Not present corresponding SetTuneN() method !" );
  }
  //G4cout << "Called G4FTFParamCollection::SetTune with tuneIndex=" << tuneIndex
  //	   << "  " <<  G4FTFTunings::Instance()->fNameOfTunes[tuneIndex] << G4endl;
}


//====================================================================
//    1st (alternative) TUNE (i.e. indexTune == 1 )
//
// This realistic, alternative tune has been determined by Julia Yarba
// and presented at the hadronic group meeting on 20-Jul-2022
//
// Note: the 0th tune, i.e. indexTune == 0, corresponds to the default
//       set of FTF parameters.
//
//====================================================================

void G4FTFParamCollection::SetTune1() {
  // An empty method implies the take the default values,
  // i.e. as set in G4FTFParamCollection::G4FTFParamCollection()
}


void G4FTFParamCollBaryonProj::SetTune1() {
  G4FTFParamCollection::SetTune1();
  /*
  // parameters of participating hadron (baryon) excitation
  // baryons projectile
  // Proc=0 --> Qexchg w/o excitation
  fProc0A1 = ???;    // FTF_BARYON_PROC0_A1   
  fProc0B1 = ???;    // FTF_BARYON_PROC0_B1
  fProc0A2 = ???;    // FTF_BARYON_PROC0_A2
  fProc0B2 = ???;    // FTF_BARYON_PROC0_B2
  fProc0A3 = ???;    // FTF_BARYON_PROC0_A3
  fProc0Atop = ???;  // FTF_BARYON_PROC0_ATOP
  fProc0Ymin = ???;  // FTF_BARYON_PROC0_YMIN
  // Proc=1 --> Qexchg w/excitation
  fProc1A1 = ???;    // FTF_BARYON_PROC1_A1
  fProc1B1 = ???;    // FTF_BARYON_PROC1_B1
  fProc1A2 = ???;    // FTF_BARYON_PROC1_A2
  fProc1B2 = ???;    // FTF_BARYON_PROC1_B2
  fProc1A3 = ???;    // FTF_BARYON_PROC1_A3
  fProc1Atop = ???;  // FTF_BARYON_PROC1_ATOP
  fProc1Ymin = ???;  // FTF_BARYON_PROC1_YMIN
  // Proc=2 & Proc=3 for the case (diffraction dissociation)
  fProjDiffDissociation = ???;  // FTF_BARYON_DIFF_DISSO_PROJ
  fTgtDiffDissociation = ???;   // FTF_BARYON_DIFF_DISSO_TGT
  // Proc=4 --> Qexchg "w/additional multiplier" in excitation 
  fProc4A1 = ???;    // FTF_BARYON_PROC4_A1
  fProc4B1 = ???;    // FTF_BARYON_PROC4_B1
  fProc4A2 = ???;    // FTF_BARYON_PROC4_A2
  fProc4B2 = ???;    // FTF_BARYON_PROC4_B2
  fProc4A3 = ???;    // FTF_BARYON_PROC4_A3
  fProc4Atop = ???;  // FTF_BARYON_PROC4_ATOP
  fProc4Ymin = ???;  // FTF_BARYON_PROC4_YMIN
  //
  fDeltaProbAtQuarkExchange = ???;  // FTF_BARYON_DELTA_PROB_QEXCHG
  fProbOfSameQuarkExchange = ???;   // FTF_BARYON_PROB_SAME_QEXCHG
  fProjMinDiffMass = ???;           // FTF_BARYON_DIFF_M_PROJ
  fProjMinNonDiffMass = ???;        // FTF_BARYON_NONDIFF_M_PROJ
  fTgtMinDiffMass = ???;            // FTF_BARYON_DIFF_M_TGT
  fTgtMinNonDiffMass = ???;         // FTF_BARYON_NONDIFF_M_TGT
  fAveragePt2 = ???;                // FTF_BARYON_AVRG_PT2
  // nuclear destruction
  fNuclearProjDestructP1 = ???;          // FTF_BARYON_NUCDESTR_P1_PROJ
  fNuclearProjDestructP1_NBRNDEP = ???;  // FTF_BARYON_NUCDESTR_P1_NBRN_PROJ
  //
  fNuclearTgtDestructP1 = ???;       // FTF_BARYON_NUCDESTR_P1_TGT
  fNuclearTgtDestructP1_ADEP = ???;  // FTF_BARYON_NUCDESTR_P1_ADEP_TGT
  fNuclearTgtDestructP2 = ???;       // FTF_BARYON_NUCDESTR_P2_TGT
  fNuclearTgtDestructP3 = ???;       // FTF_BARYON_NUCDESTR_P3_TGT
  //
  fPt2NuclearDestructP1 = ???;  // FTF_BARYON_PT2_NUCDESTR_P1
  fPt2NuclearDestructP2 = ???;  // FTF_BARYON_PT2_NUCDESTR_P2
  fPt2NuclearDestructP3 = ???;  // FTF_BARYON_PT2_NUCDESTR_P3
  fPt2NuclearDestructP4 = ???;  // FTF_BARYON_PT2_NUCDESTR_P4
  //
  fR2ofNuclearDestruct = ???;          // FTF_BARYON_NUCDESTR_R2
  fExciEnergyPerWoundedNucleon = ???;  // FTF_BARYON_EXCI_E_PER_WNDNUCLN 
  //
  fDofNuclearDestruct = ???;  // FTF_BARYON_NUCDESTR_DISP
  */
  // Values below from Julia Yarba's slides at the hadronic group meeting on 20-Jul-2022
  fExciEnergyPerWoundedNucleon = 26.1;  // +/- 0.4      // FTF_BARYON_EXCI_E_PER_WNDNUCLN 
  fNuclearTgtDestructP1 = 0.00173;      // +/- 0.00004  // FTF_BARYON_NUCDESTR_P1_TGT
  fNuclearTgtDestructP1_ADEP = true;                    // FTF_BARYON_NUCDESTR_P1_ADEP_TGT
  fProc1A1 = 23.6;                      // +/- 0.8      // FTF_BARYON_PROC1_A1
  fProc1A2 = -99.3;                     // +/- 0.4      // FTF_BARYON_PROC1_A2
  fProc1B1 = 0.815;                     // +/- 0.007    // FTF_BARYON_PROC1_B1
  fProc1B2 = 1.98;                      // +/- 0.03     // FTF_BARYON_PROC1_B2
}


void G4FTFParamCollMesonProj::SetTune1() {
  G4FTFParamCollection::SetTune1();
  /*
  // nuclear destruction
  fNuclearTgtDestructP1 = ???;         // FTF_MESON_NUCDESTR_P1_TGT 
  fNuclearTgtDestructP1_ADEP = ???;    // FTF_MESON_NUCDESTR_P1_ADEP_TGT
  fNuclearTgtDestructP2 = ???;         // FTF_MESON_NUCDESTR_P2_TGT
  fNuclearTgtDestructP3 = ???;         // FTF_MESON_NUCDESTR_P3_TGT 
  //
  fPt2NuclearDestructP1 = ???;         // FTF_MESON_PT2_NUCDESTR_P1
  fPt2NuclearDestructP2 = ???;         // FTF_MESON_PT2_NUCDESTR_P2
  fPt2NuclearDestructP3 = ???;         // FTF_MESON_PT2_NUCDESTR_P3
  fPt2NuclearDestructP4 = ???;         // FTF_MESON_PT2_NUCDESTR_P4
  //
  fR2ofNuclearDestruct = ???;          // FTF_MESON_NUCDESTR_R2
  fExciEnergyPerWoundedNucleon = ???;  // FTF_MESON_EXCI_E_PER_WNDNUCLN
  fDofNuclearDestruct = ???;           // FTF_MESON_NUCDESTR_DISP
  */
}


void G4FTFParamCollPionProj::SetTune1( ) {
  G4FTFParamCollMesonProj::SetTune1();
  /*
  // parameters of participating pion excitation (pi+/- or pi0)
  // Proc=0 --> Qexchg w/o excitation
  fProc0A1 = ???;    // FTF_PION_PROC0_A1
  fProc0B1 = ???;    // FTF_PION_PROC0_B1
  fProc0A2 = ???;    // FTF_PION_PROC0_A2
  fProc0B2 = ???;    // FTF_PION_PROC0_B2   
  fProc0A3 = ???;    // FTF_PION_PROC0_A3
  fProc0Atop = ???;  // FTF_PION_PROC0_ATOP 
  fProc0Ymin = ???;  // FTF_PION_PROC0_YMIN
  // Proc=1 --> Qexchg w/excitation
  fProc1A1 = ???;    // FTF_PION_PROC1_A1
  fProc1B1 = ???;    // FTF_PION_PROC1_B1
  fProc1A2 = ???;    // FTF_PION_PROC1_A2
  fProc1B2 = ???;    // FTF_PION_PROC1_B2
  fProc1A3 = ???;    // FTF_PION_PROC1_A3
  fProc1Atop = ???;  // FTF_PION_PROC1_ATOP
  fProc1Ymin = ???;  // FTF_PION_PROC1_YMIN
  // Proc=2 --> Projectile diffraction : keep these fixed 
  //Fixed  fProc2A1 = ???;    // FTF_PION_PROC2_A1   
  //Fixed  fProc2B1 = ???;    // FTF_PION_PROC2_B1
  //Fixed  fProc2A2 = ???;    // FTF_PION_PROC2_A2   
  //Fixed  fProc2B2 = ???;    // FTF_PION_PROC2_B2
  //Fixed  fProc2A3 = ???;    // FTF_PION_PROC2_A3
  //Fixed  fProc2Atop = ???;  // FTF_PION_PROC2_ATOP 
  //Fixed  fProc2Ymin = ???;  // FTF_PION_PROC2_YMIN
  // Proc=3 --> Target diffraction
  fProc3A1 = ???;    // FTF_PION_PROC3_A1
  fProc3B1 = ???;    // FTF_PION_PROC3_B1   
  fProc3A2 = ???;    // FTF_PION_PROC3_A2
  fProc3B2 = ???;    // FTF_PION_PROC3_B2
  fProc3A3 = ???;    // FTF_PION_PROC3_A3   
  fProc3Atop = ???;  // FTF_PION_PROC3_ATOP 
  fProc3Ymin = ???;  // FTF_PION_PROC3_YMIN 
  //
  fProjDiffDissociation = ???;  // FTF_PION_DIFF_DISSO_PROJ 
  fTgtDiffDissociation = ???;   // FTF_PION_DIFF_DISSO_TGT
  //
  // Proc=4 --> Qexchg "w/additional multiplier" in excitation 
  fProc4A1 = ???;    // FTF_PION_PROC4_A1   
  fProc4B1 = ???;    // FTF_PION_PROC4_B1   
  fProc4A2 = ???;    // FTF_PION_PROC4_A2   
  fProc4B2 = ???;    // FTF_PION_PROC4_B2   
  fProc4A3 = ???;    // FTF_PION_PROC4_A3   
  fProc4Atop = ???;  // FTF_PION_PROC4_ATOP 
  fProc4Ymin = ???;  // FTF_PION_PROC4_YMIN
  //
  fDeltaProbAtQuarkExchange = ???;  // FTF_PION_DELTA_PROB_QEXCHG
  fProjMinDiffMass = ???;           // FTF_PION_DIFF_M_PROJ 
  fProjMinNonDiffMass = ???;        // FTF_PION_NONDIFF_M_PROJ 
  fTgtMinDiffMass = ???;            // FTF_PION_DIFF_M_TGT
  fTgtMinNonDiffMass = ???;         // FTF_PION_NONDIFF_M_TGT
  fAveragePt2 = ???;                // FTF_PION_AVRG_PT2
  */
}

//====================================================================
//    2nd (alternative) TUNE (i.e. indexTune == 2 )
//
// This is work-in-progress, very preliminary alternative tune that 
// has been outlined by Julia Yarba presented at the 27th CM on 26-JSept-2022
// The nuclear destruction parameters that are common for *all mesons*
// (see G4FTFParamCollMesonProj::SetTune2)
// while the Qexchg with excitationn of participants parameters are
// for *pion" projectile only (see G4FTFParamCollPionProj::SetTune2)
//
//====================================================================

void G4FTFParamCollection::SetTune2() {
  // An empty method implies the take the default values,
  // i.e. as set in G4FTFParamCollection::G4FTFParamCollection()
}


void G4FTFParamCollBaryonProj::SetTune2() {
  G4FTFParamCollection::SetTune2();
  //...
}


void G4FTFParamCollMesonProj::SetTune2() {

  G4FTFParamCollection::SetTune2();

  // nuclear detsruction
  //
  // NOTE: These values are the same for all mesons 
  //       (although bear in mind that they've been obtained for the pion projectile
  //        via fits against experimaental data for the pion beam)
  //
  fExciEnergyPerWoundedNucleon = 58.1;  // +/- 0.7      // FTF_MESON_EXCI_E_PER_WNDNUCLN 
  fNuclearTgtDestructP1 = 0.001026;     // +/- 0.00003  // FTF_MESON_NUCDESTR_P1_TGT 
  fNuclearTgtDestructP1_ADEP = true;                    // FTF_MESON_NUCDESTR_P1_ADEP_TGT
  
  return;
  
}


void G4FTFParamCollPionProj::SetTune2( ) {

  G4FTFParamCollMesonProj::SetTune2();

  // Proc=1 --> Qexchg w/excitation
  fProc1A1 = 5.84;       // +/- 0.12     // FTF_PION_PROC1_A1
  fProc1B1 = 0.337;      // +/- 0.006    // FTF_PION_PROC1_B1
  fProc1A2 = -7.57;      // +/- 0.08     // FTF_PION_PROC1_A2
  fProc1B2 = 0.44;       // +/- 0.008    // FTF_PION_PROC1_B2

  return;
}

//====================================================================
//    3rd (alternative) TUNE (i.e. indexTune == 3 )
//
//    Combination of the 1st (baryon) and 2nd (pion) tunes
//
//====================================================================

void G4FTFParamCollection::SetTune3() {
  // An empty method implies the take the default values,
  // i.e. as set in G4FTFParamCollection::G4FTFParamCollection()
}


void G4FTFParamCollBaryonProj::SetTune3() {

  G4FTFParamCollection::SetTune3();

  // Values below from Julia Yarba's slides at the hadronic group meeting on 20-Jul-2022
  //
  fExciEnergyPerWoundedNucleon = 26.1;  // +/- 0.4      // FTF_BARYON_EXCI_E_PER_WNDNUCLN 
  fNuclearTgtDestructP1 = 0.00173;      // +/- 0.00004  // FTF_BARYON_NUCDESTR_P1_TGT
  fNuclearTgtDestructP1_ADEP = true;                    // FTF_BARYON_NUCDESTR_P1_ADEP_TGT
  fProc1A1 = 23.6;                      // +/- 0.8      // FTF_BARYON_PROC1_A1
  fProc1A2 = -99.3;                     // +/- 0.4      // FTF_BARYON_PROC1_A2
  fProc1B1 = 0.815;                     // +/- 0.007    // FTF_BARYON_PROC1_B1
  fProc1B2 = 1.98;                      // +/- 0.03     // FTF_BARYON_PROC1_B2
  
  return;
  
}


void G4FTFParamCollMesonProj::SetTune3() {

  G4FTFParamCollection::SetTune3();

  // nuclear detsruction
  //
  // NOTE: These values are the same for all mesons 
  //       (although bear in mind that they've been obtained for the pion projectile
  //        via fits against experimaental data for the pion beam)
  //
  fExciEnergyPerWoundedNucleon = 58.1;  // +/- 0.7      // FTF_MESON_EXCI_E_PER_WNDNUCLN 
  fNuclearTgtDestructP1 = 0.001026;     // +/- 0.00003  // FTF_MESON_NUCDESTR_P1_TGT 
  fNuclearTgtDestructP1_ADEP = true;                    // FTF_MESON_NUCDESTR_P1_ADEP_TGT
  
  return;
  
}


void G4FTFParamCollPionProj::SetTune3( ) {

  G4FTFParamCollMesonProj::SetTune3();

  // Proc=1 --> Qexchg w/excitation
  fProc1A1 = 5.84;       // +/- 0.12     // FTF_PION_PROC1_A1
  fProc1B1 = 0.337;      // +/- 0.006    // FTF_PION_PROC1_B1
  fProc1A2 = -7.57;      // +/- 0.08     // FTF_PION_PROC1_A2
  fProc1B2 = 0.44;       // +/- 0.008    // FTF_PION_PROC1_B2
  
  return;
  
}

//====================================================================
//    4th (alternative) TUNE (i.e. indexTune == 4 )
//
//    DUMMY tune: identical to the default set of parameters.
//                You can replace it with a "real" tune by specifying
//                only the non-default parameters.
//====================================================================

void G4FTFParamCollection::SetTune4() {
  // An empty method implies the take the default values,
  // i.e. as set in G4FTFParamCollection::G4FTFParamCollection()
}


void G4FTFParamCollBaryonProj::SetTune4() {
  G4FTFParamCollection::SetTune4();
  //...
}


void G4FTFParamCollMesonProj::SetTune4() {
  G4FTFParamCollection::SetTune4();
  //...
}


void G4FTFParamCollPionProj::SetTune4( ) {
  G4FTFParamCollMesonProj::SetTune4();
  //...
}

//====================================================================
//    5th (alternative) TUNE (i.e. indexTune == 5 )
//
//    DUMMY tune: identical to the default set of parameters.
//                You can replace it with a "real" tune by specifying
//                only the non-default parameters.
//====================================================================

void G4FTFParamCollection::SetTune5() {
  // An empty method implies the take the default values,
  // i.e. as set in G4FTFParamCollection::G4FTFParamCollection()
}


void G4FTFParamCollBaryonProj::SetTune5() {
  G4FTFParamCollection::SetTune5();
  //...
}


void G4FTFParamCollMesonProj::SetTune5() {
  G4FTFParamCollection::SetTune5();
  //...
}


void G4FTFParamCollPionProj::SetTune5( ) {
  G4FTFParamCollMesonProj::SetTune5();
  //...
}

//====================================================================
//    6th (alternative) TUNE (i.e. indexTune == 6 )
//
//    DUMMY tune: identical to the default set of parameters.
//                You can replace it with a "real" tune by specifying
//                only the non-default parameters.
//====================================================================

void G4FTFParamCollection::SetTune6() {
  // An empty method implies the take the default values,
  // i.e. as set in G4FTFParamCollection::G4FTFParamCollection()
}


void G4FTFParamCollBaryonProj::SetTune6() {
  G4FTFParamCollection::SetTune6();
  //...
}


void G4FTFParamCollMesonProj::SetTune6() {
  G4FTFParamCollection::SetTune6();
  //...
}


void G4FTFParamCollPionProj::SetTune6( ) {
  G4FTFParamCollMesonProj::SetTune6();
  //...
}

//====================================================================
//    7th (alternative) TUNE (i.e. indexTune == 7 )
//
//    DUMMY tune: identical to the default set of parameters.
//                You can replace it with a "real" tune by specifying
//                only the non-default parameters.
//====================================================================

void G4FTFParamCollection::SetTune7() {
  // An empty method implies the take the default values,
  // i.e. as set in G4FTFParamCollection::G4FTFParamCollection()
}


void G4FTFParamCollBaryonProj::SetTune7() {
  G4FTFParamCollection::SetTune7();
  //...
}


void G4FTFParamCollMesonProj::SetTune7() {
  G4FTFParamCollection::SetTune7();
  //...
}


void G4FTFParamCollPionProj::SetTune7( ) {
  G4FTFParamCollMesonProj::SetTune7();
  //...
}

//====================================================================
//    8th (alternative) TUNE (i.e. indexTune == 8 )
//
//    DUMMY tune: identical to the default set of parameters.
//                You can replace it with a "real" tune by specifying
//                only the non-default parameters.
//====================================================================

void G4FTFParamCollection::SetTune8() {
  // An empty method implies the take the default values,
  // i.e. as set in G4FTFParamCollection::G4FTFParamCollection()
}


void G4FTFParamCollBaryonProj::SetTune8() {
  G4FTFParamCollection::SetTune8();
  //...
}


void G4FTFParamCollMesonProj::SetTune8() {
  G4FTFParamCollection::SetTune8();
  //...
}


void G4FTFParamCollPionProj::SetTune8( ) {
  G4FTFParamCollMesonProj::SetTune8();
  //...
}

//====================================================================
//    9th (alternative) TUNE (i.e. indexTune == 9 )
//
//    DUMMY tune: identical to the default set of parameters.
//                You can replace it with a "real" tune by specifying
//                only the non-default parameters.
//====================================================================

void G4FTFParamCollection::SetTune9() {
  // An empty method implies the take the default values,
  // i.e. as set in G4FTFParamCollection::G4FTFParamCollection()
}


void G4FTFParamCollBaryonProj::SetTune9() {
  G4FTFParamCollection::SetTune9();
  //...
}


void G4FTFParamCollMesonProj::SetTune9() {
  G4FTFParamCollection::SetTune9();
  //...
}


void G4FTFParamCollPionProj::SetTune9( ) {
  G4FTFParamCollMesonProj::SetTune9();
  //...
}

//...
