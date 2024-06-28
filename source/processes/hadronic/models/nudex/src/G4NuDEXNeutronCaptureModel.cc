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
// -------------------------------------------------------------------
//
//      GEANT4 source file 
//
//      File name:     G4NuDEXNeutronCaptureModel
//
//      Author:        E.Mendoza & A.Ribon
// 
//      Creation date: 29 May 2024
//
//      Description:   This class (a proxy of the class G4NuDEX) uses
//                     the NuDEX model to produce gammas and internal
//                     conversion electrons from neutron capture.
//                     Whenever NuDEX is not applicable, G4PhotonEvaporation
//                     is used.
//                     The implementation of this class follows the code
//                     of the class G4NeutronRadCapture.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//

#include "G4NuDEXNeutronCaptureModel.hh"
#include "G4NuDEXStatisticalNucleus.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4NucleiProperties.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4HadronicParameters.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"
#include "G4PhotonEvaporation.hh"
#include "G4PhysicsModelCatalog.hh"


G4NuDEXNeutronCaptureModel::G4NuDEXNeutronCaptureModel() : G4HadronicInteraction( "nuDEX_neutronCapture" ) {
  for ( G4int i = 0; i < G4NUDEX_MAXZA; i++ ) {
    theStatisticalNucleus[i] = nullptr;
    HasData[i] = 0;
  }
  BrOption = -1;
  BandWidth = 0;
  NuDEXLibDirectory = "";
  photonEvaporation = nullptr;
  auto ch = G4FindDataDir( "G4NUDEXLIBDATA" );
  if ( ch == nullptr ) {
    G4Exception( "G4NuDEXNeutronCaptureModel()", "had0707", FatalException, "Environment variable G4NUDEXLIBDATA is not defined" );
  } else {
    NuDEXLibDirectory = G4String(ch);
  }
}


void G4NuDEXNeutronCaptureModel::InitialiseModel() {
  if ( photonEvaporation != nullptr ) return;
  G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
  minExcitation = param->GetMinExcitation();
  photonEvaporation = new G4PhotonEvaporation;
  photonEvaporation->Initialise();
  photonEvaporation->SetICM( true );
  secID = G4PhysicsModelCatalog::GetModelID( "model_" + GetModelName() );
  lowestEnergyLimit = 10.0*CLHEP::eV;
  minExcitation = 0.1*CLHEP::keV;
}


G4NuDEXNeutronCaptureModel::~G4NuDEXNeutronCaptureModel(){
  for ( G4int i = 0; i < G4NUDEX_MAXZA; i++ ) {
    if ( theStatisticalNucleus[i] ) delete theStatisticalNucleus[i];
  }
}


G4HadFinalState* G4NuDEXNeutronCaptureModel::ApplyYourself( const G4HadProjectile &aTrack, G4Nucleus &theNucleus ) {
  theParticleChange.Clear();
  theParticleChange.SetStatusChange( stopAndKill );
  G4int A = theNucleus.GetA_asInt();
  G4int Z = theNucleus.GetZ_asInt();
  G4double time = aTrack.GetGlobalTime();  // Time in the lab frame
  // Create initial state
  G4LorentzVector lab4mom( 0.0, 0.0, 0.0, G4NucleiProperties::GetNuclearMass(A, Z) );
  lab4mom += aTrack.Get4Momentum();
  G4double systemMass = lab4mom.mag();
  ++A;  // Compound nucleus: target nucleus + neutron
  G4double compoundMass = G4NucleiProperties::GetNuclearMass(A, Z);
  // If the energy available is to small to do anything interesting, gives up
  if ( systemMass - compoundMass  <= lowestEnergyLimit ) return &theParticleChange;
  G4ThreeVector boostFromCMtoLAB = lab4mom.boostVector();
  G4double neutronEnergy = aTrack.GetKineticEnergy();

  // Try to apply NuDEX
  //G4int lspin = 0;     // l-spin = 0, 1, 2 --> s-wave, p-wave, d-wave ... 
  //G4int jspinx2 = -1;  // A negative value of jspinx2 means that is sampled according to the 2J+1 rule.
  //G4int initialLevel = SelectInitialLevel( Z, A, neutronEnergy, lspin, jspinx2 );
  G4int initialLevel = -1;  // thermal neutron capture
  std::vector< char > pType;
  std::vector< G4double > pEnergy, pTime;
  G4int npar = GenerateNeutronCaptureCascade( Z, A, neutronEnergy, initialLevel, pType, pEnergy, pTime );
  if ( npar > 0 ) {  // NuDEX can be applied
    
    G4LorentzVector remainingLab4mom = lab4mom;
    G4double latestEmission = time;
    // Loop over the EM particles produced by 'GenerateNeutronCaptureCascade' and add them to the
    // theParticleChange as secondaries. These particles are produced by NuDEX in the nucleus' rest-frame.
    for ( G4int i = 0; i < npar; i++ ) {
      G4ParticleDefinition* particleDef = nullptr;
      if ( pType.at(i) == 'g' ) {
        particleDef = G4Gamma::Definition();
      } else if  (pType.at(i) == 'e' ) {
        particleDef = G4Electron::Definition();
      } else if ( pType.at(i) == 'p' ) {
        particleDef = G4Positron::Definition();
      } else {
        G4Exception( "G4NUDEXNeutronCaptureModel::ApplyYourself()", "had0707", FatalException, "Unknown particle type" );
      }
      G4double phi = G4UniformRand()*twopi;
      G4double costheta = 2.0*G4UniformRand() - 1.0;
      G4double sintheta = std::sqrt( 1.0 - costheta*costheta );
      G4ThreeVector direction( sintheta*std::cos(phi), sintheta*std::sin(phi), costheta );
      G4double mass = particleDef->GetPDGMass();
      G4double particle3momMod = std::sqrt( pEnergy.at(i) * ( pEnergy.at(i) + 2.0*mass ) );
      G4LorentzVector particle4mom( particle3momMod*direction, mass + pEnergy.at(i) );  // In the center-of-mass frame
      particle4mom.boost( boostFromCMtoLAB );  // Now in the Lab frame
      G4HadSecondary* secondary = new G4HadSecondary( new G4DynamicParticle( particleDef, particle4mom ) );
      remainingLab4mom -= particle4mom;
      // For simplicity, we neglect below the different frames of time (Lab) and pTime (center-of-mass) 
      secondary->SetTime( time + pTime.at(i) );
      if ( latestEmission < time + pTime.at(i) ) latestEmission = time + pTime.at(i);
      secondary->SetCreatorModelID( secID );
      theParticleChange.AddSecondary( *secondary );
      delete secondary;
    }
    // Treat now the residual nucleus (which is neglected by NuDEX)
    const G4ParticleDefinition* resNuclDef = nullptr;
    if      ( Z == 1  &&  A == 2 ) resNuclDef = G4Deuteron::Definition();
    else if ( Z == 1  &&  A == 3 ) resNuclDef = G4Triton::Definition();
    else if ( Z == 2  &&  A == 3 ) resNuclDef = G4He3::Definition();
    else if ( Z == 2  &&  A == 4 ) resNuclDef = G4Alpha::Alpha();
    else resNuclDef = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon( Z, A, 0.0, noFloat, 0 );
    if ( resNuclDef ) {
      // To conserve energy-momentum, remainingLab4mom should be the Lorentz 4-momentum of the residual nucleus.
      // Imposing the mass 'compoundMass' to the residual nucleus, and trying to conserve the total energy
      // while keeping as low as possible the violation of the 3-momentum; in the case that the total energy
      // cannot be conserved, the residual nucleus is produced at rest (in the Lab frame).
      G4double resNuclLabEkin = std::max( remainingLab4mom.e() - compoundMass, 0.0 );
      G4double resNuclLab3momMod = 0.0;
      G4ThreeVector resNuclLabDir( 0.0, 0.0, 0.0 );
      if ( resNuclLabEkin > 0.0 ) {
        resNuclLab3momMod = std::sqrt( resNuclLabEkin * ( resNuclLabEkin + 2.0*compoundMass ) );
        resNuclLabDir = remainingLab4mom.vect().unit();
      }
      G4LorentzVector resNuclLab4mom( resNuclLab3momMod*resNuclLabDir, resNuclLabEkin + compoundMass );
      G4HadSecondary* secondary = new G4HadSecondary( new G4DynamicParticle( resNuclDef, resNuclLab4mom ) );
      secondary->SetTime( latestEmission );
      secondary->SetCreatorModelID( secID );
      theParticleChange.AddSecondary( *secondary );
      delete secondary;
    }

  } else {  // NuDEX cannot be applied: use G4PhotonEvaporation

    // Code taken from G4NeutronRadCapture

    G4Fragment* aFragment = new G4Fragment( A, Z, lab4mom );
    G4FragmentVector* fv = photonEvaporation->BreakUpFragment( aFragment );
    if ( fv == nullptr ) fv = new G4FragmentVector;
    fv->push_back( aFragment );
    size_t n = fv->size();
    for ( size_t i = 0; i < n; ++i ) {
      G4Fragment* f = (*fv)[i];    
      G4double etot = f->GetMomentum().e();
      Z = f->GetZ_asInt();
      A = f->GetA_asInt();
      const G4ParticleDefinition* theDef = nullptr;
      if      ( Z == 0  && A == 0 ) { theDef = f->GetParticleDefinition(); }
      else if ( Z == 1  && A == 2 ) { theDef = G4Deuteron::Definition(); }
      else if ( Z == 1  && A == 3 ) { theDef = G4Triton::Definition(); }
      else if ( Z == 2  && A == 3 ) { theDef = G4He3::Definition(); }
      else if ( Z == 2  && A == 4 ) { theDef = G4Alpha::Definition(); }
      else {
        G4double eexc = f->GetExcitationEnergy();
        if ( eexc <= minExcitation ) eexc = 0.0;
        theDef = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon( Z, A, eexc, noFloat, 0 );
      }
      G4double ekin = std::max( 0.0, etot - theDef->GetPDGMass() );
      G4HadSecondary* news = new G4HadSecondary( new G4DynamicParticle( theDef, f->GetMomentum().vect().unit(), ekin ) );
      G4double timeF = f->GetCreationTime();
      if ( timeF < 0.0 ) timeF = 0.0;
      news->SetTime( time + timeF );
      news->SetCreatorModelID( secID );
      theParticleChange.AddSecondary( *news );
      delete news;
      delete f;
    }
    delete fv;
  }
  
  return &theParticleChange;
}


G4int G4NuDEXNeutronCaptureModel::GenerateNeutronCaptureCascade( G4int theZ, G4int theA, G4double NeutronEnergy, G4int InitialLevel,
					                         std::vector< char >& pType, std::vector< G4double >& pEnergy,
					                         std::vector< G4double >& pTime ) {
  // Returns the number of emitted particles. Returns -1 if the nucleus is not in the database.
  G4int theZA = 1000*theZ + theA;
  G4int check = Init( theZA );
  if ( check < 0 ) return -1;
  G4double Sn, I0;
  theStatisticalNucleus[theZA]->GetSnAndI0( Sn, I0 ); Sn *= MeV;  // I0 is the spin of the A-1 nucleus in the g.s.
  G4double ExcitationEnergy = Sn + (theA-1.0)/(G4double)theA*NeutronEnergy;
  G4int nPar = theStatisticalNucleus[theZA]->GenerateCascade( InitialLevel, ExcitationEnergy/MeV, pType, pEnergy, pTime );
  for ( G4int i = 0; i < nPar; i++ ) {
    pEnergy.at(i) *= MeV;
    pTime.at(i) *= s;
  }
  return nPar;
}


G4int G4NuDEXNeutronCaptureModel::Init( G4int theZA, unsigned int seed1, unsigned int seed2, unsigned int seed3 ) {
  if ( HasData[theZA] == -1 ) return -1;
  if ( HasData[theZA] ==  1 ) return 0;
  if ( theStatisticalNucleus[theZA] == 0 ) {
    G4int theZ = theZA/1000;
    G4int theA = theZA-1000*theZ;
    theStatisticalNucleus[theZA] = new G4NuDEXStatisticalNucleus( theZ, theA );
    if ( BandWidth != 0 ) theStatisticalNucleus[theZA]->SetBandWidth( BandWidth );
    theStatisticalNucleus[theZA]->SetBrOption( BrOption );
    if ( seed1 > 0 ) theStatisticalNucleus[theZA]->SetRandom1Seed( seed1 );
    if ( seed2 > 0 ) theStatisticalNucleus[theZA]->SetRandom1Seed( seed2 );
    if ( seed3 > 0 ) theStatisticalNucleus[theZA]->SetRandom1Seed( seed3 );
    G4int check = theStatisticalNucleus[theZA]->Init( NuDEXLibDirectory.c_str() );
    if ( check < 0 ) {
      HasData[theZA] = -1;
      return -1;
    } else {
      HasData[theZA] = 1;
    }
  }
  return 0;
}


G4int G4NuDEXNeutronCaptureModel::SelectInitialLevel( G4int theCompoundZ, G4int theCompoundA, G4double NeutronEnergy, G4int lspin, G4int jspinx2 ) {
  // Initial level for neutron capture. If jspinx2 < 0 it is sampled according to the 2J+1 rule.
  // l-spin = 0, 1, 2 --> s-wave, p-wave, d-wave ...
  G4int theZ = theCompoundZ;
  G4int theA = theCompoundA;
  G4int theZA = 1000*theZ + theA;
  G4int check = Init( theZA );
  if ( check < 0 ) return -1;
  G4double Sn, I0;
  theStatisticalNucleus[theZA]->GetSnAndI0( Sn, I0 ); Sn *= MeV;  // I0 is the spin of the A-1 nucleus in the g.s.
  G4double ExcitationEnergy = Sn + (theA-1.0)/(G4double)theA*NeutronEnergy;
  if ( lspin < 0 ) lspin = 0;
  if ( jspinx2 < 0 ) jspinx2 = SampleJ( theZ, theA, lspin );
  G4bool parity = false;
  if ( ( I0 >= 0  &&  (lspin%2) == 0 ) || ( I0 < 0  &&  (lspin%2) == 1 ) ) parity = true;
  G4int InitialLevel = theStatisticalNucleus[theZA]->GetClosestLevel( ExcitationEnergy/MeV, jspinx2, parity );
  return InitialLevel;
}


G4int G4NuDEXNeutronCaptureModel ::SampleJ( G4int theCompoundZ, G4int theCompoundA, G4int lspin ) {
  // Samples J for this l-spin (l-spin = 0, 1, 2 --> s-wave, p-wave, d-wave ...)
  // The probability will be proportional to 2J+1
  // Returns Jx2
  G4int AllowedJx2[100];
  G4int NAllowedJvals = GetAllowedJx2values( theCompoundZ, theCompoundA, lspin, AllowedJx2 );
  G4double AllowedJx2CumulProb[100], TotalCumul = 0.0;
  for ( G4int i = 0; i < NAllowedJvals; i++ ) {
    AllowedJx2CumulProb[i] = AllowedJx2[i] + 1.0;
    TotalCumul += AllowedJx2CumulProb[i];
  }
  for ( G4int i = 0; i < NAllowedJvals; i++ ) {
    AllowedJx2CumulProb[i] /= TotalCumul;
    if ( i > 0 ) AllowedJx2CumulProb[i] += AllowedJx2CumulProb[i-1];
  }
  G4double rand = G4UniformRand();
  G4int i_result = -1;
  for ( G4int i = 0; i < NAllowedJvals; i++ ) {
    if ( rand < AllowedJx2CumulProb[i] ) {
      i_result = i;  break;
    }
  }
  if ( i_result < 0 ) {
    G4cerr << " ############ Error in " << __FILE__ << ", line " << __LINE__ << " ############"<< G4endl;
    exit(1);
  }
  G4int jspinx2 = AllowedJx2[i_result];
  return jspinx2;
}


G4int G4NuDEXNeutronCaptureModel::GetAllowedJx2values( G4int theCompoundZ, G4int theCompoundA, G4int lspin, G4int* jx2vals ) {
  // Provides the allowed jx2 values in neutron capture for a certain l-spin (l-spin = 0, 1, 2 --> s-wave, p-wave, d-wave ...)
  G4int theZA = 1000*theCompoundZ + theCompoundA;
  G4int check = Init( theZA );
  if ( check < 0 ) return -1;
  G4double Sn, I0;
  theStatisticalNucleus[theZA]->GetSnAndI0( Sn, I0 ); Sn *= MeV;  // I0 is the spin of the A-1 nucleus in the g.s.
  G4int Ix2 = (G4int)( ( std::fabs(I0) + 0.1 )*2.0 );
  G4int Jx2min = std::min( std::abs( Ix2-1-2*lspin ), std::abs( Ix2+1-2*lspin ) );
  G4int Jx2max = Ix2 + 1 + 2*lspin;
  G4int NAllowedJvals = 0;
  for ( G4int Jx2 = Jx2min; Jx2 <= Jx2max; Jx2 += 2 ) {
    if ( Jx2 >= 0 ) {
      jx2vals[NAllowedJvals] = Jx2; 
      NAllowedJvals++;
    }
  }
  return NAllowedJvals;
}
