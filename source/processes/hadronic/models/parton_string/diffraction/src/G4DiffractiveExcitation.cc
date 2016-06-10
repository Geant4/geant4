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
// $Id: G4DiffractiveExcitation.cc 87254 2014-11-28 07:49:23Z gcosmo $
//

// ------------------------------------------------------------
//      GEANT 4 class implemetation file
//
//      ---------------- G4DiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//        diffractive Excitation used by strings models
//             Take a projectile and a target
//             excite the projectile and target
//  Essential changed by V. Uzhinsky in November - December 2006
//  in order to put it in a correspondence with original FRITIOF
//  model. Variant of FRITIOF with nucleon de-excitation is implemented.
//  Other changes by V.Uzhinsky in May 2007 were introduced to fit
//  meson-nucleon interactions. Additional changes by V. Uzhinsky
//  were introduced in December 2006. They treat diffraction dissociation
//  processes more exactly.
//  Correct treatment of the diffraction dissociation - 2012, V. Uzhinsky
//  Mass distributions for resonances and uu-diquark suppression in protons,
//  and dd-diquarks suppression in neutrons were introduced by V. Uzhinsky, 2014
// ---------------------------------------------------------------------

#include "globals.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4DiffractiveExcitation.hh"
#include "G4FTFParameters.hh"
#include "G4ElasticHNScattering.hh"

#include "G4LorentzRotation.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh" 
#include "G4ParticleTable.hh"
#include "G4SampleResonance.hh"
#include "G4VSplitableHadron.hh"
#include "G4ExcitedString.hh"
#include "G4Neutron.hh"


//#include "G4ios.hh"
//#include "UZHI_diffraction.hh"


//============================================================================

//#define debugFTFexictation


//============================================================================

G4DiffractiveExcitation::G4DiffractiveExcitation() {}


//============================================================================

G4DiffractiveExcitation::~G4DiffractiveExcitation() {}


//============================================================================

G4bool G4DiffractiveExcitation::ExciteParticipants( G4VSplitableHadron*    projectile, 
                                                    G4VSplitableHadron*    target,
                                                    G4FTFParameters*       theParameters,
                                                    G4ElasticHNScattering* theElastic ) const {

  #ifdef debugFTFexictation
  G4cout << G4endl << "FTF ExciteParticipants --------------" << G4endl;
  #endif

  // Projectile parameters
  G4LorentzVector Pprojectile = projectile->Get4Momentum();
  if ( Pprojectile.z() < 0.0 ) return false;

  G4int ProjectilePDGcode = projectile->GetDefinition()->GetPDGEncoding();
  G4int absProjectilePDGcode = std::abs( ProjectilePDGcode );
  G4double M0projectile = Pprojectile.mag();
  G4double ProjectileRapidity = Pprojectile.rapidity();

  // Target parameters
  G4LorentzVector Ptarget = target->Get4Momentum();
  G4int TargetPDGcode = target->GetDefinition()->GetPDGEncoding();
  G4int absTargetPDGcode = std::abs( TargetPDGcode );
  G4double M0target = Ptarget.mag();
  //G4double TargetRapidity = Ptarget.rapidity();

  // Kinematical properties of the interactions
  G4LorentzVector Psum = Pprojectile + Ptarget;  // 4-momentum in CMS
  G4double S = Psum.mag2();
  G4double SqrtS = std::sqrt( S ); 
  //Uzhi_SqrtS = std::sqrt( S );

  // Check off-shellness of the participants
  G4SampleResonance BrW;                                               // Uzhi Oct. 2014

  G4bool PutOnMassShell( false );

  G4double MminProjectile(0.);                                         // Uzhi Oct. 2014
  MminProjectile = BrW.GetMinimumMass(projectile->GetDefinition());    // Uzhi Oct. 2014
//G4double M0projectile = MminProjectile;                              // With de-excitation    Uzhi Oct. 2014
//G4double M0projectile = Pprojectile.mag();                           // Without de-excitation, see above

  if ( M0projectile < MminProjectile ) 
  {
    PutOnMassShell = true;
    M0projectile = BrW.SampleMass(projectile->GetDefinition(),projectile->GetDefinition()->GetPDGMass() + 5.0*projectile->GetDefinition()->GetPDGWidth());
  }

  G4double M0projectile2 = M0projectile * M0projectile;
  G4double ProjectileDiffStateMinMass    = theParameters->GetProjMinDiffMass();    // Uzhi Oct 2014
  G4double ProjectileNonDiffStateMinMass = theParameters->GetProjMinNonDiffMass(); // Uzhi Oct 2014
  if ( M0projectile > ProjectileDiffStateMinMass ) {                               // Uzhi Oct 2014
    ProjectileDiffStateMinMass    = M0projectile + 220.0*MeV;
    ProjectileNonDiffStateMinMass = M0projectile + 220.0*MeV;
    if(absProjectilePDGcode > 3000) {                          // Strange baryon   // Uzhi Nov. 2014              
      ProjectileDiffStateMinMass    += 140.0*MeV;                                  // Uzhi Nov. 2014
      ProjectileNonDiffStateMinMass += 140.0*MeV;                                  // Uzhi Nov. 2014
    }                                                                              // Uzhi Nov. 2014
  }

  G4double MminTarget(0.);                                                         // Uzhi Oct. 2014
  MminTarget = BrW.GetMinimumMass(target->GetDefinition());                        // Uzhi Oct. 2014
  if ( M0target < MminTarget ) 
  {
    PutOnMassShell = true;
    M0target = BrW.SampleMass(target->GetDefinition(),target->GetDefinition()->GetPDGMass() + 5.0*target->GetDefinition()->GetPDGWidth());
  }

  G4double M0target2 = M0target * M0target; 
  G4double TargetDiffStateMinMass    = theParameters->GetTarMinDiffMass();         // Uzhi Oct 2014
  G4double TargetNonDiffStateMinMass = theParameters->GetTarMinNonDiffMass();      // Uzhi Oct 2014
  if ( M0target > TargetDiffStateMinMass ) {                                       // Uzhi Oct 2014
    TargetDiffStateMinMass    = M0target + 220.0*MeV;
    TargetNonDiffStateMinMass = M0target + 220.0*MeV;
    if(absTargetPDGcode > 3000) {                          // Strange baryon       // Uzhi Nov. 2014              
      TargetDiffStateMinMass    += 140.0*MeV;                                      // Uzhi Nov. 2014
      TargetNonDiffStateMinMass += 140.0*MeV;                                      // Uzhi Nov. 2014
    }                                                                              // Uzhi Nov. 2014
  };

  #ifdef debugFTFexictation
  G4cout << "Proj Targ PDGcodes " << ProjectilePDGcode << " " << TargetPDGcode << G4endl
         << "M0projectile Y " << M0projectile << " " << ProjectileRapidity << G4endl;
  //G4cout << "M0target     Y " << M0target << " " << TargetRapidity << G4endl;
  G4cout << "Pproj " << Pprojectile << G4endl << "Ptarget " << Ptarget << G4endl;
  #endif

  G4double AveragePt2 = theParameters->GetAveragePt2();
//  G4double ProbLogDistrPrD = theParameters->GetProbLogDistrPrD();                  // Uzhi Oct 2014
  G4double ProbLogDistr = theParameters->GetProbLogDistr();
  G4double SumMasses = M0projectile + M0target; // + 220.0*MeV;                      // Uzhi Nov. 2014

  // Transform momenta to cms and then rotate parallel to z axis;
  G4LorentzRotation toCms( -1 * Psum.boostVector() );

  G4LorentzVector Ptmp = toCms * Pprojectile;
  if ( Ptmp.pz() <= 0.0 ) return false;        // "String" moving backwards in  CMS, abort collision!

  toCms.rotateZ( -1*Ptmp.phi() );
  toCms.rotateY( -1*Ptmp.theta() );
  G4LorentzRotation toLab(toCms.inverse());
  Pprojectile.transform( toCms );
  Ptarget.transform( toCms );

  G4double PZcms2, PZcms;

  #ifdef debugFTFexictation
  G4cout << "SqrtS     " << SqrtS << G4endl << "M0pr M0tr SumM+220 " << M0projectile << " "
         << M0target << " " << SumMasses << G4endl;
  #endif

  if ( SqrtS < M0projectile + M0target ) return false;
  if ( SqrtS < SumMasses ) return false;
  // The model cannot work at low energy

  PZcms2 = ( S*S + M0projectile2*M0projectile2 + M0target2*M0target2
             - 2*S*M0projectile2 - 2*S*M0target2 - 2*M0projectile2*M0target2 ) / 4.0 / S;

  #ifdef debugFTFexictation
  G4cout << "PZcms2 after PutOnMassShell " << PZcms2 << G4endl;
  #endif

  if ( PZcms2 < 0 ) return false;
  // It can be in an interaction with off-shell nuclear nucleon

  PZcms = std::sqrt( PZcms2 );
  if ( PutOnMassShell ) {
    if ( Pprojectile.z() > 0.0 ) {
      Pprojectile.setPz(  PZcms );
      Ptarget.setPz(     -PZcms );
    } else {
      Pprojectile.setPz( -PZcms );
      Ptarget.setPz(      PZcms );
    };
    Pprojectile.setE( std::sqrt( M0projectile2 +
                                 Pprojectile.x()*Pprojectile.x() +
                                 Pprojectile.y()*Pprojectile.y() +
                                 PZcms2 ) );
    Ptarget.setE( std::sqrt( M0target2 +
                             Ptarget.x()*Ptarget.x() +
                             Ptarget.y()*Ptarget.y() +
                             PZcms2 ) );
  }

  G4double maxPtSquare(0.);  // = PZcms2;

  //Uzhi_QEnex = 0;
  //Uzhi_QEexc = 0;
  //Uzhi_targetdiffraction     = 0;
  //Uzhi_projectilediffraction = 0;
  //Uzhi_nondiffraction        = 0;
  //G4int UzhiPrD( 0 ), UzhiTrD( 0 ), UzhiND( 0 );

  #ifdef debugFTFexictation
  G4cout << "Start --------------------" << G4endl << "Proj M0 Mdif Mndif " << M0projectile 
         << " " << ProjectileDiffStateMinMass << "  " << ProjectileNonDiffStateMinMass << G4endl
         << "Targ M0 Mdif Mndif " << M0target << " " << TargetDiffStateMinMass << " " 
         << TargetNonDiffStateMinMass << G4endl << "SqrtS " << SqrtS << G4endl
         << "Proj CMS " << Pprojectile << G4endl << "Targ CMS " << Ptarget << G4endl;
  #endif

  // Charge exchange can be possible
  // Getting the values needed for exchange
  // Check for possible quark exchange
  G4double QeNoExc = theParameters->GetProcProb( 0, ProjectileRapidity );
  G4double QeExc   = theParameters->GetProcProb( 1, ProjectileRapidity )*theParameters->GetProcProb( 4, ProjectileRapidity );
  G4double ProbProjectileDiffraction = theParameters->GetProcProb( 2, ProjectileRapidity );
  G4double ProbTargetDiffraction     = theParameters->GetProcProb( 3, ProjectileRapidity );

  if(QeNoExc+QeExc+ProbProjectileDiffraction+ProbTargetDiffraction > 1.)   // Uzhi Nov. 2014 
    {QeNoExc=1.0-QeExc-ProbProjectileDiffraction-ProbTargetDiffraction;}

  G4double ProbExc( 0.0 ); 
  if ( QeExc + QeNoExc != 0.0 ) ProbExc = QeExc/(QeExc + QeNoExc);
  G4double DeltaProbAtQuarkExchange = theParameters->GetDeltaProbAtQuarkExchange();
  G4double DeltaMass = G4ParticleTable::GetParticleTable()->FindParticle( 2224 )->GetPDGMass();

  G4double ProbOfDiffraction = ProbProjectileDiffraction + ProbTargetDiffraction;

  #ifdef debugFTFexictation
  G4cout << "Proc Probs " << QeNoExc << " " << QeExc << " " << ProbProjectileDiffraction
         << " " << ProbTargetDiffraction << G4endl 
         << "ProjectileRapidity " << ProjectileRapidity << G4endl;
//  G4int Uzhi; G4cin >> Uzhi;
  #endif
 
  G4ParticleDefinition* TestParticle(0);                                   // Uzhi Oct. 2014
  G4double MtestPr(0.), MtestTr(0.);                                       // Uzhi Oct. 2014

  if ( 1.0 - QeExc - QeNoExc > 0.0 ) { 
    ProbProjectileDiffraction /= ( 1.0 - QeExc - QeNoExc );
    ProbTargetDiffraction     /= ( 1.0 - QeExc - QeNoExc );
  }

  if ( G4UniformRand() < QeExc + QeNoExc ) {    

    #ifdef debugFTFexictation
    G4cout << "Q exchange --------------------------" << G4endl;
    #endif

    G4int NewProjCode( 0 ), NewTargCode( 0 );
    G4int ProjQ1( 0 ), ProjQ2( 0 ), ProjQ3( 0 );

    //  Projectile unpacking
    if ( absProjectilePDGcode < 1000 ) {  // projectile is meson 
      UnpackMeson( ProjectilePDGcode, ProjQ1, ProjQ2 );  
    } else {  // projectile is baryon
      UnpackBaryon( ProjectilePDGcode, ProjQ1, ProjQ2, ProjQ3 );
    }

    //  Target unpacking
    G4int TargQ1( 0 ), TargQ2( 0 ), TargQ3( 0 );
    UnpackBaryon( TargetPDGcode, TargQ1, TargQ2, TargQ3 ); 

    #ifdef debugFTFexictation
    G4cout << "Proj Quarks " << ProjQ1 << " " << ProjQ2 << " " << ProjQ3 << G4endl
           << "Targ Quarks " << TargQ1 << " " << TargQ2 << " " << TargQ3 << G4endl;
    #endif

    // Sampling of exchanged quarks
    G4int ProjExchangeQ( 0 );
    G4int TargExchangeQ( 0 );

    if ( absProjectilePDGcode < 1000 ) 
    {                                          // projectile is meson 

      if ( ProjQ1 > 0 ) 
      {                                        // ProjQ1 is quark
        ProjExchangeQ = ProjQ1;
        //------------------------------- Exchange of non-identical quarks is allowed
        G4int NpossibleStates=3; // =====================================================
//
        NpossibleStates=0;       // =====================================================
        if(ProjQ1 != TargQ1) NpossibleStates++;
        if(ProjQ1 != TargQ2) NpossibleStates++;
        if(ProjQ1 != TargQ3) NpossibleStates++;  
//
        G4int Nsampled = G4RandFlat::shootInt( G4long( NpossibleStates ) ) + 1;
//G4cout<<"NpossibleStates Nsampled "<<NpossibleStates<<" "<<Nsampled<<G4endl;
//
        NpossibleStates=0;
        if(ProjQ1 != TargQ1) 
        {
         NpossibleStates++;
         if(NpossibleStates == Nsampled) {TargExchangeQ = TargQ1; TargQ1 = ProjExchangeQ; ProjQ1 = TargExchangeQ;}
        }
        if(ProjQ1 != TargQ2) 
        {
         NpossibleStates++;
         if(NpossibleStates == Nsampled) {TargExchangeQ = TargQ2; TargQ2 = ProjExchangeQ; ProjQ1 = TargExchangeQ;}
        }
        if(ProjQ1 != TargQ3) 
        {
         NpossibleStates++;
         if(NpossibleStates == Nsampled) {TargExchangeQ = TargQ3; TargQ3 = ProjExchangeQ; ProjQ1 = TargExchangeQ;}
        }
//

//if(Nsampled == 1)      {TargExchangeQ = TargQ1; TargQ1 = ProjExchangeQ; ProjQ1 = TargExchangeQ;}
//else if(Nsampled == 2) {TargExchangeQ = TargQ2; TargQ2 = ProjExchangeQ; ProjQ1 = TargExchangeQ;}
//else                   {TargExchangeQ = TargQ3; TargQ3 = ProjExchangeQ; ProjQ1 = TargExchangeQ;}
      } 
      else 
      {                                        // ProjQ2 is quark
        ProjExchangeQ = ProjQ2;
        //------------------------------- Exchange of non-identical quarks is allowed
        G4int NpossibleStates=3; 
//
        NpossibleStates=0;
        if(ProjQ2 != TargQ1) NpossibleStates++;
        if(ProjQ2 != TargQ2) NpossibleStates++;
        if(ProjQ2 != TargQ3) NpossibleStates++;  
//
        G4int Nsampled = G4RandFlat::shootInt( G4long( NpossibleStates ) ) + 1;
//
        NpossibleStates=0;
        if(ProjQ2 != TargQ1) 
        {
         NpossibleStates++;
         if(NpossibleStates == Nsampled) {TargExchangeQ = TargQ1; TargQ1 = ProjExchangeQ; ProjQ2 = TargExchangeQ;}
        }
        if(ProjQ2 != TargQ2) 
        {
         NpossibleStates++;
         if(NpossibleStates == Nsampled) {TargExchangeQ = TargQ2; TargQ2 = ProjExchangeQ; ProjQ2 = TargExchangeQ;}
        }
        if(ProjQ2 != TargQ3) 
        {
         NpossibleStates++;
         if(NpossibleStates == Nsampled) {TargExchangeQ = TargQ3; TargQ3 = ProjExchangeQ; ProjQ2 = TargExchangeQ;}
        }
//

//if(Nsampled == 1)      {TargExchangeQ = TargQ1; TargQ1 = ProjExchangeQ; ProjQ2 = TargExchangeQ;}
//else if(Nsampled == 2) {TargExchangeQ = TargQ2; TargQ2 = ProjExchangeQ; ProjQ2 = TargExchangeQ;}
//else                   {TargExchangeQ = TargQ3; TargQ3 = ProjExchangeQ; ProjQ2 = TargExchangeQ;}
      }        // End of if ( ProjQ1 > 0 )

      #ifdef debugFTFexictation
      G4cout << "Exchanged Qs in Pr Tr " << ProjExchangeQ << " " << TargExchangeQ << G4endl;
      #endif

      G4int aProjQ1 = std::abs( ProjQ1 );
      G4int aProjQ2 = std::abs( ProjQ2 );

      G4bool ProjExcited = false;                                          // Uzhi Oct 2014

      G4int attempts=0;                                                    // Uzhi Oct 2014 start
      while(attempts < 50)
      {// Determination of a new projectile ID which garanty energy-momentum conservation
       attempts++;

       G4double Ksi = G4UniformRand();

       if ( aProjQ1 == aProjQ2 ) 
       {
         if ( aProjQ1 != 3 ) 
         {
           NewProjCode = 111;                       // Pi0-meson
           if ( Ksi < 0.5 ) 
           {
             NewProjCode = 221;                     // Eta -meson
             if ( Ksi < 0.25 ) {NewProjCode = 331;} // Eta'-meson
           }
         } else 
         {
           NewProjCode = 221;                      // Eta -meson
           if( Ksi < 0.5 ) {NewProjCode = 331;}    // Eta'-meson
         }
       } else 
       {
         if ( aProjQ1 > aProjQ2 ) 
         {
           NewProjCode = aProjQ1*100 + aProjQ2*10 + 1;
         } else 
         {
           NewProjCode = aProjQ2*100 + aProjQ1*10 + 1;
         }
       }

       #ifdef debugFTFexictation
       G4cout << "NewProjCode " << NewProjCode << G4endl;
       #endif

       ProjExcited = false;
       if ( G4UniformRand() < 0.5 ) 
       {
         NewProjCode += 2;                         // Excited meson 
         ProjExcited = true;
       }
//       if ( aProjQ1 != aProjQ2 ) NewProjCode *= ( ProjectilePDGcode / absProjectilePDGcode ); // Uzhi 27 Nov. 2014
// Uzhi 27 Nov. 2014
       G4int Qquarks=0;
       if     ( aProjQ1 == 1 ) {Qquarks -= ProjQ1;}
       else if( aProjQ1 == 2 ) {Qquarks += ProjQ1;}
       else                    {Qquarks -= ProjQ1/aProjQ1;}

       if     ( aProjQ2 == 1 ) {Qquarks -= ProjQ2;}
       else if( aProjQ2 == 2 ) {Qquarks += ProjQ2;}
       else                    {Qquarks -= ProjQ2/aProjQ2;}

       if( Qquarks < 0 ) NewProjCode *=(-1);
// Uzhi 27 Nov. 2014

       #ifdef debugFTFexictation
       G4cout << "NewProjCode +2 or 0 " << NewProjCode << G4endl;
       G4cout<<"+++++++++++++++++++++++++++++++++++++++"<<G4endl;
       G4cout<<ProjQ1<<" "<<ProjQ2<<" "<<Qquarks<<G4endl;
       G4cout<<"+++++++++++++++++++++++++++++++++++++++"<<G4endl;
       #endif

//     --------------------------------------------------------------------------------- Proj 
       TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewProjCode );
       if(!TestParticle) continue;

//MminProjectile=TestParticle->GetPDGMass();            // ??????????????????????
       MminProjectile=BrW.GetMinimumMass(TestParticle);

       if(SqrtS-M0target < MminProjectile) continue;

       MtestPr = BrW.SampleMass(TestParticle, TestParticle->GetPDGMass() + 5.0*TestParticle->GetPDGWidth());
//            G4ParticleTable::GetParticleTable()->FindParticle( NewProjCode )->GetPDGMass();  // Uzhi 2014 

       #ifdef debugFTFexictation
       G4cout << "TestParticle Name " << NewProjCode << " " << TestParticle->GetParticleName()<< G4endl;
       G4cout << "MtestPart MtestPart0 "<<MtestPr<<" "<<TestParticle->GetPDGMass()<<G4endl;
       G4cout << "M0projectile projectile PDGMass " << M0projectile << " " 
              << projectile->GetDefinition()->GetPDGMass() << G4endl;
       #endif

//     --------------------------------------------------------------------------------- Targ
       NewTargCode = NewNucleonId( TargQ1, TargQ2, TargQ3 );

       #ifdef debugFTFexictation
       G4cout << "New TrQ " << TargQ1 << " " << TargQ2 << " " << TargQ3 << G4endl
              << "NewTargCode " << NewTargCode << G4endl;
       #endif

       if( TargQ1 != TargQ2  &&  TargQ1 != TargQ3  &&  TargQ2 != TargQ3 ) 
       {                                                               // Lambda or Sigma0 ???
         if   ( G4UniformRand() < 0.5 ) {NewTargCode+=2;}
         else { if ( G4UniformRand() < 0.75 ) NewTargCode=3122;}
       } 
       else if( TargQ1 == TargQ2  &&  TargQ1 == TargQ3 ) 
       {
               NewTargCode += 2; ProjExcited = true;                   //Create Delta isobar
       } else if ( target->GetDefinition()->GetPDGiIsospin() == 3 ) {  // Delta was the target
         if ( G4UniformRand() > DeltaProbAtQuarkExchange ) 
         {NewTargCode += 2; ProjExcited = true;}                       // Save Delta isobar
         else {}                                                       // De-excite initial Delta isobar
       } else if ( ! ProjExcited  &&
                   G4UniformRand() < DeltaProbAtQuarkExchange  &&      // Nucleon was the target
                   SqrtS > M0projectile + DeltaMass ) {                // Create Delta isobar
         NewTargCode +=2;                                              // Save initial nucleon
       } else {}

       TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewTargCode );

       if(!TestParticle) continue;
       
       #ifdef debugFTFexictation
       G4cout << "New targ " << NewTargCode << " " << TestParticle->GetParticleName() << G4endl;
       #endif

       MminTarget=BrW.GetMinimumMass(TestParticle);

       if(SqrtS-MtestPr < MminTarget) continue;
      
       MtestTr = BrW.SampleMass(TestParticle,TestParticle->GetPDGMass() + 5.0*TestParticle->GetPDGWidth());

       if(SqrtS > MtestPr+MtestTr) break;
      }  // End of while(attempts < 50)//===============================

      if(attempts >= 50) return false; // ==============================
/*
      if ( MtestPr > Pprojectile.mag() ) {M0projectile = MtestPr;}
      else 
      {
        if ( std::abs( MtestPr - M0projectile )  //projectile->GetDefinition()->GetPDGMass() ) // Uzhi Oct. 2014
             < 140.0*MeV ) {
            M0projectile = MtestPr;
        }
      }
*/
      if ( MtestPr >= Pprojectile.mag() ) {M0projectile = MtestPr;}                            // Uzhi 18 Nov. 2014
      else if (projectile->GetStatus() != 0  ) {M0projectile = MtestPr;}                       // Uzhi 18 Nov. 2014


      #ifdef debugFTFexictation
      G4cout << "M0projectile After check " << M0projectile << G4endl;
      #endif

      M0projectile2 = M0projectile * M0projectile;
      ProjectileDiffStateMinMass    = M0projectile + 220.0*MeV;  //220 MeV=m_pi+80 MeV 
      ProjectileNonDiffStateMinMass = M0projectile + 220.0*MeV;  //220 MeV=m_pi+80 MeV

      if ( MtestTr >= Ptarget.mag()     ) {M0target = MtestTr;}                                // Uzhi 18 Nov. 2014
      else if (target->GetStatus() != 0 ) {M0target = MtestTr;}                                // Uzhi 18 Nov. 2014
/*
M0target = MtestTr;                                                                            // Uzhi 18 Nov. 2014
      if ( MtestTr > Ptarget.mag() ) {M0target = MtestTr;}
      else 
      {
        if ( std::abs( MtestTr - M0target )             // target->GetDefinition()->GetPDGMass() ) Uzhi Oct. 2014
             < 140.0*MeV ) {
            M0target = MtestTr;
        }
      }
*/
      M0target2 = M0target * M0target;

      #ifdef debugFTFexictation
      G4cout << "New targ M0 M0^2 " << M0target << " " << M0target2 << G4endl;
      #endif

      TargetDiffStateMinMass    = M0target + 220.0*MeV;          // 220 MeV=m_pi+80 MeV;    
      TargetNonDiffStateMinMass = M0target + 220.0*MeV;          // 220 MeV=m_pi+80 MeV; 

    } else {  // of the if ( absProjectilePDGcode < 1000 ) ;  
              // The projectile is baryon now

      G4double Same = theParameters->GetProbOfSameQuarkExchange();  //0.3; //0.5; 0.
//      G4bool ProjDeltaHasCreated( false );                     // Uzhi Oct. 2014
//      G4bool TargDeltaHasCreated( false );                     // Uzhi Oct. 2014
 
      G4double Ksi = G4UniformRand();
      if ( G4UniformRand() < 0.5 ) {  // Sampling exchange quark from proj. or targ.
        // Sampling exchanged quark from the projectile

        if ( Ksi < 0.333333 ) {
          ProjExchangeQ = ProjQ1;
        } else if ( 0.333333 <= Ksi  &&  Ksi < 0.666667 ) {
          ProjExchangeQ = ProjQ2;
        } else {
          ProjExchangeQ = ProjQ3;
        }

        if ( ProjExchangeQ != TargQ1  ||  G4UniformRand() < Same ) {
          TargExchangeQ = TargQ1; TargQ1 = ProjExchangeQ; ProjExchangeQ = TargExchangeQ;
        } else {
          if ( ProjExchangeQ != TargQ2  ||  G4UniformRand() < Same ) {
            TargExchangeQ = TargQ2; TargQ2 = ProjExchangeQ; ProjExchangeQ = TargExchangeQ;
          } else {
            TargExchangeQ = TargQ3; TargQ3 = ProjExchangeQ; ProjExchangeQ = TargExchangeQ;
          }
        }

        #ifdef debugFTFexictation
        G4cout << "Exchange Qs Pr  Tr " << ProjExchangeQ << " " << TargExchangeQ << G4endl;
        #endif

        if ( Ksi < 0.333333 ) {
          ProjQ1 = ProjExchangeQ;
        } else if ( 0.333333 <= Ksi  &&  Ksi < 0.666667 ) {
          ProjQ2 = ProjExchangeQ;
        } else {
          ProjQ3 = ProjExchangeQ;
        }

      } else {  // Sampling exchanged quark from the target

        if ( Ksi < 0.333333 ) {
          TargExchangeQ = TargQ1;
        } else if ( 0.333333 <= Ksi  &&  Ksi < 0.666667 ) {
          TargExchangeQ = TargQ2;
        } else {
          TargExchangeQ = TargQ3;
        }
        if ( TargExchangeQ != ProjQ1  ||  G4UniformRand() < Same ) {
          ProjExchangeQ = ProjQ1; ProjQ1 = TargExchangeQ; TargExchangeQ = ProjExchangeQ;
        } else {
          if ( TargExchangeQ != ProjQ2  ||  G4UniformRand() < Same ) {
            ProjExchangeQ = ProjQ2; ProjQ2 = TargExchangeQ; TargExchangeQ = ProjExchangeQ;
          } else {
            ProjExchangeQ = ProjQ3; ProjQ3 = TargExchangeQ; TargExchangeQ = ProjExchangeQ;
          }
        }

        if ( Ksi < 0.333333 ) {
          TargQ1 = TargExchangeQ;
        } else if ( 0.333333 <= Ksi  &&  Ksi < 0.666667 )  {
          TargQ2 = TargExchangeQ;
        } else {
          TargQ3 = TargExchangeQ;
        }

      } // End of quark sampling for the baryons

      NewProjCode = NewNucleonId( ProjQ1, ProjQ2, ProjQ3 );
      NewTargCode = NewNucleonId( TargQ1, TargQ2, TargQ3 );

      G4int attempts=0;                                                    // Uzhi Oct 2014 start
      while(attempts < 50)
      {// Determination of a new projectile ID which garanty energy-momentum conservation
       attempts++;

       if ( ProjQ1 == ProjQ2  &&  ProjQ1 == ProjQ3 ) {
         NewProjCode += 2; // ProjDeltaHasCreated = true;        // Uzhi Oct. 2014
       } else if ( projectile->GetDefinition()->GetPDGiIsospin() == 3 ) {  // Projectile was Delta
         if ( G4UniformRand() > DeltaProbAtQuarkExchange ) { 
           NewProjCode += 2; //ProjDeltaHasCreated = true;       // Uzhi Oct. 2014
         } else {
           NewProjCode += 0; //ProjDeltaHasCreated = false;      // Uzhi Oct. 2014
         }
       } else {  // Projectile was Nucleon
         if ( G4UniformRand() < DeltaProbAtQuarkExchange  &&  SqrtS > DeltaMass + M0target ) {
           NewProjCode += 2; //ProjDeltaHasCreated = true;       // Uzhi Oct. 2014
         } else { 
           NewProjCode += 0; //ProjDeltaHasCreated = false;      // Uzhi Oct. 2014
         }
       } 

       if ( TargQ1 == TargQ2  &&  TargQ1 == TargQ3 ) { 
         NewTargCode += 2; //TargDeltaHasCreated = true;         // Uzhi Oct. 2014
       } else if ( target->GetDefinition()->GetPDGiIsospin() == 3 ) {  // Target was Delta
         if ( G4UniformRand() > DeltaProbAtQuarkExchange ) {
           NewTargCode += 2; //TargDeltaHasCreated = true;       // Uzhi Oct. 2014
         } else {
           NewTargCode += 0; //TargDeltaHasCreated = false;      // Uzhi Oct. 2014
         }
       } else {  // Target was Nucleon
         if ( G4UniformRand() < DeltaProbAtQuarkExchange  &&  SqrtS > M0projectile + DeltaMass ) {
           NewTargCode += 2; //TargDeltaHasCreated = true;       // Uzhi Oct. 2014
         } else {
           NewTargCode += 0; //TargDeltaHasCreated = false;      // Uzhi Oct. 2014
         }
       }         

       #ifdef debugFTFexictation
       G4cout << "NewProjCode NewTargCode " << NewProjCode << " " << NewTargCode << G4endl;
//       G4int Uzhi; G4cin >> Uzhi;
       #endif

       if ( absProjectilePDGcode == NewProjCode  &&  absTargetPDGcode == NewTargCode ) { 
       }  // Nothing was changed! It is not right!?
       
       // Forming baryons
/*
       if ( G4UniformRand() > 0.5 ) {                                      // Uzhi Oct. 2014
         ProbProjectileDiffraction = 0.0; ProbTargetDiffraction = 1.0;
       } else {
         ProbProjectileDiffraction = 1.0; ProbTargetDiffraction = 0.0;
       }
*/
/*
       if ( ProjDeltaHasCreated ) {
         if ( G4UniformRand() > 0.5 ) {
           ProbProjectileDiffraction = 0.0; ProbTargetDiffraction = 1.0;
         } else {
           ProbProjectileDiffraction = 1.0; ProbTargetDiffraction = 0.0;
         }
       }

       if ( TargDeltaHasCreated ) {
         if ( G4UniformRand() > 0.5 ) {
           ProbProjectileDiffraction = 1.0; ProbTargetDiffraction = 0.0;
         } else {
           ProbProjectileDiffraction = 0.0; ProbTargetDiffraction = 1.0;
         }
       }
*/
//     --------------------------------------------------------------------------------- Proj 
       TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewProjCode );
       if(!TestParticle) continue;

       MminProjectile=BrW.GetMinimumMass(TestParticle);

       if(SqrtS-M0target < MminProjectile) continue;

       MtestPr = BrW.SampleMass(TestParticle,TestParticle->GetPDGMass() + 5.0*TestParticle->GetPDGWidth()); 

//     --------------------------------------------------------------------------------- Targ 
       TestParticle = G4ParticleTable::GetParticleTable()->FindParticle( NewTargCode );
       if(!TestParticle) continue;

       MminTarget=BrW.GetMinimumMass(TestParticle);

       if(SqrtS-MtestPr < MminTarget) continue;

       MtestTr = BrW.SampleMass(TestParticle,TestParticle->GetPDGMass() + 5.0*TestParticle->GetPDGWidth());

       if(SqrtS > MtestPr+MtestTr) break;
      }  // End of while(attempts < 50)//===============================

      if(attempts >= 50) return false; // ==============================

      if ( MtestPr >= Pprojectile.mag() ) {M0projectile = MtestPr;}
      else if (projectile->GetStatus() != 0  ) {M0projectile = MtestPr;}        // Uzhi 18 Nov. 2014
      M0projectile2 = M0projectile * M0projectile;
      ProjectileDiffStateMinMass =    M0projectile + 220.0*MeV;  //220 MeV=m_pi+80 MeV
      ProjectileNonDiffStateMinMass = M0projectile + 220.0*MeV;  //220 MeV=m_pi+80 MeV

      if ( MtestTr >= Ptarget.mag()     ) {M0target = MtestTr;}
      else if (target->GetStatus() != 0 ) {M0target = MtestTr;}                 // Uzhi 18 Nov. 2014
      M0target2 = M0target * M0target;
      TargetDiffStateMinMass    = M0target + 220.0*MeV;          //220 MeV=m_pi+80 MeV;    
      TargetNonDiffStateMinMass = M0target + 220.0*MeV;          //220 MeV=m_pi+80 MeV; 

    }  // End of if ( absProjectilePDGcode < 1000 )
//--------------------------------------------------------------------------------------

    // If we assume that final state hadrons after the charge exchange will be
    // in the ground states, we have to put 
    if ( SqrtS < M0projectile + M0target ) return false;

    PZcms2 = ( S*S + M0projectile2*M0projectile2 + M0target2*M0target2
               - 2*S*M0projectile2 - 2*S*M0target2 - 2*M0projectile2*M0target2 ) / 4.0 / S;

    #ifdef debugFTFexictation
    G4cout << "At the end// NewProjCode " << NewProjCode << G4endl 
           << "At the end// NewTargCode " << NewTargCode << G4endl
           << "M0pr  M0tr  SqS " << M0projectile << " " << M0target << " " << SqrtS << G4endl
           << "M0pr2 M0tr2 SqS " << M0projectile2 << " " << M0target2 << " " << SqrtS << G4endl
           << "PZcms2 after the change " << PZcms2 << G4endl << G4endl;
    #endif

    if ( PZcms2 < 0 ) return false;  // It can be if energy is not sufficient for Delta

    projectile->SetDefinition( G4ParticleTable::GetParticleTable()->FindParticle( NewProjCode ) ); 
    target->SetDefinition( G4ParticleTable::GetParticleTable()->FindParticle( NewTargCode ) ); 

    PZcms = std::sqrt( PZcms2 );
    Pprojectile.setPz( PZcms );
    Pprojectile.setE( std::sqrt( M0projectile2 + PZcms2 ) );
    Ptarget.setPz(    -PZcms );
    Ptarget.setE( std::sqrt( M0target2 + PZcms2 ) );

    if(projectile->GetStatus() != 0 ) projectile->SetStatus(2);                                    // Uzhi 18 Nov. 2014
    if(target->GetStatus()     != 0 ) target->SetStatus(2);                                        // Uzhi 18 Nov. 2014

    #ifdef debugFTFexictation
    G4cout << "Proj Targ and Proj+Targ in CMS" << G4endl << Pprojectile << G4endl << Ptarget
           << G4endl << Pprojectile + Ptarget << G4endl;
    #endif

//--------------------- Check for possible excitation of the participants -------------------
    if((SqrtS < M0projectile + TargetDiffStateMinMass) ||             // Uzhi Oct 2014
       (SqrtS < ProjectileDiffStateMinMass + M0target) ||
       (ProbOfDiffraction == 0.)                         ) ProbExc=0.;// Uzhi Oct 2014

    if ( G4UniformRand() > ProbExc ) {  // Make elastic scattering

      #ifdef debugFTFexictation
      G4cout << "Make elastic scattering of new hadrons" << G4endl;
      #endif
 
      Pprojectile.transform( toLab );
      Ptarget.transform( toLab );

      projectile->Set4Momentum( Pprojectile );
      target->Set4Momentum( Ptarget );

      G4bool Result = theElastic->ElasticScattering( projectile, target, theParameters );

      #ifdef debugFTFexictation
      G4cout << "Result of el. scatt " << Result << G4endl << "Proj Targ and Proj+Targ in Lab"
             << G4endl << projectile->Get4Momentum() << G4endl << target->Get4Momentum() << G4endl
             << projectile->Get4Momentum() + target->Get4Momentum() << " " << (projectile->Get4Momentum() + target->Get4Momentum()).mag() << G4endl;
      #endif

//Uzhi_QEnex++;
      return Result;
    }
//Uzhi_QEexc++;

    #ifdef debugFTFexictation
    G4cout << "Make excitation of new hadrons" << G4endl;
    #endif

// Redefinition of ProbOfDiffraction because the probabilities are changed due to quark exchange

    ProbOfDiffraction = ProbProjectileDiffraction + ProbTargetDiffraction; // Uzhi Oct. 2014
    if ( ProbOfDiffraction != 0.0 ) {
      ProbProjectileDiffraction /= ProbOfDiffraction;
      ProbTargetDiffraction     /= ProbOfDiffraction;
      ProbOfDiffraction=1.0;
    }
//Uzhi_QEnex++;
  }  // End of if ( G4UniformRand() < QeExc + QeNoExc ) , i.e. of the charge exchange part

  ProbOfDiffraction = ProbProjectileDiffraction + ProbTargetDiffraction;

  #ifdef debugFTFexictation
  G4cout << "Excitation --------------------" << G4endl
         << "Proj M0 MdMin MndMin " << M0projectile << " " << ProjectileDiffStateMinMass << "  "
         << ProjectileNonDiffStateMinMass << G4endl
         << "Targ M0 MdMin MndMin " << M0target << " " << TargetDiffStateMinMass << " " 
         << TargetNonDiffStateMinMass << G4endl << "SqrtS " << SqrtS << G4endl
         << "Prob: ProjDiff TargDiff + Sum " << ProbProjectileDiffraction << " " 
         << ProbTargetDiffraction << " " << ProbOfDiffraction << G4endl;
  #endif

  if ( ProbOfDiffraction != 0.0 ) {
    ProbProjectileDiffraction /= ProbOfDiffraction;
  } else {
    ProbProjectileDiffraction = 0.0;
  }

  #ifdef debugFTFexictation
  G4cout << "Prob: ProjDiff TargDiff + Sum " << ProbProjectileDiffraction << " " 
         << ProbTargetDiffraction << " " << ProbOfDiffraction << G4endl;
  #endif

  G4double ProjectileDiffStateMinMass2    = sqr( ProjectileDiffStateMinMass );
  G4double ProjectileNonDiffStateMinMass2 = sqr( ProjectileNonDiffStateMinMass );
  G4double TargetDiffStateMinMass2        = sqr( TargetDiffStateMinMass );
  G4double TargetNonDiffStateMinMass2     = sqr( TargetNonDiffStateMinMass );

  G4double Pt2;
  G4double ProjMassT2, ProjMassT;
  G4double TargMassT2, TargMassT;
  G4double PMinusMin, PMinusMax;
  //G4double PPlusMin , PPlusMax;
  G4double TPlusMin, TPlusMax;
  G4double PMinusNew, PPlusNew, TPlusNew, TMinusNew;
  G4LorentzVector Qmomentum;
  G4double Qminus, Qplus;
  G4int whilecount = 0;  

  // Choose a process
  if ( G4UniformRand() < ProbOfDiffraction ) {

    if ( G4UniformRand() < ProbProjectileDiffraction ) {  // projectile diffraction

      #ifdef debugFTFexictation
      G4cout << "projectile diffraction" << G4endl;
      #endif

      //UzhiPrD++;

      do {  // while ( ( Pprojectile + Qmomentum ).mag2() <  ProjectileDiffStateMinMass2 )

        //Uzhi_projectilediffraction = 1;
        //Uzhi_targetdiffraction = 0;
        //Uzhi_Mx2 = 1.0;

        // Generate pt and mass of projectile

        whilecount++;
        if ( whilecount > 1000 ) {
          Qmomentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
          return false;  //  Ignore this interaction
        };

        // Check that the interaction is possible
        ProjMassT2 = ProjectileDiffStateMinMass2;
        ProjMassT  = ProjectileDiffStateMinMass;
        TargMassT2 = M0target2;
        TargMassT  = M0target;
        if ( SqrtS < ProjMassT + TargMassT ) return false;

        PZcms2 =( S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2
                  - 2.0*S*ProjMassT2 - 2.0*S*TargMassT2 - 2.0*ProjMassT2*TargMassT2 ) / 4.0 / S;

        if ( PZcms2 < 0 ) return false; 

        maxPtSquare = PZcms2;

        Qmomentum = G4LorentzVector( GaussianPt( AveragePt2, maxPtSquare ), 0 );

        Pt2 = G4ThreeVector( Qmomentum.vect() ).mag2();
        ProjMassT2 = ProjectileDiffStateMinMass2 + Pt2;
        ProjMassT = std::sqrt( ProjMassT2 );
        TargMassT2 = M0target2 + Pt2;
        TargMassT = std::sqrt( TargMassT2 );
        if ( SqrtS < ProjMassT + TargMassT ) continue;

        PZcms2 = ( S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2
                   - 2.0*S*ProjMassT2 - 2.0*S*TargMassT2 - 2.0*ProjMassT2*TargMassT2 ) / 4.0 / S;

        if ( PZcms2 < 0 ) continue;

        PZcms = std::sqrt( PZcms2 );
        PMinusMin = std::sqrt( ProjMassT2 + PZcms2 ) - PZcms;
        PMinusMax = SqrtS - TargMassT;

        PMinusNew = ChooseP( PMinusMin, PMinusMax );

        TMinusNew = SqrtS - PMinusNew;
        Qminus = Ptarget.minus() - TMinusNew;
        TPlusNew = TargMassT2 / TMinusNew;
        Qplus = Ptarget.plus() - TPlusNew;
        Qmomentum.setPz( (Qplus - Qminus)/2 );
        Qmomentum.setE(  (Qplus + Qminus)/2 );

      } while ( ( Pprojectile + Qmomentum ).mag2() <  ProjectileDiffStateMinMass2 ); 
                // Repeat the sampling because there was not any excitation

//      projectile->SetStatus( 1*projectile->GetStatus() );                    // Uzhi Oct 2014

      if(projectile->GetStatus() == 2) projectile->SetStatus(1);               // Uzhi Oct 2014
      if((target->GetStatus() == 1) && (target->GetSoftCollisionCount() == 0)) // Uzhi Oct 2014
          target->SetStatus(2);                                                // Uzhi Oct 2014

    } else {  // Target diffraction

      #ifdef debugFTFexictation
      G4cout << "Target diffraction" << G4endl;
      #endif

      //UzhiTrD++;

      do {  // while ( ( Ptarget - Qmomentum ).mag2() <  TargetDiffStateMinMass2 )

        //Uzhi_projectilediffraction = 0;
        //Uzhi_targetdiffraction = 1;
        //Uzhi_Mx2 = 1.0;

        // Generate pt and target mass

        whilecount++;
        if ( whilecount > 1000 ) {
          Qmomentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
          return false;  //  Ignore this interaction
        };

        // Check that the interaction is possible
        ProjMassT2 = M0projectile2;
        ProjMassT  = M0projectile;
 
        TargMassT2 = TargetDiffStateMinMass2;
        TargMassT  = TargetDiffStateMinMass;

        if ( SqrtS < ProjMassT + TargMassT ) return false;

        PZcms2 = ( S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2
                   - 2.0*S*ProjMassT2 - 2.0*S*TargMassT2 - 2.0*ProjMassT2*TargMassT2 ) / 4.0 / S;

        if ( PZcms2 < 0 ) return false; 

        maxPtSquare = PZcms2;

        Qmomentum = G4LorentzVector( GaussianPt( AveragePt2, maxPtSquare ), 0 );

        Pt2 = G4ThreeVector( Qmomentum.vect() ).mag2();
        ProjMassT2 = M0projectile2 + Pt2;
        ProjMassT  = std::sqrt( ProjMassT2 );
        TargMassT2 = TargetDiffStateMinMass2 + Pt2;
        TargMassT  = std::sqrt( TargMassT2 );
        if ( SqrtS < ProjMassT + TargMassT ) continue;

        PZcms2 = ( S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2
                   - 2.0*S*ProjMassT2 - 2.0*S*TargMassT2 - 2.0*ProjMassT2*TargMassT2 ) / 4.0 / S;

        if ( PZcms2 < 0 ) continue;

        PZcms = std::sqrt( PZcms2 );
        TPlusMin = std::sqrt( TargMassT2 + PZcms2 ) - PZcms;
//TPlusMax = std::sqrt( TargMassT2 + PZcms2 ) + PZcms;
        TPlusMax = SqrtS - ProjMassT;

        TPlusNew = ChooseP( TPlusMin, TPlusMax );
//TPlusNew = TPlusMin;

        PPlusNew = SqrtS - TPlusNew;
        Qplus = PPlusNew - Pprojectile.plus();
        PMinusNew = ProjMassT2 / PPlusNew;
        Qminus = PMinusNew - Pprojectile.minus();
        Qmomentum.setPz( (Qplus - Qminus)/2 );
        Qmomentum.setE(  (Qplus + Qminus)/2 );

      } while ( ( Ptarget - Qmomentum ).mag2() <  TargetDiffStateMinMass2 );
                 // Repeat the sampling because there was not any excitation

//      target->SetStatus( 1*target->GetStatus() );                                    // Uzhi Oct 2014

      if((projectile->GetStatus() == 1) && (projectile->GetSoftCollisionCount() == 0)) // Uzhi Oct 2014
          projectile->SetStatus(2);                                                    // Uzhi Oct 2014
      if(target->GetStatus() == 2) target->SetStatus(1);                               // Uzhi Oct 2014
    
    } // End of if ( G4UniformRand() < ProbProjectileDiffraction )
  
  } else {  // Non-diffraction process

    #ifdef debugFTFexictation
    G4cout << "Non-diffraction process" << G4endl;
    #endif

//UzhiND++;
//Uzhi_QEnex++;
    do {  // while ( ( Pprojectile + Qmomentum ).mag2() <  ProjectileNonDiffStateMinMass2 || ...

      //Uzhi_projectilediffraction = 0;
      //Uzhi_targetdiffraction = 0;
      //Uzhi_Mx2 = 1.0;

      // Generate pt and masses

      whilecount++;
      if ( whilecount > 1000 ) {
        Qmomentum = G4LorentzVector( 0.0, 0.0, 0.0, 0.0 );
        return false;  // Ignore this interaction
      };

      // Check that the interaction is possible
      ProjMassT2 = ProjectileNonDiffStateMinMass2;
      ProjMassT  = ProjectileNonDiffStateMinMass;
      TargMassT2 = TargetNonDiffStateMinMass2;
      TargMassT  = TargetNonDiffStateMinMass;
      if ( SqrtS < ProjMassT + TargMassT ) return false;

      PZcms2 = ( S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2
                 - 2.0*S*ProjMassT2 - 2.0*S*TargMassT2 - 2.0*ProjMassT2*TargMassT2 ) / 4.0 / S;

      if ( PZcms2 < 0 ) return false; 

      maxPtSquare = PZcms2;

      Qmomentum = G4LorentzVector( GaussianPt( AveragePt2, maxPtSquare ), 0 );

      Pt2 = G4ThreeVector( Qmomentum.vect() ).mag2();
      ProjMassT2 = ProjectileNonDiffStateMinMass2 + Pt2;
      ProjMassT  = std::sqrt( ProjMassT2 );
      TargMassT2 = TargetNonDiffStateMinMass2 + Pt2;
      TargMassT  = std::sqrt( TargMassT2 );
      if ( SqrtS < ProjMassT + TargMassT ) continue;

      PZcms2 =( S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2
                -2.0*S*ProjMassT2 - 2.0*S*TargMassT2 -2.0*ProjMassT2*TargMassT2 ) / 4.0 / S;

      if ( PZcms2 < 0 ) continue;

      PZcms = std::sqrt( PZcms2 );
      PMinusMin = std::sqrt( ProjMassT2 + PZcms2 ) - PZcms;
//PMinusMax = std::sqrt( ProjMassT2 + PZcms2 ) + PZcms;
      PMinusMax = SqrtS - TargMassT;

      TPlusMin = std::sqrt( TargMassT2 + PZcms2 ) - PZcms;
//TPlusMax = std::sqrt( TargMassT2 + PZcms2 ) + PZcms;
      TPlusMax = SqrtS - ProjMassT;                              // Uzhi 18 Sept. 2014

/*
      if ( G4UniformRand() < ProbLogDistrPrD ) {                 // Uzhi Oct 2014
        PMinusNew = ChooseP( PMinusMin, PMinusMax );
      } else {
        PMinusNew = ( PMinusMax - PMinusMin )*G4UniformRand() + PMinusMin;
      }
//      Qminus = PMinusNew - Pprojectile.minus();
*/

//      TPlusMax = SqrtS - PMinusNew;                   

      if ( G4UniformRand() < ProbLogDistr ) {
        PMinusNew = ChooseP( PMinusMin, PMinusMax );
        TPlusNew = ChooseP( TPlusMin, TPlusMax );
      } else {
        PMinusNew = ( PMinusMax - PMinusMin )*G4UniformRand() + PMinusMin;
        TPlusNew = ( TPlusMax - TPlusMin )*G4UniformRand() + TPlusMin;
      } 

      Qminus = PMinusNew - Pprojectile.minus();

      Qplus = -( TPlusNew - Ptarget.plus() );
      Qmomentum.setPz( (Qplus - Qminus)/2 );
      Qmomentum.setE(  (Qplus + Qminus)/2 );

      #ifdef debugFTFexictation
      G4cout << ( Pprojectile + Qmomentum ).mag2() << " " << ProjectileNonDiffStateMinMass2
             << G4endl << ( Ptarget - Qmomentum ).mag2() << " "
             << TargetNonDiffStateMinMass2 << G4endl;
      G4cout<<"To continue - enter any integer"<<G4endl;
//      G4int Uzhi; G4cin >> Uzhi;
      #endif

    } while ( ( Pprojectile + Qmomentum ).mag2() <  ProjectileNonDiffStateMinMass2  ||  //No double Diffraction
              ( Ptarget     - Qmomentum ).mag2() <  TargetNonDiffStateMinMass2      ||  // ); //
              ( Pprojectile + Qmomentum ).pz()   <  0.);

    projectile->SetStatus( 0*projectile->GetStatus() );
    target->SetStatus( 0*target->GetStatus() );

  }  // End of if ( G4UniformRand() < ProbOfDiffraction )

  Pprojectile += Qmomentum;
  Ptarget     -= Qmomentum;

  // Transform back and update SplitableHadron Participant.
  Pprojectile.transform( toLab );
  Ptarget.transform( toLab );

  // Calculation of the creation time
  //Uzhi 9.11     projectile->SetTimeOfCreation( target->GetTimeOfCreation() );
  //Uzhi 9.11     projectile->SetPosition( target->GetPosition() );
  // Creation time and position of target nucleon were determined in
  // ReggeonCascade() of G4FTFModel
  // 
  //if ( Uzhi_projectilediffraction != 0 ) { 
  //  Uzhi_Mx2 = Pprojectile.mag2(); Uzhi_modT = ( target->Get4Momentum() - Ptarget ).mag2();
  //}
  //if ( Uzhi_targetdiffraction != 0 ) { 
  //  Uzhi_Mx2 = Ptarget.mag2(); Uzhi_modT = ( projectile->Get4Momentum() - Pprojectile ).mag2();
  //}
  //if ( Uzhi_QE != 0 ) {
  //  Uzhi_projectilediffraction = 0;
  //  Uzhi_targetdiffraction     = 0;
  //  Uzhi_Mx2                   = 1.0;
  //}

  #ifdef debugFTFexictation
  G4cout << "Mproj " << Pprojectile.mag() << G4endl << "Mtarg " << Ptarget.mag() << G4endl;
  #endif

  projectile->Set4Momentum( Pprojectile );
  target->Set4Momentum( Ptarget );
  projectile->IncrementCollisionCount( 1 );
  target->IncrementCollisionCount( 1 );

  //Uzhi_projectilediffraction = UzhiPrD;
  //Uzhi_targetdiffraction = UzhiTrD;
  //Uzhi_nondiffraction = UzhiND;
  //G4cout << Uzhi_projectilediffraction << " " << Uzhi_targetdiffraction << " "
  //       << Uzhi_nondiffraction << G4endl;

  return true;
}


//============================================================================

void G4DiffractiveExcitation::CreateStrings( G4VSplitableHadron* hadron, 
                                             G4bool isProjectile,
                                             G4ExcitedString*& FirstString, 
                                             G4ExcitedString*& SecondString,
                                             G4FTFParameters* theParameters ) const {

  //G4cout << "Create Strings SplitUp " << hadron << G4endl
  //       << "Defin " << hadron->GetDefinition() << G4endl
  //       << "Defin " << hadron->GetDefinition()->GetPDGEncoding() << G4endl;

  hadron->SplitUp();

  G4Parton* start = hadron->GetNextParton();
  if ( start == NULL ) { 
    G4cout << " G4FTFModel::String() Error: No start parton found" << G4endl;
    FirstString = 0; SecondString = 0;
    return;
  }

  G4Parton* end = hadron->GetNextParton();
  if ( end == NULL ) { 
    G4cout << " G4FTFModel::String() Error: No end parton found" <<  G4endl;
    FirstString = 0; SecondString = 0;
    return;
  }

  //G4cout << start << " " << start->GetPDGcode() << " " << end << " " << end->GetPDGcode()
  //       << G4endl
  //       << "Create string " << start->GetPDGcode() << " " << end->GetPDGcode() << G4endl;

  G4LorentzVector Phadron = hadron->Get4Momentum();
  //G4cout << "String mom " << Phadron << G4endl;
  G4LorentzVector Pstart( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector Pend( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector Pkink( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PkinkQ1( 0.0, 0.0, 0.0, 0.0 );
  G4LorentzVector PkinkQ2( 0.0, 0.0, 0.0, 0.0 );

  G4int PDGcode_startQ = std::abs( start->GetDefinition()->GetPDGEncoding() );
  G4int PDGcode_endQ   = std::abs(   end->GetDefinition()->GetPDGEncoding() );
  //G4cout << "PDGcode_startQ " << PDGcode_startQ << " PDGcode_endQ   " << PDGcode_endQ << G4endl;

  G4double Wmin( 0.0 );
  if ( isProjectile ) {
    Wmin = theParameters->GetProjMinDiffMass();
  } else {
    Wmin = theParameters->GetTarMinDiffMass();
  }

  G4double W = hadron->Get4Momentum().mag();
  //G4cout << "Wmin W " << Wmin << " " << W << G4endl;
  //G4int Uzhi; G4cin >> Uzhi;
  G4double W2 = W*W;
  G4double Pt( 0.0 ), x1( 0.0 ), x3( 0.0 );  // x2( 0.0 ) 
  G4bool Kink = false;

  if ( ! ( ( start->GetDefinition()->GetParticleSubType() == "di_quark"  &&
               end->GetDefinition()->GetParticleSubType() == "di_quark"  )  ||
           ( start->GetDefinition()->GetParticleSubType() == "quark"     &&
               end->GetDefinition()->GetParticleSubType() == "quark"     ) ) ) { 
    // Kinky strings are allowed only for qq-q strings;
    // Kinky strings are impossible for other systems (qq-qqbar, q-qbar)
    // according to the analysis of Pbar P interactions

    if ( W > Wmin ) {  // Kink is possible
      if ( hadron->GetStatus() == 0 ) {  // VU 10.04.2012
        G4double Pt2kink = theParameters->GetPt2Kink(); // For non-diffractive
//        Pt = std::sqrt( Pt2kink * ( std::pow( W2/16.0/Pt2kink + 1.0, G4UniformRand() ) - 1.0 ) ); // Uzhi 18 Sept. 2014
        if(Pt2kink)                                                                                 // Uzhi 18 Sept. 2014
        {Pt = std::sqrt( Pt2kink * ( std::pow( W2/16.0/Pt2kink + 1.0, G4UniformRand() ) - 1.0 ) );} // Uzhi 18 Sept. 2014
        else {Pt=0.;}                                                                               // Uzhi 18 Sept. 2014
      } else {
        Pt = 0.0;
      }

      if ( Pt > 500.0*MeV ) {
        G4double Ymax = std::log( W/2.0/Pt + std::sqrt( W2/4.0/Pt/Pt - 1.0 ) );
        G4double Y = Ymax*( 1.0 - 2.0*G4UniformRand() );
        x1 = 1.0 - Pt/W * std::exp( Y );
        x3 = 1.0 - Pt/W * std::exp(-Y );
        //x2 = 2.0 - x1 - x3;

        G4double Mass_startQ = 650.0*MeV;
        if ( PDGcode_startQ <  3 ) Mass_startQ =  325.0*MeV;
        if ( PDGcode_startQ == 3 ) Mass_startQ =  500.0*MeV;
        if ( PDGcode_startQ == 4 ) Mass_startQ = 1600.0*MeV;
        G4double Mass_endQ = 650.0*MeV;
        if ( PDGcode_endQ <  3 ) Mass_endQ =  325.0*MeV;
        if ( PDGcode_endQ == 3 ) Mass_endQ =  500.0*MeV;
        if ( PDGcode_endQ == 4 ) Mass_endQ = 1600.0*MeV;

        G4double P2_1 = W2*x1*x1/4.0 - Mass_endQ*Mass_endQ;
        G4double P2_3 = W2*x3*x3/4.0 - Mass_startQ*Mass_startQ;
        G4double P2_2 = sqr( (2.0 - x1 - x3)*W/2.0 );
        if ( P2_1 <= 0.0  ||  P2_3 <= 0.0 ) { 
          Kink = false;
        } else {
          G4double P_1 = std::sqrt( P2_1 );
          G4double P_2 = std::sqrt( P2_2 );
          G4double P_3 = std::sqrt( P2_3 );
          G4double CosT12 = ( P2_3 - P2_1 - P2_2 ) / (2.0*P_1*P_2);
          G4double CosT13 = ( P2_2 - P2_1 - P2_3 ) / (2.0*P_1*P_3);
          //Pt = P_2 * std::sqrt( 1.0 - CosT12*CosT12 );  // because system was rotated 11.12.09

          if ( std::abs( CosT12 ) > 1.0  ||  std::abs( CosT13 ) > 1.0 ) {
            Kink = false;
          } else { 
            Kink = true; 
            Pt = P_2 * std::sqrt( 1.0 - CosT12*CosT12 );  // because system was rotated 11.12.09
            Pstart.setPx( -Pt ); Pstart.setPy( 0.0 ); Pstart.setPz( P_3*CosT13 ); 
            Pend.setPx(   0.0 ); Pend.setPy(   0.0 ); Pend.setPz(          P_1 ); 
            Pkink.setPx(   Pt ); Pkink.setPy(  0.0 ); Pkink.setPz(  P_2*CosT12 );
            Pstart.setE( x3*W/2.0 );                
            Pkink.setE( Pkink.vect().mag() );
            Pend.setE( x1*W/2.0 );

            G4double XkQ = GetQuarkFractionOfKink( 0.0, 1.0 );
            if ( Pkink.getZ() > 0.0 ) {
              if ( XkQ > 0.5 ) {
                PkinkQ1 = XkQ*Pkink;
              } else {
                PkinkQ1 = (1.0 - XkQ)*Pkink;
              }
            } else {
              if ( XkQ > 0.5 ) {
                PkinkQ1 = (1.0 - XkQ)*Pkink;
              } else {
                PkinkQ1 = XkQ*Pkink;
              }
            }

            PkinkQ2 = Pkink - PkinkQ1;
            // Minimizing Pt1^2+Pt3^2
            G4double Cos2Psi = ( sqr(x1) - sqr(x3) + 2.0*sqr( x3*CosT13 ) ) /
                               std::sqrt( sqr( sqr(x1) - sqr(x3) ) + sqr( 2.0*x1*x3*CosT13 ) );
            G4double Psi = std::acos( Cos2Psi );

            G4LorentzRotation Rotate;
            if ( isProjectile ) {
              Rotate.rotateY( Psi );
            } else {
              Rotate.rotateY( pi + Psi );                  // Uzhi 18 Sept. 2014
//            Rotate.rotateY( pi - Psi );                  // Uzhi 18 Sept. 2014
            }                   
            Rotate.rotateZ( twopi * G4UniformRand() );
            Pstart *= Rotate;
            Pkink *= Rotate;
            PkinkQ1 *= Rotate;
            PkinkQ2 *= Rotate;
            Pend *= Rotate;
          }
        }  // End of if ( P2_1 <= 0.0  ||  P2_3 <= 0.0 )
      }  // End of if ( Pt > 500.0*MeV )
    } // End of if ( W > Wmin ) : check for a kink
  }  // end of qq-q string selection

  //G4cout << "Kink " << Kink << " " << start->GetDefinition()->GetParticleSubType() << " "
  //       << end->GetDefinition()->GetParticleSubType() << G4endl;
  //G4cout << "Kink " << Kink << " " << start->GetDefinition()->GetPDGEncoding() << " "
  //       << end->GetDefinition()->GetPDGEncoding() << G4endl;
  //G4int Uzhi; G4cin >> Uzhi;

  if ( Kink ) {  // Kink is possible

    //G4cout << "Kink is sampled!" << G4endl;
    std::vector< G4double > QuarkProbabilitiesAtGluonSplitUp = 
        theParameters->GetQuarkProbabilitiesAtGluonSplitUp();

    G4int QuarkInGluon( 1 ); G4double Ksi = G4UniformRand();
    for ( unsigned int Iq = 0; Iq < 3; Iq++ ) {
      //G4cout << "Iq " << Iq << G4endl;
      if ( Ksi > QuarkProbabilitiesAtGluonSplitUp[Iq] ) QuarkInGluon++;
    }
    //G4cout << "Last Iq " << QuarkInGluon << G4endl;
    G4Parton* Gquark = new G4Parton( QuarkInGluon );
    G4Parton* Ganti_quark = new G4Parton( -QuarkInGluon );
    //G4cout << "Lorentz " << G4endl;

    G4LorentzRotation toCMS( -1 * Phadron.boostVector() );
    G4LorentzRotation toLab( toCMS.inverse() );
    //G4cout << "Pstart " << Pstart << G4endl;
//G4cout << "Pend   " << Pend << G4endl;
//G4cout << "Kink1  " <<PkinkQ1<<G4endl;
//G4cout << "Kink2  " <<PkinkQ2<<G4endl;
//G4cout << "Pstart " << Pstart << G4endl<<G4endl;

    Pstart.transform( toLab );  start->Set4Momentum( Pstart );
    PkinkQ1.transform( toLab );
    PkinkQ2.transform( toLab );
    Pend.transform( toLab );    end->Set4Momentum( Pend );
    //G4cout << "Pstart " << Pstart << G4endl;
    //G4cout << "Pend   " << Pend << G4endl;
    //G4cout << "Defin " << hadron->GetDefinition()<< G4endl;
    //G4cout << "Defin " << hadron->GetDefinition()->GetPDGEncoding()<< G4endl;

    //G4int absPDGcode = std::abs( hadron->GetDefinition()->GetPDGEncoding() );
    G4int absPDGcode = 1500;  // 23 Dec
    if ( start->GetDefinition()->GetParticleSubType() == "quark"  &&
         end->GetDefinition()->GetParticleSubType() == "quark" ) {
      absPDGcode = 110;
    }
    //G4cout << "absPDGcode " << absPDGcode << G4endl;

    if ( absPDGcode < 1000 ) {  // meson
      if ( isProjectile ) { // Projectile
        if ( end->GetDefinition()->GetPDGEncoding() > 0 ) {  // A quark on the end
          
FirstString  = new G4ExcitedString( end   , Ganti_quark, +1 );
          SecondString = new G4ExcitedString( Gquark, start      , +1 );
          Ganti_quark->Set4Momentum( PkinkQ1 );
          Gquark->Set4Momentum( PkinkQ2 );
        } else {  // Anti_Quark on the end
          FirstString  = new G4ExcitedString( end        , Gquark, +1 );
          SecondString = new G4ExcitedString( Ganti_quark, start , +1 );
          Gquark->Set4Momentum( PkinkQ1 );
          Ganti_quark->Set4Momentum( PkinkQ2 );
        }
      } else {  // Target
        if ( end->GetDefinition()->GetPDGEncoding() > 0 ) { // A quark on the end
          FirstString  = new G4ExcitedString( Ganti_quark, end   , -1 );
          SecondString = new G4ExcitedString( start      , Gquark, -1 );
          Ganti_quark->Set4Momentum( PkinkQ2 );
          Gquark->Set4Momentum( PkinkQ1 );
        } else {  // Anti_Quark on the end
          FirstString  = new G4ExcitedString( Gquark, end        , -1 );
          SecondString = new G4ExcitedString( start , Ganti_quark, -1 );
          Gquark->Set4Momentum( PkinkQ2 );
          Ganti_quark->Set4Momentum( PkinkQ1 );
        }
      }
    } else {  // Baryon/AntiBaryon

//G4cout<<"isProjectile "<<isProjectile<<G4endl;
//G4cout<<"end "<<end->GetDefinition()->GetPDGEncoding()<<" "<<end->Get4Momentum()<<G4endl;
//G4cout<<PkinkQ1<<G4endl;
//G4cout<<PkinkQ2<<G4endl;
//G4cout<<"start "<<start->GetDefinition()->GetPDGEncoding()<<" "<<start->Get4Momentum()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
      if ( isProjectile ) {  // Projectile
        if ( end->GetDefinition()->GetParticleType() == "diquarks"  &&
             end->GetDefinition()->GetPDGEncoding() > 0 ) {  // DiQuark on the end
          FirstString  = new G4ExcitedString( end        , Gquark, +1 );  // Open Uzhi
          SecondString = new G4ExcitedString( Ganti_quark, start , +1 );  // Open Uzhi
          Gquark->Set4Momentum( PkinkQ1 );
          Ganti_quark->Set4Momentum( PkinkQ2 );

//          FirstString  = new G4ExcitedString( Gquark,        end, +1 );   // Uzhi 18 Sept. 2014
//          SecondString = new G4ExcitedString(  start, Ganti_quark, +1 );  // Uzhi 18 Sept. 2014

        } else {                            // Anti_DiQuark on the end or quark
          FirstString  = new G4ExcitedString( end   , Ganti_quark, +1 );
          SecondString = new G4ExcitedString( Gquark, start      , +1 );
          Ganti_quark->Set4Momentum( PkinkQ1 );
          Gquark->Set4Momentum( PkinkQ2 );
        }
      } else {  // Target
        if ( end->GetDefinition()->GetParticleType() == "diquarks"  &&
             end->GetDefinition()->GetPDGEncoding() > 0 ) {  // DiQuark on the end
//          FirstString  = new G4ExcitedString( Gquark, end        , -1 );
//          SecondString = new G4ExcitedString( start , Ganti_quark, -1 );
          Gquark->Set4Momentum( PkinkQ1 );
          Ganti_quark->Set4Momentum( PkinkQ2 );

          FirstString  = new G4ExcitedString(         end, Gquark, -1 );
          SecondString = new G4ExcitedString( Ganti_quark,  start, -1 );

        } else {  // Anti_DiQuark on the end or Q
          FirstString  = new G4ExcitedString( Ganti_quark, end   , -1 );  // Uzhi ?????
          SecondString = new G4ExcitedString( start      , Gquark, -1 );
          Gquark->Set4Momentum( PkinkQ2 );
          Ganti_quark->Set4Momentum( PkinkQ1 );
        }
      }
    }

    FirstString->SetTimeOfCreation( hadron->GetTimeOfCreation() );
    FirstString->SetPosition( hadron->GetPosition() );
    SecondString->SetTimeOfCreation( hadron->GetTimeOfCreation() );
    SecondString->SetPosition( hadron->GetPosition() );
  
  } else { // End of kink is possible: Kink is impossible

    //G4cout << start << " " << start->GetPDGcode() << " " << end << " " << end->GetPDGcode() 
    //       << G4endl;
/*                                                         // Uzhi 18 Sept. 2014
    if ( isProjectile ) {
      FirstString = new G4ExcitedString( end, start, +1 );
    } else {
      FirstString = new G4ExcitedString( start, end, -1 );
    }
*/                                                         // Uzhi 18 Sept. 2014

    FirstString = new G4ExcitedString( end, start, +1 );   // Uzhi 18 Sept. 2014

    FirstString->SetTimeOfCreation( hadron->GetTimeOfCreation() );
    FirstString->SetPosition( hadron->GetPosition() );
    SecondString = 0;

    // momenta of string ends
    G4double Momentum = hadron->Get4Momentum().vect().mag();
    G4double Plus  = hadron->Get4Momentum().e() + Momentum;
    G4double Minus = hadron->Get4Momentum().e() - Momentum;
    G4ThreeVector tmp;
    if ( Momentum > 0.0 ) {
      tmp.set( hadron->Get4Momentum().px(), 
               hadron->Get4Momentum().py(), 
               hadron->Get4Momentum().pz() );
      tmp /= Momentum;
    } else {
      tmp.set( 0.0, 0.0, 1.0 );
    }
    G4LorentzVector Pstart1( tmp, 0.0 );
    G4LorentzVector   Pend1( tmp, 0.0 );
    if ( isProjectile ) {
      Pstart1 *= (-1.0)*Minus/2.0;
      Pend1   *= (+1.0)*Plus /2.0;
    } else {
      Pstart1 *= (+1.0)*Plus/ 2.0;
      Pend1   *= (-1.0)*Minus/2.0;
    }
    Momentum = -Pstart1.mag();
    Pstart1.setT( Momentum );  // It is assumed that quark has m=0.
    Momentum = -Pend1.mag();
    Pend1.setT( Momentum );    // It is assumed that di-quark has m=0.
    start->Set4Momentum( Pstart1 );
    end->Set4Momentum( Pend1 );
    SecondString = 0;
      
  } // End of kink is impossible 

  //G4cout << "Quarks in the string at creation" << FirstString->GetRightParton()->GetPDGcode()
  //       << " " << FirstString->GetLeftParton()->GetPDGcode() << G4endl
  //       << FirstString << " " << SecondString << G4endl;

  #ifdef G4_FTFDEBUG
  G4cout << " generated string flavors          " << start->GetPDGcode() << " / " 
         << end->GetPDGcode() << G4endl << " generated string momenta:   quark " 
         << start->Get4Momentum() << "mass : " << start->Get4Momentum().mag() << G4endl
         << " generated string momenta: Diquark " << end->Get4Momentum() << "mass : " 
         << end->Get4Momentum().mag() << G4endl << " sum of ends                       "
         << Pstart + Pend << G4endl << " Original                          " 
         << hadron->Get4Momentum() << G4endl;
  #endif

  return;
}


//============================================================================

G4double G4DiffractiveExcitation::ChooseP( G4double Pmin, G4double Pmax ) const {
  // Choose an x between Xmin and Xmax with P(x) ~ 1/x . 
  // To be improved...
  G4double range = Pmax - Pmin;                    
  if ( Pmin <= 0.0 || range <= 0.0 ) {
    G4cout << " Pmin, range : " << Pmin << " , " << range << G4endl;
    throw G4HadronicException( __FILE__, __LINE__,
                               "G4DiffractiveExcitation::ChooseP : Invalid arguments " );
  }
  G4double P = Pmin * std::pow( Pmax/Pmin, G4UniformRand() ); 
  //G4double P = (Pmax - Pmin) * G4UniformRand() + Pmin;
  return P;
}


//============================================================================

G4ThreeVector G4DiffractiveExcitation::GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const {
  //  @@ this method is used in FTFModel as well. Should go somewhere common!
  G4double Pt2( 0.0 );
  if ( AveragePt2 <= 0.0 ) {
    Pt2 = 0.0;
  } else {
    Pt2 = -AveragePt2 * std::log( 1.0 + G4UniformRand() * 
                                        ( std::exp( -maxPtSquare/AveragePt2 ) - 1.0 ) );
  }
  G4double Pt = std::sqrt( Pt2 );
  G4double phi = G4UniformRand() * twopi;
  return G4ThreeVector( Pt * std::cos( phi ), Pt * std::sin( phi ), 0.0 );
}


//============================================================================

G4double G4DiffractiveExcitation::GetQuarkFractionOfKink( G4double zmin, G4double zmax ) const {
  G4double z, yf;
  do {
    z = zmin + G4UniformRand() * (zmax - zmin);
    yf = z*z + sqr(1.0 - z);
  } while ( G4UniformRand() > yf ); 
  return z;
}


//============================================================================

void G4DiffractiveExcitation::UnpackMeson( const G4int IdPDG, G4int& Q1, G4int& Q2 ) const {
  G4int absIdPDG = std::abs( IdPDG );
  if(!((absIdPDG == 111)||(absIdPDG == 221)||(absIdPDG == 331)))         // Uzhi Oct. 2014
  {                          // Ordinary mesons =======================
   Q1 =  absIdPDG / 100;
   Q2 = (absIdPDG % 100) / 10;
   G4int anti = 1 - 2 * ( std::max( Q1, Q2 ) % 2 );
   if ( IdPDG < 0 ) anti *= -1;
   Q1 *= anti;
   Q2 *= -1 * anti;
  }
  else 
  {                         // Pi0, Eta, Eta' =======================
   if( G4UniformRand() < 0.5 ) {Q1 = 1; Q2 = -1;}
   else                        {Q1 = 2; Q2 = -2;}
  }
  return;
}


//============================================================================

void G4DiffractiveExcitation::UnpackBaryon( G4int IdPDG, 
                                            G4int& Q1, G4int& Q2, G4int& Q3 ) const {
  Q1 = IdPDG          / 1000;
  Q2 = (IdPDG % 1000) / 100;
  Q3 = (IdPDG % 100)  / 10;
  return;
}


//============================================================================

G4int G4DiffractiveExcitation::NewNucleonId( G4int Q1, G4int Q2, G4int Q3 ) const {
  G4int TmpQ( 0 );
  if ( Q3 > Q2 ) {
    TmpQ = Q2;
    Q2 = Q3;
    Q3 = TmpQ;
  } else if ( Q3 > Q1 ) {
    TmpQ = Q1;
    Q1 = Q3;
    Q3 = TmpQ;
  }
  if ( Q2 > Q1 ) {
    TmpQ = Q1;
    Q1 = Q2;
    Q2 = TmpQ;
  }
  G4int NewCode = Q1*1000 + Q2*100 + Q3*10 + 2; 
  return NewCode;
}


//============================================================================

G4DiffractiveExcitation::G4DiffractiveExcitation( const G4DiffractiveExcitation& ) {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation copy contructor not meant to be called" );
}


//============================================================================

const G4DiffractiveExcitation & G4DiffractiveExcitation::operator=( const G4DiffractiveExcitation& ) {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation = operator not meant to be called" );
  return *this;
}


//============================================================================

int G4DiffractiveExcitation::operator==( const G4DiffractiveExcitation& ) const {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation == operator not meant to be called" );
}


//============================================================================

int G4DiffractiveExcitation::operator!= ( const G4DiffractiveExcitation& ) const {
  throw G4HadronicException( __FILE__, __LINE__, 
                             "G4DiffractiveExcitation != operator not meant to be called" );
}
