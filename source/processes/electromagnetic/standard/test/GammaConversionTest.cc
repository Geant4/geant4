// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaConversionTest.cc,v 1.1 1999-01-08 16:32:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "G4ios.hh"
#include <fstream.h>
#include <iomanip.h>

#include "g4templates.hh"

#include "G4Material.hh"

#include "G4ProcessManager.hh"
#include "G4GammaConversion.hh"

#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4GRSVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"


//    it tests the G4GammaConversion process ----------------
//    ---------- M.Maire on 22/05/96 ---------------------
//
// Updates : 24-05-96, 14-01-97, 14-03-97, 22-09-97

int main()
{
  //-------- set output format-------
   G4cout.setf( ios::scientific, ios::floatfield );
  //-------- write results onto a file --------
   ofstream outFile( "GammaConversion.out", ios::out);
   outFile.setf( ios::scientific, ios::floatfield );

  //
  //--------- Materials definition ---------
  //
  G4Material* Al = new G4Material("Aluminium", 13.,  26.98*g/mole,  2.7 *g/cm3 );
  G4Material* Fe = new G4Material("Iron",      26.,  55.85*g/mole,  7.87*g/cm3 );
  G4Material* Pb = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3 );

  G4Element*   H = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Material* wa = new G4Material ("Water" , 1.*g/cm3, 2);
  wa->AddElement(H,2);
  wa->AddElement(O,1);

  G4Element::DumpInfo(); 
  G4Material::DumpInfo();

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // fix a material
  G4int imat = 0;

  // Geometry definitions
  //
    G4Box* theFrame = new G4Box ("Frame",10*m, 10*m, 10*m);

    G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
                                                  (*theMaterialTable)(imat),
						  "LFrame", 0, 0, 0);

    G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
                                          "PFrame",LogicalFrame,0,false,0);

  //--------- Particle definition ---------
  //
  G4ParticleDefinition* Gamma    = G4Gamma::GammaDefinition();
  G4ParticleDefinition* Electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* Positron = G4Positron::PositronDefinition();

  //--------- Processes definition ---------
  //
  G4ProcessManager* theGammaProcessManager = Gamma->GetProcessManager();

  G4GammaConversion theGammaConversion;
  theGammaProcessManager->AddProcess(&theGammaConversion,-1,-1,0);

  G4ForceCondition* condition;

  // ------- set cut and Build CrossSection Tables -------
  //

  Gamma->SetCuts(1.*mm);
  Electron->SetCuts(1.*mm);  
  Positron->SetCuts(1.*mm);

  // -------- create 1 Dynamic Particle ----

  G4double photonenergy = 10.*MeV;
  G4ParticleMomentum photondirection(0.,0.,1.);
  G4DynamicParticle aPhoton(G4Gamma::Gamma(),photondirection,photonenergy);

  //--------- track definition (for this test ONLY!)------------
  //
    G4ThreeVector aPosition(0.,0.,0.);
    G4double aTime = 0. ;

    G4Track* ptrack = new G4Track(&aPhoton,aTime,aPosition) ;
    G4Track& aTrack = (*ptrack) ;
    //ptrack->SetVolume(PhysicalFrame);
    G4GRSVolume* touche = new G4GRSVolume(PhysicalFrame, NULL, aPosition);   
    ptrack->SetTouchable(touche);

  // -------- create 1 Step (for this test only)----  

    G4Step* Step = new G4Step();  G4Step& aStep = (*Step);
    Step->SetTrack(ptrack);
    ptrack->SetStep(Step);
      
  // ---------- Print the tables

  G4Material* apttoMaterial ;
  G4String MaterialName ;
  G4double meanFreePath ;
  G4int nkin = 15; 
  G4double Tkin[15];
  Tkin[ 0] =  0.5*MeV; Tkin[ 1] =  1.*MeV; Tkin[ 2] =  5.*MeV; Tkin[ 3] =  10.*MeV; 
  Tkin[ 4] =  20.*MeV; Tkin[ 5] = 30.*MeV; Tkin[ 6] = 50.*MeV; Tkin[ 7] = 100.*MeV;
  Tkin[ 8] = 500.*MeV; Tkin[ 9] =  1.*GeV; Tkin[10] = 10.*GeV; Tkin[11] =  90.*GeV;
  Tkin[12] = 500.*GeV; Tkin[13] =  1.*TeV; Tkin[14] = 10.*TeV;

  for ( G4int J = 0 ; J < theMaterialTable->length() ; J++ )
  {
    apttoMaterial = (*theMaterialTable)[ J ] ;
    MaterialName  = apttoMaterial->GetName() ; 

    LogicalFrame->SetMaterial(apttoMaterial); 

    outFile << " " << endl;   // Print table of mean free path
    outFile <<"  " << MaterialName  << "    Gamma Conversion Mean Free Path in cm" << endl;
    outFile << "----------------------------------------------------" << endl;
    outFile << " " << endl;
    outFile << "kinetic energy (MeV)     mean free path (cm)" << endl ;
    outFile << " " << endl;

    for ( G4int i=0 ; i<nkin ; i++)
      {
       aPhoton.SetKineticEnergy(Tkin[i]);

       meanFreePath = theGammaConversion.GetMeanFreePath(aTrack, 1., condition);

       outFile <<"  " <<  Tkin[i]/MeV << "            " << meanFreePath/cm << endl ;
      }
    outFile << " " << endl;
   }
  outFile << " -----------------------------------------------------------" << endl;
  //
  // --------- Test the DoIt for the Gamma conversion process
  //

    apttoMaterial = (*theMaterialTable)(imat) ;
    LogicalFrame->SetMaterial(apttoMaterial); 
    G4double EPhot = 1.*GeV;
    aPhoton.SetKineticEnergy(EPhot);
    aPhoton.SetMomentumDirection(0., 0., 1.);
    G4VParticleChange* adummy;
    G4ParticleChange* aFinalPhoton;
    G4Track* aCreatedParticle;
 
    G4int iteration = 0;   
    while ((EPhot > 2*electron_mass_c2)&&(iteration < 10))
       {
        meanFreePath = theGammaConversion.GetMeanFreePath(aTrack, 1., condition);
        adummy       = theGammaConversion.PostStepDoIt(aTrack, aStep);
        aFinalPhoton = (G4ParticleChange*)adummy;

        // check the kinematic
        //
        G4double Tkin0 = aPhoton.GetKineticEnergy();
        G4double Px0   = aPhoton.GetMomentum().x() ,
                 Py0   = aPhoton.GetMomentum().y() ,
                 Pz0   = aPhoton.GetMomentum().z() ;
        outFile << " Initial Photon : Tkin= " << Tkin0/MeV << "   Px= " << Px0/MeV 
                << "   Py= " << Py0/MeV << "   Pz= " << Pz0/MeV << endl << endl ;

        G4String  Name;
        G4double  Px2 , Py2 , Pz2 , Etot2 ;
        G4double sEtot = Tkin0 , sPx = Px0 , sPy = Py0 , sPz = Pz0 ;

        G4int NumberOfSecondary = aFinalPhoton->GetNumberOfSecondaries();

        for (G4int i=0; i < NumberOfSecondary; i++) {
           aCreatedParticle = aFinalPhoton->GetSecondary(i) ;
           Name  =  aCreatedParticle->GetDefinition()->GetParticleName();
           Etot2 =  aCreatedParticle->GetTotalEnergy() ;
           Px2   = (aCreatedParticle->GetMomentum()).x() ;
           Py2   = (aCreatedParticle->GetMomentum()).y() ;
           Pz2   = (aCreatedParticle->GetMomentum()).z() ;
           outFile << " No: " << i << " (" << Name << ")     : Etot= " << Etot2/MeV   
                   << "   Px= " << Px2/MeV << "   Py= " << Py2/MeV << "   Pz= " << Pz2/MeV 
                   << endl;

          // NOTE - Secondaries are normally deleted by the track destructor !
          delete aFinalPhoton->GetSecondary(i);

          // Balance
             sEtot -= Etot2;  sPx -= Px2;  sPy -= Py2;  sPz -= Pz2;
        }

        G4double Edep  =  aFinalPhoton->GetLocalEnergyDeposit();
        outFile << " Local Energ deposit= " << Edep/MeV << endl ;
        sEtot -= Edep;

        outFile << endl;
        outFile << " Balance        : Ener= " << sEtot/MeV
             << "   Px= " << sPx/MeV << "   Py= " << sPy/MeV << "   Pz= " << sPz/MeV << endl; 
        outFile << " -----------------------------------------------------------" << endl;

        // 'Build' a new initial photon
        //
        EPhot = Etot2;
        aPhoton.SetKineticEnergy(EPhot);
        aPhoton.SetMomentumDirection(Px2/EPhot, Py2/EPhot, Pz2/EPhot);
        iteration++; 

   } ;

  return EXIT_SUCCESS;
}
