// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PhotoElectricTest.cc,v 1.1 1999-01-08 16:32:26 gunter Exp $
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
#include "G4PhotoElectricEffect.hh"

#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

#include "G4GRSVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"

//    it tests the G4PhotoElectric effect ----------------
//    ---------- M.Maire on 29/04/96 ---------------------
//
//  Modifs:
//  04-07-96, revised by H.Kurashige
//  13-03-97, GetMeanFreePath and PostStepDoIt, by M.Maire
//  22-09-97, geometry definition modified for the touchable
//  20-11-97, test out of limits

int main()
{
  //-------- set output format-------
   G4cout.setf( ios::scientific, ios::floatfield );
  //-------- write results onto a file --------
   ofstream outFile( "PhotoElectric.out", ios::out);
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
  G4ParticleDefinition* Gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* Electron = G4Electron::ElectronDefinition();

  //--------- Process definition ---------
  //
  G4ProcessManager* theGammaProcessManager = Gamma->GetProcessManager();
  G4PhotoElectricEffect thePhotoElectricEffect;
  theGammaProcessManager->AddProcess(&thePhotoElectricEffect,-1,-1,0);

  G4ForceCondition* condition;

  // ------- set cut and Build CrossSection Tables -------
  //
  Gamma->SetCuts(1.*mm);
  Electron->SetCuts(1.*mm);  

  // -------- create 1 Dynamic Particle  ----
  //
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
  G4int nkin = 16; 
  G4double Tkin[16];
  Tkin[ 0] =   5.*keV; Tkin[ 1] =  10.*keV; Tkin[ 2] =  20.*keV; Tkin[ 3] =  30.*keV; 
  Tkin[ 4] =  40.*keV; Tkin[ 5] =  50.*keV; Tkin[ 6] = 100.*keV; Tkin[ 7] = 200.*keV;
  Tkin[ 8] = 300.*keV; Tkin[ 9] = 500.*keV; Tkin[10] = 700.*keV; Tkin[11] =   1.*MeV;
  Tkin[12] =   2.*MeV; Tkin[13] =   5.*MeV; Tkin[14] =  10.*MeV; Tkin[15] =  60.*MeV;

  for ( G4int J = 0 ; J < theMaterialTable->length() ; J++ )
  {
    apttoMaterial = (*theMaterialTable)[ J ] ;
    MaterialName  = apttoMaterial->GetName() ; 
    LogicalFrame->SetMaterial(apttoMaterial); 

    outFile << " " << endl;   // Print table of mean free path
    outFile <<"  " << MaterialName  << "    Photo Electric Mean Free Path in cm" << endl;
    outFile << "----------------------------------------------------" << endl;
    outFile << " " << endl;
    outFile << "kinetic energy (MeV)     mean free path (cm)" << endl ;
    outFile << " " << endl;

    for ( G4int i=0 ; i<nkin ; i++)
      {
       aPhoton.SetKineticEnergy(Tkin[i]);

       meanFreePath = thePhotoElectricEffect.GetMeanFreePath(aTrack, 1., condition);

       outFile <<"  " <<  Tkin[i]/MeV << "            " << meanFreePath/cm << endl ;
      }
    outFile << " " << endl;
   }

  //
  // --------- Test the PostStepDoIt for the Photo Electric effect
  //

    apttoMaterial = (*theMaterialTable)(imat);
    LogicalFrame->SetMaterial(apttoMaterial); 

    G4double Tkin2 = 1.*MeV;
    aPhoton.SetKineticEnergy(Tkin2);
    aPhoton.SetMomentumDirection(0., 0., 1.);
    G4VParticleChange* adummy;
    G4ParticleChange* aFinalPhoton;
    G4Track* aFinalElectr;
 
    G4int iteration = 0;   
    while ((Tkin2 > 0.)&&(iteration < 10))
       {
        meanFreePath = thePhotoElectricEffect.GetMeanFreePath(aTrack, 1., condition);
        adummy       = thePhotoElectricEffect.PostStepDoIt(aTrack, aStep);
        aFinalPhoton = (G4ParticleChange*)adummy;

         outFile << " -----------------------------------------------------------" << endl;   
        // check the kinematic
        //
        G4double Tkin0 = aPhoton.GetKineticEnergy();
        G4double Px0   = aPhoton.GetMomentum().x() ,
                 Py0   = aPhoton.GetMomentum().y() ,
                 Pz0   = aPhoton.GetMomentum().z() ;
        outFile << " Initial Photon : Tkin= " << Tkin0/MeV
                << "  Px= " << Px0/MeV << "   Py= " << Py0/MeV << "   Pz= " << Pz0/MeV << endl; 

        G4double  Px2 = 0., Py2 = 0., Pz2 = 0. ; Tkin2 = 0. ;
        if (aFinalPhoton->GetNumberOfSecondaries()) {
            aFinalElectr = aFinalPhoton->GetSecondary(0) ;
            Tkin2 =  aFinalElectr->GetKineticEnergy() ;
            Px2   = (aFinalElectr->GetMomentum()).x() ;
            Py2   = (aFinalElectr->GetMomentum()).y() ;
            Pz2   = (aFinalElectr->GetMomentum()).z() ;
            outFile << "   final Electr : Tkin= " << Tkin2/MeV  
                    << "  Px= " << Px2/MeV << "   Py= " << Py2/MeV << "   Pz= " << Pz2/MeV << endl;
           // NOTE - Secondaries are normally deleted by the track destructor !
           delete aFinalPhoton->GetSecondary(0);
           }

        G4double Edep  =  aFinalPhoton->GetLocalEnergyDeposit();
        outFile << " Energy deposit = " << Edep/MeV << endl;

        outFile << endl;
        outFile << " Balance  = " << (Tkin0-Tkin2-Edep)/MeV << endl ;

        // 'Build' a new initial photon
        //
        aPhoton.SetKineticEnergy(Tkin2);
        aPhoton.SetMomentumDirection(Px2/Tkin2, Py2/Tkin2, Pz2/Tkin2);
        iteration++; 

   } ;

  return EXIT_SUCCESS;
}
