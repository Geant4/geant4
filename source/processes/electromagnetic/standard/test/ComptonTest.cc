// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ComptonTest.cc,v 1.1 1999-01-08 16:32:24 gunter Exp $
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
#include "G4ComptonScattering.hh"

#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

#include "G4GRSVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"

//    it tests the G4ComptonScattering process ----------------
//    ---------- M.Maire on 19/04/96 --------------------------
//
//  Modifs:
//  19-02-97, Modifs in DoIt for new physic scheme
//  25-02-97, 'simplified' version for the process manager
//  12-03-97, GetMeanFreePath and PostStepDoIt
//  21-03-97, calling sequence of AddProcess modified
//  22-09-97, geometry definition modified for the touchable

int main()
{
  //-------- set output format-------
   G4cout.setf( ios::scientific, ios::floatfield );
  //-------- write results onto a file --------
   ofstream outFile( "Compton.out", ios::out);
   outFile.setf( ios::scientific, ios::floatfield );

  //
  //--------- Materials definition ---------
  //
  G4Material* Al = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3 );
  G4Material* Fe = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3 );
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

  //--------- Processes definition ---------
  //
  G4ProcessManager* theGammaProcessManager = Gamma->GetProcessManager();

  G4ComptonScattering theComptonScatteringProcess;
  theGammaProcessManager->AddProcess(&theComptonScatteringProcess,-1,-1,0);

  G4ForceCondition* condition;

  // ------- set cut and Build CrossSection Tables -------
  //
  Gamma->SetCuts(1.*mm);
  Electron->SetCuts(1.*mm);  

  // -------- create 1 Dynamic Particle  ----

  G4double gammaenergy = 10.*MeV;
  G4ParticleMomentum gammadirection(0.,0.,1.);
  G4DynamicParticle aGamma(G4Gamma::Gamma(),gammadirection,gammaenergy);

  //--------- track definition (for this test ONLY!)------------
  //
    G4ThreeVector aPosition(0.,0.,0.);
    G4double aTime = 0. ;

    G4Track* ptrack = new G4Track(&aGamma,aTime,aPosition) ;
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
  Tkin[ 0] =   5.*keV; Tkin[ 1] =  10.*keV; Tkin[ 2] =  50.*keV; Tkin[ 3] = 100.*keV; 
  Tkin[ 4] = 200.*keV; Tkin[ 5] = 500.*keV; Tkin[ 6] =   1.*MeV; Tkin[ 7] =   2.*MeV;
  Tkin[ 8] =   5.*MeV; Tkin[ 9] =  10.*MeV; Tkin[10] =  20.*MeV; Tkin[11] =  50.*MeV;
  Tkin[12] = 100.*MeV; Tkin[13] = 200.*MeV; Tkin[14] = 500.*MeV; Tkin[15] = 200.*GeV;

  for ( G4int J = 0 ; J < theMaterialTable->length() ; J++ )
  {
    apttoMaterial = (*theMaterialTable)[ J ] ;
    MaterialName  = apttoMaterial->GetName() ;

    LogicalFrame->SetMaterial(apttoMaterial); 

    outFile << " " << endl;   // Print table of mean free path
    outFile <<"  " << MaterialName  << "    Compton Scattering Mean Free Path in cm" << endl;
    outFile << "----------------------------------------------------" << endl;
    outFile << " " << endl;
    outFile << "kinetic energy (MeV)     mean free path (cm)" << endl ;
    outFile << " " << endl;

    for ( G4int i=0 ; i<nkin ; i++)
      {
       aGamma.SetKineticEnergy(Tkin[i]);

       meanFreePath = theComptonScatteringProcess.GetMeanFreePath(aTrack, 1., condition);

       outFile <<"  " <<  Tkin[i]/MeV << "            " << meanFreePath/cm << endl ;
      }
    outFile << " " << endl;
   }

  //
  // --------- Test the DoIt for the Compton Scattering
  //

    apttoMaterial = (*theMaterialTable)(imat) ;
    LogicalFrame->SetMaterial(apttoMaterial); 
    aGamma.SetKineticEnergy(10.*MeV);
    aGamma.SetMomentumDirection(0., 0., 1.);
    G4VParticleChange* adummy;
    G4ParticleChange* aFinalGamma;
    G4Track* aFinalElectr;
 
    G4int iteration = 0;   
    do {
        adummy      = theComptonScatteringProcess.PostStepDoIt(aTrack, aStep);
        aFinalGamma = (G4ParticleChange*)adummy;

         outFile << " -----------------------------------------------------------" << endl;   
        // check the kinematic
        //
        G4double Tkin0 = aGamma.GetKineticEnergy();
        G4double Px0   = aGamma.GetMomentum().x() ,
                 Py0   = aGamma.GetMomentum().y() ,
                 Pz0   = aGamma.GetMomentum().z() ;
        outFile << " Initial Gamma :  Px= " << Px0/MeV << "   Py= " << Py0/MeV 
                << "   Pz= " << Pz0/MeV << "   Tkin= " << Tkin0/MeV << endl ;

        G4double Tkin1 =  aFinalGamma->GetEnergyChange();
        G4double Px1   = (aFinalGamma->GetMomentumChange()->x())*Tkin1 ,
                 Py1   = (aFinalGamma->GetMomentumChange()->y())*Tkin1 ,
                 Pz1   = (aFinalGamma->GetMomentumChange()->z())*Tkin1 ;
        G4double Edep  =  aFinalGamma->GetLocalEnergyDeposit();
        outFile << "   final Gamma :  Px= " << Px1/MeV << "   Py= " << Py1/MeV 
                << "   Pz= " << Pz1/MeV << "   Tkin= " << Tkin1/MeV 
                << "   Edep= " << Edep/MeV << endl ;

        G4double Tkin2 = 0., Px2 = 0., Py2 = 0., Pz2 = 0. ;
        if (aFinalGamma->GetNumberOfSecondaries()) {
            aFinalElectr = aFinalGamma->GetSecondary(0) ;
            Tkin2 =  aFinalElectr->GetKineticEnergy() ;
            Px2   = (aFinalElectr->GetMomentum()).x() ;
            Py2   = (aFinalElectr->GetMomentum()).y() ;
            Pz2   = (aFinalElectr->GetMomentum()).z() ;
            outFile << "   final Electr :  Px= " << Px2/MeV << "   Py= " << Py2/MeV 
                    << "   Pz= " << Pz2/MeV << "   Tkin= " << Tkin2/MeV << endl ;
           // NOTE - Secondaries are normally deleted by the track destructor !
           delete aFinalGamma->GetSecondary(0);
           }

        G4double sTkin = Tkin0 - Tkin1 - Tkin2 - Edep,    
                 sPx = Px0 - Px1 - Px2, sPy = Py0 - Py1 - Py2, sPz = Pz0 - Pz1 - Pz2;
        outFile << endl;
        outFile << " Balance        :  Px= " << sPx/MeV << "   Py= " << sPy/MeV << "   Pz= " << sPz/MeV 
                           << "   Tkin= " << sTkin/MeV << endl ;

        // update the initial gamma
        //
        aGamma.SetKineticEnergy(Tkin1);
        aGamma.SetMomentumDirection(Px1/Tkin1, Py1/Tkin1, Pz1/Tkin1);
        iteration++; 

   }  while ((!aFinalGamma->GetStatusChange())&&(iteration < 10)) ;

  return EXIT_SUCCESS;
}
