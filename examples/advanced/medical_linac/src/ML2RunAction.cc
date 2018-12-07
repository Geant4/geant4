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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "G4ios.hh"
#include "ML2RunAction.hh"
#include "ML2Run.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"   

CML2RunAction::CML2RunAction(CML2Convergence *conv, G4int nB, G4bool bOV, G4int voxelX, G4int voxelY, G4int voxelZ): fNx(voxelX), fNy(voxelY), fNz(voxelZ)
{
    bRotationTranslationFileNames = true;
    convergence = conv;
    nBeam = nB;
    bOnlyVisio = bOV;
    nLoop = 0;


  fSDName.push_back(G4String("PhantomSD"));
}

CML2RunAction::~CML2RunAction(void)
{
 fSDName.clear();
}

G4Run* CML2RunAction::GenerateRun()
{
  // SUSANNA
  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  return new ML2Run(fSDName);
}

void CML2RunAction::BeginOfRunAction(const G4Run * aRun)
{
 G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

    G4String fullName;
    if (bRotationTranslationFileNames)
    {
    	fullName = CML2AcceleratorConstruction::GetInstance()->getCurrentRotationString()+
    			CML2PhantomConstruction::GetInstance()->getCurrentTranslationString();
    }
    else
    {
    	fullName = "";
    }
//    CML2PhantomConstruction::GetInstance()->setNewName(fullName);

    CML2AcceleratorConstruction::GetInstance()->writeInfo();
    CML2PhantomConstruction::GetInstance()->writeInfo();


    G4cout << "*********************************************" << G4endl;
    if (convergence -> getNMaxLoops() < 0 || bOnlyVisio)
    {
        G4cout << "loop n. " << ++nLoop << G4endl;
        G4cout << "Launched " << nBeam << " random primary particles" << G4endl;
    }
    else
    {
    	G4cout << "loop n. " << ++nLoop << "/" << convergence->getNMaxLoops() << G4endl;
    	G4cout << "Launched " << nBeam << " random primary particles" << G4endl;
    }
    if (!bOnlyVisio)
    {
    	G4cout <<"Launched " << nBeam << " random primary particles" << G4endl;
    }
    G4cout<<"*********************************************"<<'\n';
    MyTime.Start();
}
void CML2RunAction::EndOfRunAction(const G4Run * aRun)
{

   if(!IsMaster()) return;

  ML2Run* ml2Run = (ML2Run*)aRun;
  //--- Dump all socred quantities involved in RE02Run.
  ml2Run->DumpAllScorer();
  //---

  //---------------------------------------------
  // Dump accumulated quantities for this RUN.
  //  (Display only central region of x-y plane)
  //---------------------------------------------
  G4THitsMap<G4double>* totDose  = ml2Run->GetHitsMap("PhantomSD/TotalDose");

  G4int ix;  
  G4int iy;
  G4int iz;

  std::ofstream  file("totDose.txt");
  for ( iz = 0; iz < fNz; iz++){   
    for ( iy = 0; iy < fNy; iy++){ 
      for (ix = 0; ix < fNx; ix++){ 
        G4double* TotD = (*totDose)[CopyNo(ix,iy,iz)];
        if ( !TotD ) TotD = new G4double(0.0);
        if (TotD!=0) file << ix << " "<<iy<<" "<<iz<<" "<< *TotD/gray << G4endl;
      }
    }
  }
  file.close();

  //  CML2WorldConstruction::GetInstance()->savePhantomData();
  //  CML2WorldConstruction::GetInstance()->savePhaseSpaceData();
    convergence->saveResults();

    MyTime.Stop();
    loopElapsedTime = MyTime.GetUserElapsed();

    G4cout << "loop elapsed time [s] : " << loopElapsedTime << '\n' << G4endl;
}
