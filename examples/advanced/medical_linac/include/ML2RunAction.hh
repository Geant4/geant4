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
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanit�, Italy
// *Istituto Superiore di Sanit� and INFN Roma, gruppo collegato Sanit�, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#ifndef CML2RunActionH
#define CML2RunActionH

#include "G4UserRunAction.hh"
#include "G4Timer.hh"
#include "ML2WorldConstruction.hh"
#include "ML2Convergence.hh"
#include "ML2SinputData.hh"
#include "ML2CInputData.hh"

class CML2RunAction : public G4UserRunAction
{
public:
	CML2RunAction(CML2Convergence *convergence, G4int nBeam, G4bool bOnlyVisio, G4int voxelX, G4int voxelY, G4int voxelz);
	~CML2RunAction(void); 
        virtual G4Run* GenerateRun();
	void BeginOfRunAction(const G4Run *aRun);
	void EndOfRunAction(const G4Run *aRun);
	void setActualLoop(G4int nL){ nLoop = nL; } 

       void ChangeOutputFileName(G4String name){filename =name;}
  
  // Utility method for converting segment number of
  // water phantom to copyNo of HitsMap.
       G4int CopyNo(G4int ix, G4int iy, G4int iz)
       {return (iy*(fNx*fNz)+ix*fNz+iz); }

private:

	G4bool bRotationTranslationFileNames;
	CML2Convergence *convergence; 
	G4Timer MyTime;
	G4double loopElapsedTime;
	G4int nBeam, nLoop;
	G4bool bOnlyVisio; 
        G4String filename;

  std::vector<G4String> fSDName; 

  // for conversion of sengment number to copyNo.
  G4int fNx, fNy, fNz; 
};

#endif
