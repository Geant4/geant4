//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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

class CML2RunAction : public G4UserRunAction
{
public:
	CML2RunAction(CML2Convergence *convergence, G4int nBeam, G4bool bOnlyVisio);
	~CML2RunAction(void);
	void BeginOfRunAction(const G4Run *aRun);
	void EndOfRunAction(const G4Run *aRun);
	void setActualLoop(G4int nLoop){this->nLoop=nLoop;};
private:

	G4bool bRotationTranslationFilesNames;
	CML2Convergence *convergence; 
	G4Timer MyTime;
	G4double loopElapsedTime;
	G4int nBeam, nLoop;
	G4bool bOnlyVisio;
};

#endif