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
// ^INAIL DIPIA - ex ISPESL  and INFN Roma, gruppo collegato SanitÃ , Italy
// *Istituto Superiore di SanitÃ  and INFN Roma, gruppo collegato SanitÃ , Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#ifndef CML2_DicomMessengerH
#define CML2_DicomMessengerH


#include "globals.hh"
#include "G4UImessenger.hh"
#include "DicomDetectorConstruction.hh"

class G4UIcmdWithAnInteger;
class G4UImessenger;
class G4UIcmdWithAString;
class DicomDetectorConstruction;

class CML2_DicomMessenger : public G4UImessenger 
{
public:
	CML2_DicomMessenger(DicomDetectorConstruction *dicomDetectorConstruction);
	~CML2_DicomMessenger(void);
	void SetNewValue(G4UIcommand* cmd, G4String newValue);
private:
	DicomDetectorConstruction *dicomDetectorConstruction;
	G4UIcmdWithAString *dicomDirectory, *fileData, *calibrationDensityFileName, *dicomColorMap;
	G4UIcmdWithAnInteger *nROGVoxelsX, *nROGVoxelsY, *nROGVoxelsZ;
};

#endif
