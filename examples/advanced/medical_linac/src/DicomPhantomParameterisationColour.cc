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
// --------------------------------------------------------------------
//
// The code is part of the DICOM extended example and it was modified by :
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


#include "DicomPhantomParameterisationColour.hh"

#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

//------------------------------------------------------------------
DicomPhantomParameterisationColour::DicomPhantomParameterisationColour(G4String dicomDirectory, G4String dicomColorMap)
{
        this->dicomDirectory=dicomDirectory;
	this->dicomColorMap=dicomColorMap;
      	ReadColourData();
}


//------------------------------------------------------------------
DicomPhantomParameterisationColour::~DicomPhantomParameterisationColour()
{
}

//------------------------------------------------------------------
void DicomPhantomParameterisationColour::ReadColourData()
{
  //----- Add a G4VisAttributes for materials not defined in file;
  G4VisAttributes* blankAtt = new G4VisAttributes;
  blankAtt->SetVisibility( FALSE );
  fColours["Default"] = blankAtt;

  //----- Read file
  G4String fullColorMapName=this->dicomDirectory+this->dicomColorMap;
  std::ifstream fin(fullColorMapName);
  G4int nMate;
  G4String mateName;
  G4double cred, cgreen, cblue, copacity;
  fin >> nMate;
  for( G4int ii = 0; ii < nMate; ii++ ){
    fin >> mateName >> cred >> cgreen >> cblue >> copacity;
    G4Colour colour( cred, cgreen, cblue, copacity );
    G4VisAttributes* visAtt = new G4VisAttributes( colour );
    fColours[mateName] = visAtt;
  }

}
       

//------------------------------------------------------------------
G4Material* DicomPhantomParameterisationColour::
ComputeMaterial(const G4int copyNo, G4VPhysicalVolume * physVol, const G4VTouchable *) 
{ 
  G4Material* mate = G4PhantomParameterisation::ComputeMaterial( copyNo, physVol, 0 );
  if( physVol ) {
    G4String mateName = mate->GetName();
    unsigned int iuu = mateName.find("__");
    if( iuu != std::string::npos ) {
      mateName = mateName.substr( 0, iuu );
    }
    std::map<G4String,G4VisAttributes*>::const_iterator ite = fColours.find(mateName);
	//static G4String matNames[10000];
	//static G4int nMateNames=0;
	//static G4bool bFirst=true;
	//G4bool bFound=false;
    if( ite != fColours.end() ){

/*std::cout << physVol->GetLogicalVolume()->GetName() <<" - phyVolName - " <<mateName <<" - mateName - " << (*ite).first<< G4endl;*/
		//if (bFirst)
		//{
		//	nMateNames=0;
		//	bFirst=false;
		//	matNames[nMateNames]=mateName;
		//	std::cout<<nMateNames<<") mateName: "<<mateName << G4endl;
		//	nMateNames++;
		//}
		//else
		//{
		//	for (int j=0;j<nMateNames;j++)
		//	{
		//		if (mateName==matNames[j])
		//		{
		//			bFound=true;
		//		}
		//		if (!bFound)
		//		{
		//			matNames[nMateNames]=mateName;
		//			std::cout<<nMateNames<<") mateName: "<<mateName << G4endl;
		//			nMateNames++;
		//			break;
		//		}
		//	}
		//}
      const G4Colour col = ((*ite).second)->GetColour();
      physVol->GetLogicalVolume()->SetVisAttributes( (*ite).second ); 
    } else {
      physVol->GetLogicalVolume()->SetVisAttributes( (*(fColours.begin()) ).second ); // set it as unseen
    }
  }  

  return mate;
}

