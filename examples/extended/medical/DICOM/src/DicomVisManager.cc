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
//
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************

#ifdef G4VIS_USE 

#include "DicomVisManager.hh"

#ifdef G4VIS_USE_OPENGLX
// Supported drivers...
#include "G4OpenGLImmediateX.hh"
#include "G4OpenGLStoredX.hh"
#endif

DicomVisManager::DicomVisManager ()
{}


void DicomVisManager::RegisterGraphicsSystems ()
{

#ifdef G4VIS_USE_OPENGLX
    RegisterGraphicsSystem (new G4OpenGLImmediateX);
    RegisterGraphicsSystem (new G4OpenGLStoredX );
#endif
     if (fVerbose > 0)
    {
        G4cout << "\nYou have successfully chosen to use the following graphics systems." << G4endl;
        PrintAvailableGraphicsSystems ();
    }
}
#endif




