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
///////////////////////////////////////////////////////////////////////////////
// File: CCalEcal.hh
// Description: CCalEcal Geometry factory class for crystal matrix
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalEcal_h
#define CCalEcal_h 1

#include "CCalDetector.hh"

class CCalEcal: public CCalDetector
{

public:
  //Constructor and Destructor
  CCalEcal(const G4String &name):
    CCalDetector(name) {}
  virtual ~CCalEcal();

  //Get Methods
  G4String getGenMat()                 const {return genMat;}
  G4double getWidBox()                 const {return widBox;}
  G4double getLengBox()                const {return lengBox;}
  G4double getXpos()                   const {return xpos;}
  G4double getYpos()                   const {return ypos;}
  G4double getZpos()                   const {return zpos;}
  G4double getThetaX()                 const {return thetaX;}
  G4double getPhiX()                   const {return phiX;}
  G4double getThetaY()                 const {return thetaY;}
  G4double getPhiY()                   const {return phiY;}
  G4double getThetaZ()                 const {return thetaZ;}
  G4double getPhiZ()                   const {return phiZ;}
  G4String getLayMat()                 const {return layMat;}
  G4int getLayNum()                    const {return layNum;}
  G4double getLayRadius()              const {return layRadius;}
  G4double getLayAngle()               const {return layAngle;}
  G4double getLengFront()              const {return lengFront;}
  G4double getLayPar(unsigned int i)   const {return layPar[i];}
  G4String getCrystMat()               const {return crystMat;}
  G4int getCrystNum()                  const {return crystNum;}
  G4double getCrystLength()            const {return crystLength;}
  G4double getCrystTol()               const {return crystTol;}
  G4double getCrystPar(unsigned int i) const {return crystPar[i];}
  G4String getSuppMat()                const {return suppMat;}
  G4double getDxSupp()                 const {return dxSupp;}
  G4double getDySupp()                 const {return dySupp;}
  G4double getDzSupp()                 const {return dzSupp;}
  G4double getDistSupp()               const {return distSupp;}


protected:
  virtual G4int readFile();
  virtual void constructDaughters();

private:
  G4String genMat;                            //General material
  G4double widBox, lengBox;                   //Box parameters
  G4double xpos, ypos, zpos;                  //Translation matrix definition
  G4double thetaX,phiX,thetaY,phiY,thetaZ,phiZ; //Rotation matrix definition
  G4String layMat;                            //Material for the layers
  G4int layNum;                               //Layer numbers
  G4double layRadius,layAngle;                //Positioning parameters
  G4double lengFront;                         //Distance from front
  G4double layPar[5];                         //Layer parameters
  G4String crystMat;                          //Material for the crystals
  G4int crystNum;                             //Crystal numbers
  G4double crystLength;                       //Crystal length
  G4double crystTol;                          //Tolerance for position
  G4double crystPar[5];                       //Crystal parameters
  G4String suppMat;                           //Material of support system
  G4double dxSupp, dySupp, dzSupp;            //Dimension of support material
  G4double distSupp;                          //Separation of support material
};

#endif
