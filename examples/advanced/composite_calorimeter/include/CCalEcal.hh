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

class CCalEcal: public CCalDetector {
public:
  //Constructor and Destructor
  CCalEcal(const G4String &name):
    CCalDetector(name) {}
  virtual ~CCalEcal();

  //Get Methods
  G4String getGenMat()               const {return genMat;}
  double getWidBox()                 const {return widBox;}
  double getLengBox()                const {return lengBox;}
  double getXpos()                   const {return xpos;}
  double getYpos()                   const {return ypos;}
  double getZpos()                   const {return zpos;}
  double getThetaX()                 const {return thetaX;}
  double getPhiX()                   const {return phiX;}
  double getThetaY()                 const {return thetaY;}
  double getPhiY()                   const {return phiY;}
  double getThetaZ()                 const {return thetaZ;}
  double getPhiZ()                   const {return phiZ;}
  G4String getLayMat()               const {return layMat;}
  int getLayNum()                    const {return layNum;}
  double getLayRadius()              const {return layRadius;}
  double getLayAngle()               const {return layAngle;}
  double getLengFront()              const {return lengFront;}
  double getLayPar(unsigned int i)   const {return layPar[i];}
  G4String getCrystMat()             const {return crystMat;}
  int getCrystNum()                  const {return crystNum;}
  double getCrystLength()            const {return crystLength;}
  double getCrystTol()               const {return crystTol;}
  double getCrystPar(unsigned int i) const {return crystPar[i];}
  G4String getSuppMat()              const {return suppMat;}
  double getDxSupp()                 const {return dxSupp;}
  double getDySupp()                 const {return dySupp;}
  double getDzSupp()                 const {return dzSupp;}
  double getDistSupp()               const {return distSupp;}


protected:
  virtual int readFile();
  virtual void constructDaughters();

private:
  G4String genMat;                            //General material
  double widBox, lengBox;                     //Box parameters
  double xpos, ypos, zpos;                    //Translation matrix definition
  double thetaX,phiX,thetaY,phiY,thetaZ,phiZ; //Rotation matrix definition
  G4String layMat;                            //Material for the layers
  int layNum;                                 //Layer numbers
  double layRadius,layAngle;                  //Positioning parameters
  double lengFront;                           //Distance from front
  double layPar[5];                           //Layer parameters
  G4String crystMat;                          //Material for the crystals
  int crystNum;                               //Crystal numbers
  double crystLength;                         //Crystal length
  double crystTol;                            //Tolerance for position
  double crystPar[5];                         //Crystal parameters
  G4String suppMat;                           //Material of support system
  double dxSupp, dySupp, dzSupp;              //Dimension of support material
  double distSupp;                            //Separation of support material
};

#endif
















