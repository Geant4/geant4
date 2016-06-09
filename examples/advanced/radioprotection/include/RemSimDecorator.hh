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
//    ***************************
//    *                         *
//    *    RemSimDecorator.hh   *
//    *                         *          
//    ***************************
//
// $Id: RemSimDecorator.hh,v 1.8 2005/05/27 14:21:42 guatelli Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#ifndef RemSimDecorator_h
#define RemSimDecorator_h 1

#include "RemSimVGeometryComponent.hh"

class G4VPhysicalVolume;
class RemSimVGeometryComponent;

class RemSimDecorator: public RemSimVGeometryComponent
{
public:
  RemSimDecorator(RemSimVGeometryComponent*);
  ~RemSimDecorator();

  virtual void ConstructComponent(G4VPhysicalVolume*);
  virtual void DestroyComponent(); 
  virtual void ChangeThickness(G4double) = 0;

private:
   RemSimVGeometryComponent* component;
};
#endif
