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
//
// File name:     RadmonPhysicsInfoList.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfoList.cc,v 1.3 2006/06/29 16:18:53 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//

#include "RadmonPhysicsInfoList.hh"

                                                RadmonPhysicsInfoList :: RadmonPhysicsInfoList(const RadmonPhysicsInfoList & copy)
:
 infoVector(copy.infoVector)
{
}





RadmonPhysicsInfoList &                         RadmonPhysicsInfoList :: operator=(const RadmonPhysicsInfoList & copy)
{
 infoVector=copy.infoVector;
 
 return (*this);
}





G4bool                                          RadmonPhysicsInfoList :: CollidesWith(const RadmonPhysicsInfoList & other) const
{
 InfoVector::const_iterator j;
 InfoVector::const_iterator i(infoVector.begin());
 const InfoVector::const_iterator endJ(other.infoVector.end());
 const InfoVector::const_iterator endI(infoVector.end());
 
 while (i!=endI)
 {
  j=other.infoVector.begin();
  
  while (j!=endJ)
  {
   if (i->CollidesWith(*j))
    return true;
    
   j++;
  }
  
  i++;
 }
 
 return false;
}





void                                            RadmonPhysicsInfoList :: InsertPhysicsInfo(const RadmonPhysicsInfo & info)
{
 infoVector.push_back(info);
}





std::ostream &                                  operator<<(std::ostream & out, const RadmonPhysicsInfoList & infoList)
{
 const G4int n(infoList.GetNPhysicsInfos());
 
 for (G4int i(0); i<n; i++)
 {
  if (i>0)
   out << ", ";
   
  out << infoList.GetPhysicsInfo(i);
 }
 
 return out;
}
