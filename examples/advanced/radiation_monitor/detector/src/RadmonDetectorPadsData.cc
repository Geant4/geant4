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
// File name:     RadmonDetectorPadsData.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorPadsData.cc,v 1.3 2006/06/29 16:14:11 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//

// Include files
#include "RadmonDetectorPadsData.hh"
#include "RadmonMessenger.hh"

                                                RadmonDetectorPadsData :: RadmonDetectorPadsData(const RadmonDetectorPadsData & copy)
:
 padWidth(copy.padWidth),
 padHeight(copy.padHeight),
 padMaterial(copy.padMaterial),
 padPositions(copy.padPositions),
 padRotations(copy.padRotations)
{
}
 
 


 
G4int                                           RadmonDetectorPadsData :: GetNPads(void) const
{
 return padPositions.size();
}



const G4ThreeVector &                           RadmonDetectorPadsData :: GetPosition(G4int index) const
{
 return padPositions[index];
}



const G4RotationMatrix &                        RadmonDetectorPadsData :: GetRotation(G4int index) const
{
 return padRotations[index];
}





void                                            RadmonDetectorPadsData :: AppendPositionAndRotation(const G4ThreeVector & position, const G4RotationMatrix & rotation)
{
 padPositions.push_back(position);
 padRotations.push_back(rotation);
}





RadmonDetectorPadsData &                        RadmonDetectorPadsData :: operator=(const RadmonDetectorPadsData & copy)
{
 padWidth=copy.padWidth;
 padHeight=copy.padHeight;
 padMaterial=copy.padMaterial;
 padPositions=copy.padPositions;
 padRotations=copy.padRotations;
 
 return (*this);
}





bool                                            RadmonDetectorPadsData :: ReadPositionsAndRotationsFromString(const G4String & positionStr)
{
 // (xx mm, yy mm, [delta deg])
 
 G4bool somethingDone(false);
 G4String positions(positionStr);
 str_size i(0);
 const str_size size(positions.length());
 
 while (i<size)
 {
  if (positions(i)==' ')
  {
   i++;
   continue;
  }
  
  if (positions(i)=='(')
  {
   str_size j(positions.index(')', i));
   
   if (j==G4String::npos)
    return false;
    
   somethingDone=true;
   G4String content(positions(i+1, j-i-1));
   
   G4double x;
   G4double y;
   G4double delta;
   if (!ProcessElement(x, y, delta, content))
    return false;
    
   G4RotationMatrix rotation(G4RotationMatrix::IDENTITY);
   rotation.setDelta(delta);
   AppendPositionAndRotation(G4ThreeVector(x, y, 0.), rotation);
   
   i=j+1;
   continue;
  }
  
  return false;
 }
 
 return somethingDone;
}



bool                                            RadmonDetectorPadsData :: ProcessElement(G4double & x, G4double & y, G4double & delta, const G4String & content)
{
 G4String cont(content);
 
 str_size i(cont.index(','));
 str_size size(cont.length());
 while (size>1)
 {
  if (cont(size-1)!=' ')
   break;

  size--;
 }
 
 if (i==G4String::npos)
  return false;
 
 if (!ReadUmis(x, cont(0, i), "Length"))
  return false;
 

 do
 {
  i++;
  
  if (i>=size)
   return false;
 }
 while (cont(i)==' ');

 str_size j(cont.index(',', i));
 
 if (j==G4String::npos)
 {
  delta=0.*deg;

  return ReadUmis(y, cont(i, size-i), "Length");
 }
 
 if (!ReadUmis(y, cont(i+1, j-i-1), "Length"))
  return false;
  
 do
 {
  j++;
  
  if (j>=size)
   return false;
 }
 while (cont(j)==' ');
  
 return ReadUmis(y, cont(j, size-j), "Angle");
}



bool                                            RadmonDetectorPadsData :: ReadUmis(G4double & value, const G4String & text, const char * category)
{
 G4String args[2];

 if (!RadmonMessenger::ProcessArguments(text, 2, args))
  return false;

 value=RadmonMessenger::GetUnit(args[1], category);
 if (value<=0.)
  return false;
  
 value*=G4UIcommand::ConvertToDouble(args[0]);
  
 return true;
}

