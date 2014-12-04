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
#ifdef WIN32
#include <unordered_map>
#define UOM "UOM"
#elif (defined(__GNUC__) && ((__GNUC__==4 && __GNUC_MINOR__>=1) || __GNUC__>4 ))
#include <tr1/unordered_map>
#define unordered_map tr1::unordered_map
#define UOM "UOM"
#else
#include <map>
#define unordered_map map
#define UOM "MAP"
#endif 

#include <iostream>

void Print(const char* p, int i )
{
  std::cout << p << " = " << i << std::endl;
}

int main()
{
  std::unordered_map<int,int> uom;

#ifdef __GNUC__
  Print("__GNUC__",__GNUC__);
  Print("__GNUC_MINOR__",__GNUC_MINOR__);
#endif

  std::cout << UOM << std::endl;
  return 0;
}
