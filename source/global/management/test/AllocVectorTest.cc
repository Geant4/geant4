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
// $Id: AllocVectorTest.cc,v 1.2 2006-06-29 19:04:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
#include "G4Timer.hh"
#include "G4Allocator.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <vector>

int main()
{
  G4Timer timer;
  G4double d = 1.;
  const size_t maxiter = 10000000;
  const size_t modulo  = 2000000;

  G4cout << "+++++  G4Allocator test ..........." << G4endl;
  timer.Start();

  std::vector< double, G4Allocator<double> > tvec1;
  std::vector< double, G4Allocator<double> > tvec2(maxiter);
  for (size_t i=0; i<maxiter; i++)
  {
    tvec1.push_back(d);
    tvec2.push_back(d);
    if (i%modulo == 0)
    {
      G4cout << "Dynamic vector - size: " << tvec1.size() << G4endl
             << "      ...and capacity: " << tvec1.capacity() << G4endl;
      G4cout << "Fixed vector   - size: " << tvec2.size() << G4endl
             << "      ...and capacity: " << tvec2.capacity() << G4endl;
    }
  }

  timer.Stop();  

  G4cout << "System time: " << timer.GetSystemElapsed() << G4endl
         << "User time  : " << timer.GetUserElapsed() << G4endl;

  G4cout << "+++++  std::allocator test ..........." << G4endl;
  timer.Start();

  std::vector< double > tvec3;
  std::vector< double > tvec4(maxiter);
  for (size_t j=0; j<maxiter; j++)
  {
    tvec3.push_back(d);
    tvec4.push_back(d);
    if (j%modulo == 0)
    {
      G4cout << "Dynamic vector - size: " << tvec3.size() << G4endl
             << "      ...and capacity: " << tvec3.capacity() << G4endl;
      G4cout << "Fixed vector   - size: " << tvec4.size() << G4endl
             << "      ...and capacity: " << tvec4.capacity() << G4endl;
    }
  }

  timer.Stop();  

  G4cout << "System time: " << timer.GetSystemElapsed() << G4endl
         << "User time  : " << timer.GetUserElapsed() << G4endl;

  return 0;
}
