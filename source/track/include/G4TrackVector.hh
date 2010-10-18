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
// $Id: G4TrackVector.hh,v 1.3 2010-10-18 23:52:04 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------
//
//  G4TrackVector.hh
//
//  class description:
//    This class keeps a List of G4Track objects. It is implemented 
//    as a RougeWave pointer ordered vector.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4TrackVector_h
#define G4TrackVector_h 1

#include <vector>
class G4Track;
class G4Step;


///////////////////////////////////////////////////
//typedef std::vector<G4Track*> G4TrackVector;
///////////////////////////////////////////////////
class  G4TrackVector : public  std::vector<G4Track*>
{
  typedef G4Track* T;

  // this class is used by G4Step!!!
  friend class G4Step;
  private:
  std::vector<const G4Track*> secondaryInCurrent;
  

  public:
  virtual ~G4TrackVector(){
             secondaryInCurrent.clear();
          }
   
  void push_back( const T& x){
                std::vector<G4Track*>::push_back(x);
                secondaryInCurrent.push_back(x);
       }

};


#endif


