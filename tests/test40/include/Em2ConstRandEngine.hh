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
// $Id: Em2ConstRandEngine.hh,v 1.1 2004-05-28 08:51:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em2ConstRandEngine_h
#define Em2ConstRandEngine_h 1

#include "G4Types.hh"
#include "Randomize.hh"

class Em2ConstRandEngine : public HepRandomEngine
{

public:

  Em2ConstRandEngine(std::istream& is);
  Em2ConstRandEngine();
  Em2ConstRandEngine(long seed);
  virtual ~Em2ConstRandEngine();
  // Constructors and destructor

  G4double flat();
  // It returns a constant random number 0.5

  void flatArray (const G4int size, G4double* vect);
  // Fills the array "vect" of specified size with flat random values.

  void setSeed(G4long seed, G4int dum=0);
  // Sets the state of the algorithm according to seed.

  void setSeeds(const long * seeds, int dum=0);
  // Sets the state of the algorithm according to the zero terminated
  // array of seeds. Only the first seed is used.

  void saveStatus( const char filename[] = "Rand.conf" ) const;
  // Saves on file Rand.conf the current engine status.

  void restoreStatus( const char filename[] = "Rand.conf" );
  // Reads from file Rand.conf the last saved engine status
  // and restores it.

  void showStatus() const;
  // Dumps the engine status on the screen.
 
  operator unsigned int(); // 32-bit flat value, quickest of all.

  friend std::ostream& operator<< (std::ostream& os,
                                   const Em2ConstRandEngine& e);
  friend std::istream& operator>> (std::istream& is,
                                         Em2ConstRandEngine& e);

private:

  Em2ConstRandEngine(const Em2ConstRandEngine &p);
  Em2ConstRandEngine & operator = (const Em2ConstRandEngine &p);
  // Private copy constructor and assignment operator.

private:

  static G4int numEngines;
};

#endif
