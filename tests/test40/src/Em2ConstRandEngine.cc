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
// $Id: Em2ConstRandEngine.cc,v 1.1 2004-05-28 08:51:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em2ConstRandEngine.hh"
#include "G4ios.hh"

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 

// number of instances with automatic seed selection
G4int Em2ConstRandEngine::numEngines = 0;

Em2ConstRandEngine::Em2ConstRandEngine(G4long) 
{
}

Em2ConstRandEngine::Em2ConstRandEngine()
{
}

Em2ConstRandEngine::Em2ConstRandEngine(std::istream& is)
{
   is >> *this;
}

Em2ConstRandEngine::~Em2ConstRandEngine() {}

Em2ConstRandEngine::Em2ConstRandEngine(const Em2ConstRandEngine &p)
  : HepRandomEngine(p)
{
  *this = p;
}

Em2ConstRandEngine&
Em2ConstRandEngine::operator = (const Em2ConstRandEngine &p)
{
  *this = p;
  return *this;
}

void Em2ConstRandEngine::setSeed(G4long, G4int)
{
}

void Em2ConstRandEngine::setSeeds(const G4long*, G4int)
{
}

void Em2ConstRandEngine::saveStatus( const char filename[] ) const
{
   std::ofstream outFile( filename, std::ios::out ) ;

   if (!outFile.bad())
   {
     outFile << theSeed << std::endl;
   }
}

void Em2ConstRandEngine::restoreStatus( const char filename[] )
{
   std::ifstream inFile( filename, std::ios::in);

   if (!inFile.bad() && !inFile.eof())
   {
     inFile >> theSeed;
     setSeed(theSeed,0);
   }
}

void Em2ConstRandEngine::showStatus() const
{
   std::cout << std::endl;
   std::cout << "------ Em2ConstRand engine status ------" << std::endl;
   std::cout << " Initial seed  = " << theSeed << std::endl;
   std::cout << "----------------------------------------" << std::endl;
}

G4double Em2ConstRandEngine::flat()
{
   return 0.555555;
}

void Em2ConstRandEngine::flatArray(const G4int size, G4double* vect)
{
   for (G4int i=0; i<size; ++i)
     vect[i]=flat();
}

Em2ConstRandEngine::operator unsigned int()
{
  return 0;
}

std::ostream & operator << ( std::ostream& os, const Em2ConstRandEngine& e ) 
{
     char beginMarker[] = "Em2ConstRandEngine-begin";
     char endMarker[]   = "Em2ConstRandEngine-end";

     os << " " << beginMarker << " ";
     os << e.theSeed << " ";
     os << endMarker << " ";
     return os;
}

std::istream & operator >> ( std::istream& is, Em2ConstRandEngine& e )
{
  char beginMarker [MarkerLen];
  char endMarker   [MarkerLen];

  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"Em2ConstRandEngine-begin"))
  {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nInput stream mispositioned or"
	       << "\nEm2ConstRandEngine state description missing or"
	       << "\nwrong engine type found." << std::endl;
     return is;
  }
  is >> e.theSeed;
  is >> std::ws;
  is.width(MarkerLen);  
  is >> endMarker;
  if (strcmp(endMarker,"Em2ConstRandEngine-end"))
  {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nEm2ConstRandEngine state description incomplete."
	       << "\nInput stream is probably mispositioned now." << std::endl;
     return is;
  }

  e.setSeed(e.theSeed,0);

  return is;
}
