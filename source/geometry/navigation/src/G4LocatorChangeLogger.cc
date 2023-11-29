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
// G4LocatorChangeLogger class implementation
//
// Author: John Apostolakis, 04.09.19 - First version
// --------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <locale>
// #include <cassert>

#include "G4LocatorChangeLogger.hh"

// --------------------------------------------------------------------
// Streaming operator dumping record 
//
std::ostream& operator<< ( std::ostream& os,
                           const G4LocatorChangeLogger& logger )
{
  return logger.StreamInfo(os);
}
      
// --------------------------------------------------------------------
// Stream object contents to an output stream
//
std::ostream& G4LocatorChangeLogger::StreamInfo(std::ostream& os) const
{
  G4long oldprc = os.precision(16);
  G4LocatorChangeRecord::ReportVector( os, this->fName, *this );
  os.precision(oldprc);
  return os;
}

// --------------------------------------------------------------------
//  Print the changes in start, end points in columns -- one event per row
//
std::ostream& G4LocatorChangeLogger::ReportEndChanges( std::ostream& os,
                                       const G4LocatorChangeLogger & startA,
                                       const G4LocatorChangeLogger & endB )
{
  using std::setw;
  G4int prec= 16;
  const G4bool confirm = true;
  G4long oldprc = os.precision(prec);

  auto itrecA= startA.cbegin();
  auto itrecB= endB.cbegin();

  os << "====================================================================="
     << G4endl;
  os << "  Size of individual change record:  startA : " << startA.size() 
     << "  endB : " <<   endB.size() << G4endl;
  os << "====================================================================="
     << G4endl;

  os << setw( 7 ) << "Change#"  << "  "
     << setw( 4 ) << "Iter"  << "  "
     << setw( 20 ) << "CodeLocation"  << "  "
     << setw( prec+9 ) << "Length-A (start)" << "  " 
     << setw( prec+9 ) << "Length-B (end)" << "  "
     << G4endl;
  os << "=====================================================================";

  auto eventA = (*itrecA).GetCount();
  auto eventB = (*itrecB).GetCount();

  G4bool isLastA= false;
  G4bool isLastB= false;

  G4int maxEvent = std::max( startA[ startA.size() - 1 ].GetCount() ,
                             endB[   endB.size() - 1 ].GetCount() );
  G4int prevA = -1;
  G4int prevB = -1;  

  G4bool advanceA= false, advanceB= false;
  do
  {
     advanceA= false;
     advanceB= false;

     if( ((G4int)eventA>prevA) && ((G4int)eventB>prevB) )
     {
        auto codeLocA= (*itrecA).GetLocation();
        
        os << G4endl;                   
        os << setw( 7 ) << eventA  << "  " 
           << setw( 4 ) << (*itrecA).GetIteration()  << "  "
           << setw( 3 ) << codeLocA << " " 
           << setw( 15 )
           << G4LocatorChangeRecord::GetNameChangeLocation( codeLocA ) << " "
           << setw( prec+9 ) << (*itrecA).GetLength() << "  " 
           << setw( prec+9 ) << (*itrecB).GetLength() << "  ";
        if( confirm )
        { 
          os << setw( 4 ) << (*itrecB).GetIteration()  << "  "
             << setw( 15 ) << (*itrecB).GetLocation();
        }
     }
     else
     {
        if ( (G4int)eventA > prevA )
        {
          auto codeLocA= (*itrecA).GetLocation();           
          os << G4endl;           
          os << setw( 7 ) << (*itrecA).GetCount() << "  " 
             << setw( 4 ) << (*itrecA).GetIteration()  << "  "
             << setw( 3 ) << codeLocA << " " 
             << setw( 15 )
             << G4LocatorChangeRecord::GetNameChangeLocation( codeLocA ) << " "
             << setw( prec+9 ) << (*itrecA).GetLength() << "  " 
             << setw( prec+9 ) << "       " << "  ";
        }
        else
        {
          // assert( (G4int)eventB > prevB );
          auto codeLocB= (*itrecB).GetLocation();
           
          os << G4endl;           
          os << setw( 7 ) << eventB  << "  " 
             << setw( 4 ) << (*itrecB).GetIteration()  << "  "
             << setw( 3 ) << codeLocB << " " 
             << setw( 15 )
             << G4LocatorChangeRecord::GetNameChangeLocation( codeLocB ) << " "
             << setw( prec+9 ) << "       " << "  "
             << setw( prec+9 ) << (*itrecB).GetLength() << "  " ;           
        }
     }

     prevA= eventA;
     prevB= eventB;
     
     auto nextA= itrecA;
     auto nextB= itrecB;
     
     G4int nextAct = maxEvent, nextBct = maxEvent;
     ++nextA;
     ++nextB;
     if ( nextA != startA.end() ) { nextAct = (*nextA).GetCount(); }
     if ( nextB !=   endB.end() ) { nextBct = (*nextB).GetCount(); }

     isLastA=  ( nextA >= startA.end() );
     isLastB=  ( nextB >=   endB.end() );
     
     advanceA= ( nextAct <= nextBct ) && !isLastA;
     advanceB= ( nextBct <= nextAct ) && !isLastB;

     if( advanceA )
     {
        ++itrecA;
        eventA = isLastA ? maxEvent : (*itrecA).GetCount();
     }
     
     if( advanceB )
     {
        ++itrecB;
        eventB = isLastB ? maxEvent : (*itrecB).GetCount();
     }

     // Checks
     if( isLastA !=  ( nextA == startA.end() ) )
     {
        os << G4endl;
        os << "  Checking isLastA= " << isLastA
           << " vs expected :  " << ( itrecA == startA.end() );
        os << " BAD --- ERROR " << G4endl;
     }
     if( isLastB !=  ( nextB == endB.end() ) )
     {
        os << G4endl;        
        os << "  Checking isLastB= " << isLastB
           << " vs expected :  " << ( itrecB ==  endB.end() );        
        os << " BAD --- ERROR " << G4endl;        
     }       
  } while ( ! ( isLastA && isLastB ) );
  
  os << G4endl;           
  os.precision(oldprc);
  return os;
}
