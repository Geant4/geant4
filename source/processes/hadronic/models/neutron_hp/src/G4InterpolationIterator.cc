#include "G4InterpolationIterator.hh"

   G4InterpolationIterator::G4InterpolationIterator() {}
   G4InterpolationIterator::G4InterpolationIterator(G4InterpolationManager * aManager)
   {
     started = false;
     theManager = aManager;
   }
   G4InterpolationIterator::~G4InterpolationIterator(){}
