//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2ContainerElement.hh,v 1.1 1999-06-17 04:41:22 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#ifndef tst2ContainerElement_h
#define tst2ContainerElement_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
class tst2ContainerElement
{
  public:
    friend class tst2ParticleContainer;
 
	tst2ContainerElement(G4ParticleDefinition * aParticle, G4int anEncoding)
    {
		particle = aParticle;
		encoding = anEncoding;
    }
	
	tst2ContainerElement(const tst2ContainerElement& right)
    {
		particle = right.particle;
		encoding = right.encoding;
    }
	
	tst2ContainerElement()
    {
		particle = 0;
		encoding = 0;
    }

    // less-than operator 
    G4int operator<(const tst2ContainerElement &right) const
	{
	  G4double massdiff 
	   =  ( this->particle->GetPDGMass() - right.particle->GetPDGMass() );
      G4int value =0;
      if ( abs(massdiff) < DBL_MIN) {
	    value = ( this->encoding < right.encoding );
      } else {
        value = ( massdiff < 0. );
      }
	  return value;
    }

    // equality operators
    G4int operator==(const tst2ContainerElement &right) const 
    {
       return (this->encoding == right.encoding);
    }

    G4int operator!=(const tst2ContainerElement &right) const 
    {
       return (this->encoding != right.encoding);
    }

  private:
	G4ParticleDefinition *   particle;
	G4int					 encoding;
};


#endif
