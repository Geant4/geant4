// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//  
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EnergyLossTables.cc,v 1.1 1999-01-07 16:11:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// first version created by P.Urban , 06/04/1998
// modifications + "precise" functions added by L.Urban , 27/05/98
// modifications , TOF functions , 26/10/98, L.Urban

#include "G4EnergyLossTables.hh"

RWTValHashDictionary<const G4ParticleDefinition*, G4EnergyLossTablesHelper>
G4EnergyLossTables::dict(G4EnergyLossTables::HashFun);

#if defined(GNU_GCC) && ((__GNUC__ * 10 + __GNUC_MINOR__) < 28)
template RWTValHashDictionary
  <G4ParticleDefinition const *, G4EnergyLossTablesHelper>
  ::RWTValHashDictionary(unsigned int (*)(G4ParticleDefinition
  const *const &), unsigned int);
#endif

void G4EnergyLossTables::Register(
  const G4ParticleDefinition* p,
  const G4PhysicsTable* tDEDX,
  const G4PhysicsTable* tRange,
  const G4PhysicsTable* tInverseRange,
  const G4PhysicsTable* tLabTime,
  const G4PhysicsTable* tProperTime,
  G4double lowestKineticEnergy,
  G4double highestKineticEnergy,
  G4double massRatio,
  G4int NumberOfBins)
{
  dict[p]= G4EnergyLossTablesHelper(tDEDX, tRange,tInverseRange,
                    tLabTime,tProperTime,lowestKineticEnergy,
		    highestKineticEnergy, massRatio,NumberOfBins);
           
}

///////////////////////////////////////////////////////////////////////

 G4double G4EnergyLossTables::GetPreciseDEDX(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    G4Material *aMaterial)
{
  G4EnergyLossTablesHelper t= GetTables(aParticle);
  const G4PhysicsTable* rangeTable= t.theRangeTable;
  const G4PhysicsTable*  dEdxTable= t.theDEDXTable;

  G4int materialIndex = aMaterial->GetIndex();
  G4double scaledKineticEnergy = KineticEnergy*t.theMassRatio;
  G4double dEdx;
  G4bool isOut;

  if (scaledKineticEnergy<t.theLowestKineticEnergy) {

     dEdx = sqrt(scaledKineticEnergy/t.theLowestKineticEnergy)*
            (*dEdxTable)(materialIndex)->GetValue(
              t.theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>t.theHighestKineticEnergy) {

     dEdx = (*dEdxTable)(materialIndex)->GetValue(
	      t.theHighestKineticEnergy,isOut);

  } else {
    
      dEdx = (*dEdxTable)(materialIndex)->GetValue(
                          scaledKineticEnergy,isOut) ;

  }

  G4double Charge = aParticle->GetPDGCharge() ;
  G4double Chargesquare = Charge*Charge ;
  return dEdx*Chargesquare;
}

 G4double G4EnergyLossTables::GetPreciseRangeFromEnergy(
    const G4ParticleDefinition *aParticle,
    G4double KineticEnergy,
    G4Material *aMaterial)
{
  G4EnergyLossTablesHelper t= GetTables(aParticle);
  const G4PhysicsTable* rangeTable= t.theRangeTable;
  const G4PhysicsTable*  dEdxTable= t.theDEDXTable;

  G4int materialIndex = aMaterial->GetIndex();

  G4double Thigh = t.theHighestKineticEnergy*t.theLowestKineticEnergy/
                   (*rangeTable)(materialIndex)->
                   GetLowEdgeEnergy(1) ;

  G4double scaledKineticEnergy = KineticEnergy*t.theMassRatio;
  G4double Range;
  G4bool isOut;

  if (scaledKineticEnergy<t.theLowestKineticEnergy) {

    Range = sqrt(scaledKineticEnergy/t.theLowestKineticEnergy)*
            (*rangeTable)(materialIndex)->GetValue(
              t.theLowestKineticEnergy,isOut);

  } else if (scaledKineticEnergy>Thigh) {

    Range = (*rangeTable)(materialIndex)->GetValue(
	      Thigh,isOut)+
            (scaledKineticEnergy-Thigh)/
            (*dEdxTable)(materialIndex)->GetValue(
              Thigh,isOut);

  } else {
    
     Range = (*rangeTable)(materialIndex)->GetValue(
                       scaledKineticEnergy,isOut) ;

  }

  G4double Charge = aParticle->GetPDGCharge() ;
  G4double Chargesquare = Charge*Charge ;
  return Range/(Chargesquare*t.theMassRatio);
}



G4double G4EnergyLossTables::GetPreciseEnergyFromRange(
                                     const G4ParticleDefinition *aParticle,
                                           G4double range,
                                           G4Material *aMaterial)
// it returns the value of the kinetic energy for a given range
{
  G4EnergyLossTablesHelper t= GetTables(aParticle);
  const G4PhysicsTable* rangeTable= t.theRangeTable;
  const G4PhysicsTable*  dEdxTable= t.theDEDXTable;
  const G4PhysicsTable*  inverseRangeTable= t.theInverseRangeTable;

  G4double rmin,rmax,scaledrange,scaledKineticEnergy ;
  G4int materialIndex ;
  G4bool isOut ;

  materialIndex = aMaterial->GetIndex() ;
 
  rmin = (*inverseRangeTable)(materialIndex)->
                   GetLowEdgeEnergy(0) ;  

  rmax = (*inverseRangeTable)(materialIndex)->
                   GetLowEdgeEnergy(t.theNumberOfBins-2) ;

  G4double Thigh = (*inverseRangeTable)(materialIndex)->
                   GetValue(rmax,isOut) ;

  G4double Charge = aParticle->GetPDGCharge() ;
  G4double Chargesquare = Charge*Charge ;
  scaledrange = range*Chargesquare*t.theMassRatio ;

  if(scaledrange < rmin)
  {
    scaledKineticEnergy = t.theLowestKineticEnergy*
                       exp(2.*log(scaledrange/rmin)) ;

  }
  else
  {
    if(scaledrange < rmax)
    {

         scaledKineticEnergy = (*inverseRangeTable)(materialIndex)->
                              GetValue( scaledrange,isOut) ;

    }
    else
    {
      scaledKineticEnergy = Thigh +
                      (scaledrange-rmax)*
                      (*dEdxTable)(materialIndex)->
                      GetValue(Thigh,isOut) ;
    }
  }

  return scaledKineticEnergy/t.theMassRatio ;
}

