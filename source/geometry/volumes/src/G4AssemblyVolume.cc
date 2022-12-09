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
// Class G4AssemblyVolume - implementation
//
// ----------------------------------------------------------------------

#include "G4AssemblyVolume.hh"
#include "G4AssemblyStore.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ReflectionFactory.hh"

#include <sstream>

G4ThreadLocal unsigned int G4AssemblyVolume::fsInstanceCounter = 0;

// Default constructor
//
G4AssemblyVolume::G4AssemblyVolume()
{
  InstanceCountPlus();
  SetAssemblyID( GetInstanceCount() );
  SetImprintsCount( 0 );
  G4AssemblyStore* aStore = G4AssemblyStore::GetInstance();
  if (aStore->GetAssembly(fAssemblyID,false))
  {
    std::ostringstream message;
    message << "The assembly has NOT been registered !" << G4endl
            << "          Assembly " << fAssemblyID
            << " already existing in store !" << G4endl;
    G4Exception("G4AssemblyVolume::G4AssemblyVolume()", "GeomVol1001",
                JustWarning, message);
  }
  else
  {
    aStore->Register(this);
  }
}

// Composing constructor
//
G4AssemblyVolume::G4AssemblyVolume( G4LogicalVolume* volume,
                                    G4ThreeVector& translation,
                                    G4RotationMatrix* rotation )
{
  InstanceCountPlus();
  SetAssemblyID( GetInstanceCount() );
  SetImprintsCount( 0 );
  AddPlacedVolume(volume, translation, rotation);
  G4AssemblyStore* aStore = G4AssemblyStore::GetInstance();
  if (aStore->GetAssembly(fAssemblyID,false))
  {
    std::ostringstream message;
    message << "The assembly has NOT been registered !" << G4endl
            << "          Assembly " << fAssemblyID
            << " already existing in store !" << G4endl;
    G4Exception("G4Assembly::G4Assembly()", "GeomVol1001",
                JustWarning, message);
  }
  else
  {
    aStore->Register(this);
  }
}

// Destructor
//
G4AssemblyVolume::~G4AssemblyVolume()
{
  std::size_t howmany = fTriplets.size();
  if( howmany != 0 )
  {
    for( std::size_t i = 0; i < howmany; ++i )
    {
      G4RotationMatrix* pRotToClean = fTriplets[i].GetRotation();
      if( pRotToClean != 0 )
      {
        delete pRotToClean;
      }
    }
  }
  fTriplets.clear();
  
  // No need to delete physical volumes, as these are deleted
  // by physical-volumes store.
  //
  fPVStore.clear();
  InstanceCountMinus();
  G4AssemblyStore::GetInstance()->DeRegister(this);
}

// Add and place the given volume according to the specified
// translation and rotation.
//
// The rotation matrix passed in can be 0 = identity or an address even of an
// object on the upper stack frame. During assembly imprint, it creates anyway
// a new matrix and keeps track of it so it can delete it later at destruction
// time.
// This policy has been adopted since user has no control on the way the
// rotations are combined.
//
void G4AssemblyVolume::AddPlacedVolume( G4LogicalVolume*  pVolume,
                                        G4ThreeVector&    translation,
                                        G4RotationMatrix* pRotation )
{
  G4RotationMatrix*  toStore  = new G4RotationMatrix;
  
  if( pRotation != nullptr )  { *toStore = *pRotation; }
  
  G4AssemblyTriplet toAdd( pVolume, translation, toStore );
  fTriplets.push_back( toAdd );
}

// Add and place the given volume according to the specified transformation
//
void G4AssemblyVolume::AddPlacedVolume( G4LogicalVolume*  pVolume,
                                        G4Transform3D&    transformation )
{
  // Decompose transformation
  G4Scale3D     scale;
  G4Rotate3D    rotation;
  G4Translate3D translation;
  transformation.getDecomposition(scale, rotation, translation);

  G4ThreeVector      v = translation.getTranslation();
  G4RotationMatrix*  r = new G4RotationMatrix;
                    *r = rotation.getRotation();
  
  G4bool isReflection = false;
  if (scale(0,0)*scale(1,1)*scale(2,2) < 0.)  { isReflection = true; }

  G4AssemblyTriplet toAdd( pVolume, v, r, isReflection );
  fTriplets.push_back( toAdd );
}

// Add and place the given assembly volume according to the specified
// translation and rotation.
//
void G4AssemblyVolume::AddPlacedAssembly( G4AssemblyVolume* pAssembly,
                                          G4ThreeVector&    translation,
                                          G4RotationMatrix* pRotation )
{
  G4RotationMatrix*  toStore  = new G4RotationMatrix;
  
  if( pRotation != nullptr )  { *toStore = *pRotation; }
  
  G4AssemblyTriplet toAdd( pAssembly, translation, toStore );
  fTriplets.push_back( toAdd );
}

// Add and place the given assembly volume according to the specified 
// transformation
//
void G4AssemblyVolume::AddPlacedAssembly( G4AssemblyVolume* pAssembly,
                                          G4Transform3D&    transformation )
{
  // Decompose transformation
  //
  G4Scale3D     scale;
  G4Rotate3D    rotation;
  G4Translate3D translation;
  transformation.getDecomposition(scale, rotation, translation);

  G4ThreeVector      v = translation.getTranslation();
  G4RotationMatrix*  r = new G4RotationMatrix;
                    *r = rotation.getRotation();
  
  G4bool isReflection = false;
  if (scale(0,0)*scale(1,1)*scale(2,2) < 0.)  { isReflection = true; }
  
  G4AssemblyTriplet toAdd( pAssembly, v, r, isReflection );
  fTriplets.push_back( toAdd );
}

// Create an instance of an assembly volume inside of the specified
// mother volume. This works analogically to making stamp imprints.
// This method makes use of the Geant4 affine transformation class.
// The algorithm is defined as follows:
//  
// Having rotation matrix Rm and translation vector Tm to be applied
// inside the mother and rotation matrix Ra and translation vector Ta
// to be applied inside the assembly itself for each of the participating
// volumes the resulting transformation is
//  
// Tfinal = Ta * Tm
//  
// where Ta and Tm are constructed as
//  
//        -1                                     -1
// Ta = Ra  * Ta           and            Tm = Rm  * Tm
//  
// which in words means that we create first the affine transformations
// by inverse rotation matrices and translations for mother and assembly.
// The resulting final transformation to be applied to each of the
// participating volumes is their product.
//
// IMPORTANT NOTE!
// The order of multiplication is reversed when comparing to CLHEP 3D
// transformation matrix(G4Transform3D class).
//  
// The rotation matrix passed in can be 0 = identity or an address even of an
// object on the upper stack frame. During assembly imprint, it creates anyway
// a new matrix and keeps track of it so it can delete it later at destruction
// time.
// This policy has been adopted since user has no control on the way the
// rotations are combined.
// 
// If the assembly volume contains assembly (a'), the function is called
// recursively with composed transformation:
//
// Tanew =  Ta * Ta'
//
void G4AssemblyVolume::MakeImprint( G4AssemblyVolume* pAssembly,
                                    G4LogicalVolume*  pMotherLV,
                                    G4Transform3D&    transformation,
                                    G4int copyNumBase,
                                    G4bool surfCheck )
{
  std::size_t numberOfDaughters;
    
  if( copyNumBase == 0 )
  {
    numberOfDaughters = pMotherLV->GetNoDaughters();
  }
  else
  {
    numberOfDaughters = copyNumBase;
  }

  // We start from the first available index
  //
  ++numberOfDaughters;

  ImprintsCountPlus();
  
  auto triplets = pAssembly->fTriplets;

  // store the transformation in a container (for GDML persistency)
  fImprintsTransf[GetImprintsCount()] = transformation;

  for( std::size_t i = 0; i < triplets.size(); ++i )
  {
    G4Transform3D Ta( *(triplets[i].GetRotation()),
                      triplets[i].GetTranslation() );
    if ( triplets[i].IsReflection() )  { Ta = Ta * G4ReflectZ3D(); }

    G4Transform3D Tfinal = transformation * Ta;
    
    if ( triplets[i].GetVolume() )
    {
      // Generate the unique name for the next PV instance
      // The name has format:
      //
      // av_WWW_impr_XXX_YYY_ZZZ
      // where the fields mean:
      // WWW - assembly volume instance number
      // XXX - assembly volume imprint number
      // YYY - the name of a log. volume we want to make a placement of
      // ZZZ - the log. volume index inside the assembly volume
      //
      std::stringstream pvName;
      pvName << "av_"
             << GetAssemblyID()
             << "_impr_"
             << GetImprintsCount()
             << "_"
             << triplets[i].GetVolume()->GetName().c_str()
             << "_pv_"
             << i
             << std::ends;

      // Generate a new physical volume instance inside a mother
      // (as we allow 3D transformation use G4ReflectionFactory to 
      //  take into account eventual reflection)
      //
      G4PhysicalVolumesPair pvPlaced
        = G4ReflectionFactory::Instance()->Place( Tfinal,
                                                  pvName.str().c_str(),
                                                  triplets[i].GetVolume(),
                                                  pMotherLV,
                                                  false,
                                                  G4int(numberOfDaughters+i),
                                                  surfCheck );

      // Register the physical volume created by us so we can delete it later
      //
      fPVStore.push_back( pvPlaced.first );
      if ( pvPlaced.second )  { fPVStore.push_back( pvPlaced.second ); }

      // Here I want to save the imprint transformation to some container, so I can retrieve the original values later
      // imprintTrans[GetImprintsCount()] = transformation;
      
    }
    else if ( triplets[i].GetAssembly() )
    {
      // Place volumes in this assembly with composed transformation
      //
      MakeImprint( triplets[i].GetAssembly(), pMotherLV,
                   Tfinal, (G4int)i*100+copyNumBase, surfCheck ); 
    }
    else
    {
      G4Exception("G4AssemblyVolume::MakeImprint(..)",
                  "GeomVol0003", FatalException,
                  "Triplet has no volume and no assembly");
    }  
  }  
}    

void G4AssemblyVolume::MakeImprint( G4LogicalVolume*  pMotherLV,
                                    G4ThreeVector&    translationInMother,
                                    G4RotationMatrix* pRotationInMother,
                                    G4int copyNumBase,
                                    G4bool surfCheck )
{
  // If needed user can specify explicitely the base count from which to start
  // off for the generation of phys. vol. copy numbers.
  // The old behaviour is preserved when copyNumBase == 0, e.g. the generated
  // copy numbers start from the count equal to current number of daughter
  // volumes before an imprint is made

  // Compose transformation
  //
  if( pRotationInMother == nullptr )
  {
    // Make it by default an indentity matrix
    //
    pRotationInMother =
      const_cast<G4RotationMatrix*>( &G4RotationMatrix::IDENTITY );
  }

  G4Transform3D transform( *pRotationInMother,
                            translationInMother );
  MakeImprint(this, pMotherLV, transform, copyNumBase, surfCheck);
}

void G4AssemblyVolume::MakeImprint( G4LogicalVolume*  pMotherLV,
                                    G4Transform3D&    transformation,
                                    G4int copyNumBase,
                                    G4bool surfCheck )
{
  // If needed user can specify explicitely the base count from which to start
  // off for the generation of phys. vol. copy numbers.
  // The old behaviour is preserved when copyNumBase == 0, e.g. the generated
  // copy numbers start from the count equal to current number of daughter
  // volumes before a imprint is made

  MakeImprint(this, pMotherLV, transformation, copyNumBase, surfCheck);
}

unsigned int G4AssemblyVolume::GetInstanceCount() const
{
  return G4AssemblyVolume::fsInstanceCounter;
}

void G4AssemblyVolume::SetInstanceCount( unsigned int value )
{
  G4AssemblyVolume::fsInstanceCounter = value;
}

void G4AssemblyVolume::InstanceCountPlus()
{
  G4AssemblyVolume::fsInstanceCounter++;
}

void G4AssemblyVolume::InstanceCountMinus()
{
  G4AssemblyVolume::fsInstanceCounter--;
}
