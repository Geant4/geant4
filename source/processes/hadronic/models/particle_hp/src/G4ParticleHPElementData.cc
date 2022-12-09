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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 02-08-06 Modified Harmonise to reslove cross section trouble at high-end. T. KOI
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPElementData.hh"

  G4ParticleHPElementData::G4ParticleHPElementData()
  {
     precision = 0.02;
     theFissionData = new G4ParticleHPVector;
     theCaptureData = new G4ParticleHPVector;
     theElasticData = new G4ParticleHPVector;
     theInelasticData = new G4ParticleHPVector;
     theIsotopeWiseData = 0;
     theBuffer = nullptr;
 }
  
  G4ParticleHPElementData::~G4ParticleHPElementData()
  {
    delete theFissionData;
    delete theCaptureData;
    delete theElasticData;
    delete theInelasticData;
    delete [] theIsotopeWiseData;
  }
  
  void G4ParticleHPElementData::Init(G4Element * theElement, G4ParticleDefinition* projectile, const char* dataDirVariable )
  {
    G4int count = (G4int)theElement->GetNumberOfIsotopes();
      if(count == 0)
        count += theStableOnes.GetNumberOfIsotopes((G4int)theElement->GetZ());
    theIsotopeWiseData = new G4ParticleHPIsoData[count];
    // filename = ein data-set je isotope.
    count = 0;
    G4int nIso = (G4int)theElement->GetNumberOfIsotopes();
    G4int Z = (G4int)theElement->GetZ();

    if(nIso!=0)
    {
      for (G4int i1=0; i1<nIso; ++i1)
      {
        G4int A = theElement->GetIsotope(i1)->GetN();
        G4int M = theElement->GetIsotope(i1)->Getm();
        G4double frac = theElement->GetRelativeAbundanceVector()[i1]/CLHEP::perCent;
        //UpdateData(A, Z, count++, frac);
        UpdateData(A, Z, M, count++, frac, projectile, dataDirVariable);
      }
    }
    else
    {
      G4int first = theStableOnes.GetFirstIsotope(Z);
      for(G4int i1=0; 
          i1<theStableOnes.GetNumberOfIsotopes((G4int)theElement->GetZ()); ++i1)
      {
        G4int A = theStableOnes.GetIsotopeNucleonCount(first+i1);
        G4double frac = theStableOnes.GetAbundance(first+i1);
        UpdateData(A, Z, count++, frac, projectile, dataDirVariable);
      }
    }
    theElasticData->ThinOut(precision);
    if( projectile == G4Neutron::Neutron() ) theInelasticData->ThinOut(precision);

    theCaptureData->ThinOut(precision);
    theFissionData->ThinOut(precision);
  }
  
  void G4ParticleHPElementData::UpdateData(G4int A, G4int Z, G4int M, G4int index, G4double abundance, G4ParticleDefinition* projectile, const char* dataDirVariable )
  {
    //Reads in the Data, using G4ParticleHPIsoData[], and its Init
    //theIsotopeWiseData[index].Init(A, Z, abundance);
    theIsotopeWiseData[index].Init(A, Z, M, abundance,projectile, dataDirVariable);

    theBuffer = theIsotopeWiseData[index].MakeElasticData();
    Harmonise(theElasticData, theBuffer);
    delete theBuffer;
    
    theBuffer = theIsotopeWiseData[index].MakeInelasticData();
    Harmonise(theInelasticData, theBuffer);
    delete theBuffer;
    
    theBuffer = theIsotopeWiseData[index].MakeCaptureData();
    Harmonise(theCaptureData, theBuffer);
    delete theBuffer;
    
    theBuffer = theIsotopeWiseData[index].MakeFissionData();
    Harmonise(theFissionData, theBuffer);
    delete theBuffer;
  }
  
  void G4ParticleHPElementData::Harmonise(G4ParticleHPVector *& theStore, G4ParticleHPVector * theNew)
  {
    if(theNew == 0) { return; }
    G4int s_tmp = 0, n=0, m_tmp=0;
    G4ParticleHPVector * theMerge = new G4ParticleHPVector(theStore->GetVectorLength());
    while ( theStore->GetEnergy(s_tmp)<theNew->GetEnergy(0)&&s_tmp<theStore->GetVectorLength() ) // Loop checking, 11.05.2015, T. Koi
    {
      theMerge->SetData(m_tmp++, theStore->GetEnergy(s_tmp), theStore->GetXsec(s_tmp));
      ++s_tmp;
    }
    G4ParticleHPVector *active = theStore;
    G4ParticleHPVector * passive = theNew;
    G4ParticleHPVector * tmp;
    G4int a = s_tmp, p = n, t;
    while (a<active->GetVectorLength()&&p<passive->GetVectorLength()) // Loop checking, 11.05.2015, T. Koi
    {
      if(active->GetEnergy(a) <= passive->GetEnergy(p))
      {
        theMerge->SetData(m_tmp, active->GetEnergy(a), active->GetXsec(a));
        G4double x  = theMerge->GetEnergy(m_tmp);
        G4double y = std::max(0., passive->GetXsec(x)); 
        theMerge->SetData(m_tmp, x, theMerge->GetXsec(m_tmp)+y);
        ++m_tmp;
        ++a;
      }
      else
      {
        tmp = active; t=a;
        active = passive; a=p;
        passive = tmp; p=t;
      }
    }
    while (a!=active->GetVectorLength()) // Loop checking, 11.05.2015, T. Koi
    {
      theMerge->SetData(m_tmp++, active->GetEnergy(a), active->GetXsec(a));
      ++a;
    }
    while (p!=passive->GetVectorLength()) // Loop checking, 11.05.2015, T. Koi
    {
      G4double x = passive->GetEnergy(p);
      G4double y = std::max(0., active->GetXsec(x));
      theMerge->SetData(m_tmp++, x, passive->GetXsec(p)+y);
      ++p;
    }
    delete theStore;
    theStore = theMerge;
  }

  G4ParticleHPVector * G4ParticleHPElementData::MakePhysicsVector(G4Element * theElement,
      						G4ParticleDefinition * projectile,
						G4ParticleHPFissionData* theSet,
  				      char* dataDirVariable)
  {
    if(projectile != G4Neutron::Neutron()) throw G4HadronicException(__FILE__, __LINE__, "not a neutron");
   Init ( theElement, projectile, dataDirVariable );
   return GetData(theSet);
  }
  G4ParticleHPVector * G4ParticleHPElementData::MakePhysicsVector(G4Element * theElement,
		                      G4ParticleDefinition * projectile,
                                      G4ParticleHPCaptureData * theSet,
  				      char* dataDirVariable)
  {
    if(projectile != G4Neutron::Neutron()) throw G4HadronicException(__FILE__, __LINE__, "not a neutron");
   Init ( theElement, projectile, dataDirVariable );
   return GetData(theSet);
  }
  G4ParticleHPVector * G4ParticleHPElementData::MakePhysicsVector(G4Element * theElement,
				      G4ParticleDefinition * projectile,
                                      G4ParticleHPElasticData * theSet,
  				      char* dataDirVariable)
  {
    if(projectile != G4Neutron::Neutron()) throw G4HadronicException(__FILE__, __LINE__, "not a neutron");
   Init ( theElement, projectile, dataDirVariable );
   return GetData(theSet);
  }
    G4ParticleHPVector * G4ParticleHPElementData::MakePhysicsVector(G4Element * theElement,
				      G4ParticleDefinition * projectile,
                                      G4ParticleHPInelasticData * theSet,
	   			      char* dataDirVariable)
  {
    if(projectile != G4Neutron::Neutron()) throw G4HadronicException(__FILE__, __LINE__, "not a neutron");
   Init ( theElement, projectile, dataDirVariable );
   return GetData(theSet);
  }
