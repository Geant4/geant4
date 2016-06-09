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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 02-08-06 Modified Harmonise to reslove cross section trouble at high-end. T. KOI
//
#include "G4NeutronHPElementData.hh"
#include "G4SystemOfUnits.hh"

  G4NeutronHPElementData::G4NeutronHPElementData()
  {
     precision = 0.02;
     theFissionData = new G4NeutronHPVector;
     theCaptureData = new G4NeutronHPVector;
     theElasticData = new G4NeutronHPVector;
     theInelasticData = new G4NeutronHPVector;
    theIsotopeWiseData = 0;
  }
  
  G4NeutronHPElementData::~G4NeutronHPElementData()
  {
    delete theFissionData;
    delete theCaptureData;
    delete theElasticData;
    delete theInelasticData;
    delete [] theIsotopeWiseData;
  }
  
  void G4NeutronHPElementData::Init(G4Element * theElement)  
  {
    G4int count = theElement->GetNumberOfIsotopes();
      if(count == 0) count +=
         theStableOnes.GetNumberOfIsotopes(static_cast<G4int>(theElement->GetZ()));
    theIsotopeWiseData = new G4NeutronHPIsoData[count];
    // filename = ein data-set je isotope.
    count = 0;
    G4int nIso = theElement->GetNumberOfIsotopes();
    G4int Z = static_cast<G4int> (theElement->GetZ());
    //G4int i1;
    if(nIso!=0)
    {
      for (G4int i1=0; i1<nIso; i1++)
      {
//        G4cout <<" Init: normal case"<<G4endl;
        G4int A = theElement->GetIsotope(i1)->GetN();
        G4int M = theElement->GetIsotope(i1)->Getm();
        G4double frac = theElement->GetRelativeAbundanceVector()[i1]/perCent;
        //UpdateData(A, Z, count++, frac);
        UpdateData(A, Z, M, count++, frac);
      }
    }else{
//      G4cout <<" Init: theStableOnes case: Z="<<Z<<G4endl;
      G4int first = theStableOnes.GetFirstIsotope(Z);
//      G4cout <<"first="<<first<<" "<<theStableOnes.GetNumberOfIsotopes(theElement->GetZ())<<G4endl;
      for(G4int i1=0; 
        i1<theStableOnes.GetNumberOfIsotopes(static_cast<G4int>(theElement->GetZ()) );
        i1++)
      {
//        G4cout <<" Init: theStableOnes in the loop"<<G4endl;
        G4int A = theStableOnes.GetIsotopeNucleonCount(first+i1);
        G4double frac = theStableOnes.GetAbundance(first+i1);
//        G4cout <<" Init: theStableOnes in the loop: "<<A<<G4endl;
        UpdateData(A, Z, count++, frac);
      }
    }
    theElasticData->ThinOut(precision);
    theInelasticData->ThinOut(precision);
    theCaptureData->ThinOut(precision);
    theFissionData->ThinOut(precision);
  }
  
  //void G4NeutronHPElementData::UpdateData(G4int A, G4int Z, G4int index, G4double abundance)
  void G4NeutronHPElementData::UpdateData(G4int A, G4int Z, G4int M, G4int index, G4double abundance)
  {
    //Reads in the Data, using G4NeutronHPIsoData[], and its Init
//    G4cout << "entered: ElementWiseData::UpdateData"<<G4endl;
    //theIsotopeWiseData[index].Init(A, Z, abundance);
    theIsotopeWiseData[index].Init(A, Z, M, abundance);
//    G4cout << "ElementWiseData::UpdateData Init finished"<<G4endl;

    theBuffer = theIsotopeWiseData[index].MakeElasticData();
//    G4cout << "ElementWiseData::UpdateData MakeElasticData finished: "
//         <<theBuffer->GetVectorLength()<<G4endl;
    Harmonise(theElasticData, theBuffer);
//    G4cout << "ElementWiseData::UpdateData Harmonise finished: "
//         <<theElasticData->GetVectorLength()<<G4endl;
    delete theBuffer;
    
    theBuffer = theIsotopeWiseData[index].MakeInelasticData();
//    G4cout << "ElementWiseData::UpdateData MakeInelasticData finished: "
//         <<theBuffer->GetVectorLength()<<G4endl;
    Harmonise(theInelasticData, theBuffer);
//    G4cout << "ElementWiseData::UpdateData Harmonise finished: "
//         <<theInelasticData->GetVectorLength()<<G4endl;
    delete theBuffer;
    
    theBuffer = theIsotopeWiseData[index].MakeCaptureData();
//    G4cout << "ElementWiseData::UpdateData MakeCaptureData finished: "
//         <<theBuffer->GetVectorLength()<<G4endl;
    Harmonise(theCaptureData, theBuffer);
//    G4cout << "ElementWiseData::UpdateData Harmonise finished: "
//         <<theCaptureData->GetVectorLength()<<G4endl;
    delete theBuffer;
    
    theBuffer = theIsotopeWiseData[index].MakeFissionData();
//    G4cout << "ElementWiseData::UpdateData MakeFissionData finished: "
//         <<theBuffer->GetVectorLength()<<G4endl;
    Harmonise(theFissionData, theBuffer);
//    G4cout << "ElementWiseData::UpdateData Harmonise finished: "
//         <<theFissionData->GetVectorLength()<<G4endl;
    delete theBuffer;
    
//    G4cout << "ElementWiseData::UpdateData finished"<endl;
  }
  
  void G4NeutronHPElementData::Harmonise(G4NeutronHPVector *& theStore, G4NeutronHPVector * theNew)
  {
    if(theNew == 0) { return; }
    G4int s_tmp = 0, n=0, m_tmp=0;
    G4NeutronHPVector * theMerge = new G4NeutronHPVector(theStore->GetVectorLength());
//    G4cout << "Harmonise 1: "<<theStore->GetEnergy(s)<<" "<<theNew->GetEnergy(0)<<G4endl;
    while ( theStore->GetEnergy(s_tmp)<theNew->GetEnergy(0)&&s_tmp<theStore->GetVectorLength() )
    {
      theMerge->SetData(m_tmp++, theStore->GetEnergy(s_tmp), theStore->GetXsec(s_tmp));
      s_tmp++;
    }
    G4NeutronHPVector *active = theStore;
    G4NeutronHPVector * passive = theNew;
    G4NeutronHPVector * tmp;
    G4int a = s_tmp, p = n, t;
//    G4cout << "Harmonise 2: "<<active->GetVectorLength()<<" "<<passive->GetVectorLength()<<G4endl;
    while (a<active->GetVectorLength()&&p<passive->GetVectorLength())
    {
      if(active->GetEnergy(a) <= passive->GetEnergy(p))
      {
        theMerge->SetData(m_tmp, active->GetEnergy(a), active->GetXsec(a));
        G4double x  = theMerge->GetEnergy(m_tmp);
        G4double y = std::max(0., passive->GetXsec(x)); 
        theMerge->SetData(m_tmp, x, theMerge->GetXsec(m_tmp)+y);
        m_tmp++;
        a++;
      } else {
//        G4cout << "swapping in Harmonise"<<G4endl;
        tmp = active; t=a;
        active = passive; a=p;
        passive = tmp; p=t;
      }
    }
//    G4cout << "Harmonise 3: "<< a <<" "<<active->GetVectorLength()<<" "<<m<<G4endl;
    while (a!=active->GetVectorLength())
    {
      theMerge->SetData(m_tmp++, active->GetEnergy(a), active->GetXsec(a));
      a++;
    }
//    G4cout << "Harmonise 4: "<< p <<" "<<passive->GetVectorLength()<<" "<<m<<G4endl;
    while (p!=passive->GetVectorLength())
    {
      // Modified by T. KOI
      //theMerge->SetData(m++, passive->GetEnergy(p), passive->GetXsec(p));
      G4double x = passive->GetEnergy(p);
      G4double y = std::max(0., active->GetXsec(x));
      theMerge->SetData(m_tmp++, x, passive->GetXsec(p)+y);
      p++;
    }
//    G4cout <<"Harmonise 5: "<< theMerge->GetVectorLength() << " " << m << G4endl;
    delete theStore;
    theStore = theMerge;
//    G4cout <<"Harmonise 6: "<< theStore->GetVectorLength() << " " << m << G4endl;
  }

  G4NeutronHPVector * G4NeutronHPElementData::MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPFissionData* theSet)
  {
   if(theP != G4Neutron::Neutron()) throw G4HadronicException(__FILE__, __LINE__, "not a neutron");
   Init ( theElement );
   return GetData(theSet);
  }
  G4NeutronHPVector * G4NeutronHPElementData::MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPCaptureData * theSet)
  {
   if(theP != G4Neutron::Neutron()) throw G4HadronicException(__FILE__, __LINE__, "not a neutron");
   Init ( theElement );
   return GetData(theSet);
  }
  G4NeutronHPVector * G4NeutronHPElementData::MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPElasticData * theSet)
  {
   if(theP != G4Neutron::Neutron()) throw G4HadronicException(__FILE__, __LINE__, "not a neutron");
   Init ( theElement );
   return GetData(theSet);
  }
    G4NeutronHPVector * G4NeutronHPElementData::MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPInelasticData * theSet)
  {
   if(theP != G4Neutron::Neutron()) throw G4HadronicException(__FILE__, __LINE__, "not a neutron");
   Init ( theElement );
   return GetData(theSet);
  }
