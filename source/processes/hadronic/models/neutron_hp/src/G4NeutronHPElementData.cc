// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPElementData.hh"

  G4NeutronHPElementData::G4NeutronHPElementData()
  {
     precision = 0.02;
     theFissionData = new G4NeutronHPVector;
     theCaptureData = new G4NeutronHPVector;
     theElasticData = new G4NeutronHPVector;
     theInelasticData = new G4NeutronHPVector;
    theIsotopeWiseData = NULL;
  }
  
  G4NeutronHPElementData::~G4NeutronHPElementData()
  {
    if(theFissionData!=NULL) delete theFissionData;
    if(theCaptureData!=NULL) delete theCaptureData;
    if(theElasticData!=NULL) delete theElasticData;
    if(theInelasticData!=NULL) delete theInelasticData;
    if(theIsotopeWiseData!=NULL) delete [] theIsotopeWiseData;
  }
  
  void G4NeutronHPElementData::Init(G4Element * theElement)  
  {
    G4int count = theElement->GetNumberOfIsotopes();
      if(count == 0) count +=
         theStableOnes.GetNumberOfIsotopes(theElement->GetZ());
    theIsotopeWiseData = new G4NeutronHPIsoData[count];
    // filename = ein data-set je isotope.
    count = 0;
    G4int nIso = theElement->GetNumberOfIsotopes();
    G4int Z = theElement->GetZ();
    G4int i1;
    if(nIso!=0)
    {
      for (i1=0; i1<nIso; i1++)
      {
//        G4cout <<" Init: normal case"<<G4endl;
        G4int A = theElement->GetIsotope(i1)->GetN();
        G4double frac = theElement->GetRelativeAbundanceVector()[i1]/perCent;
        UpdateData(A, Z, count++, frac);
      }
    }else{
//      G4cout <<" Init: theStableOnes case: Z="<<Z<<G4endl;
      G4int first = theStableOnes.GetFirstIsotope(Z);
//      G4cout <<"first="<<first<<" "<<theStableOnes.GetNumberOfIsotopes(theElement->GetZ())<<G4endl;
      for(G4int i1=0; 
        i1<theStableOnes.GetNumberOfIsotopes(theElement->GetZ());
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
  
  void G4NeutronHPElementData::UpdateData(G4int A, G4int Z, G4int index, G4double abundance)
  {
    //Reads in the Data, using G4NeutronHPIsoData[], and its Init
//    G4cout << "entered: ElementWiseData::UpdateData"<<G4endl;
    theIsotopeWiseData[index].Init(A, Z, abundance);
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
    if(theNew == NULL) return;
    G4int s = 0, n=0, i=0, m=0;
    G4NeutronHPVector * theMerge = new G4NeutronHPVector(theStore->GetVectorLength());
    G4bool flag;
//    G4cout << "Harmonise 1: "<<theStore->GetEnergy(s)<<" "<<theNew->GetEnergy(0)<<G4endl;
    while ( theStore->GetEnergy(s)<theNew->GetEnergy(0)&&s<theStore->GetVectorLength() )
    {
      theMerge->SetData(m++, theStore->GetEnergy(s), theStore->GetXsec(s));
      s++;
    }
    G4NeutronHPVector *active = theStore;
    G4NeutronHPVector * passive = theNew;
    G4NeutronHPVector * tmp;
    G4int a = s, p = n, t;
//    G4cout << "Harmonise 2: "<<active->GetVectorLength()<<" "<<passive->GetVectorLength()<<G4endl;
    while (a<active->GetVectorLength()&&p<passive->GetVectorLength())
    {
      if(active->GetEnergy(a) <= passive->GetEnergy(p))
      {
        theMerge->SetData(m, active->GetEnergy(a), active->GetXsec(a));
        G4double x  = theMerge->GetEnergy(m);
        G4double y = passive->GetXsec(x); 
        theMerge->SetData(m, x, theMerge->GetXsec(m)+y);
        m++;
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
      theMerge->SetData(m++, active->GetEnergy(a), active->GetXsec(a));
      a++;
    }
//    G4cout << "Harmonise 4: "<< p <<" "<<passive->GetVectorLength()<<" "<<m<<G4endl;
    while (p!=passive->GetVectorLength())
    {
      theMerge->SetData(m++, passive->GetEnergy(p), passive->GetXsec(p));
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
   if(theP != G4Neutron::Neutron()) G4Exception();
   Init ( theElement );
   return GetData(theSet);
  }
  G4NeutronHPVector * G4NeutronHPElementData::MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPCaptureData * theSet)
  {
   if(theP != G4Neutron::Neutron()) G4Exception();
   Init ( theElement );
   return GetData(theSet);
  }
  G4NeutronHPVector * G4NeutronHPElementData::MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPElasticData * theSet)
  {
   if(theP != G4Neutron::Neutron()) G4Exception();
   Init ( theElement );
   return GetData(theSet);
  }
    G4NeutronHPVector * G4NeutronHPElementData::MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPInelasticData * theSet)
  {
   if(theP != G4Neutron::Neutron()) G4Exception();
   Init ( theElement );
   return GetData(theSet);
  }
