// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronCaptureTest.cc,v 1.1 1999-01-08 16:33:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Test of Low-Energy Capture class
//
// F.W. Jones, TRIUMF, 03-DEC-96
//
// 27-MAR-97 F.W.Jones: first version for Alpha release
//


#include "G4ios.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4ProcessManager.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4LCapture.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include "G4PionPlus.hh"
#include "G4Neutron.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"

#ifdef G4_SOLVE_TEMPLATES
#include "g4templates.hh"
#endif

// The following can be used to produce a Fortran-readable file
// of the events:
// extern "C" fwrite_(float*, float*, float*, float*);
// where:
//      subroutine fwrite(ek,px,py,pz)
//      write(44)ek,px,py,pz
//      return
//      end

// Declarations for HBOOK (for future use?)
// extern "C" hlimit_(int*);
// extern "C" hropen_(int*, char*, char*, char*, int*, int*,
//                    int, int, int);
// extern "C" hbookn_(int*, char*, int*, char*, int*, char*,
//                   int, int, int);
// extern "C" hfn_(int*, float*);
// extern "C" hrout_(int*, int*, char*,
//                   int);
// extern "C" hrend_(char*,
//                  int);


int main()
{
// Material definition

// G4Material Uranium("uranium 238", 92, 238.03*g/mole, 18.95*g/cm3,
//                    kStateSolid);
   G4Material Uranium("uranium 235", 92, 235.*g/mole, 18.95*g/cm3,
                      kStateSolid);

//   G4cout << "Dumping material table info:" << endl;
//   Uranium.DumpInfo();

// Particle definition

   G4Neutron* neutron = G4Neutron::NeutronDefinition();
// neutron->DumpTable();
   G4ProcessManager* theNeutronProcessManager = neutron->GetProcessManager();

// Process definition

   G4HadronCaptureProcess theLCaptureProcess;
   G4LCapture the;
   theLCaptureProcess.RegisterMe(&the);
   //   theLCaptureProcess.SetVerboseLevel(2);
   //   theLCaptureProcess.setCrossSectionDataVerboseLevel(2);
   theNeutronProcessManager->AddDiscreteProcess(&theLCaptureProcess);
   theLCaptureProcess.BuildPhysicsTable(*neutron);
   theLCaptureProcess.DumpPhysicsTable(*neutron);

// Dynamic Particle Definition
   G4DynamicParticle aNeutron(neutron, G4ParticleMomentum(1.,0.,0.), 
                              1.*GeV);

// Physics

   G4double T;
   G4cout << "Enter neutron kinetic energy in GeV: ";
   cin >> T;
   aNeutron.SetKineticEnergy(T*GeV);
   G4double p = aNeutron.GetTotalMomentum()/GeV;
   G4cout << "neutron momentum=" << p << endl;
   G4double px, py, pz;
   G4cout << "Enter momentum direction px/p, py/p, pz/p: ";
   cin >> px >> py >> pz;
   aNeutron.SetMomentumDirection(G4ParticleMomentum(px, py, pz));
   aNeutron.DumpInfo();

   G4double meanFreePath;
   meanFreePath = theLCaptureProcess.GetMeanFreePathBasic(&aNeutron, &Uranium);
   G4cout << endl << "Neutron mean free path = " << meanFreePath << endl;

   G4int nevent;
   G4cout << endl << "Enter number of events: ";
   cin >> nevent;

   G4Step aStep;

   G4double truePathLength = 1.;
   G4VParticleChange* pch;
   G4Track* theSecondary;
   G4double eknew;
   G4double ekn, ekn2, ekg, ekg2;
   G4double pnew, pxnew, pynew, pznew;
   G4double pxg, pxg2, pyg, pyg2, pzg, pzg2;
   G4int ngnew;
   G4int ng, ng2;
   float eknew4, pxnew4, pynew4, pznew4;
   ekn = ekn2 = ekg = ekg2 = 0.;
   pxg = pxg2 = pyg = pyg2 = pzg = pzg2 = 0.;
   ng = ng2 = 0;

   if (nevent > 10) G4cout << "Printing first 10 events only" << endl;
   G4cout << endl << "Event  sec#  type  Ek(MeV)  Px/P  Py/P  Pz/P"
        << endl << endl;

// Open HBOOK RZ file (not commissioned)
   //   int NWPAW = 10000;
   //   hlimit_(&NWPAW);
   //   int LUN = 1;
   //   char* CHTOP = "RZTOP";
   //   char* CHFILE = "CAPTURE";
   //   char* CHOPT = "N";
   //   int LREC = 1024;
   //   int ISTAT;
   //   hropen_(&LUN, CHTOP, CHFILE, CHOPT, &LREC, &ISTAT,
   //           strlen(CHTOP), strlen(CHFILE), strlen(CHOPT));
   //   G4cout << "HROPEN ISTAT=" << ISTAT << endl;

// Book ntuple
   //   int ID = 10;
   //   char* CHTITL = "G4LCapture test";
   //   int NVAR = 4;
   //   char* CHRZPA = "RZTOP";
   //   int NWBUFF = NVAR*nevent;
   //   char* TAGS = "EkPxPyPz";
   //   hbookn_(&ID, CHTITL, &NVAR, CHRZPA, &NWBUFF, TAGS,
   //           strlen(CHTITL), strlen(CHRZPA), 2);
   //   float X[4];
   
   G4Box* theFrame = new G4Box ("Frame",10*m, 10*m, 10*m);
   G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
                                                       &Uranium,
                                                      "LFrame", 0, 0, 0);
   G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
                                                    "PFrame",LogicalFrame,0,false,0);
  G4ThreeVector aPosition(0.,0.,0.);
   G4double aTime = 0.0;
   G4Track* aTrack = new G4Track( &aNeutron, aTime, aPosition );
   aTrack->SetVolume( PhysicalFrame );
   aStep.SetTrack( aTrack );


   for (G4int ievent = 1; ievent <= nevent; ievent++) {
      pch = theLCaptureProcess.PostStepDoIt( *aTrack, aStep );
//      pch = theLCaptureProcess.DoIt(&aNeutron, &Uranium, aStep);
      ngnew = 0;
      for (G4int i = 1; i <= pch->GetNumberOfSecondaries(); i++) {
         theSecondary = pch->GetSecondary(i - 1);
         if (theSecondary->GetDefinition() == 
                  G4Gamma::GammaDefinition()) {
            ngnew++;
            eknew = theSecondary->GetKineticEnergy()/MeV;
            //      G4cout << "eknew=" << eknew << endl;
            ekg += eknew;
            ekg2 += eknew*eknew;
            pnew = theSecondary->GetDynamicParticle()->GetTotalMomentum()/MeV;
            pxnew = pnew*theSecondary->GetMomentumDirection().x();
            pxg += pxnew;
            pxg2 += pxnew*pxnew;
            pynew = pnew*theSecondary->GetMomentumDirection().y();
            pyg += pynew;
            pyg2 += pynew*pynew;
            pznew = pnew*theSecondary->GetMomentumDirection().z();
            pzg += pznew;
            pzg2 += pznew*pznew;
            // Fill ntuple
            //            X[0] = eknew;
            //            X[1] = pxnew;
            //            X[2] = pynew;
            //            X[3] = pznew;
            //            hfn_(&ID, X);
            eknew4 = eknew;
            //            G4cout << "eknew4=" << eknew4 << endl;
            pxnew4 = pxnew;
            pynew4 = pynew;
            pznew4 = pznew;
            //            fwrite_(&eknew4, &pxnew4, &pynew4, &pznew4);
         }
         if (ievent <= 10) {
            printf("%d  %d  %s  %f  %f  %f  %f\n",
                   ievent, i,
                   theSecondary->GetDefinition()->GetParticleName(),
                   eknew, pxnew, pynew, pznew);
         }
      }
      ng += ngnew;
      ng2 += ngnew*ngnew;
   }

// Close HBOOK RZ file
   //   int IDZERO = 0;
   //   static int ICYCLE;
   //   hrout_(&IDZERO, &ICYCLE, " ",
   //          1);
   //   hrend_(CHTOP,
   //          strlen(CHTOP));

   G4cout << endl << "GAMMA:" << endl;
   G4cout << "  Total number     " << ng << endl;
   G4cout << "  Number per event " << (G4double)ng/nevent << endl;
   G4cout << "  Ek " << ekg/ng << endl;
   G4cout << "  Px " << pxg/ng << endl;
   G4cout << "  Py " << pyg/ng << endl;
   G4cout << "  Pz " << pzg/ng << endl;

   //   return EXIT_SUCCESS;
}
