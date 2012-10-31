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
//
#include "BeamTestEventAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

////////////////////////////////////////////////////////////////////////

BeamTestEventAction::BeamTestEventAction(/*Parameters* parameter,*/ BeamTestRunAction* run) 
: runAct(run),fHitsCollectionID(0), fHitsCollectionID_monitor(0),p(0),pT(0),angle(0),pz(0), I(0),E(0),PD(0),b(0),
 numberOfChambers(0),abortEvent(false)
{ }

void BeamTestEventAction::Initialize( G4ThreeVector aP, G4ThreeVector aPosition ) {
    p = aP.mag();
    pT = aP.perp();
    I = aPosition;
    //BeamTestConversion A(ap,apt);
    //p = A.getP();
    //pT = A.getPT();
    //angle = xzAngle(A);
    //pz = pZ(A);
    //G4double pz = pZ(A);
    //particleGun->SetParticleEnergy(10.*GeV);
    //G4double dx = 6.*mm;
    //G4double dz = dx/tan(angle);
    //G4cout << "pT " << pT << "  and tan(theta_xz) = " << tan(xzAngle(A)) << G4endl;
    //G4ThreeVector Position(-dx,0.,-dz);
    //I = Position;	
}

BeamTestEventAction::~BeamTestEventAction() 
{
	//delete f;
}

void BeamTestEventAction::BeginOfEventAction(const G4Event*) 
{
	G4SDManager * SDman = G4SDManager::GetSDMpointer();
	////////////////////////////////////////////////////////////////////////
	// Getting code for HitsCollection of Silicon Monitor
	fHitsCollectionID_monitor = SDman->GetCollectionID("MonitorCollection");

	//opens existing file for updating if no file exists it creates it
	//f = new TFile("ntuple.root","UPDATE");
    abortEvent = false;
}


void BeamTestEventAction::EndOfEventAction(const G4Event* event)
{
    if ( abortEvent ) {
        //G4cout<<"Event to be aborted..."<<G4endl;
        return;
    }
    //G4cout<<"Event";

	G4HCofThisEvent *hitsCollectionOfThisEvent = event->GetHCofThisEvent();
	G4ThreeVector M;
	G4ThreeVector P;
	// Too many variables here but...ipT is the main one used for the inverseof the transverse momentum 
	//G4double En(0.),x(0.),y(0.),z(0.),px(0.),py(0.),pz(0.);//,ipT(0.);
	
	////////////////////////////////////////////////////////////////////////
	// Output code of Silicon Monitor Hits
	if ( fHitsCollectionID_monitor >= 0 )
	{
		BeamTestSiliconMonitorHitsCollection* hitsCollection_monitor = 
		dynamic_cast<BeamTestSiliconMonitorHitsCollection*>(hitsCollectionOfThisEvent->GetHC(fHitsCollectionID_monitor));
        assert(hitsCollection_monitor!=0);
		BeamTestSiliconMonitorHit* aHit = 0;
		//tree->Branch("aHit", "aHit", &aHit, bsize, split );

		G4int numberOfHits = hitsCollection_monitor->GetSize();
        //G4cout<<numberOfHits<<G4endl;
        G4double totE = 0; //Accumulate total energy released by Primary
        G4double totL = 0; //Accumulate total length by Primary
        for ( G4int idx = 0 ; idx < numberOfHits ; ++idx ) 
            {
                aHit = (*hitsCollection_monitor)[idx]; 
                //Check that hit is produced by the primary
                totE += aHit->GetEDep();
                totL += aHit->GetStepLength();
                if ( aHit->GetTrackId() == 1 ) 
                {
                    //Check that hit is placed in last chamber,
                    //and that PostStep was limited by boundary
                    if (aHit->GetChamberNumber() == numberOfChambers-1 &&
                        aHit->GetStepStatus() == fGeomBoundary ) 
                    {
                        /*En = aHit->GetExitKineticEnergy()/MeV;
                        P = aHit->GetExitMomentum()/mm ;
                        M = aHit->GetExitMomentumDirection();
                        x = P.x();
                        y = P.y();
                        z = P.z();
                        px = M.x();
                        py = M.y();
                        pz = M.z();*/
                
                        E = aHit->GetExitPosition()/mm;
                        //G4cout << "Position of Particle Gun: " << I << G4endl;
                        //G4cout << "Exiting Position: " << E << G4endl;
                    
                        // We want the vector tracking back to the PV so ive direction
                        PD = -1.0*aHit->GetExitMomentumDirection();
                        //G4cout << "Exiting Direction: " << PD << G4endl;
                        G4ThreeVector diff = I-E;
                        //G4ThreeVector diff = E-I;
                        //G4cout << "DIFF " << diff << G4endl;
                    
                        G4double t = diff.dot(PD);
                        //G4cout << "dot product  " << t << G4endl;
                    
                        G4ThreeVector X = E + t*PD;
                        //G4cout << "X " << X << G4endl;
                    
                        // check as this could be I-X in this backward case but want to got from IX^-> vector hence should be ok
                        b = X-I;
                        //G4cout << "AND FINALLY B: " << b << G4endl;
                    
                    
                    
                        //G4cout << "Information of " << i+1 << " th Silicon Monitor Hit of this event." << G4endl; 
                    
                        //G4cout << "Exiting Particle Name and Kinetic Energy " << aHit->GetExitDefinition()->GetParticleName() 	<< " " << E << " MeV" << G4endl; 
                        //G4cout << "Exiting position in silicon monitor "<< P << " in mm" << G4endl; 
                        //G4cout << "Exiting Momentum  " << M << G4endl;
                    
                        //G4double theta_x = ScatteredAngle(x, z);
                        //G4double theta_y = ScatteredAngle(y, z);
                    
                        //G4cout << "SCATTERED ANGLE IN X IS: " << theta_x << G4endl;
                        //G4cout << "SCATTERED ANGLE IN Y IS: " << theta_y << G4endl;
                        //G4double dz = 10.*cm;
                        //G4double dx = dz*tan(angle);
                        //G4double dc = sqrt(dx*dx +dz*dz);
                        // Due to rotation about the y-axis to get intial pT we have an extra angle
                        //G4double phi_x = ScatteredAngle(dx, dz);
                    
                        //impactParameter_x = sin(phi_x - theta_x)*dc;
                        //impactParameter_x = tan(phi_x - theta_x)*dc;
                        //impactParameter_y = sin(theta_y)*dc;
                        //impactParameter_z = impactParameter_x*tan(phi_x - theta_x);
                    
                        runAct->FillEvents(b.x()/CLHEP::mm,b.y()/CLHEP::mm,b.z()/CLHEP::mm, b.mag()/mm);
                    
                    } //endif Hit is in last chamber and boundary limited
                } //endif Hit is from primary
        }//End loop on hits
        static G4double density = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al")->GetDensity();
        G4double dedx = (totE/totL/density);
        runAct->SetEnergyDeposit(dedx/(MeV*cm2/g));
        
    }
}

G4double BeamTestEventAction::ScatteredAngle(G4double x_i, G4double z_FinalExitPos)
{
	return std::atan(x_i/z_FinalExitPos);
}

std::string makeFilename( const std::string& basename, double index )
{
	std::ostringstream result;
	result << basename << index << ".dat";
	return result.str();
}
