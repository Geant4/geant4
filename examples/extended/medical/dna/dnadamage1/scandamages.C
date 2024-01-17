// -------------------------------------------------------------------
// $Id: scandamages.C 
// -------------------------------------------------------------------
// 13/10/2023 - Le Tuan Anh renamed and improved
// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT v6.xx (usually type 'root' at your machine's prompt)
//   2 - type '.X analysis.C' at the ROOT session prompt
// *********************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

struct Damage;
struct Cluster;
void CreateSDDfile(std::map<int,std::vector<Damage>>,long int);
void ExtractDirectAndIndirectDamages(std::map<int,std::vector<Damage>>&,long int&,double,double);
void ClassifyDamages(std::map<int,std::vector<Damage>>&,int, 
                    std::map<int,vector<Cluster>>&,int verbose =0);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// play with scandamages() by changing firect/indirect criterion selection or Clusterdistance

void scandamages()
{
    ifstream fs("output_t0.root");
    if ( fs.good() ) {
      system ("rm -rf output.root");
      system ("hadd output.root output_*.root");
    }
    
    double fEnergyThresholdForDirectSelection = 17.5; // eV
    double fProbabilityForIndirectionSelection = 0.40; // 40 %
    int fClusterdistance = 10; // bp, minium distance that separate two clusters

    std::map<int,std::vector<Damage>> fDamagesEventMap;
    std::map<int,vector<Cluster>> fClustersMap;
    long int totalDamages = 0;

    ExtractDirectAndIndirectDamages(fDamagesEventMap,totalDamages,
                                    fEnergyThresholdForDirectSelection, 
                                    fProbabilityForIndirectionSelection);

    ClassifyDamages(fDamagesEventMap,fClusterdistance,fClustersMap);
    // print results in each event for checking:
    //ClassifyDamages(fDamagesEventMap,fClusterdistance,fClustersMap,1);
    
    CreateSDDfile(fDamagesEventMap,totalDamages);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

struct Damage
{   
    Damage(int a, int b, int c=-1): strand(a),copyNb(b), cause(c){}
    int strand;
    int copyNb; 
    int cause; // Unknown = -1, Direct = 0, Indirect = 1
    bool operator<(const Damage& r) const {return copyNb<r.copyNb;} // for sorting
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

struct Cluster // to score SSB or DSB
{   
    std::vector<Damage> fDamagesVector;
    bool HasLeftStrand=false, HasRightStrand=false;
    int firstCopyNb;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
ClassifyDamages(std::map<int,std::vector<Damage>> &fDamages,int fClusterditance, 
                std::map<int,vector<Cluster>>& fClustersMap,int verbose)
{
    long int dirSB = 0, indirSB = 0;
    long int nDSB = 0, nsDSB = 0, ncDSB=0, nSSB = 0;
    std::cout<<"-----------------------------------------------"<<std::endl;
    if (verbose>0)std::cout<<" ------> Print few events for checking <------"<<std::endl;
    if (verbose>0)std::cout<<"\t\tCopyNb\tStrand\tCause"<<std::endl;
    for (auto &[eventID, dmgVector] : fDamages) {
        std::sort(dmgVector.begin(),dmgVector.end());
        if (verbose>0) std::cout<<"Event #"<<eventID<<":  "<<std::endl; 
        unsigned long fSBcountThisEvent = 0;
        long copyNBofPrevDamage = 0;
        long firstCopyNb;
        std::map<int,Cluster> clustersInthisEvent;
        for (auto dmg : dmgVector) {
            if (dmg.cause == 0) dirSB++;
            if (dmg.cause == 1) indirSB++;
            if (verbose>0) std::cout<<"\t\t"<<dmg.copyNb<<"\t"
                                    <<dmg.strand<<"\t"<<dmg.cause<<"\t"<<std::endl;
            if (clustersInthisEvent.size()==0 || 
                ((dmg.copyNb - copyNBofPrevDamage) >=  fClusterditance)) {
                Cluster newCluster;
                newCluster.fDamagesVector.push_back(dmg);
                copyNBofPrevDamage = dmg.copyNb;
                firstCopyNb = dmg.copyNb;
                if (dmg.strand == 0) newCluster.HasLeftStrand = true;
                if (dmg.strand == 1) newCluster.HasRightStrand = true;
                clustersInthisEvent.insert({firstCopyNb,newCluster});
            } else {
                clustersInthisEvent[firstCopyNb].fDamagesVector.push_back(dmg);
                copyNBofPrevDamage = dmg.copyNb;
                if (dmg.strand == 0 ) clustersInthisEvent[firstCopyNb].HasLeftStrand = true;
                if (dmg.strand == 1 ) clustersInthisEvent[firstCopyNb].HasRightStrand = true;
            }
        }

        for (auto [firstCpNb,acluster] : clustersInthisEvent) {
            if (fClustersMap.find(eventID) == fClustersMap.end()) {
                std::vector<Cluster> sclv{acluster};
                fClustersMap.insert({eventID,sclv});
            } else {
                fClustersMap[eventID].push_back(acluster);
            }
        }
        clustersInthisEvent.clear();
    }
    if (verbose>0)std::cout<<" \n------> Checking SSB and DSB:"<<std::endl;
    if (verbose>0)std::cout<<"\t\t#SSB\t#DSB\t#cDSB\t#sDSB"<<std::endl;

    for (auto [evt,clustervector] : fClustersMap) {
        long int nDSBinThisEvt = 0, nSSBinThisEvt = 0;
        long int nsDSBinThisEvt = 0, ncDSBinThisEvt = 0;
        if (verbose>0) std::cout<<"Event #"<<evt<<":  "<<std::endl; // maximum 5 event for cheking
        for (auto sb : clustervector) {
            if (sb.HasLeftStrand && sb.HasRightStrand) {
                nDSBinThisEvt++;
                if (sb.fDamagesVector.size()>2) ncDSBinThisEvt++;
                else nsDSBinThisEvt++;
            }
            else nSSBinThisEvt++;
        }
        nSSB += nSSBinThisEvt;
        nDSB += nDSBinThisEvt;
        ncDSB += ncDSBinThisEvt;
        nsDSB += nsDSBinThisEvt;
        if (verbose>0) std::cout<<"\t\t"<<nSSBinThisEvt<<"\t"<<nDSBinThisEvt
                                <<"\t"<<ncDSBinThisEvt<<"\t"<<nsDSBinThisEvt<<std::endl;
    }
    std::cout<<"-----------------------------------------------"<<std::endl;
    std::cout<<"\tSummary of results:"  
                <<"\n #Total SBs = "<<(dirSB+indirSB)
                <<"\n #dirSB     = "<<dirSB <<" \t;\t #indirSB = "<<indirSB
                <<"\n #SSB       = "<<nSSB  <<" \t;\t #DSB     = "<<nDSB
                <<"\n #sDSB      = "<<nsDSB <<" \t;\t #cDSB    = "<<ncDSB<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExtractDirectAndIndirectDamages(std::map<int,std::vector<Damage>> &fDamagesEventMap,
long int &totalDamages,double Ethresh,double ChemPro)
{
    TRandom rdomN;
           
    vector<Molecule> fMolecules = molecule();
    using Table = map<int,array<map<int,double>, 2>>;
    double xphy, yphy, zphy;
    double xChe, yChe, zChe;
    double xrad, yrad, zrad;

    map<int,int> CopyNumTable;
    Table physTable; 

    double EnergyDeposit;
    double kineticEnergyDifference;
    int flagVolume;
    double copyN;
    int EventIDp;
    int EventIDc;

    TFile f("output.root"); 

// Load the file, the directory and the ntuple
    TDirectoryFile* d = dynamic_cast<TDirectoryFile*>(f.Get("ntuple") );   
//analyse Physical stage 
    TTree* tPhys = dynamic_cast<TTree*> (d->Get("ntuple_1") );
    tPhys->SetBranchAddress("x",&xphy);
    tPhys->SetBranchAddress("y",&yphy);
    tPhys->SetBranchAddress("z",&zphy);
    tPhys->SetBranchAddress("edep",&EnergyDeposit);
    tPhys->SetBranchAddress("diffKin",&kineticEnergyDifference);
    tPhys->SetBranchAddress("volumeName",&flagVolume);
    tPhys->SetBranchAddress("CopyNumber",&copyN);
    tPhys->SetBranchAddress("EventID",&EventIDp);
    
    auto entryPhyNumber = tPhys->GetEntries();
    
    bool addCopyN = true;
    for(decltype(entryPhyNumber) i = 0; i < entryPhyNumber; ++i)
    {
        tPhys->GetEntry(i);
//avoid histone
        if(flagVolume == DNAVolumeType::histone)
        {
            continue;
        }
//Determine the strand
        int strand = -1;
        bool isLeftStrand = DNAVolumeType::deoxyribose1 == flagVolume || 
                            DNAVolumeType::phosphate1 == flagVolume || 
                            DNAVolumeType::deoxyribose1_water == flagVolume || 
                            DNAVolumeType::phosphate1_water == flagVolume;
                            
        bool isRightStrand = DNAVolumeType::deoxyribose2 == flagVolume || 
                             DNAVolumeType::phosphate2 == flagVolume || 
                             DNAVolumeType::deoxyribose2_water== flagVolume || 
                             DNAVolumeType::phosphate2_water == flagVolume;
                            
            
        if( isLeftStrand )
        {   
            strand = 0;
        }
        else if( isRightStrand )
        {
            strand = 1;
        }
        else 
        {
            //in case molecules are not assigned "strand" 
            continue;
        }

//Determine the name
        bool isDeoxyribode = flagVolume == DNAVolumeType::deoxyribose1 || 
                             flagVolume == DNAVolumeType::deoxyribose2 || 
                             flagVolume == DNAVolumeType::deoxyribose1_water || 
                             flagVolume == DNAVolumeType::deoxyribose2_water;
                             
        bool  isPhosphate =  flagVolume == DNAVolumeType::phosphate1 || 
                             flagVolume == DNAVolumeType::phosphate2 || 
                             flagVolume == DNAVolumeType::phosphate1_water|| 
                             flagVolume == DNAVolumeType::phosphate2_water; 

        if(isDeoxyribode || isPhosphate)
        {
            physTable[EventIDp][strand][copyN] += EnergyDeposit;
        }
        
        if(physTable[EventIDp][strand][copyN] < Ethresh)
        {
            continue;
        }
        
        if(CopyNumTable.empty())
        {
            CopyNumTable.insert(pair<int,int>(copyN,strand));
            if (fDamagesEventMap.find(EventIDp) == fDamagesEventMap.end()) {
                std::vector<Damage> admg{Damage(strand,copyN,0)};
                fDamagesEventMap.insert({EventIDp,admg});
            } else {
                fDamagesEventMap[EventIDp].push_back(Damage(strand,copyN,0));
            }
            totalDamages++;
        }
        else                           
        {
            addCopyN = true;
        }
        
        auto itCopyNumTable = CopyNumTable.find(copyN);
        if (itCopyNumTable != CopyNumTable.end())
        {
            if (itCopyNumTable->second == strand)
            {
                addCopyN = false;
            }
        }
        
        if(addCopyN)
        {
            CopyNumTable.insert(pair<int,int>(copyN,strand));
            if (fDamagesEventMap.find(EventIDp) == fDamagesEventMap.end()) {
                std::vector<Damage> admg{Damage(strand,copyN,0)};
                fDamagesEventMap.insert({EventIDp,admg});
            } else {
                fDamagesEventMap[EventIDp].push_back(Damage(strand,copyN,0));
            }
            totalDamages++;
        }
    }    
//Chemistry analyse
    TTree* tChem = dynamic_cast<TTree*> (d->Get("ntuple_2") );
    tChem->SetBranchAddress("x",&xrad);
    tChem->SetBranchAddress("y",&yrad);
    tChem->SetBranchAddress("z",&zrad);
    tChem->SetBranchAddress("EventID",&EventIDc);
    
    auto entryNChem = tChem->GetEntries() ;
    
    for(decltype(entryNChem) i=0; i < entryNChem; ++i)
    {
        tChem->GetEntry(i);
        ThreeVector<double> DNAElement(xrad,yrad,zrad);
        
        for(const auto& moleculeIt : fMolecules)
        {
            if(moleculeIt.fPosition == DNAElement)
            {
                int strand = -1;
                if(moleculeIt.fName == "deoxyribose1") 
                {
                    strand = 0;
                }
                else if(moleculeIt.fName == "deoxyribose2") 
                {
                    strand = 1;
                }
                else 
                {   
                    string msg = "Unknown DNA component";
                    throw std::runtime_error(msg);
                }
                
                if(rdomN.Rndm() > ChemPro)//probability of reaction make damages
                {
                    continue;
                }
                
                if (!CopyNumTable.empty())
                {
                    addCopyN = true;
                }

                auto itCopyNumTable = CopyNumTable.find(moleculeIt.fCopyNumber);
                
                if (itCopyNumTable != CopyNumTable.end())
                {
                    if (itCopyNumTable->second == strand)
                    {
                        addCopyN = false;
                    }
                }
                if(addCopyN)
                {
                    if (fDamagesEventMap.find(EventIDc) == fDamagesEventMap.end()) {
                        std::vector<Damage> admg{Damage(strand,moleculeIt.fCopyNumber,1)};
                        fDamagesEventMap.insert({EventIDc,admg});
                    } else {
                        fDamagesEventMap[EventIDc].push_back(
                                                Damage(strand,moleculeIt.fCopyNumber,1));
                    }
                    totalDamages++;
                }
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CreateSDDfile(std::map<int,std::vector<Damage>> fDamagesEventMap,long int totalDamages)
{
    ofstream ofile("SDDFormat.txt");
    if(ofile.is_open())
    {
        ofile<<'\n'
            <<"SDD version, SDDv1.0;\n"
            <<"Software, chromatin fiber Model;\n"
            <<"Author contact, Carmen Villagrasa,carmen.villagrasa@irsn.fr, "
              "02/04/2018, Melyn et al.,Sc. Rep. 7 (2017)11923;\n"
            <<"Simulation Details, DNA damages from direct and "
              "indirect effects;\n"
            <<"Source, Monoenergetic cylinder-parallel proton beam uniformly " 
              "in 5 nm radius from left side of the cube"
            <<" exposing chromatin fiber. Energy: 5 keV;\n"
            <<"Source type, 1;\n"
            <<"Incident particles, 2212;\n"
            <<"Mean particle energy, 5 keV;\n"
            <<"Energy distribution, M, 0;\n"
            <<"Particle fraction, 1.0;\n"
            <<"Dose or fluence, 0, 0;\n"
            <<"Dose rate, 0.0;\n"
            <<"Irradiation target, Simple spherical chromatin fiber model"
              " in a voxel 40 nm;\n"
            <<"Volumes, 0,20,20,20,-20,-20,-20,2,15,15,20,-15,-15,20;\n"
            <<"Chromosome sizes, ;\n"
            <<"DNA Density, ;"
            <<"Cell Cycle Phase, G0/G1;\n"
            <<"DNA Structure, 4, 1;\n"
            <<"In vitro / in vivo, 0;\n"
            <<"Proliferation status, 1;\n"
            <<"Microenvironment, 27, 0.0;\n"
            <<"Damage definition, 1, 1, 10, -1, 17.5;\n"
            <<"Time, 2.5ns;\n"
            <<"Damage and primary count, "<<totalDamages<<", 0;\n"
            <<"Data entries, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;\n"
            <<"Data field explaination, Field 1: [2]-eventID, Field 3: [4]-strand, " 
            <<"Field 4:damage position (copynb), Field 5: Cause (direct: [1]=0) "
            <<" (indirect: [1]=1), Field 6: Damage types (Base:[1]>0) (Backbone: [2]>0);\n"
            <<"\n"
            <<"***EndOfHeader***;\n"
            <<endl;
//For detail, please refer to Schuemann, J., et al. (2019). A New Standard DNA 
//Damage (SDD) Data Format. Radiation Research, 191(1), 76. 

    bool Primaryflag = true;
    for (auto const [eventID,dmgvector] : fDamagesEventMap)
    {
        Primaryflag = true; // update Primaryflag for new primary
        for (auto const dmg:dmgvector) {
            //Field 1 Calassification
            int newEvtFlag = 0; // = 2 if new event;
            if (Primaryflag) {
                newEvtFlag = 2;
                Primaryflag = false;
            }
            ofile<<newEvtFlag<<", "<<eventID<<"; ";
            //Field 3 Chromosome IDs and strand
            ofile<<0<<", "<<0<<", "<<0<<", "<<dmg.strand<<"; ";
            //Field 4, Chromosome position 
            ofile<<dmg.copyNb<<"; ";
            //Field 5, Cause: Unknown = -1, Direct = 0, Indirect = 1
            ofile<<dmg.cause<<", "<<0<<", "<<0<<"; ";
            //Field 6, Damage types: only Backbone is considered here
            ofile<<0<<", "<<1<<", "<<0<<"; ";
            ofile<<"\n";
        }
    }
    
    ofile<<endl;	
    }
    
    ofile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....