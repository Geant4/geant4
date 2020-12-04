// -------------------------------------------------------------------
// $Id: plot.C 
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT v6.xx (usually type 'root' at your machine's prompt)
//   2 - type '.X analysis.C' at the ROOT session prompt
// *********************************************************************

{
    gROOT->Reset();
    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");
    TRandom rdomN;
           
    vector<Molecule> fMolecules = molecule();
    using Table = map<int,array<map<int,double>, 2>>;
    using ResultTable = map<int,array<vector<int>, 2>>;

    system ("rm -rf output.root");
    system ("hadd output.root output_*.root");
    
    double xphy;
    double yphy;
    double zphy;
    
    double xChe;
    double yChe;
    double zChe;

    double xrad;
    double yrad;
    double zrad;

    map<int,int> CopyNumTable;
    Table physTable; 
    ResultTable resultPhysTable;
    ResultTable resultChemTable;
    ResultTable MergedResultTable;

    double EnergyDeposit;
    double kineticEnergyDifference;
    int flagVolume;
    double copyN;
    int EventIDp;
    int EventIDc;

    TFile f("output.root"); 

// Load the file, the directory and the ntuple
    TDirectoryFile* d = dynamic_cast<TDirectoryFile*>(f.Get("ntuple") );   
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
        
        if(physTable[EventIDp][strand][copyN] < 17.5)
        {
            continue;
        }
        
        if(CopyNumTable.empty())
        {
            CopyNumTable.insert(pair<int,int>(copyN,strand));
            resultPhysTable[EventIDp][strand].push_back(copyN);
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
            resultPhysTable[EventIDp][strand].push_back(copyN);
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
                
                if(rdomN.Rndm() > 0.42)//42% of reaction make damages
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
                    resultChemTable[EventIDc][strand].push_back(
                    moleculeIt.fCopyNumber);
                }
            }
        }
    }
//SDD Format
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
            <<"Incident particles, 200;\n"
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
            <<"Cell Cycle Phase, G1/G2;\n"
            <<"DNA Structure, 4, 1;\n"
            <<"In vitro / in vivo, 0;\n"
            <<"Proliferation status, 1;\n"
            <<"Microenvironment, 27, 0.0;\n"
            <<"Damage definition, 1, 1, 10, -1, 0, 17.5;\n"
            <<"Time, 2.5ns;\n"
            <<"Damage and primary count, 638, 200;\n"
            <<"Data entries, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;\n"
            <<"Data was generated in minimal output format and used as"
              " an example. Please modify the primary particle number"
              "for this file;\n"
            <<"\n"
            <<"***EndOfHeader***;\n"
            <<endl;
//For detail, please refer to Schuemann, J., et al. (2019). A New Standard DNA 
//Damage (SDD) Data Format. Radiation Research, 191(1), 76. 

        bool Primaryflag = true;
    for (int i = 0; i < 200; i++)//for 200 primary particles
    {
        for(int j = 0; j < 2; j++)//for strand
        {
            const std::string positionAndCopyN = "; 0, 0, 0; 3, 0, 0, ";
            const std::string SSBType = " 0, 1, 0";
            if (!resultPhysTable[i][j].empty())
            {
                for(int e = 0; e < resultPhysTable[i][j].size(); ++e)
                {
                    if(Primaryflag)
                    {
                        ofile<<"2, "<<i<<positionAndCopyN
                             <<resultPhysTable[i][j][e]
                             <<"; 0;"<<SSBType<<"\n" ;//physics
                        Primaryflag = false;
                    }else
                    {
                        ofile<<"0, "<<i<<positionAndCopyN
                            <<resultPhysTable[i][j][e]
                            <<"; 0;"<<SSBType<<"\n" ;
                    }
                }
            }
            if(resultChemTable[i][j].empty())
            {
                continue;
            }
            
            for(int ee = 0; ee < resultChemTable[i][j].size(); ++ee)
            {
                if(Primaryflag)
                {
                    ofile<<"2, "<<i<<positionAndCopyN
                    <<resultChemTable[i][j][ee]
                    <<"; 1;"<<SSBType<<"\n" ;//chemistry
                    Primaryflag = false;
                }else
                {
                    ofile<<"0, "<<i<<positionAndCopyN
                    <<resultChemTable[i][j][ee]
                    <<"; 1;"<<SSBType<<"\n" ;
                }
            }
        }
    }
    
    ofile<<endl;	
    }
    
    ofile.close();
}
