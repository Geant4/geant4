#include <string>
#include <sstream>
#include <iomanip>
#include "gzstream.h"

#include "TROOT.h" 
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "precoreaction.h"

#include <unistd.h>

double ToDouble(const std::string& s)
{
  std::istringstream tonum(s.c_str());
  double result;
  tonum >> result;
  return result;
}

int ToInt(const std::string& s)
{
  std::istringstream tonum(s.c_str());
  int result;
  tonum >> result;
  return result;
}

int ToInt(const char * t)
{
  std::istringstream tonum(t);
  int result;
  tonum >> result;
  return result;
}

void usage(const std::string& progname)
{
  std::cout << "Usage: " << progname << " filename>\n" 
            << "       " << progname << " <input filename> <output filename>\n"
            << "       Options:    -h        prints this help\n"
            << "                   -l num    conver only the first num events"
            << std::endl;
  return;
}


int main(int argc, char **argv)
{
  TROOT rconversion("rconversion","G4PreCompoundTest to root data conversion");

  std::string progname(argv[0]);
  int limit_event = -1;
  bool hflag = false;
  int st;
  
  while ((st = getopt(argc,argv,"hl:")) != -1)
    {
      switch(st)
        {
        case 'l':
          limit_event = ToInt(optarg);
          break;
        case 'h':
          hflag = true;
          break;
        default:
          break;
        }
    }
  
  if (hflag)
    {
      usage(progname);
      return 0;
    }

  std::string filename,basename;
  std::string fname_root;
  std::string fname_ascii;

  if (argc-optind == 1)
    {
      filename = argv[optind];
      fname_ascii = filename;
      // Is the path included in the filename?
      std::string::size_type idx = filename.rfind('/');
      if (idx == std::string::npos) 
	{
	  basename = filename;
	}
      else 
	{
	  basename = filename.substr(idx+1);
	}

      // Clear the extesion if it exists
      idx = basename.find_last_of('.');
      if (idx == std::string::npos) 
	{
	  // No extension nor dot
	  fname_root = basename + ".root";
	} 
      else 
	{
	  // There is a dot
	  std::string extension = basename.substr(idx+1);
	  if (extension.empty() ) 
	    {
	      // There is no extension but filename is finished with dot
	      fname_root = basename + "root";
	    } 
	  else if ( extension.size() > 3 ) 
	    {
	      // There is not extension and dot belongs to filename
	      fname_root = basename + ".root";
	    } 
	  else 
	    {
	      // remove the extension
	      fname_root = basename.substr(0,idx);
	      fname_root += ".root";
	    }
	}
    } 
  else if (argc-optind == 2) 
    {
      fname_ascii = argv[optind];
      fname_root = argv[optind+1];
    }
  else 
    {
      usage(progname);
      return 1;
    }
  

  // Create root file
  TFile * file = new TFile(fname_root.c_str(),"RECREATE",
			   "G4PreCompoundTest program output",1);
  
  
  if (!file) 
    {
      std::cerr << "Can not create " << fname_root << '\n';
      return 2;
    }
  file->SetCompressionLevel(1);
  // Create tree
  TTree * tree = new TTree("PCT","G4PreCompound Tree");
  if (tree) tree->SetAutoSave(1000000); // Auto save when 1 Mb written 
  else 
    {
      std::cerr << "Can not create tree" << '\n';
      return 3;
    }
    
  Int_t bufsize = 256000;
  Int_t split = 99; // 1
  if (split != 0) bufsize /= 4; 
    
  // Create branch
  precoreaction * reaction = new precoreaction();
  TBranch * branch = tree->Branch("reaction","precoreaction",&reaction,bufsize,split);
  if (!branch) 
    {
      std::cerr << "Can't create branch\n";
      return 4;
    }
  //  branch->SetAutoDelete(kFALSE);

  
  // Open the data file 
  igzstream is(fname_ascii.c_str());
  if (!is) 
    {
      std::cerr << "Can not open " << fname_ascii << '\n';
      return 5;
    }

  char c;
  // Discard first lines
  while ( (c = is.peek()) != '#') {
    is.ignore(1000,'\n');
  }
  
  // Counting variables
  int REACTN = 0, TOTALFRAGS = 0;

  // Variables to be read from file
  //-------------------------------- 
  // Misc:
  double ProjectileXS = 0.0;
  // Projectile:
  TLorentzVector ProjectileP;
  double ProjectileM(0.0);
  int ProjectileA(-1);
  int ProjectileZ(-1);
  // Target:
  int TargetA(-1);
  int TargetZ(-1);
  double TargetM(0.0);
  // Compound:
  int CompoundA(-1);
  double CompoundU = 0.0;
  TLorentzVector CompoundP;
  double CompoundM = 0.0;
  // Fragments:
  std::string FragmentName;
  int FragmentA;
  int FragmentZ;
  double FragmentMass;
  TLorentzVector FragmentMomentum;
  std::string FragmentCreator;


  std::cout << "Reading reaction number " << std::setw(8) << 0;
  
  // Start reading reactions

  // Read the header
  while ( (c = is.peek()) != '*' ) 
    {
      // Forget about comment lines (#) and empty lines. 
      if ( c == '#' || c == '\n' || c == ' ') is.ignore(10000,'\n');
      // Lines starting with % contains reaction info
      else if (c == '%')
	{
	  std::string Line;
	  char tmp[1001];
	  is.getline(tmp,1000);
	  Line = tmp;
	  if (Line.find("XS =",0) != std::string::npos) 
	    {
	      std::string::size_type idx = Line.find('=');
	      ProjectileXS =  ToDouble(Line.substr(idx+1,std::string::npos));
	    }
	  else if (Line.find("PA =",0) != std::string::npos) 
	    {
	      // Get the Projectile baryonic number
	      std::string::size_type idx = Line.find('=');
	      ProjectileA =  ToInt(Line.substr(idx+1,std::string::npos));
	    }
	  else if (Line.find("PZ =",0) != std::string::npos) 
	    {
	      // Get the Projectile charge
	      std::string::size_type idx = Line.find('=');
	      ProjectileZ =  ToInt(Line.substr(idx+1,std::string::npos));
	    }
	  else if (Line.find("PP =",0) != std::string::npos) 
	    {
	      // Get the Projectile Momentum
	      std::string::size_type idx1 = Line.find('(');
	      std::string::size_type idx2 = Line.find(',',idx1+1);
	      double p(ToDouble(Line.substr(idx1+1,idx2-idx1)));
	      ProjectileP.SetX(p);
	      idx1 = idx2;
	      idx2 = Line.find(',',idx1+1);
	      p = ToDouble(Line.substr(idx1+1,idx2-idx1));
	      ProjectileP.SetY(p);
	      idx1 = idx2;
	      idx2 = Line.find(';',idx1+1);
	      p = ToDouble(Line.substr(idx1+1,idx2-idx1));
	      ProjectileP.SetZ(p);
	      idx1 = idx2;
	      idx2 = Line.find(')',idx1+1);
	      p = ToDouble(Line.substr(idx1+1,idx2-idx1));
	      ProjectileP.SetE(p);
	    }
	  else if (Line.find("PM =",0) != std::string::npos) 
	    {
	      std::string::size_type idx = Line.find('=');
	      ProjectileM = ToDouble(Line.substr(idx+1,std::string::npos));
	    }
	  else if (Line.find("TA =",0) != std::string::npos) 
	    {
	      // Get the Target baryonic number
	      std::string::size_type idx = Line.find('=');
	      TargetA =  ToInt(Line.substr(idx+1,std::string::npos));
	    }
	  else if (Line.find("TZ =",0) != std::string::npos) 
	    {
	      // Get the Target charge
	      std::string::size_type idx = Line.find('=');
	      TargetZ =  ToInt(Line.substr(idx+1,std::string::npos));
	    }
	  else if (Line.find("TM =",0) != std::string::npos) 
	    {
	      // Get the Target mass
	      std::string::size_type idx = Line.find('=');
	      TargetM =  ToDouble(Line.substr(idx+1,std::string::npos));
	    }
	  //  	else if (Line.find("CU =",0) != std::string::npos)
	  //  	{
	  //  	   // Get the Compound Excitation Energy
	  //  	   std::string::size_type idx = Line.find('=');
	  //  	   CompoundU = ToDouble(Line.substr(idx+1,std::string::npos));
	  //  	}
	  //  	else if (Line.find("CP =",0) != std::string::npos)
	  //  	{
	  //  	   // Get the Compound Momentum
	  //  	   std::string::size_type idx1 = Line.find('(');
	  //  	   std::string::size_type idx2 = Line.find(',',idx1+1);
	  //  	   double p(ToDouble(Line.substr(idx1+1,idx2-idx1)));
	  //  	   CompoundP.SetX(p);
	  //  	   idx1 = idx2;
	  //  	   idx2 = Line.find(',',idx1+1);
	  //  	   p = ToDouble(Line.substr(idx1+1,idx2-idx1));
	  //  	   CompoundP.SetY(p);
	  //  	   idx1 = idx2;
	  //  	   idx2 = Line.find(';',idx1+1);
	  //  	   p = ToDouble(Line.substr(idx1+1,idx2-idx1));
	  //  	   CompoundP.SetZ(p);
	  //  	   idx1 = idx2;
	  //  	   idx2 = Line.find(')',idx1+1);
	  //  	   p = ToDouble(Line.substr(idx1+1,idx2-idx1));
	  //  	   CompoundP.SetE(p);
	  //  	}
	  //  	else if (Line.find("CM =",0) != std::string::npos)
	  //  	{
	  //  	   // Get the Compound Mass
	  //  	   std::string::size_type idx = Line.find('=');
	  //  	   CompoundM = ToDouble(Line.substr(idx+1,std::string::npos));
	  //  	}
	}
    }
  
  while ( (c = is.peek()) != EOF ) 
    {
      if (c == '*')
	{
	  double px,py,pz,e;
	  is.ignore(1000000,'\n');
       
	  std::cout << "\b\b\b\b\b\b\b\b" << std::setw(8) << ++REACTN;
	  std::cout.flush();

	  is >> CompoundA
	     >> px
	     >> py
	     >> pz
	     >> e
	     >> CompoundU
	     >> CompoundM;

	  CompoundP.SetX(px);
	  CompoundP.SetY(py);
	  CompoundP.SetZ(pz);
	  CompoundP.SetE(e);
	    
	  reaction->SetCrossSection(ProjectileXS);
	  reaction->SetProjectileA(ProjectileA);
	  reaction->SetProjectileZ(ProjectileZ);
	  reaction->SetProjectileMomentum(ProjectileP);
	  reaction->SetProjectileMass(ProjectileM);
	  reaction->SetTargetA(CompoundA-ProjectileA);
	  reaction->SetTargetZ(TargetZ);
	  reaction->SetTargetMass(TargetM);
	  reaction->SetCompoundU(CompoundU);
	  reaction->SetCompoundMomentum(CompoundP);
	  reaction->SetCompoundMass(CompoundM);
	  reaction->SetReactionNumber(REACTN);
	
	  do 
	    {
	      is >> FragmentName
		 >> FragmentA
		 >> FragmentZ
		 >> FragmentMass
		 >> px
		 >> py
		 >> pz
		 >> e;
	   
	      FragmentMomentum.SetPxPyPzE(px,py,pz,e);

	      is >> FragmentCreator;
	      if (FragmentCreator == "No") 
		{
		  FragmentCreator += " Name";
		}
	      is.ignore(1000000,'\n');
	   
	      reaction->AddFragment(FragmentName,FragmentA,FragmentZ,
				    FragmentMass,FragmentMomentum,
				    FragmentCreator);
	      TOTALFRAGS++;
	    }
	  while ( (c = is.peek()) != '*');
	  tree->Fill();
	  reaction->Clear(); 
	  is.ignore(1000000,'\n');
	}
      if (limit_event > 0 && REACTN == limit_event) break;
    }
  
  std::cout
    << "\n\n" 
    << "Total number of reactions read: "
    << REACTN
    << '\n'
    << "Total number of fragments read: "
    << TOTALFRAGS
    << '\n'
    << "Mean number of fragments per reaction: " 
    << double(TOTALFRAGS)/double(REACTN)
    << "\n\n";  
  
  file->Write();
  file->Purge(1);
  file->Close();
  delete file;
  
  file = new TFile(fname_root.c_str());
  std:: cout << file->GetEND() << " bytes written\n";
  delete file;
  
  return 0;
}
