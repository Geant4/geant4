#include <string>

#include "TApplication.h"

#include "mkplfilemanager.h"
#include "mkplsimmanager.h"
#include "mkplexpdatamanager.h"

//#include "cmdargs.h"

#include <unistd.h>

using namespace std;

void usage(const string& progname)
{
	cout << "Usage: " << progname << " [-c] [-t] [-p] simulated_file experimental_file\n\n"
	     <<  "      simfile:        simulated data file\n"
	     <<  "      expfile:        experimental data file\n"
	     <<  "      -c:             draw model components\n"
	     <<  "      -t:             enable basic tests\n"
	     <<  "      -p:             generate postscript file for plots\n";
	return;
}

int main(int argc, char **argv)
{  
//     // Declare arguments
//     CmdArgStr   siminput("sim-file","simulated data file");
//     CmdArgStr   expinput("exp-file","experimental data file");
//     CmdArgBool  testflag('t',"test","turn on simple tests.");
//     CmdArgBool  psflag('p',"postscript","generate postscrip files of the plots");
//     
//     // Declare commnad object and its argument-iterator
//     CmdLine cmd(*argv,&testflag,&psflag,&siminput,&expinput,NULL);
//     CmdArgvIter arg_iter(--argc,++argv);
//     
//     // Initialize optional arguments to appropriate default values
//     testflag = false;
//     psflag = false;
//     
//     // Parse arguments
//     cmd.parse(arg_iter);
// 
// 
//     string simfilename(siminput);
//     string expfilename(expinput);
//     bool opt_ps(psflag);
//     bool opt_test(testflag);

	string progname(argv[0]);
	int compflag(0);
	int psflag(0);
	int testflag(0);
	int c;
	opterr = 0;
	while ((c = getopt(argc,argv,"ctp")) != -1)
	{
		switch(c)
			{
				case 'c':
					compflag = 1;
					break;
				case 't':
					testflag = 1;
					break;
				case 'p':
					psflag = 1;
					break;
				case '?':
					cout << "Unknown option " << optopt << "\n\n";
					usage(progname);
					return 1;
				default:
					break;
			}		
	}

	if (argc-optind != 2) 
	{
		usage(progname);
		return 2;
	}
	
	bool opt_comp(false);
	bool opt_ps(false);
	bool opt_test(false);
	if (compflag != 0) opt_comp = true;
	if (psflag != 0) opt_ps = true;
	if (testflag != 0) opt_test = true;
	
	string simfilename(argv[optind]);
	string expfilename(argv[optind+1]);

    TApplication makeplotsApp("makeplots", &argc, argv);
  
    mkplfilemanager * fm = new mkplfilemanager();
    fm->OpenExpFile(expfilename);
    fm->OpenSimFile(simfilename);

    mkplexpdatamanager * expm = new mkplexpdatamanager(fm);
    mkplsimmanager * simm = new mkplsimmanager(fm);
    simm->PrintInfo();
    
    
    //+++++++++++++++++++++
    //+ Create histograms +
    //+++++++++++++++++++++

    if (opt_test) simm->PrepareTestHistograms();
  
    simm->PrepareComparisonHistograms(expm,opt_comp);

    
    //++++++++++++++++++++
    //+ Fill histograms) +
    //++++++++++++++++++++
	  
    simm->FillHistograms();

    //+++++++++++++++++++++++
    //+ Plot the histograms +
    //+++++++++++++++++++++++

    simm->Draw(&makeplotsApp,opt_ps);

    delete simm;
    delete expm;
    delete fm;

    return 0;
}

