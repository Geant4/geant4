//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/dir.h>
#include <dirent.h>
#include <unistd.h>
#include <string.h>

#include <fstream>
#include <strstream>
#include <iostream>
#include <iostream>
#include <sstream>
#include "globals.hh"


#include "G4HadFileFinder.hh"


////////////////////////////////////////////////////////////////
//
//

G4HadFileFinder::G4HadFileFinder(G4String& currdir, G4String& G4filetofind)
{

  char filetofindone;
  char* pfiletofindone=&filetofindone;

  char directory;
  char* pdirectory=&directory;

  char nameChar[100] = {""};

  std::ostrstream ost(nameChar, 100, std::ios::out);
  ost << G4filetofind ;

  pfiletofindone=nameChar;


  char nameDir[100] = {""};

  std::ostrstream ostdir(nameDir, 100, std::ios::out);
  ostdir << currdir;

  pdirectory=nameDir;
  

  // G4cout<<currdir<<nameChar<<G4endl;


  char* workingpath = getenv("G4HADWORKING");
  if ( !workingpath )
  {
    G4String excep = "";
    excep += "G4HADWORKING environment variable not set";
    G4Exception(excep);
  }

  int errnum=chdir(workingpath);

  G4cout << "File to find " << G4filetofind << G4endl;

  crawldir(pdirectory,pfiletofindone);

  errnum=chdir(workingpath);
}



////////////////////////////////////////////////////////////////
//
//

void G4HadFileFinder::crawldir(char* currdir, char* filetofind)
{
  /* Open the directory */
  DIR *dir;
  struct dirent *de;
  struct stat st;
  char* totaldir;
  char tempname[MAXNAMLEN];

  char* helpme="HELPME";

  
    if (strcmp(filetofind,helpme)!=0) 
      {
	chdir (currdir);
	totaldir=getwd(tempname);
	G4cout << totaldir << G4endl;
	dir = opendir(".");
	while (de = readdir(dir)) {
	  lstat(de->d_name, &st);
	  if (S_ISDIR(st.st_mode)) {
	    if ((strcmp(de->d_name, ".") != 0) && (strcmp(de->d_name, "..") != 0)) {
	      crawldir(de->d_name, filetofind);
	    }
	  }
	  else {
	    StripFile(de->d_name, totaldir, filetofind);
	  }
	}
	closedir(dir);
	//  G4cout<<'Switching back one dir level'<<G4endl;
	chdir("..");
      }
    else
      {
	 char* workingpath = getenv("G4HADATASET");
	 if ( !workingpath )
	   {
	     G4String excep = "G4HADATASET ";
	     excep += "environment variable not set";
	     G4Exception(excep);
	   }

	 int errnum=chdir(workingpath);

	std::ifstream fin("G4HDSFileFinderREADME.txt");
	std::filebuf* finbuf = fin.rdbuf();
	
	if ( !( finbuf->is_open() ) )
	  {
	    G4String excep = "G4HadDataSet -  G4HDSFileFinderREADME.txt not found. Did you delete something?";
	    excep +="The help file G4HDSFileFinderREADME.txt is no more where it should be";
	    G4Exception(excep);
	  }
	else
	  {
	    G4String tmpString;
	    while (!fin.eof())
	      {
		getline(fin,tmpString);
		G4cout << tmpString << G4endl;
	      }
	  }
      }
}

///////////////////////////////////////////////////////////////////
//
//

void G4HadFileFinder::StripFile(char* filename, char* dirposition, char* findmefile) 
{
  /*
  if (strcmp(filename,findmefile)==0) 
    {
      G4cout<< filename<<G4endl;
      G4cout<< dirposition<<G4endl;
      }
  */
  
  char* controlA="A*";
  char* controlZ="Z*";
  char* control="*";
  
  char *condition;
  char *pcontrolchar="<>";
  char numbers[]="1234567890.";
  char letters[]="qwertyuiopasdfghjklzxcvbnm";
  int numberposition;
  int letterposition;

  int length=strlen(findmefile);
  char tempfilename[length];
  char tempfilenameone[length];
  char tempfilenametwo[length];
  char copyfindmefile[length];
  

  int lengthtwo=strlen(filename);
  char copymainfile[lengthtwo];
  char *ptempfilename=tempfilename;
  char *ptempfilenameone=tempfilenameone;
  char *ptempfilenametwo=tempfilenametwo;
  char *pcopymainfile=copymainfile;
  char *pcopyfindmefile=copyfindmefile;

  int pch=-1;
  char* wherechar;
  int position;
    
  double Zfindmefile;
  double Afindmefile;
  double Zfilename;
  double Afilename;
  char *Acharfilename;
  char *Acharfindmefile;
    
  
  //--------------------------------------------------------------------------
  //
  // Case when is looking for a specified file or for all the files
  //
  //
  //--------------------------------------------------------------------------

  if ((strcmp(filename,findmefile)==0)||(strcmp(findmefile,control)==0))
    {
      G4cout << "File found in" << G4endl;
      G4cout<< dirposition << "/"<< filename << G4endl;
    }

  else
    pch=strcspn(findmefile,control);
    strcpy(copymainfile,filename);
 

  //--------------------------------------------------------------------------
  //
  // Special case if there is a * in the findmefile for the process
  //
  //
  //--------------------------------------------------------------------------   
    
    if (pch==0) 
      {
	ptempfilename=strtok(pcopymainfile,"Z");
	ptempfilename=strtok(NULL,"\t");
	strcpy(ptempfilenametwo,"*Z");
	strcat(ptempfilenametwo,ptempfilename);
	strcpy(copymainfile,ptempfilenametwo);

	wherechar=strstr(findmefile,controlZ);
	if (wherechar==NULL)
	  {
	    wherechar=strstr(findmefile,controlA);
	    if (wherechar!=NULL)
	      {
		strcpy(ptempfilenametwo,findmefile);
		//strcpy(pcopyfindmefile,findmefile);
		position=strcspn(copymainfile,"A");
		strncpy(ptempfilenametwo,copymainfile,position);
		// G4cout<< " help 1 " << ptempfilenametwo << " " << copymainfile << G4endl;
		ptempfilenametwo=strtok(ptempfilenametwo,"A");		  
		//G4cout<< " help 2" << ptempfilenametwo << G4endl;
		strcat(ptempfilenametwo,"A*");
		//G4cout<< " help 3" << ptempfilenametwo << G4endl;
		
		ptempfilenameone=strtok(copymainfile,"A");
		ptempfilenameone=strtok(NULL,"\t");
		letterposition=strspn(ptempfilenameone,numbers);
		//G4cout<< " help 4 " << ptempfilenameone << " " << letterposition <<G4endl;
		Acharfilename=strtok(ptempfilenameone,&ptempfilenameone[letterposition-1]);
		//G4cout<< " help 4.5" << Acharfilename << G4endl;
		Acharfilename=strtok(NULL,"\t");
		//G4cout<< " help 5" << Acharfilename << G4endl;
		
		strcat(ptempfilenametwo,Acharfilename);
		//G4cout<< " help 6" << ptempfilenametwo << G4endl;
		
		if (strcmp(ptempfilenametwo,findmefile)==0) 
		  {
		    G4cout << "File found in" << G4endl;
		    G4cout<< dirposition << "/"<< filename << G4endl;
		  }
	      }
	    else 
	      if (strcmp(ptempfilenametwo,findmefile)==0) 
		{
		  G4cout << "File found in" << G4endl;
		  G4cout<< dirposition << "/"<< filename << G4endl;  
		} 
	    
	  }
	else 
	  {
	    wherechar=strstr(findmefile,controlA);
	    if (wherechar==NULL)
	      {
		strcpy(ptempfilenametwo,findmefile);
		position=strcspn(copymainfile,"Z");
		strncpy(ptempfilenametwo,copymainfile,position);
		//	  G4cout<< " help " << ptempfilenametwo<< G4endl;
		
		ptempfilenametwo=strtok(ptempfilenametwo,"Z");
		strcat(ptempfilenametwo,"Z*");
		
		//G4cout<< " help 1 " << ptempfilenametwo << " " << copymainfile << G4endl;
		
		ptempfilenameone=strtok(copymainfile,"Z");
		ptempfilenameone=strtok(NULL,"\t");
		
		//  G4cout<< " help 4 " << ptempfilenameone << " " << letterposition <<G4endl;
		
		Acharfilename=strtok(ptempfilenameone,"A");
		
		//	  G4cout<< " help 4.5 " << Acharfilename << G4endl;
		Acharfilename=strtok(NULL,"\t");
		//	  G4cout<< " help 5" << Acharfilename << G4endl;
		
		strcat(ptempfilenametwo,"A");
		strcat(ptempfilenametwo,Acharfilename);
		
		//G4cout<< " help " << ptempfilenametwo << G4endl;
		if (strcmp(ptempfilenametwo,findmefile)==0) 
		  {
		    G4cout << "File found in" << G4endl;
		    G4cout<< dirposition << "/"<< filename << G4endl;
		  }
	      }
	      
	  } 
      }
    
  //--------------------------------------------------------------------------
  //
  // Special case if there is * in the findmefile for the for the Z or A
  //
  //--------------------------------------------------------------------------   

    else
      {
      if (pch<=length)
	{
	  wherechar=strstr(findmefile,controlZ);
	  if (wherechar==NULL)
	    {
	      wherechar=strstr(findmefile,controlA);
	      if (wherechar!=NULL)
		{
		  strcpy(ptempfilenametwo,findmefile);
		  //strcpy(pcopyfindmefile,findmefile);
		  position=strcspn(copymainfile,"A");
		  strncpy(ptempfilenametwo,copymainfile,position);
		  // G4cout<< " help 1 " << ptempfilenametwo << " " << copymainfile << G4endl;
		  ptempfilenametwo=strtok(ptempfilenametwo,"A");		  
		  //G4cout<< " help 2" << ptempfilenametwo << G4endl;
		  strcat(ptempfilenametwo,"A*");
		  //G4cout<< " help 3" << ptempfilenametwo << G4endl;
		
		  ptempfilenameone=strtok(copymainfile,"A");
		  ptempfilenameone=strtok(NULL,"\t");
		  letterposition=strspn(ptempfilenameone,numbers);
		  //G4cout<< " help 4 " << ptempfilenameone << " " << letterposition <<G4endl;
		  Acharfilename=strtok(ptempfilenameone,&ptempfilenameone[letterposition-1]);
		  //G4cout<< " help 4.5" << Acharfilename << G4endl;
		  Acharfilename=strtok(NULL,"\t");
		  //G4cout<< " help 5" << Acharfilename << G4endl;
		  
		  strcat(ptempfilenametwo,Acharfilename);
		  //G4cout<< " help 6" << ptempfilenametwo << G4endl;
		
		  if (strcmp(ptempfilenametwo,findmefile)==0) 
		    {
		      G4cout << "File found in" << G4endl;
		      G4cout<< dirposition << "/"<< filename << G4endl;
		    }
		}
	    }
	  else 
	    {
	      wherechar=strstr(findmefile,controlA);
	      if (wherechar==NULL)
		{
		  strcpy(ptempfilenametwo,findmefile);
		  position=strcspn(copymainfile,"Z");
		  strncpy(ptempfilenametwo,copymainfile,position);
		  //	  G4cout<< " help " << ptempfilenametwo<< G4endl;
		  
		  ptempfilenametwo=strtok(ptempfilenametwo,"Z");
		  strcat(ptempfilenametwo,"Z*");
		  
		  //G4cout<< " help 1 " << ptempfilenametwo << " " << copymainfile << G4endl;
	  
		  ptempfilenameone=strtok(copymainfile,"Z");
		  ptempfilenameone=strtok(NULL,"\t");
		 
		  //  G4cout<< " help 4 " << ptempfilenameone << " " << letterposition <<G4endl;
		 
		  Acharfilename=strtok(ptempfilenameone,"A");
	
		  //	  G4cout<< " help 4.5 " << Acharfilename << G4endl;
		  Acharfilename=strtok(NULL,"\t");
		  //	  G4cout<< " help 5" << Acharfilename << G4endl;
	
		  strcat(ptempfilenametwo,"A");
		  strcat(ptempfilenametwo,Acharfilename);
		  	 
		  //G4cout<< " help " << ptempfilenametwo << G4endl;
		  if (strcmp(ptempfilenametwo,findmefile)==0) 
		    {
		      G4cout << "File found in" << G4endl;
		      G4cout<< dirposition << "/"<< filename << G4endl;
		    }
		}
	      else if (wherechar!=NULL)
		{
		  strcpy(ptempfilenametwo,findmefile);
		  position=strcspn(copymainfile,"Z");
		  strncpy(ptempfilenametwo,copymainfile,position);
		  
		  ptempfilenametwo=strtok(ptempfilenametwo,"Z");
		  strcat(ptempfilenametwo,"Z*A*");

		  ptempfilenameone=strtok(copymainfile,"A");
		  ptempfilenameone=strtok(NULL,"\t");
		  letterposition=strspn(ptempfilenameone,numbers);
		  G4cout<< " help 4 " << ptempfilenameone << " " << letterposition <<G4endl;
		  Acharfilename=strtok(ptempfilenameone,&ptempfilenameone[letterposition-1]);
		  G4cout<< " help 4.5" << Acharfilename << G4endl;
		  Acharfilename=strtok(NULL,"\t");
		  G4cout<< " help 5" << Acharfilename << G4endl;
		  
		  strcat(ptempfilenametwo,Acharfilename);
		  G4cout<< " help 6" << ptempfilenametwo << G4endl;
		

		  if (strcmp(ptempfilenametwo,findmefile)==0) 
		    {
		      G4cout << "File found in" << G4endl;
		      G4cout<< dirposition << "/"<< filename << G4endl;
		    }
		}
	    }
	}
      // end special case
      
  //--------------------------------------------------------------------------
  //
  // Special case if there is a condition <> in the findmefile for the for the Z or A
  //
  //--------------------------------------------------------------------------   

	      position=strcspn(findmefile,pcontrolchar);
	      // G4cout << "I am in this position " << position << G4endl;
	      if (position<length)
		{
		  strcpy(ptempfilenametwo,findmefile);
		  
		  // Looking for condition on Z
		  if (ptempfilenametwo[position-1]=='Z') 
		    {
		      //      G4cout << "I am here " << G4endl;
		      ptempfilenametwo=strtok(ptempfilenametwo,"Z");
		      ptempfilenametwo=strtok(NULL,"A");
		      
		      ptempfilename=strtok(copymainfile,"Z");
		      ptempfilename=strtok(NULL,"A");
		      // G4cout << "Find this Z " << ptempfilename << G4endl;
		      Zfilename=atof(ptempfilename);
		      // G4cout << "Find this Z number " << Zfilename<< G4endl;
		      
		      // G4cout << "I am probably here " << ptempfilenametwo << G4endl;
		      
		      if (ptempfilenametwo[0]=='<')
			{
			  ptempfilenametwo=strtok(ptempfilenametwo,"<");
			  Zfindmefile=atof(ptempfilenametwo);
			  //   G4cout << "I am in this position <" << ptempfilenametwo << G4endl;
			  // G4cout << " my Z " << Zfindmefile << G4endl;
			  if (Zfilename<Zfindmefile)
			    {
			      G4cout << "File found in" << G4endl;
			      G4cout<< dirposition << "/"<< filename << G4endl;
			    }
			}
		      else 
			if (ptempfilenametwo[0]=='>')
			  {
			    ptempfilenametwo=strtok(ptempfilenametwo,">");
			    Zfindmefile=atof(ptempfilenametwo);
			    //G4cout << "I am in this position >" << ptempfilenametwo << G4endl;
			    //G4cout << " my Z " << Zfindmefile << G4endl;
			    if (Zfilename>Zfindmefile)
			      {
				G4cout << "File found in" << G4endl;
				G4cout<< dirposition << "/"<< filename << G4endl;
			      }
			  }
		      // G4cout << "I am in this position " << ptempfilenametwo << G4endl;
		    }
		  else 
		    
		    
		    // Looking for condition on A
		    if (ptempfilenametwo[position-1]=='A') 
		      {
			//G4cout << "I am here " << G4endl;
			ptempfilenametwo=strtok(ptempfilenametwo,"A");
			ptempfilenametwo=strtok(NULL,"\t");
			
			ptempfilename=strtok(copymainfile,"A");
			ptempfilename=strtok(NULL,"\t");			
			//G4cout << "I don t like that " << ptempfilename << G4endl;
			numberposition = strspn(ptempfilename,numbers);
			//G4cout << "I don t like this " << numberposition << G4endl;
			Acharfilename=strtok(ptempfilename,&ptempfilename[numberposition]);
			//G4cout << "I don t like this " << Acharfilename<< G4endl;
			Afilename=atof(Acharfilename);
			
			//G4cout << "Find this A number " << Afilename<< G4endl;
			
			//G4cout << "I am probably here " << ptempfilenametwo << G4endl;
			
			if (ptempfilenametwo[0]=='<')
			  {
			    ptempfilenametwo=strtok(ptempfilenametwo,"<");
			    numberposition = strspn(ptempfilenametwo,numbers);
			    Acharfindmefile=strtok(ptempfilenametwo,&ptempfilenametwo[numberposition]);
			    Afindmefile=atof(Acharfindmefile);
			    //  G4cout << " my A " << Afindmefile << G4endl;
			    
			  if (Afilename<Afindmefile)
			    {
			      G4cout << "File found in" << G4endl;
			      G4cout<< dirposition << "/"<< filename << G4endl;
			    }
			  }
			else 
			  if (ptempfilenametwo[0]=='>')
			    {
			      ptempfilenametwo=strtok(ptempfilenametwo,">");
			      //G4cout << "I don t like that " << ptempfilenametwo << G4endl;
			      numberposition = strspn(ptempfilenametwo,numbers);
			      Acharfindmefile=strtok(ptempfilenametwo,&ptempfilenametwo[numberposition]);
			      //G4cout << "I don t like that too  " << Acharfindmefile << G4endl;
			      Afindmefile=atof(Acharfindmefile);
			      //G4cout << " my A " << Afindmefile << G4endl;
			      if (Afilename>Afindmefile)
				{
				  G4cout << "File found in" << G4endl;
				  G4cout<< dirposition << "/"<< filename << G4endl;
				}
			    }
			//G4cout << "I am in this position " << ptempfilenametwo << G4endl;
		      }
		  
		}
      }
}


