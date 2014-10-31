## @package BRICC
#  Define functions to run BRICC 
#
#
#author:  L.Desorgher
#
#History: 
#---------
#      26/07/2014     Creation


import os
import subprocess
import numpy as np
from utilities import *

Bricc_dir="%s/../extra_executable/Bricc_V23" %(os.path.dirname(__file__))

## Run of bricc program to compute internal conversion coefficient
#
#
def RunLocalBriccWithEnsdfFile(data_card_ensdf,verbose=False,parallel_run_nb=None):
        #Go on Bricc local directory 
        ###########################
        #print data_card_ensdf
        if verbose:
            print "In Bricc", data_card_ensdf
        Bricc_run_dir=Bricc_dir
        if parallel_run_nb is not None and type(parallel_run_nb) == type(1): 
            Bricc_run_dir="%s/run%i" %(Bricc_dir,parallel_run_nb)
            Bricc_run_dir=Bricc_run_dir.replace("//","")
        loc_dir=os.getcwd()
        os.chdir(Bricc_run_dir)
        
        #Make the input ensdf file
        #######################
        file_object=open("input.ensdf",'w')
        name_nuc=data_card_ensdf.split()
        file_object.write(" 21NE    ADOPTED LEVELS, GAMMAS                                  04NDS    200412\n")
        print data_card_ensdf
        file_object.write(data_card_ensdf)
        file_object.write("\n")
        file_object.close()
        
        #Run the Bricc program
        ######################
       
        cmd ="env  BrIccHome=%s  ./bricc input.ensdf" %(Bricc_dir)
        
        p=subprocess.Popen(cmd, shell=True,
                                stderr=subprocess.PIPE,
                                           stdin=subprocess.PIPE,
                                           stdout=subprocess.PIPE)
        
        p.stdin.write("output.txt\n" )
        p.stdin.write("\n")
        p.stdin.write("\n")
        p.stdin.write("Y\n")
        p.stdin.write("N\n")
        p.stdin.write("\n")
        p.stdin.write("\n")
        p.wait()
        """
        p.wait()
        while(p.poll() == None):
            i=0
        """  
        
        
        #Get the results from the output.txt file
        #######################################
        file_object=open("output.txt",'r')
        lines=file_object.readlines()
        file_object.close()
        file_object=open("res_table.txt",'w')
        write_line=False
        table_finished=False
        id_Icc=-1
        table_is_valid=False
        for line in lines:
            if verbose:
                print line
            words=line.split()
            id_dIcc=line.find("dIcc")
            if  id_Icc==-1 and id_dIcc !=-1 and not table_finished:
                write_line=True
                id_Icc=line.find("Icc")
               
            if (write_line and not table_finished and line!='\n' and len(words)>1 and is_number(words[1])):
                line_to_write=line[0:8]+line[id_Icc-1:id_Icc+9]
                file_object.write(line_to_write+"\n")
                if verbose:
                    print line_to_write
                    print line
            if(write_line and line.find("Tot")!=-1):
                table_finished=True
                if len(words)>1 and is_number(words[1]):
                    table_is_valid=True
        
        file_object.close()
        icc_vec=np.zeros(10)
        ipf=0.

        if (table_is_valid):
            shell_id_vec,ICC_vec=np.loadtxt("res_table.txt", 
                                                      usecols=(0,1),
                                                     unpack=True,converters={0:shell_label_to_id,
                                                                             1:lambda s:float(s or 0.)})
            if verbose:
                print shell_id_vec,ICC_vec
            for i in range(len(shell_id_vec)):
                shell_id=int(shell_id_vec[i])
                if (shell_id>=0):
                    if (shell_id<10):
                        icc_vec[shell_id]+=ICC_vec[i]
                    else:
                        ipf+=ICC_vec[i]
        
        os.chdir(loc_dir)
        #print icc_vec,sum(icc_vec)
        return icc_vec,ipf

     