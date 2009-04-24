import os
import sys
import urllib
from qm.test.classes.command import ShellCommandTest

class Geant4Test (ShellCommandTest):

  def Run(self, context, result):
    installdir = os.environ['G4INSTALL']
    bindir = os.environ['G4BIN']
    system = os.environ['G4SYSTEM']
   
    if sys.platform == 'win32':
        self.command = self.command.replace('$G4INSTALL',installdir)
        self.command = self.command.replace('$G4BIN',bindir)
        self.command = self.command.replace('$G4SYSTEM',system)
        self.command = self.command.replace(';',' && ')
       
        self.command = self.command.replace('/','\\')
       
        folder = self.command.split('&&')
    else:
        folder = self.command.split(';')

    folder = folder[0].split(' ')
    folder = folder[1]
   
    print self.command
    ShellCommandTest.Run(self, context, result)
   
    os.chdir(folder)
    stdout_file = open(result.GetId()+'_stdout.txt', 'w')
    stdout_file.write(result["ExecTest.stdout"])
    stdout_file.close()
   
    stderr_file = open(result.GetId()+'_stderr.txt', 'w')
    stderr_file.write(result["ExecTest.stderr"])
    stderr_file.close()
