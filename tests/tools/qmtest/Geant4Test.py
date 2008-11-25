import os
import sys
from qm.test.classes.command import ShellCommandTest

class Geant4Test (ShellCommandTest):

  def Run(self, context, result):
    installdir = os.environ['G4INSTALL']
    workdir = os.environ['G4WORKDIR']
    system = os.environ['G4SYSTEM']
    
    if sys.platform == 'win32':
        self.command = self.command.replace('$G4INSTALL',installdir)
        self.command = self.command.replace('$G4BIN',installdir)
        self.command = self.command.replace('$G4SYSTEM',system)
        self.command = self.command.replace(';',' && ')
        if self.command.count('make') > 0:
            # We are starting a build test
            self.command = self.command.replace('\\','/')
        else:
            # We are starting a run test
            self.command = self.command.replace('/','\\')
            
    #TEMPORAL: EVERYTHING TO BE TESTED ON CYGWIN
    #TODO: JUST BUILD TESTS ON CYGWIN; RUN ON WINDOWS CMD
    self.command = self.command.replace('\\','/')
    
    print self.command
    ShellCommandTest.Run(self, context, result)