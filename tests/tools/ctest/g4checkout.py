#!/usr/bin/env python


import sys, os, shutil

if sys.platform == 'win32':
  cern_svn_repos = 'https://svn.cern.ch/reps'
else :
  cern_svn_repos = 'svn+ssh://svn.cern.ch/reps'

def svn_geant4(config_version, destdir, quiet, update, proposed):
    if not os.path.exists(destdir):
        try :
            os.makedirs(destdir)
        except Exception, e:
            print self.argv0, ': ERROR: problem creating directory %s in %s: ' % (destdir, os.getcwd()), e
            sys.exit(1)

    os.chdir(destdir)
    options = ''
    if(quiet) : options += '--quiet '    

    if  config_version.endswith("_branch") :
        command = "svn checkout %s svn+ssh://svn.cern.ch/reps/geant4/branches/geant4/_symbols/%s ." % (options, config_version)
        print command
        sys.stdout.flush()
        os.system(command)
        command="svn checkout %s svn+ssh://svn.cern.ch/reps/g4tests/tags/benchmarks/_symbols/%s benchmarks" % (options, os.environ['geant4_benchmarks'])
        print command
        sys.stdout.flush()
        os.system(command)
        return

    if update:
      for l in os.popen('svn status').readlines():
        if l[0] == '?' : 
          path = l.split()[1]
          print 'removing file/directory: %s' % path
          if os.path.isdir(path): shutil.rmtree(path)
          elif os.path.isfile(path): os.remove(path)

    import urllib
    if proposed: tags_url = "http://lcgapp.cern.ch/spi/cgi-bin/g4tags.py?devline=%s" % config_version
    else:        tags_url = "http://sftweb.cern.ch/geant4/geant4tags/gettaglist/%s" % config_version
    taglist = urllib.urlopen(tags_url).read()
    used_tags = open('gettags.txt', 'w')
    used_tags.write(taglist)
    used_tags.close()

    replace_dot_with_geant4 = lambda str:str.replace(".", "geant4",1)
    remove_dot_slach = lambda str:str.replace("./", "",1)

    class modes:
        NONE = 0
        CHECKOUT = 1 
        ROOT_UPDATE = 2
        SPECIAL_CHECKOUT = 3
        CATEGORIES = 4

    def checkout(line):
        repository, global_name = line.split()
        if update:
          command = "svn switch %s %s/%s/tags/geant4/_symbols/%s ." % (options,cern_svn_repos, repository, global_name)          
        else:
          command = "svn checkout %s %s/%s/tags/geant4/_symbols/%s . -N" % (options, cern_svn_repos, repository, global_name)
        return command

    def root_update(line):
        what, = line.split()
        what = remove_dot_slach(what)
        if update and os.path.exists(what):
          command = ""
        else:
          command = "svn update %s %s"% (options, what)
        return command

    def special_checkout(line):
        repository, name, path, origin_path = line.split()
        command = "svn checkout %s %s/%s/tags/%s/_symbols/%s %s"%(options, cern_svn_repos, repository, origin_path, name, path)
        return command

    def switch(line):
        temp = line.split()
        file = None;
        if len(temp) == 4: status, repository, name, path = temp
        elif len(temp) == 5: status, repository, name, path, file = temp

        origin_path = replace_dot_with_geant4(path)
        if file:
            if os.path.exists(path+"/"+file):
                command = "svn switch %s %s/%s/tags/%s/_symbols/%s/%s %s"%(options, cern_svn_repos, repository, origin_path, name, file, path+"/"+file)
            else:
                if not os.path.exists(path): os.makedirs(path)
                command = "svn cp %s %s/%s/tags/%s/_symbols/%s/%s %s"%(options, cern_svn_repos, repository, origin_path, name, file, path+"/"+file)
        else:
            if os.path.exists(path) :
                command = "svn switch %s %s/%s/tags/%s/_symbols/%s %s"%(options, cern_svn_repos, repository, origin_path, name, path)
            else:
                os.makedirs(path)
                command = "svn co %s %s/%s/tags/%s/_symbols/%s %s"%(options, cern_svn_repos, repository, origin_path, name, path)
        return command

    mode = modes.NONE;
    lines = taglist.split('\n')
    for line in lines:
        if line.startswith('#') or line == '': continue
        elif line == "CHECKOUT SECTION": mode = modes.CHECKOUT
        elif line == "ROOT UPDATE SECTION": mode = modes.ROOT_UPDATE
        elif line == "SPECIAL CHECKOUT SECTION": mode = modes.SPECIAL_CHECKOUT
        elif line == "CATEGORIES SECTION": mode = modes.CATEGORIES
        else:
            if mode == modes.CHECKOUT: command = checkout(line)
            elif mode == modes.ROOT_UPDATE: command = root_update(line) 
            elif mode == modes.SPECIAL_CHECKOUT: command = special_checkout(line)
            elif mode == modes.CATEGORIES: command = switch(line)
            
            if command: 
              print command
              sys.stdout.flush()
              if os.system(command) == 1: return 1
        #command = "rm -rf `find . -type d -name .svn`"
        #os.system(command) 
        
if __name__ == "__main__":
  import optparse
  parser = optparse.OptionParser(description='Geant4 checkout utility using the tag collector')
  parser.add_option("-c", "--config", dest="config", default='g4tags-dev',
                      help="Geant4 checkout configuration")
  parser.add_option("-d", "--dest", dest="dest", default=None,
                      help="checkout directory destination")
  parser.add_option("-q", "--quiet", action='store_true', dest="quiet", default=False,
                      help="don't print status messages to stdout")
  parser.add_option("-u", "--update", action='store_true', dest="update", default=False,
                      help="don't do a full checkout, only an update")
  parser.add_option("-p", "--proposed", action='store_true', dest="proposed", default=False,
                      help="add the proposed tags in addtion")

  (options, args) = parser.parse_args()
  if(not options.dest) : options.dest = options.config 
  sys.exit(svn_geant4(options.config, options.dest, options.quiet, options.update, options.proposed))
