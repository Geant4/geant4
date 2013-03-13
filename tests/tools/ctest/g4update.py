#!/usr/bin/env python
import sys, os, shutil,re

cern_svn_repos = 'svn+ssh://svn.cern.ch/reps'

replace_dot_with_geant4 = lambda str:str.replace(".", "geant4",1)
remove_dot_slash = lambda str:str.replace("./", "",1)
get_tag = lambda path:re.search('^URL: +.*_symbols/([^ \n/]+)', os.popen('svn info %s' % path).read(), 8).group(1)

def g4svn_update(devline, destdir, quiet, proposed):
  if not os.path.exists(destdir):
    try :
      os.makedirs(destdir)
    except Exception, e:
      print sys.argv[0], ': ERROR: problem creating directory %s in %s: ' % (destdir, os.getcwd()), e
      sys.exit(1)
  
  os.chdir(destdir)
  options = ''
  if(quiet) : options += ' --quiet'    
  
  #---In case the development line is a branch...---------------------------------------------------
  if devline.endswith("_branch") :
    command = "svn switch %s svn+ssh://svn.cern.ch/reps/geant4/branches/geant4/_symbols/%s ." % (options, devline)
    print '# ', command
    sys.stdout.flush()
    os.system(command)
    command="svn switch %s svn+ssh://svn.cern.ch/reps/g4tests/tags/benchmarks/_symbols/%s benchmarks" % (options, os.environ['geant4_benchmarks'])
    print '# ', command
    sys.stdout.flush()
    os.system(command)
    return
  #---Analyse the new and switched directories and files -------------------------------------------
  newpaths = []
  swipaths =[]
  for l in os.popen('svn status').readlines():
    stat, path = l.split()
    if stat == '?' and os.path.isdir(path) : newpaths.append(path)
    elif stat == 'S' : swipaths.append(path)
    
  #---Get the current list of tags ------------------------------------------------------------------
  import urllib
  tags_url = "http://lcgapp.cern.ch/spi/cgi-bin/g4tags.py?devline=%s" % devline
  taglist = urllib.urlopen(tags_url).read()
  used_tags = open('gettags.txt', 'w')
  used_tags.write(taglist)
  used_tags.close()

  class modes:
    NONE = 0
    CHECKOUT = 1 
    ROOT_UPDATE = 2
    SPECIAL_CHECKOUT = 3
    CATEGORIES = 4

  def checkout(line):
    repository, global_tag = line.split()
    if not os.path.exists('.svn'):
      command = "svn co %s %s/%s/tags/geant4/_symbols/%s ." % (options, cern_svn_repos, repository, global_tag)
    else :
      if get_tag('.') == global_tag:
        #command = "svn update %s ." % (options)
        command = None
      else:
        command = "svn switch %s %s/%s/tags/geant4/_symbols/%s ." % (options, cern_svn_repos, repository, global_tag)
    return command

  def root_update(line):
    what, = line.split()
    what = remove_dot_slash(what)
    #command = "svn update %s %s"% (options, what)
    command = None
    return command

  def special_checkout(line):
    repository, name, path, origin_path = line.split()
    nakedpath = remove_dot_slash(path) 
    if os.path.exists(path) :
      if get_tag(path) == name :
        command = None
        #command = "svn update %s %s" % (options, path)
      else :
        command = "svn switch %s %s/%s/tags/%s/_symbols/%s %s"%(options, cern_svn_repos, repository, origin_path, name, path)
    else:
      command = "svn co %s %s/%s/tags/%s/_symbols/%s %s"%(options, cern_svn_repos, repository, origin_path, name, path)
    if nakedpath in newpaths: newpaths.remove(nakedpath)
    return command

  def switch(line):
    temp = line.split()
    file = None;
    if len(temp) == 4: status, repository, name, path = temp
    elif len(temp) == 5: status, repository, name, path, file = temp
    if not proposed and status == 'P':
      command = None
      return
    origin_path = replace_dot_with_geant4(path)
    nakedpath = remove_dot_slash(path) 
    if file:
      if os.path.exists(path+"/"+file):
        command = "svn switch %s %s/%s/tags/%s/_symbols/%s/%s %s"%(options, cern_svn_repos, repository, origin_path, name, file, path+"/"+file)
      else:
        if not os.path.exists(path): os.makedirs(path)
        command = "svn cp %s %s/%s/tags/%s/_symbols/%s/%s %s"%(options, cern_svn_repos, repository, origin_path, name, file, path+"/"+file)
      if nakedpath+"/"+file in swipaths: swipaths.remove(nakedpath+"/"+file)
      if nakedpath+"/"+file in newpaths: newpaths.remove(nakedpath+"/"+file)
    else:
      if os.path.exists(path) :
        if  get_tag(path) == name :
          command = None
          #command = "svn update %s %s" % (options, path)  
        else:
          command = "svn switch %s %s/%s/tags/%s/_symbols/%s %s"%(options, cern_svn_repos, repository, origin_path, name, path)
      else:
        os.makedirs(path)
        command = "svn co %s %s/%s/tags/%s/_symbols/%s %s"%(options, cern_svn_repos, repository, origin_path, name, path)
      if nakedpath in swipaths: swipaths.remove(nakedpath)
      if nakedpath in newpaths: newpaths.remove(nakedpath)
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
        print '# ',command
        sys.stdout.flush()
        if os.system(command) == 1: return 1
  #---Check if anything needs to be deleted----------------------------------
  for path in newpaths:
    print '#  rmdir ', path
    if os.path.isdir(path): shutil.rmtree(path)
    elif os.path.isfile(path): os.remove(path)
  for path in swipaths:
    #print '# Need to do something about ', path, '  switching to parent URL'
    parent, base = os.path.split(path)
    url = re.search('^URL: +([^ \n]+)', os.popen('svn info %s' % parent).read(), 8).group(1)
    command = "svn switch %s %s/%s %s"%(options, url, base, path)
    print '# ', command
    os.system(command) 

if __name__ == "__main__":
  import optparse
  parser = optparse.OptionParser(description='Geant4 checkout utility using the tag collector')
  parser.add_option("-c", "--config", dest="config", default='g4tags-dev',
                    help="Geant4 checkout configuration")
  parser.add_option("-d", "--dest", dest="dest", default=None,
                    help="checkout directory destination")
  parser.add_option("-q", "--quiet", action='store_true', dest="quiet", default=False,
                    help="don't print status messages to stdout")
  parser.add_option("-p", "--proposed", action='store_true', dest="proposed", default=False,
                    help="add the proposed tags in addtion")
  parser.add_option("", "--non-interactive", action='store_true', dest="interactive", default=False,
                    help="to please CTest!!!")
  
  (options, args) = parser.parse_args()
  if(not options.dest) : options.dest = options.config 
  sys.exit(g4svn_update(options.config, options.dest, options.quiet, options.proposed))
