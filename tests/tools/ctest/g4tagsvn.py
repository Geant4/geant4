#!/usr/bin/env python
import sys, os, shutil,re
import subprocess as sb
import string

if sys.platform == 'win32':
  cern_svn_repos = 'https://svn.cern.ch/reps'
elif sys.platform == 'darwin':
  cern_svn_repos = 'svn+ssh://svn.cern.ch/reps'
else :
  cern_svn_repos = 'svn+ssh://svn.cern.ch/reps'

#--ignore-ancestry must be used with svn version 1.7 or higher
svn_version=sb.Popen(["svn", "--version","--quiet"],stdout=sb.PIPE ).communicate()[0]
svnver=string.split(svn_version,'.')
if int(svnver[0]) > 1 or int(svnver[1]) > 6 :
  extra_options = '--ignore-ancestry'
else:
  extra_options = ''
  
replace_dot_with_geant4 = lambda str:str.replace(".", "geant4",1)
remove_dot_slash = lambda str:str.replace("./", "",1)
def get_tag(path):
  r = re.search('^URL: +.*_symbols/([^ \n/]+)', os.popen('svn info %s' % path).read(), 8)
  if r : return r.group(1)
  else : return ''

#--------------------------------------------------------------------------------------------------
def g4svn_update(devline, destdir, quiet, proposed , accepted):
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
  if devline.endswith("_branch") or devline.startswith("geant4-"):
    svnbranch="tags"
    if devline.endswith("_branch"):
       svnbranch="branches"
    if not os.path.exists('.svn'):
      command = "svn co %s %s/geant4/%s/geant4/_symbols/%s ." % (options, cern_svn_repos, svnbranch, devline)
    else:
      command = "svn switch %s %s/geant4/%s/geant4/_symbols/%s ." % (options, cern_svn_repos, svnbranch, devline)
    print '# ', command
    sys.stdout.flush()
    os.system(command)
    #---Get now the benchmarks
    benchmarks_tag=""
    if 'VERSION_BENCHMARKS' in os.environ:
      benchmarks_tag = os.environ['VERSION_BENCHMARKS']
    if len(benchmarks_tag) == 0:
      benchmarks_tag= os.popen('svn ls %s/g4tests/tags/benchmarks/_symbols' % cern_svn_repos).readlines()[-1].strip()
    if os.path.exists('benchmarks') :
      command="svn switch %s %s/g4tests/tags/benchmarks/_symbols/%s benchmarks" % (options, cern_svn_repos, benchmarks_tag)
    else:
      command="svn co %s %s/g4tests/tags/benchmarks/_symbols/%s benchmarks" % (options, cern_svn_repos, benchmarks_tag)
    print '# ', command
    sys.stdout.flush()
    os.system(command)
    return
  #---Analyse the new and switched directories and files -------------------------------------------
  newpaths = []
  swipaths =[]
  for l in os.popen('svn status').readlines():
    # ' ' no modifications
    # 'A' Added
    # 'C' Conflicted
    # 'D' Deleted
    # 'I' Ignored
    # 'M' Modified
    # 'R' Replaced
    # 'X' an unversioned directory created by an externals definition
    # '?' item is not under version control
    # '!' item is missing (removed by non-svn command) or incomplete
    # '~' versioned item obstructed by some item of a different kind
    if l[0] not in ' ?' : continue
    stat, path = l.split()
    if stat == '?' and os.path.isdir(path) : newpaths.append(path.replace('\\','/'))
    elif stat == 'S' : swipaths.append(path.replace('\\','/'))
    
  #---Get the current list of tags ------------------------------------------------------------------
  import urllib
  if proposed : tags_url = "http://lcgapp.cern.ch/spi/cgi-bin/g4tags.py?devline=%s;proposed=true" % devline
  elif accepted :  tags_url = "http://lcgapp.cern.ch/spi/cgi-bin/g4tags.py?devline=%s;accepted=true" % devline
  else        : tags_url = "http://lcgapp.cern.ch/spi/cgi-bin/g4tags.py?devline=%s" % devline
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
        if get_tag(path+"/"+file) == name :
          command = None
        else: 
          command = "svn switch %s %s %s/%s/tags/%s/_symbols/%s/%s %s"%(options, extra_options, cern_svn_repos, repository, origin_path, name, file, path+"/"+file)
      else:
        if not os.path.exists(path): os.makedirs(path)
        command = "svn cp %s %s/%s/tags/%s/_symbols/%s/%s %s"%(options, cern_svn_repos, repository, origin_path, name, file, path+"/"+file)
      nakedfile = remove_dot_slash(nakedpath+"/"+file)       
      if nakedfile in swipaths: swipaths.remove(nakedfile)
      if nakedfile in newpaths: newpaths.remove(nakedfile)
    else:
      if os.path.exists(path) :
        if  get_tag(path) == name :
          command = None
        else:
          command = "svn switch %s %s %s/%s/tags/%s/_symbols/%s %s"%(options, extra_options, cern_svn_repos, repository, origin_path, name, path)
      else:
        os.makedirs(path)
        command = "svn co %s %s/%s/tags/%s/_symbols/%s %s"%(options, cern_svn_repos, repository, origin_path, name, path)
      if nakedpath in swipaths: swipaths.remove(nakedpath)      
      #if nakedpath in newpaths: newpaths.remove(nakedpath)
      p = nakedpath
      while(p) : 
        if p in newpaths: newpaths.remove(p)
        p = os.path.dirname(p)
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
        rc = os.system(command)
        if rc != 0: return rc
  #---Check if anything needs to be deleted----------------------------------
  for path in newpaths:
    print '#  rmdir ', path
    if os.path.isdir(path): shutil.rmtree(path)
    elif os.path.isfile(path): os.remove(path)
  for path in swipaths:
    #print '# Need to do something about ', path, '  switching to parent URL'
    parent, base = os.path.split(path)
    url = re.search('^URL: +([^ \n]+)', os.popen('svn info %s' % parent).read(), 8).group(1)
    command = "svn switch %s %s %s/%s %s"%(options, extra_options, url, base, path)
    print '# ', command
    os.system(command)

"""
  Implementation of a simple cross-platform file locking mechanism.
  This is a modified version of code retrieved on 2013-01-01 from http://www.evanfosmark.com/2009/01/cross-platform-file-locking-support-in-python.
  The original code was released under the BSD License, as is this modified version.
  """
import os
import sys
import time
import errno

class FileLock(object):
  """ A file locking mechanism that has context-manager support so
    you can use it in a ``with`` statement. This should be relatively cross
    compatible as it doesn't rely on ``msvcrt`` or ``fcntl`` for the locking.
    """
  
  class FileLockException(Exception):
    pass
  
  def __init__(self, protected_file_path, timeout=None, delay=1, lock_file_contents=None):
    """ Prepare the file locker. Specify the file to lock and optionally
      the maximum timeout and the delay between each attempt to lock.
      """
    self.is_locked = False
    self.lockfile = os.path.abspath(protected_file_path + ".lock")
    self.timeout = timeout
    self.delay = delay
    self._lock_file_contents = lock_file_contents
    if self._lock_file_contents is None:
      self._lock_file_contents = "Owning process args:\n"
      for arg in sys.argv:
        self._lock_file_contents += arg + "\n"
  
  def locked(self):
    """
      Returns True iff the file is owned by THIS FileLock instance.
      (Even if this returns false, the file could be owned by another FileLock instance, possibly in a different thread or process).
      """
    return self.is_locked
  
  def available(self):
    """
      Returns True iff the file is currently available to be locked.
      """
    return not os.path.exists(self.lockfile)
  
  def acquire(self, blocking=True):
    """ Acquire the lock, if possible. If the lock is in use, and `blocking` is False, return False.
      Otherwise, check again every `self.delay` seconds until it either gets the lock or
      exceeds `timeout` number of seconds, in which case it raises an exception.
      """
    start_time = time.time()
    while True:
      try:
        # Attempt to create the lockfile.
        # These flags cause os.open to raise an OSError if the file already exists.
        fd = os.open( self.lockfile, os.O_CREAT | os.O_EXCL | os.O_RDWR )
        f = os.fdopen( fd, 'a' )
        try:
          # Print some info about the current process as debug info for anyone who bothers to look.
          f.write( self._lock_file_contents )
        finally:
          f.close()
        break;
      except OSError, e:
        if e.errno != errno.EEXIST:
          raise
        if self.timeout is not None and (time.time() - start_time) >= self.timeout:
          raise FileLock.FileLockException("Timeout occurred.")
        if not blocking:
          return False
        time.sleep(self.delay)
    self.is_locked = True
    return True

  def release(self):
    """ Get rid of the lock by deleting the lockfile.
      When working in a `with` statement, this gets automatically
      called at the end.
      """
    self.is_locked = False
    os.unlink(self.lockfile)
  
  
  def __enter__(self):
    """ Activated when used in the with statement.
      Should automatically acquire a lock to be used in the with block.
      """
    self.acquire()
    return self
  
  
  def __exit__(self, type, value, traceback):
    """ Activated at the end of the with statement.
      It automatically releases the lock if it isn't locked.
      """
    self.release()
  
  
  def __del__(self):
    """ Make sure this ``FileLock`` instance doesn't leave a .lock file
      lying around.
      """
    if self.is_locked:
      self.release()
  
  def purge(self):
    """
      For debug purposes only.  Removes the lock file from the hard disk.
      """
    if os.path.exists(self.lockfile):
      self.release()
      return True
    return False


if __name__ == "__main__":
  import optparse, string
  if sys.argv[1] in ('update','up', 'checkout', 'co') :
    parser = optparse.OptionParser(description='Geant4 checkout utility using the tag collector')
    parser.add_option("-c", "--config", dest="config", default='g4tags-dev',
                      help="Geant4 checkout configuration")
    parser.add_option("-d", "--dest", dest="dest", default=None,
                      help="checkout directory destination")
    parser.add_option("-q", "--quiet", action='store_true', dest="quiet", default=False,
                      help="don't print status messages to stdout")
    parser.add_option("-p", "--proposed", action='store_true', dest="proposed", default=False,
                      help="add the proposed tags in addtion")
    parser.add_option("-a", "--accepted", action='store_true', dest="accepted", default=False,
                      help="add *ONLY* accepted tags")
    parser.add_option("", "--non-interactive", action='store_true', dest="interactive", default=False,
                      help="to please CTest!!!")
    parser.add_option("-r", "--revision", dest="revision", default=None,
                      help="to please CTest!!!")
    (options, args) = parser.parse_args()
    if(not options.dest) : options.dest = options.config
    #---We need to prevent several processed updating the same repository. For this we use the FileLock class
    if not os.path.exists(options.dest): os.makedirs(options.dest)
    fl = FileLock(os.path.join(options.dest,'tagsvn_update'), timeout=1000)
    try:
      fl.acquire()
      g4svn_update(options.config, options.dest, options.quiet, options.proposed, options.accepted)
    finally:
      fl.release()
  else:
    os.system("svn %s" % string.join(sys.argv[1:]))
