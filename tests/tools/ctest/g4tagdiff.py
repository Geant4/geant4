import urllib, difflib, sys, os

def g4tagdiff(args, proposed, devline) :

  if len(args) == 0 or not os.path.exists(args[0]) : o_taglist = []
  else : o_taglist = open(args[0], 'r').readlines()
  
  if proposed : tags_url = "http://lcgapp.cern.ch/spi/cgi-bin/g4tags.py?devline=%s;proposed=true" % devline
  else        : tags_url = "http://lcgapp.cern.ch/spi/cgi-bin/g4tags.py?devline=%s" % devline
  n_taglist = urllib.urlopen(tags_url).readlines()

  diff = False
  d = difflib.Differ()
  diffs = list(d.compare(o_taglist,n_taglist))
  for l in diffs:
    if len(l) > 3 and l[0:2] in ('+ ', '- ') and l[2] != '#' :
       print l[:-1]
       diff = True
  
  if diff : exit(1)
  else : exit(0)

if __name__ == "__main__":
  import optparse, string
  parser = optparse.OptionParser(description='Geant4 Tag Diff utility')
  parser.add_option('-q', '--quiet', action='store_true', dest='quiet', default=False,
                      help='don\'t print status messages to stdout')
  parser.add_option('-p', '--proposed', action='store_true', dest='proposed', default=False,
                      help='add the proposed tags in addtion')
  parser.add_option('-d', '--devline', dest='devline', default='g4tags-dev',
                      help='development line or version')
  (options, args) = parser.parse_args()
  g4tagdiff(args, options.proposed, options.devline)

