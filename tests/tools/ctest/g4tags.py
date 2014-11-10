#!/usr/bin/env python

# This script is to be installed in /afs/cern.ch/sw/lcg/app/nightlies/scripts
# The URL is http://lcgapp.cern.ch/spi/cgi-bin/g4tags.py?devline=g4tags-dev

import MySQLdb
import time, sys, re, os

host     = 'dbod-g4-tags.cern.ch'
port     = 5500
user     = 'g4tagsro'
passwd   = 'read'
dbname   = 'geant4tags'
slot     = 'g4tags-dev'

specialcases = ({'regex':'.*-gmk-.*', 'switch':'GNUmakefile', 'path':''},
                {'regex':'.*-cmk-.*', 'switch':'CMakeLists.txt', 'path':''},
                {'regex':'benchmarks-V.*', 'switch':'', 'path':'geant4/benchmarks'},
                {'regex':'Configure-.*', 'switch':'Configure', 'path':''})
smallstatusProposed = {'selected':'S', 'proposed':'P', 'accepted':'A'}
smallstatusSelected = {'selected':'S', 'accepted':'A'}
smallstatusAccepted = {'accepted':'A'}
updatesections = ('Configure', 'config', 'examples', 'source', 'tests','environments', 'cmake')

#-----------------------------------------------------------------------------------------------
def gettaglist(slot, proposed=False):
  db = MySQLdb.connect( host=host, port=port, user=user, passwd=passwd, db=dbname)
  cursor = db.cursor(MySQLdb.cursors.DictCursor)
  cursor.execute("select * from g4t_devlines where slot='%s' and status=1" % slot )
  devline = cursor.fetchone()
  if not devline:
     print "Error: Development line '%s' not found. Stopping" % slot
     return ''
  s =  "# Date: %s %s\n" % (time.asctime(), time.tzname[0])
  s += "# Development line #%d: %s on %s slot %s.\n" %(devline['dlid'],devline['name'],dbname, slot)

  s += "\n\nCHECKOUT SECTION\n"
  s += "geant4 %s\n" % devline['name']

  s += "\n\nROOT UPDATE SECTION\n"
  for sec in updatesections: 
    s +=  "./%s\n" % sec

  s += "\n\nSPECIAL CHECKOUT SECTION\n"
  cursor.execute("""
    SELECT t.tid, t.name, r.repository, p.path  
    FROM g4t_tags t 
    JOIN g4t_paths p ON t.path_id=p.pid
    JOIN g4t_repositories r ON t.repository_id=r.rid
    WHERE t.name LIKE 'benchmarks-%'""")
  btag = cursor.fetchall()[-1]
  s += "%-10s %-40s ./%-50s %s\n"%( btag['repository'],btag['name'],btag['path'], btag['path'])

  s += "\n\nCATEGORIES SECTION\n"
  #---Select all tags associated to a develment line (slot)------------------------------------
  cursor.execute("""
    SELECT t.tid, tl.tlid, t.name, a.author, r.repository, p.path, t.date, ts.status, t.sentence, t.bugfix, t.description 
    FROM g4t_tags t 
    JOIN g4t_paths p ON t.path_id = p.pid
    JOIN g4t_repositories r ON t.repository_id = r.rid
    JOIN g4t_tag_statuses ts ON t.status_id = ts.tsid
    JOIN g4t_authors a ON t.author_id = a.aid 
    JOIN g4t_taglist_tags tlt ON t.tid=tlt.tag_id 
    JOIN g4t_taglists tl ON tlt.taglist_id=tl.tlid 
    JOIN g4t_devlines dl ON tl.devline_id=dl.dlid 
    WHERE dl.slot='%s' and dl.status=1 and t.name NOT LIKE 'benchmarks-%%'
    ORDER BY t.date ASC """ % slot)
  tags = cursor.fetchall()
  #---Select all the tags that are in state 'proposed'-----------------------------------------
  if proposed :
    cursor.execute("""
      SELECT t.tid, t.name, a.author, r.repository, p.path, t.date, ts.status, t.sentence, t.bugfix, t.description 
      FROM g4t_tags t 
      JOIN g4t_paths p ON t.path_id = p.pid
      JOIN g4t_repositories r ON t.repository_id = r.rid
      JOIN g4t_tag_statuses ts ON t.status_id = ts.tsid
      JOIN g4t_authors a ON t.author_id = a.aid 
      WHERE ts.status='proposed'
      ORDER BY t.date ASC """)
    tags += cursor.fetchall()
  #---Fill the file attribute -----------------------------------------------------------------
  for tag in tags:
    tag['file'] = ''
    for case in specialcases:
      if re.match(case['regex'],tag['name']) : tag['file'] = case['switch']
  #---Remove path duplicates and sort by date-------------------------------------------------- 
  pathes = {}
  for tag in tags : 
    if tag['status'] not in smallstatus.keys(): continue
    key = os.path.join(tag['path'],tag['file'])
    #---Go back and remove any tag that has been superseeded----------------------------------
    for oldkey in pathes.keys() :
      if key in oldkey : del pathes[oldkey]
    #---Add current tag with key----------------------------------------------------------------
    pathes[key] = tag
  testing_tags = pathes.values()
  testing_tags.sort(key=lambda t:t['date'])
  #---Loop over tags and format correctly the outpout------------------------------------------
  for tag in testing_tags:
    s +=  "%s %-8s %-40s %-50s %s\n"%(smallstatus[tag['status']], tag['repository'], tag['name'], tag['path'].replace(tag['repository'],'.'), tag['file'])
  #---Additional details of the tag lists------------------------------------------------------
  s +=  "\n\n# Details about the taglists assigned to tested development line: \n\n"
  s +=  "# %-8s %-8s %-20s %-10s %-40s %-10s %s\n" % ('Tag id','TL id','Date','Author','Name','Status','Path') 
  lasttlid = 0
  for tag in tags:
    if tag.get('tlid',-1) != lasttlid : s +=  "#  \n"; lasttlid = tag.get('tlid',-1)
    s +=  "# %-8d %-8d %-20s %-10s %-40s %-10s %s\n" % (tag['tid'],tag.get('tlid',-1),tag['date'],tag['author'],tag['name'],tag['status'],tag['path']) 
  return s

if __name__ == '__main__' :
  try :
    import cgi
    form = cgi.FieldStorage()
    print "Content-type: text/plain\n\n"
    if form.has_key('devline') :  devline = form['devline'].value
    else:                         devline = 'g4tags-dev'
    if form.has_key('proposed') : 
       proposed = form['proposed'].value.lower() in ('true','ok','1','true','yes')
       smallstatus = smallstatusProposed
    elif form.has_key('accepted'):
       accepted = form['accepted'].value.lower() in ('true','ok','1','true','yes')
       smallstatus = smallstatusAccepted
    else:
       proposed = False
       smallstatus = smallstatusSelected
    print gettaglist(devline, proposed)  
  except :
    print gettaglist(len(sys.argv) > 1 and sys.argv[1] or 'g4tags-dev')






