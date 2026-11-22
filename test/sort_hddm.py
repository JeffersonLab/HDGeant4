#!/usr/bin/env python3

import sys
import hddm_s
import hddm_r

max_eventNo = 10000000000
max_readbuf = 1000
eventNo = 1
xmlout = 0
clean = 0
rest = 0

readbuf = {}

narg = 1
while narg < len(sys.argv):
   if sys.argv[narg] == "-n":
      max_eventNo = int(sys.argv[narg + 1])
      narg += 2
      print(f"stopping after {max_eventNo} events")
   elif sys.argv[narg] == "-x":
      xmlout = 1
      narg += 1
   elif sys.argv[narg] == "-c":
      clean = 1
      narg += 1
   elif len(sys.argv[narg]) > 2 and sys.argv[narg][:2] == "-n":
      max_eventNo = int(sys.argv[narg][2:])
      print(f"stopping after {max_eventNo} events")
      narg += 1
   else:
      break

def getPhysicsEvents(rec):
   if rest == 0:
      return rec.getPhysicsEvents()
   else:
      return rec.getReconstructedPhysicsEvents()

for fin in sys.argv[narg:]:
   finhead = open(fin, 'r').readline()
   if '<HDDM class="s" ' in finhead:
      rest = 0
      fini = hddm_s.istream(fin)
      fout = hddm_s.ostream("sort_hddm.hddm")
      if xmlout > 0:
         xout = open("sort_hddm.xml", 'w')
   elif '<HDDM class="r" ' in finhead:
      rest = 1
      fini = hddm_r.istream(fin)
      fout = hddm_r.ostream("sort_hddm.hddm")
      if xmlout > 0:
         xout = open("sort_hddm.xml", 'w')
   else:
      print(f"input file {fin} has unknown data format, skipping...")
      continue
   for rec in fini:
      for pev in getPhysicsEvents(rec):
         if clean > 0:
            pev.deleteHitViews()
         if pev.eventNo == eventNo:
            if clean > 0:
               for rea in pev.getReactions():
                  rea.deleteVertices(-1,1)
            fout.write(rec)
            if xmlout > 0:
               xout.write(rec.toXML())
            if eventNo % 100 == 0:
               print("wrote event", eventNo)
            eventNo += 1
            while eventNo in readbuf and eventNo < max_eventNo:
               fout.write(readbuf[eventNo])
               if xmlout > 0:
                  xout.write(readbuf[eventNo].toXML())
               if eventNo % 100 == 0:
                  print("wrote event", eventNo)
               del readbuf[eventNo]
               eventNo += 1
            if eventNo == max_eventNo:
               sys.exit(0)
         else:
            readbuf[pev.eventNo] = rec
            if len(readbuf) > max_readbuf:
               print("readbuf backed up with", len(readbuf), "waiting events")
               print("giving up, sorry.")
               sys.exit(1)
