MARSHALGEN_PATH = /afs/cern.ch/user/c/cooperma/vietha/Marshalgen
MARSHAL=$(MARSHALGEN_PATH)/marshalgen

default: all

all:	G4String G4THits G4VHits ExN02TrackerHit G4HCofThisEvent
	(cp -p ../AnnotatedFiles/Marshaled* .; cp -p $(MARSHALGEN_PATH)/MarshaledObj.h .)


ExN02TrackerHit: ExN02TrackerHit.hh
	$(MARSHAL) $?

G4String: ../AnnotatedFiles/G4String.hh
        $(MARSHAL) $?

G4HCofThisEvent: ../AnnotatedFiles/G4HCofThisEvent.hh
	$(MARSHAL) $?

G4THits: ../AnnotatedFiles/G4THitsCollection.hh
	$(MARSHAL) $?

G4VHits: ../AnnotatedFiles/G4VHitsCollection.hh
	$(MARSHAL) $?

clean:
	rm *.msh Marshaled* ../AnnotatedFiles/*.msh ../AnnotatedFiles/Marshaled*
