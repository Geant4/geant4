INSTALL_DIR = "/usr/local/javaclass"

JAVAC       = javac

g4vrmlview.class: g4vrmlview.java
	$(JAVAC) g4vrmlview.java
	$(JAVAC) g4mini.java

install: 
	@if [ ! -d $(INSTALL_DIR) ] ; then mkdir $(INSTALL_DIR) ;fi	
	cp g4vrmlview.class        $(INSTALL_DIR)
	cp g4vrmlviewThread.class  $(INSTALL_DIR)

clean:
	rm -f *.class
	rm -f \#*  *~ *.bak

wrlclean:
	rm -f g4.wrl  g4_*.wrl

