# SEQAN Development Build System
# Platform Adapter Makefile for SUN
#
# 2003/2004 by Andreas Doering
# doering@inf.fu-berlin.de
#________________________________________________________________________
#

extEXE := 

# Print compiler version
define cmdVersion
	@/opt/SUNWspro/bin/CC -V
endef

# Create timestamp file $(1)
define cmdTouch
	-@touch $(1)
endef

# Make directory $(1)
define cmdMkdir
	-@mkdir $(1) 2> /dev/null
endef

# Create forwards $(1), set optional argument $(2) to "all" for rebuilding all forwards
define cmdCreateForwards
    @python misc/build_forwards.py "projects/library/seqan/$(1)" $(2)|| true
endef


# Compile $(1) to output file $(2) including files $(3)
define cmdCompile
	@echo $(Indent)compile $(1)
	@/opt/SUNWspro/bin/CC -c $(1) -fast -o "$(2)" $(addprefix -I,$(3))
endef

# Link $(1) to output file $(2)
define cmdLink
	@echo $(Indent)link $(2)
	@/opt/SUNWspro/bin/CC $(1) -o "$(2)" 
endef

# Delete directories $(1) and files $(2)
define cmdClean 
	@rm -f $(2)
	@rm -r -f $(1)
endef

# Runs the application $(1)
define cmdRun 
	@$(1)
endef

# Create documentation
define cmdDoc
	@cd docs;\
	python main.py ../projects/library;\
	cd ..
endef

# Create snapshot
define cmdSnapshot
	FILENAME=Seqan_snap_`date +%Y%m%d`.zip;\
	rm -f $$FILENAME;\
	mkdir -p misc/snapshot;\
	cd projects/library/demos;\
	make clean;\
	cd ../../../misc/vs_demo_projectfile;\
	python vs_demo_projectfile.py;\
	cd ../vs_apps_projectfile;\
    python vs_apps_projectfile.py;\
	cd ../snapshot;\
	ln -s ../../projects/library/demos;\
	ln -s ../../projects/library/seqan;\
	ln -s ../../projects/library/apps;\
	ln -s ../../projects/library/cmake;\
	ln -s ../../docs/html doc;\
	echo Create $$FILENAME ...;\
	zip -qr ../../$$FILENAME * -x \*.svn\*;\
	rm demos seqan apps cmake doc;\
	cd ../..
endef
