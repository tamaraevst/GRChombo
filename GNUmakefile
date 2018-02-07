TestDirs := $(wildcard Tests/*/.)
ExampleDirs := $(wildcard Examples/*/.)
CleanTestDirs := $(TestDirs:%=clean-%)
CleanExampleDirs := $(ExampleDirs:%=clean-%)
RealCleanTestDirs := $(TestDirs:%=realclean-%)
RealCleanExampleDirs := $(ExampleDirs:%=realclean-%)

.PHONY: all run $(TestDirs) $(ExampleDirs)

test: $(TestDirs)

examples: $(ExampleDirs)

all: $(TestDirs) $(ExampleDirs)

clean: $(CleanTestDirs) $(CleanExampleDirs)

realclean: $(RealCleanTestDirs) $(RealCleanExampleDirs)

$(TestDirs):
	$(info ################# Making test $@ #################)
	$(MAKE) -C $@ all
	$(info ################# Running test $@ #################)
	$(MAKE) -C $@ run

$(ExampleDirs):
	$(info ################# Making example $@ #################)
	$(MAKE) -C $@ all

$(CleanTestDirs):
	$(MAKE) -C $(@:clean-%=%) clean

$(CleanExampleDirs):
	$(MAKE) -C $(@:clean-%=%) clean

$(RealCleanTestDirs):
	$(MAKE) -C $(@:realclean-%=%) clean

$(RealCleanExampleDirs):
	$(MAKE) -C $(@:realclean-%=%) clean
