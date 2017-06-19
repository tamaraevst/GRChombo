TestDirs := $(wildcard Tests/*/.)
ExampleDirs := $(wildcard Examples/*/.)

.PHONY: all run $(TestDirs) $(ExampleDirs)

test: $(TestDirs)

all: $(TestDirs) $(ExampleDirs)

$(TestDirs):
	$(info ################# Making test $@ #################)
	$(MAKE) -C $@ all
	$(info ################# Running test $@ #################)
	$(MAKE) -C $@ run

$(ExampleDirs):
	$(info ################# Making example $@ #################)
	$(MAKE) -C $@ all
