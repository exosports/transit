# Makefile - created for transit.
#
# This Makefile calls subsidiary Makefiles located in the various modules that
# make up Transit. Those Makefiles should be maintained as each module evolves.
#
# Revision    July 6, 2015    AJ Foster
#             Created; naively calls module make commands.


# Get the location of this Makefile.
mkfile_dir := $(dir $(lastword $(MAKEFILE_LIST)))

# `make [clean]` should run `make [clean]` on all of the modules.
all: make_pu make_transit make_tips
clean: clean_pu clean_transit clean_tips

# The following builds the various modules using their respective `make`
# tasks.
#

make_pu:
	@echo "\nCompiling module pu..."
	@cd $(mkfile_dir)/pu/ && make
	@echo "Finished compiling pu."

make_transit: make_pu
	@echo "\nCompiling module transit..."
	@cd $(mkfile_dir)/transit/ && make
	@echo "Finished compiling transit."

make_tips:
	@echo "\nCompiling PYTIPS for pylineread..."
	@cd $(mkfile_dir)/pylineread/src/pytips/ && make
	@echo "Finished compiling PYTIPS."


# The following tasks clean the various modules using their respective
# `make clean` commands.
#

clean_pu:
	@echo "\nCleaning module pu..."
	@cd $(mkfile_dir)/pu/ && make clean
	@echo "Finished cleaning pu."

clean_transit:
	@echo "\nCleaning module transit..."
	@cd $(mkfile_dir)/transit/ && make clean
	@echo "Finished cleaning transit."

clean_tips:
	@echo "\nCleaning PYTIPS for pylineread..."
	@cd $(mkfile_dir)/pylineread/src/pytips/ && make clean
	@echo "Finished cleaning PYTIPS."
