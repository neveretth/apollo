UNIT_TEST_SOURCE_DIR=tests

unit-test: builddir sample-test test-lib
	$(CC) $(BUILD_DIR)/tests/*.o -o $(BUILD_DIR)/unit-test $(CLIBS)
	
include tests/lib/Makefile

include tests/sample/Makefile
