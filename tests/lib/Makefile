test-lib: $(BUILD_DIR)/tests/main.c.o $(BUILD_DIR)/tests/lib.c.o
	
$(BUILD_DIR)/tests/main.c.o: $(UNIT_TEST_SOURCE_DIR)/lib/main.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/tests/lib.c.o: $(UNIT_TEST_SOURCE_DIR)/lib/lib.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
