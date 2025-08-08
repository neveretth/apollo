sample-test:  $(BUILD_DIR)/tests/logic.c.o $(BUILD_DIR)/tests/kernel.c.o
	
$(BUILD_DIR)/tests/logic.c.o: $(UNIT_TEST_SOURCE_DIR)/sample/logic.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
# This will not exist in normal test dirs.
$(BUILD_DIR)/tests/kernel.c.o: $(UNIT_TEST_SOURCE_DIR)/sample/kernel.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
