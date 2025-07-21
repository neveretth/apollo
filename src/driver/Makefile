driver: $(BUILD_DIR)/src/kernel-driver.c.o $(BUILD_DIR)/src/unified-driver.c.o
	@
	
driver-rocm: $(BUILD_DIR)/src/rocm/kernel-driver.c.o

$(BUILD_DIR)/src/unified-driver.c.o: $(SOURCE_DIR)/driver/unified-driver.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/kernel-driver.c.o: $(SOURCE_DIR)/driver/kernel-driver.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)

$(BUILD_DIR)/src/rocm/kernel-driver.c.o: $(SOURCE_DIR)/driver/rocm/kernel-driver.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
