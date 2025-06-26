kernel: $(BUILD_DIR)/src/kernel.c.o $(BUILD_DIR)/src/kernel.hip.o
	@

$(BUILD_DIR)/src/kernel.c.o: $(SOURCE_DIR)/kernel/kernel.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/kernel.hip.o: $(SOURCE_DIR)/kernel/kernel.hip
	$(HIPCC) -o $@ -c $^
