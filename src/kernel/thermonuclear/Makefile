kernel-thermonuclear: $(BUILD_DIR)/src/kernel-thermonuclear.c.o
	@
	
kernel-rocm-thermonuclear: kernel-thermonuclear $(BUILD_DIR)/src/rocm/kernel-thermonuclear.hip.o
	@

$(BUILD_DIR)/src/kernel-thermonuclear.c.o: $(SOURCE_DIR)/kernel/thermonuclear/kernel.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/rocm/kernel-thermonuclear.hip.o: $(SOURCE_DIR)/kernel/thermonuclear/kernel.hip
	$(HIPCC) -o $@ -c $^
