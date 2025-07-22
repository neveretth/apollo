kernel-neutrino: $(BUILD_DIR)/src/kernel-neutrino.c.o
	@
	
kernel-rocm-neutrino: kernel-neutrino $(BUILD_DIR)/src/rocm/kernel-neutrino.hip.o
	@

$(BUILD_DIR)/src/kernel-neutrino.c.o: $(SOURCE_DIR)/kernel/neutrino/kernel.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/rocm/kernel-neutrino.hip.o: $(SOURCE_DIR)/kernel/neutrino/kernel.hip
	$(HIPCC) -o $@ -c $^
