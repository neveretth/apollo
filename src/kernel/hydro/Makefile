kernel-hydro: $(BUILD_DIR)/src/kernel-hydro.c.o
	@
	
kernel-rocm-hydro: kernel-hydro $(BUILD_DIR)/src/kernel-hydro.hip.o
	@

$(BUILD_DIR)/src/kernel-hydro.c.o: $(SOURCE_DIR)/kernel/hydro/kernel.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/kernel-hydro.hip.o: $(SOURCE_DIR)/kernel/hydro/kernel.hip
	$(HIPCC) -o $@ -c $^
