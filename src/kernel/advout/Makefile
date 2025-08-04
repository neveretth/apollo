kernel-advout: $(BUILD_DIR)/src/kernel-advout.c.o
	@
	
# kernel-rocm-advout: kernel-advout $(BUILD_DIR)/src/kernel-advout.hip.o
# 	@

$(BUILD_DIR)/src/kernel-advout.c.o: $(SOURCE_DIR)/kernel/advout/kernel.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
# $(BUILD_DIR)/src/kernel-advout.hip.o: $(SOURCE_DIR)/kernel/advout/kernel.hip
# 	$(HIPCC) -o $@ -c $^
