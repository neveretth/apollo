types: $(BUILD_DIR)/src/hydro-types.c.o $(BUILD_DIR)/src/neutrino-types.c.o $(BUILD_DIR)/src/tnn-types.c.o
	@

$(BUILD_DIR)/src/hydro-types.c.o: $(SOURCE_DIR)/types/hydro-types.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/neutrino-types.c.o: $(SOURCE_DIR)/types/neutrino-types.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/tnn-types.c.o: $(SOURCE_DIR)/types/tnn-types.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
