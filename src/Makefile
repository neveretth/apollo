source: $(BUILD_DIR)/src/types.c.o $(BUILD_DIR)/src/hip-util.c.o $(BUILD_DIR)/src/main.c.o $(BUILD_DIR)/src/parameters.c.o $(BUILD_DIR)/src/parse-data.c.o $(BUILD_DIR)/src/args.c.o 
	
$(BUILD_DIR)/src/types.c.o: $(SOURCE_DIR)/types.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/hip-util.c.o: $(SOURCE_DIR)/hip-util.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/main.c.o: $(SOURCE_DIR)/main.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/parameters.c.o: $(SOURCE_DIR)/parameters.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/parse-data.c.o: $(SOURCE_DIR)/parse-data.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/args.c.o: $(SOURCE_DIR)/args.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)

include src/kernel/Makefile
include src/driver/Makefile
include src/types/Makefile
