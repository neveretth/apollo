source: $(BUILD_DIR)/src/types.c.o $(BUILD_DIR)/src/main.c.o $(BUILD_DIR)/src/args.c.o $(BUILD_DIR)/src/display.c.o

source-rocm: source $(BUILD_DIR)/src/rocm/hip-util.c.o 
	
$(BUILD_DIR)/src/types.c.o: $(SOURCE_DIR)/types.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/main.c.o: $(SOURCE_DIR)/main.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/args.c.o: $(SOURCE_DIR)/args.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)

$(BUILD_DIR)/src/display.c.o: $(SOURCE_DIR)/display.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
	
$(BUILD_DIR)/src/rocm/hip-util.c.o: $(SOURCE_DIR)/hip-util.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)

include src/kernel/Makefile
include src/driver/Makefile
include src/types/Makefile
include src/parser/Makefile
