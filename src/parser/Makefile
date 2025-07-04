parser: $(BUILD_DIR)/src/rate-library-parser.c.o $(BUILD_DIR)/src/tnn-parameters.c.o $(BUILD_DIR)/src/tnn-parser.c.o $(BUILD_DIR)/src/neunet-parser.c.o $(BUILD_DIR)/src/validate.c.o
	@
	
$(BUILD_DIR)/src/neunet-parser.c.o: $(SOURCE_DIR)/parser/neunet-parser.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)

$(BUILD_DIR)/src/rate-library-parser.c.o: $(SOURCE_DIR)/parser/rate-library-parser.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)

$(BUILD_DIR)/src/tnn-parameters.c.o: $(SOURCE_DIR)/parser/tnn-parameters.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)

$(BUILD_DIR)/src/tnn-parser.c.o: $(SOURCE_DIR)/parser/tnn-parser.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)

$(BUILD_DIR)/src/validate.c.o: $(SOURCE_DIR)/parser/validate.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(CLIBS)
