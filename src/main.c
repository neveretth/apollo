#include "args.h"
#include "parse-data.h"

#include <stdlib.h>

int main(int argc, char** argv) {
    struct option_values options = parse_args(argc, argv);

    struct rate_library* rates = rate_library_create(options);

    struct network* network = network_create(options);

    rate_library_destroy(rates);
    network_destroy(network);
    fclose(options.rate_library_file);
    fclose(options.network_file);
    
    return 0;
}
