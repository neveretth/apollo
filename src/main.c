#include "args.h"

int main(int argc, char** argv) {
    struct option_values options = parse_args(argc, argv);
    
    return 0;
}
