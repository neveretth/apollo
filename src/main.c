#include "args.h"
#include "driver.h"
#include "toml/tomlc17.h"

#include <stdlib.h>

int main(int argc, char** argv) {
    
    // May exit here if certain conditions are not met...
    struct option_values options = parse_args(argc, argv);

    toml_result_t config_toml = toml_parse_file_ex(options.config_file);
    toml_result_t simulation_toml = toml_parse_file_ex(options.simulation_file);

    bool fail = false;

    if (!config_toml.ok) {
        printf("==apollo== error: config file \"%s\" is invalid.\n",
               options.config_file);
        printf("[TOML ERROR]: %s\n", config_toml.errmsg);
        fail = true;
        goto exit;
    }
    if (!simulation_toml.ok) {
        printf("==apollo== error: config file \"%s\" is invalid.\n",
               options.simulation_file);
        printf("[TOML ERROR]: %s\n", simulation_toml.errmsg);
        fail = true;
        goto exit;
    }

    struct simulation_properties sim_prop =
        simulation_properties_create(simulation_toml, config_toml, &options);

#ifdef __MP_ROCM
    if (options.rocm_accel) {
        printf("==apollo== USING ROCM\n");
    }
#else
    if (options.rocm_accel) {
        printf("==apollo== NOT USING ROCM (running non-rocm apollo)\n");
    }
#endif

    printf("==apollo== Simulation config read.\n");
    printf("==apollo== Beginning simulation...\n");

    if (unified_driver(sim_prop, options) == EXIT_FAILURE) {
        printf("==apollo== error: driver failed.\n");
        fail = true;
        goto exit;
    }

exit:
    simulation_properties_clean(sim_prop);
    options_clean(options);
    toml_free(config_toml);
    toml_free(simulation_toml);
    if (fail) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
