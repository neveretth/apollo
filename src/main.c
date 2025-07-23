#include "args.h"
#include "driver.h"
#include "toml/tomlc17.h"

#include <stdlib.h>

int main(int argc, char** argv) {
    // May exit here if certain conditions are not met...
    struct option_values options = parse_args(argc, argv);

    toml_result_t config_toml = toml_parse_file_ex(options.config_file);
    toml_result_t simulation_toml = toml_parse_file_ex(options.simulation_file);

    if (!config_toml.ok) {
        printf("==apollo== error: config file \"%s\" is invalid.\n",
               options.config_file);
        printf("[TOML ERROR]: %s\n", config_toml.errmsg);
        goto exit_fail;
    }
    if (!simulation_toml.ok) {
        printf("==apollo== error: config file \"%s\" is invalid.\n",
               options.simulation_file);
        printf("[TOML ERROR]: %s\n", simulation_toml.errmsg);
        goto exit_fail;
    }

    struct simulation_properties sim_prop =
        simulation_properties_create(simulation_toml, config_toml, &options);

#ifdef __MP_ROCM
    if (options.rocm_accel) {
        printf("==apollo== USING ROCM");
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
        goto exit_fail;
    }

    // Need to destory simulation_properties.
    // But it's not a big deal since it's only made once and has like 10 bytes
    // to it's name.

exit:
    options_clean(options);
    toml_free(config_toml);
    toml_free(simulation_toml);
    return EXIT_SUCCESS;

exit_fail:
    options_clean(options);
    toml_free(config_toml);
    toml_free(simulation_toml);
    return EXIT_FAILURE;
}
