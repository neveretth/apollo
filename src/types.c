#include "types.h"

#include "numeffect.h"

#include <H5Include.h>
#include <stdlib.h>
#include <string.h>

int freenptr(void** ptr, int len) {
    if (ptr == NULL) {
        return EXIT_SUCCESS;
    }
    for (int i = 0; i < len; i++) {
        free(ptr[i]);
    }
    return EXIT_SUCCESS;
}

int print_abundances(const struct tnn* network) {
    printf("\n\nFINAL ABUNDANCES:\n");
    printf("\nIndex  Isotope   Abundance Y   Mass Frac X");
    float sum_x = 0.0f;
    for (int i = 0; i < network->info->number_species; i++) {
        float x = network->fptr->y[i] * network->fptr->aa[i];
        sum_x += x;
        printf("\n %4d     %4s    %8.4e    %8.4e", i,
               network->info->iso_label[i], network->fptr->y[i], x);
    }
    printf("\n\nsum X = %6.4f\n", sum_x);
    return EXIT_SUCCESS;
}

int print_results(const struct rate_library* rates, const struct tnn* network,
                  const struct problem_parameters* params) {
    printf("\n\nFINAL F+ VALUES:\n");
    for (int i = 0; i < params->f_plus_total; i++) {
        printf("\nF+[%d] = %7.4e  Increases Y[%s] through %s  MapIndex=%d  "
               "FplusFac=%3.1f",
               i, params->f_plus[i],
               network->info->iso_label[params->f_plus_isotope_idx[i]],
               rates->reaction_label[params->f_plus_map[i]],
               params->f_plus_map[i], params->f_plus_factor[i]);
    }
    printf("\n\n\nFINAL F- VALUES:\n");
    for (int i = 0; i < params->f_minus_total; i++) {
        printf("\nF-[%d] = %7.4e  Decreases Y[%s] through %s  MapIndex=%d  "
               "FminusFac=%3.1f",
               i, params->f_minus[i],
               network->info->iso_label[params->f_minus_isotope_idx[i]],
               rates->reaction_label[params->f_minus_map[i]],
               params->f_minus_map[i], params->f_minus_factor[i]);
    }

    printf("\n\n\nF+ and F- MIN AND MAX FOR EACH ISOTOPE:\n");
    for (int i = 0; i < network->info->number_species; i++) {
        printf("\n%3d %5s F+min=%3d F+max=%3d F-min=%3d F-max=%d", i,
               network->info->iso_label[i], params->f_plus_min[i],
               params->f_plus_max[i], params->f_minus_min[i],
               params->f_minus_max[i]);
    }

    printf("\n\n\nSUM OF FLUXES FOR EACH ISOTOPE:\n");

    float f_plus_total = 0.0f;
    float f_minus_total = 0.0f;
    for (int i = 0; i < network->info->number_species; i++) {
        printf("\n%3d %5s  sumF+=%10.4e  sumF-=%10.4e Fnet=%10.4e Y=%10.4e", i,
               network->info->iso_label[i], params->f_plus_sum[i],
               params->f_minus_sum[i],
               params->f_plus_sum[i] - params->f_minus_sum[i],
               network->fptr->y[i]);
        f_plus_total += params->f_plus_sum[i];
        f_minus_total += params->f_minus_sum[i];
    }

    printf("\n\ntotalF+ = %7.4e  totalF- = %7.4e", f_plus_total, f_minus_total);
    return EXIT_SUCCESS;
}

int options_clean(struct option_values options) {
    if (options.thermo_debug) {
        fclose(options.network_file);
        fclose(options.rate_library_file);
    }
    if (options.neutrino_debug) {
        H5Fclose(options.neutrino_file);
    }
    if (options.hydro_debug) {
    }
    return EXIT_SUCCESS;
}

// This is largely garbage... but I need it to know it will _work_ before I
// refine it.
struct simulation_properties
simulation_properties_create(toml_result_t simulation_toml,
                             toml_result_t config_toml) {
    struct simulation_properties sim_prop;
    sim_prop.hydro_out_file = NULL;
    sim_prop.hydro_temp_effect = NULL;
    sim_prop.hydro_density_effect = NULL;

    toml_datum_t data;

    // CONFIG INFO
    // NONE AT THIS TIME.

    // OUTPUT
    data = toml_seek(simulation_toml.toptab, "simulation.output.output");
    if (data.type != TOML_BOOLEAN) {
        printf("==apollo== error: simulation.output.output is invalid\n");
        goto exit_fail;
    }
    sim_prop.output = data.u.boolean;

    // 256 is a reasonable expectation here.
    char outputdir[256];
    if (sim_prop.output) {
        data = toml_seek(simulation_toml.toptab, "simulation.output.outputdir");
        if (data.type != TOML_STRING) {
            printf(
                "==apollo== error: simulation.output.outputdir is invalid\n");
            goto exit_fail;
        }
        strcpy(outputdir, data.u.s);
    }

    data = toml_seek(simulation_toml.toptab, "simulation.output.tres");
    if (data.type != TOML_INT64) {
        printf("==apollo== error: simulation.output.tres is invalid\n");
        goto exit_fail;
    }
    sim_prop.output_tres = data.u.int64;

    // TIME
    data = toml_seek(simulation_toml.toptab, "simulation.time.endtime");
    if (data.type != TOML_FP64) {
        printf("==apollo== error: simulation.time.endtime is invalid\n");
        goto exit_fail;
    }
    sim_prop.t_end = data.u.fp64;

    // RESOLUTION
    data = toml_seek(simulation_toml.toptab, "simulation.resolution.x");
    if (data.type != TOML_INT64) {
        printf("==apollo== error: simulation.output.x is invalid\n");
        goto exit_fail;
    }
    sim_prop.resolution[0] = data.u.int64;

    data = toml_seek(simulation_toml.toptab, "simulation.resolution.y");
    if (data.type != TOML_INT64) {
        printf("==apollo== error: simulation.output.y is invalid\n");
        goto exit_fail;
    }
    sim_prop.resolution[1] = data.u.int64;

    data = toml_seek(simulation_toml.toptab, "simulation.resolution.z");
    if (data.type != TOML_INT64) {
        printf("==apollo== error: simulation.output.z is invalid\n");
        goto exit_fail;
    }
    sim_prop.resolution[2] = data.u.int64;

    // HYDRO
    data = toml_seek(simulation_toml.toptab, "simulation.hydro.use");
    if (data.type != TOML_BOOLEAN) {
        printf("==apollo== error: simulation.hydro.use is invalid\n");
        goto exit_fail;
    }
    sim_prop.hydro = data.u.boolean;

    if (sim_prop.hydro && sim_prop.output) {
        data = toml_seek(simulation_toml.toptab, "simulation.hydro.outputfile");
        if (data.type != TOML_STRING) {
            printf(
                "==apollo== error: simulation.hydro.outputfile is invalid\n");
            goto exit_fail;
        }
        char tmp[256];
        snprintf(tmp, 256, "%s/%s", outputdir, data.u.s);
        sim_prop.hydro_out_file = fopen(tmp, "wa");
        if (sim_prop.hydro_out_file == NULL) {
            printf("==apollo== error: could not open file: %s\n", tmp);
            goto exit_fail;
        }
    }

    data = toml_seek(simulation_toml.toptab, "simulation.hydro.effect");
    if (data.type != TOML_BOOLEAN) {
        printf("==apollo== error: simulation.hydro.effect is invalid\n");
        goto exit_fail;
    }

    // If we are using a hydro effect:
    // NOTE: I'm not including validation, as I need to make that a func...
    if (data.u.boolean) {
        data = toml_seek(simulation_toml.toptab,
                         "simulation.hydroeffect.temp.effect");
        if (strcmp(data.u.s, "random") == 0) {
            sim_prop.hydro_temp_effect = effect_rand;
        } else if (strcmp(data.u.s, "radial") == 0) {
            sim_prop.hydro_temp_effect = effect_radial;
        } else if (strcmp(data.u.s, "gradient") == 0) {
            sim_prop.hydro_temp_effect = effect_gradient;
        }
        data = toml_seek(simulation_toml.toptab,
                         "simulation.hydroeffect.density.effect");
        if (strcmp(data.u.s, "random") == 0) {
            sim_prop.hydro_density_effect = effect_rand;
        } else if (strcmp(data.u.s, "radial") == 0) {
            sim_prop.hydro_density_effect = effect_radial;
        } else if (strcmp(data.u.s, "gradient") == 0) {
            sim_prop.hydro_density_effect = effect_gradient;
        }
    }

    // THERMO
    data = toml_seek(simulation_toml.toptab, "simulation.thermo.use");
    if (data.type != TOML_BOOLEAN) {
        printf("==apollo== error: simulation.thermo.use is invalid\n");
        goto exit_fail;
    }
    sim_prop.thermo = data.u.boolean;

    // NEUTRINO
    data = toml_seek(simulation_toml.toptab, "simulation.neutrino.use");
    if (data.type != TOML_BOOLEAN) {
        printf("==apollo== error: simulation.neutrino.use is invalid\n");
        goto exit_fail;
    }
    sim_prop.neutrino = data.u.boolean;

    return sim_prop;
exit_fail:
    toml_free(config_toml);
    toml_free(simulation_toml);
    exit(1);
}
