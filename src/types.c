#include "types.h"

#include "numeffect.h"

#include <H5Include.h>
#include <linux/limits.h>
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
    real_t sum_x = 0.0f;
    for (int i = 0; i < network->info->number_species; i++) {
        real_t x = network->fptr->y[i] * network->fptr->aa[i];
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

    real_t f_plus_total = 0.0f;
    real_t f_minus_total = 0.0f;
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
    free(options.root_dir);
    // free(options.config_file);
    // free(options.simulation_file);
    return EXIT_SUCCESS;
}

bool toml_bool(toml_result_t toml, char* request) {
    toml_datum_t data = toml_seek(toml.toptab, request);
    if (data.type != TOML_BOOLEAN) {
        printf("==apollo== error: %s is invalid\n", request);
        exit(1);
    }
    return data.u.boolean;
}

int toml_int(toml_result_t toml, char* request) {
    toml_datum_t data = toml_seek(toml.toptab, request);
    if (data.type != TOML_INT64) {
        printf("==apollo== error: %s is invalid\n", request);
        exit(1);
    }
    return data.u.int64;
}

real_t toml_real_t(toml_result_t toml, char* request) {
    toml_datum_t data = toml_seek(toml.toptab, request);
    if (data.type != TOML_FP64) {
        printf("==apollo== error: %s is invalid\n", request);
        exit(1);
    }
    return data.u.fp64;
}

char* toml_string(toml_result_t toml, char* request) {
    toml_datum_t data = toml_seek(toml.toptab, request);
    if (data.type != TOML_STRING) {
        printf("==apollo== error: %s is invalid\n", request);
        exit(1);
    }
    char* str = malloc((strlen(data.u.s) + 1) * sizeof(char));
    strcpy(str, data.u.s);
    return str;
}

// I should keep track of a failstate here somewhere for clean exit on fail...
// NOTE: should the parsing of config_toml be in a separate file?
struct simulation_properties
simulation_properties_create(toml_result_t simulation_toml,
                             toml_result_t config_toml,
                             struct option_values* opts) {
    struct simulation_properties sim_prop;
    sim_prop.hydro_out_file = NULL;
    sim_prop.hydro_temp_effect = NULL;
    sim_prop.hydro_density_effect = NULL;
    sim_prop.print_kernel_time = false;

    // CONFIG INFO
    opts->rocm_accel = toml_bool(config_toml, "base.enablerocm");

    // OUTPUT
    sim_prop.output = toml_bool(simulation_toml, "simulation.output.output");

    char* outputdir;
    if (sim_prop.output) {
        outputdir = toml_string(simulation_toml, "simulation.output.outputdir");
    }

    // TIME
    sim_prop.t_end = toml_real_t(simulation_toml, "simulation.time.endtime");
    sim_prop.output_tres = toml_int(simulation_toml, "simulation.time.tres");
    sim_prop.print_kernel_time =
        toml_bool(simulation_toml, "simulation.time.printkerneltime");
    sim_prop.dt_init = toml_real_t(simulation_toml, "simulation.time.initdt");

    // RESOLUTION
    sim_prop.resolution[0] =
        toml_int(simulation_toml, "simulation.resolution.x");
    sim_prop.resolution[1] =
        toml_int(simulation_toml, "simulation.resolution.y");
    sim_prop.resolution[2] =
        toml_int(simulation_toml, "simulation.resolution.z");

    // HYDRO
    sim_prop.hydro = toml_bool(simulation_toml, "simulation.hydro.use");

    if (sim_prop.hydro && sim_prop.output) {
        char* tmp = toml_string(simulation_toml, "simulation.hydro.outputfile");
        char tmpp[256];
        snprintf(tmpp, 256, "%s/%s/%s", opts->root_dir, outputdir, tmp);
        sim_prop.hydro_out_file = fopen(tmpp, "wa");
        if (sim_prop.hydro_out_file == NULL) {
            printf("==apollo== error: could not open file: %s\n", tmpp);
            goto exit_fail;
        }
        free(tmp);
    }

    bool usehydroeffect = toml_bool(simulation_toml, "simulation.hydro.effect");

    // If we are using a hydro effect:
    // NOTE: I'm not including validation, as I need to make that a func...
    if (usehydroeffect) {
        char* tmp =
            toml_string(simulation_toml, "simulation.hydroeffect.temp.effect");
        if (strcmp(tmp, "random") == 0) {
            sim_prop.hydro_temp_effect = effect_rand;
        } else if (strcmp(tmp, "radial") == 0) {
            sim_prop.hydro_temp_effect = effect_radial;
        } else if (strcmp(tmp, "gradient") == 0) {
            sim_prop.hydro_temp_effect = effect_gradient;
        }
        tmp = toml_string(simulation_toml,
                          "simulation.hydroeffect.density.effect");
        if (strcmp(tmp, "random") == 0) {
            sim_prop.hydro_density_effect = effect_rand;
        } else if (strcmp(tmp, "radial") == 0) {
            sim_prop.hydro_density_effect = effect_radial;
        } else if (strcmp(tmp, "gradient") == 0) {
            sim_prop.hydro_density_effect = effect_gradient;
        }
        free(tmp);

        sim_prop.hydro_temp_base =
            toml_real_t(simulation_toml, "simulation.hydroeffect.temp.base");
        sim_prop.hydro_density_base =
            toml_real_t(simulation_toml, "simulation.hydroeffect.density.base");
    }

    // THERMO
    sim_prop.thermo = toml_bool(simulation_toml, "simulation.thermo.use");
    if (sim_prop.thermo) {
        char* tmp =
            toml_string(simulation_toml, "simulation.thermo.networkfile");
        char file[PATH_MAX];
        snprintf(file, PATH_MAX, "%s/%s", opts->root_dir, tmp);
        sim_prop.network_file = fopen(file, "r");
        if (sim_prop.network_file == NULL) {
            printf("==apollo== error: cannot open thermo network \"%s\"\n",
                   tmp);
            exit(1);
        }
        tmp = toml_string(simulation_toml, "simulation.thermo.ratefile");
        snprintf(file, PATH_MAX, "%s/%s", opts->root_dir, tmp);
        sim_prop.rate_library_file = fopen(file, "r");
        if (sim_prop.rate_library_file == NULL) {
            printf("==apollo== error: cannot open thermo ratelib \"%s\"\n",
                   tmp);
            exit(1);
        }
        free(tmp);
    }

    // NEUTRINO
    sim_prop.neutrino = toml_bool(simulation_toml, "simulation.neutrino.use");
    if (sim_prop.neutrino) {
        char* tmp =
            toml_string(simulation_toml, "simulation.neutrino.opacityfile");
        hid_t fapl;
        if ((fapl = H5Pcreate(H5P_FILE_ACCESS)) == H5I_INVALID_HID) {
            printf("==apollo== error: cannot open neutrino file \"%s\"\n", tmp);
            exit(123);
        }
        char file[PATH_MAX];
        snprintf(file, PATH_MAX, "%s/%s", opts->root_dir, tmp);
        if ((sim_prop.neutrino_file = H5Fopen(file, H5F_ACC_RDONLY, fapl)) ==
            H5I_INVALID_HID) {
            printf("==apollo== error: cannot open neutrino file \"%s\"\n", tmp);
            exit(123);
        }
        free(tmp);
    }

    if (simulation_properties_validate(&sim_prop) == EXIT_FAILURE) {
        goto exit_fail;
    }

    return sim_prop;
exit_fail:
    toml_free(config_toml);
    toml_free(simulation_toml);
    exit(1);
}

int simulation_properties_validate(struct simulation_properties* sim_prop) {
    if (sim_prop->hydro == false && sim_prop->neutrino == false &&
        sim_prop->thermo == false) {
        printf(
            "==apollo== No simulation kernels are being run, is this debug?\n");
        printf("==apollo== (no output will be produced)\n");
        sim_prop->output = false;
    }

    return EXIT_SUCCESS;
}
