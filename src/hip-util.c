#include "hip-util.h"
#include "kernel/kernel.h"
#include <hip/hip_runtime.h>

struct hipDeviceProp_t* get_hip_device() {
    int num_devices = 0;
    if (hipGetDeviceCount(&num_devices)) {
        printf("==apollo== ERROR: getting number of devices\n");
        return NULL;
    }

    printf("==apollo== devices on system: %i\n", num_devices);

    int best_device_id = 0;
    float best_device_tflops = 0;
    for (int i = 0; i < num_devices; i++) {
        struct hipDeviceProp_t temp_device;
        if (hipGetDevicePropertiesR0600(&temp_device, i)) {
            printf("==apollo== ERROR: cannot get device properties\n");
            exit(1);
        }
        // clockRate * 1000 since HIP reports it in kHz
        float tflops = ((long)temp_device.clockRate * 1000) *
                       temp_device.multiProcessorCount *
                       temp_device.maxBlocksPerMultiProcessor;
        tflops *= 64;    // Estimated number of operations per cycle
        tflops *= 1e-12; // Convert to tera
        printf("==apollo== device: %i (%f TFlops/s)\n", i, tflops);
        if (tflops > best_device_tflops) {
            best_device_tflops = tflops;
            best_device_id = i;
        }
    }

    struct hipDeviceProp_t* device =
        (struct hipDeviceProp_t*)malloc(sizeof(struct hipDeviceProp_t));

    if (hipGetDevicePropertiesR0600(device, best_device_id)) {
        printf("==apollo== ERROR: cannot get device properties\n");
        return NULL;
    }

    printf("==apollo==\n");
    printf("==apollo== using device: %s\n", device->name);
    printf("==apollo== arch: %s\n", device->gcnArchName);
    printf("==apollo== dedicated gpu: %s\n", device->integrated ? "no" : "yes");
    printf("==apollo== clock: %i\n", device->clockRate);
    printf("==apollo== compute units: %i\n",
           device->multiProcessorCount * device->maxBlocksPerMultiProcessor);
    printf("==apollo== estimated TFlops/s: %f\n", best_device_tflops);

    hipSetDevice(best_device_id);

    return device;
}

int benchmark_device(struct hipDeviceProp_t* device) {
    int error;

    printf("==apollo==\n");
    printf("==apollo== beginning benchmark (%s)\n", device->name);

    float* h_A = malloc(3 * sizeof(float));
    h_A[0] = 2.0;
    h_A[1] = 4.0;
    h_A[2] = 6.0;
    float* h_B = malloc(3 * sizeof(float));
    h_B[0] = 0.0;
    h_B[1] = 1.0;
    h_B[2] = 2.0;
    float* h_C = malloc(3 * sizeof(float));

    void* d_A;
    void* d_B;
    void* d_C;
    hipMalloc(&d_A, sizeof(float) * 3);
    hipMalloc(&d_B, sizeof(float) * 3);
    if ((error = hipMalloc(&d_C, sizeof(float) * 3)) != hipSuccess) {
        printf(
            "==apollo== ERROR: encountered error allocating device mem: %i\n",
            error);
        return EXIT_FAILURE;
    }

    hipMemcpy(d_A, h_A, 3 * sizeof(float), hipMemcpyHostToDevice);
    hipMemcpy(d_C, h_C, 3 * sizeof(float), hipMemcpyHostToDevice);
    if ((error = hipMemcpy(d_B, h_B, 3 * sizeof(float),
                           hipMemcpyHostToDevice)) != hipSuccess) {
        printf("==apollo== ERROR: encountered error hip memcpy: %i\n", error);
        return EXIT_FAILURE;
    }

    int sharedmem_allocation = 0 * sizeof(float);

    struct dim3 blockdim = {3, 1, 1};
    struct dim3 griddim = {1, 1, 1};
    if ((error = hipConfigureCall(griddim, blockdim, sharedmem_allocation,
                                  hipStreamDefault)) != hipSuccess) {
        printf("==apollo== ERROR: encountered error configuring kernel: %i\n",
               error);
        return EXIT_FAILURE;
    }

    void** args = malloc(sizeof(void*) * 3);
    args[0] = &d_A;
    args[1] = &d_B;
    args[2] = &d_C;

    if ((error = hipLaunchKernel(vector_mult_kernel, griddim, blockdim, args,
                                 sharedmem_allocation, hipStreamDefault)) !=
        hipSuccess) {
        printf("==apollo== ERROR: encountered error launching kernel: %i\n",
               error);
        return 0;
    }

    hipMemcpy(h_C, d_C, 3 * sizeof(float), hipMemcpyDeviceToHost);

    printf("==apollo== Output: {  ");
    for (int i = 0; i < 3; i++) {
        printf("%f  ", h_C[i]);
    }
    printf("}\n");

    return EXIT_SUCCESS;
}

int devbuf_create(void** devptr, void* hostptr, int size) {
    hipMalloc(devptr, size);
    hipMemcpy(*devptr, hostptr, size, hipMemcpyHostToDevice);
    return EXIT_SUCCESS;
}

int devbuf_read(void* hostptr, void** devptr, int size) {
    hipMemcpy(hostptr, *devptr, size, hipMemcpyDeviceToHost);
    return EXIT_SUCCESS;
}
