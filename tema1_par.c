// Author: APD team, except where source was noted

#define _XOPEN_SOURCE 600

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
unsigned char **sample_grid(ppm_image *image, int step_x, int step_y, unsigned char sigma,
                            unsigned char **grid, int thread_id, int thread_num) {
    int p = image->x / step_x;
    int q = image->y / step_y;
    int start = thread_id * p / thread_num;
    int end = min((thread_id + 1) * p / thread_num, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    grid[p][q] = 0;

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }
    start = thread_id * q / thread_num;
    end = min((thread_id + 1) * q / thread_num, q);
    for (int j = start; j < end; j++) {
        ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }

    return grid;
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(ppm_image *image, unsigned char **grid, ppm_image **contour_map, int step_x, int step_y, int thread_id, int thread_num) {
    int p = image->x / step_x;
    int q = image->y / step_y;
    int start = thread_id * p / thread_num;
    int end = min((thread_id + 1) * p / thread_num, p);
    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

ppm_image *rescale_image(ppm_image *image, ppm_image *new_image, int thread_id, int thread_num) {
    uint8_t sample[3];

    int start = thread_id * RESCALE_X / thread_num;
    int end = min((thread_id + 1) * RESCALE_X / thread_num, RESCALE_X);

    // use bicubic interpolation for scaling
    for (int i = start; i < end; i++) {
        for (int j = 0; j < RESCALE_Y; j++) {
            float u = (float)i / (float)(RESCALE_X - 1);
            float v = (float)j / (float)(RESCALE_Y - 1);
            sample_bicubic(image, u, v, sample);

            new_image->data[i * RESCALE_Y + j].red = sample[0];
            new_image->data[i * RESCALE_Y + j].green = sample[1];
            new_image->data[i * RESCALE_Y + j].blue = sample[2];
        }
    }

    return new_image;
}
//thrad data struct
typedef struct thread_data {
    int step_x;
    int step_y;
    int sigma;
    ppm_image *image;
    ppm_image *rescale_image;
    unsigned char **grid;
    ppm_image **contour_map;
    pthread_barrier_t *barrier;
    int thread_num;
} thread_data;

//thread only thread struct
typedef struct thread {
    int thread_id;
    thread_data *args;
} thread;

int min(int a, int b) {
    if (a < b)
        return a;
    return b;
}

void* thread_function(void* args) {
    thread *t = (thread*)args;
    thread_data *args_data = t->args;

    // 1. Rescale the image
    if (args_data->rescale_image != NULL) {
        //rescale
        args_data->rescale_image = rescale_image(args_data->image, args_data->rescale_image,
                                                 t->thread_id, args_data->thread_num); 
        pthread_barrier_wait(args_data->barrier);
        //grid
        args_data->grid = sample_grid(args_data->rescale_image, args_data->step_x, args_data->step_y,
                                      args_data->sigma, args_data->grid, t->thread_id, args_data->thread_num);
        pthread_barrier_wait(args_data->barrier);
        march(args_data->rescale_image, args_data->grid, args_data->contour_map, args_data->step_x,
              args_data->step_y, t->thread_id, args_data->thread_num);
    }
    else {
        //grid
        args_data->grid = sample_grid(args_data->image, args_data->step_x, args_data->step_y,
                                      args_data->sigma, args_data->grid, t->thread_id, args_data->thread_num);
        pthread_barrier_wait(args_data->barrier);
        march(args_data->image, args_data->grid, args_data->contour_map, args_data->step_x,
              args_data->step_y, t->thread_id, args_data->thread_num);
    }
    
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    ppm_image *image = read_ppm(argv[1]);
    int step_x = STEP;
    int step_y = STEP;
    int q = image->y / step_y;
    int p = image->x / step_x;
    // 0. Initialize contour map
    ppm_image **contour_map = init_contour_map();

    // memory allocation for threads
    int threads_num = atoi(argv[3]);
    pthread_t *threads = (pthread_t*)malloc(threads_num * sizeof(pthread_t));
    if (!threads) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

        thread_data *args = (thread_data*)malloc(sizeof(thread_data));
        if (!args) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
        args->step_x = step_x;
        args->step_y = step_y;
        args->sigma = SIGMA;
        args->image = image;
        args->rescale_image = NULL;
        args->grid = NULL;
        args->contour_map = contour_map;
        args->thread_num = threads_num;
        args->barrier = (pthread_barrier_t*)malloc(sizeof(pthread_barrier_t));
        if (!args->barrier) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
        pthread_barrier_init(args->barrier, NULL, threads_num);
    
    // we only rescale downwards
    if (!(image->x <= RESCALE_X && image->y <= RESCALE_Y)) {
        args->rescale_image = (ppm_image *)malloc(sizeof(ppm_image));
        if (!args->rescale_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
        args->rescale_image->x = RESCALE_X;
        args->rescale_image->y = RESCALE_Y;

        args->rescale_image->data = (ppm_pixel*)malloc(RESCALE_X * RESCALE_Y * sizeof(ppm_pixel));
        if (!args->rescale_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }
    // memory allocation for grid
    args->grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!args->grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        args->grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!args->grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    for (int i = 0; i < threads_num; i++) {
        thread *t = (thread*)malloc(sizeof(thread));
        if (!t) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
        t->thread_id = i;
        t->args = args;
        pthread_create(&threads[i], NULL, thread_function, t);
    }

    // write output
    for (int i = 0; i < threads_num; i++) {
        pthread_join(threads[i], NULL);
    }
    if (args->rescale_image != NULL)
        write_ppm(args->rescale_image, argv[2]);
    else
        write_ppm(args->image, argv[2]);

    free(threads);
    free(args->barrier);
    free(args);

    return 0;
}
