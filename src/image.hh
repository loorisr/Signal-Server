#ifndef _IMAGE_HH_
#define _IMAGE_HH_

#include <stdint.h>
#include <stdio.h>

#define RGB_SIZE 3

typedef struct {
    size_t   width;
    size_t   height;
    uint8_t *canvas;
    uint8_t *next_pixel;
} image_ctx_t;

int  image_init(image_ctx_t *ctx, size_t width, size_t height);
int  image_add_pixel(image_ctx_t *ctx, uint8_t r, uint8_t g, uint8_t b);
void image_free(image_ctx_t *ctx);

#endif
