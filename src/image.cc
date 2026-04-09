#include <stdlib.h>
#include <errno.h>
#include "image.hh"

int image_init(image_ctx_t *ctx, size_t width, size_t height)
{
    if (ctx == NULL || width == 0 || height == 0)
        return EINVAL;

    ctx->width  = width;
    ctx->height = height;
    ctx->canvas = (uint8_t *)calloc(width * height * RGB_SIZE, 1);
    if (ctx->canvas == NULL)
        return ENOMEM;
    ctx->next_pixel = ctx->canvas;
    return 0;
}

int image_add_pixel(image_ctx_t *ctx, uint8_t r, uint8_t g, uint8_t b)
{
    ctx->next_pixel[0] = r;
    ctx->next_pixel[1] = g;
    ctx->next_pixel[2] = b;
    ctx->next_pixel += RGB_SIZE;
    return 0;
}

void image_free(image_ctx_t *ctx)
{
    free(ctx->canvas);
    ctx->canvas = NULL;
}
