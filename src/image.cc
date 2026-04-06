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

/*int image_write(image_ctx_t *ctx, FILE *fd)
{
    size_t count   = ctx->width * ctx->height * RGB_SIZE;
    fprintf(fd, "P6\n%zu %zu\n255\n", ctx->width, ctx->height);
    size_t written = fwrite(ctx->canvas, 1, count, fd);
    return written < count ? EPIPE : 0;
}*/

void image_free(image_ctx_t *ctx)
{
    free(ctx->canvas);
    ctx->canvas = NULL;
}
