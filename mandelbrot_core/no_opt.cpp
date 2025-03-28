/*============================================================================*/
#include <stdio.h>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <stdarg.h>
#include <math.h>
/*----------------------------------------------------------------------------*/
#include "colors.h"
#include "mandelbrot.h"
#include "parse_flags.h"
/*============================================================================*/

#define _EXIT_IF_ERROR(...) {                                                  \
    err_state_t _error_code = (__VA_ARGS__);                                   \
    if(_error_code != STATE_SUCCESS) {                                         \
        prog_dtor(&ctx);                                                       \
        return (int)_error_code;                                               \
    }                                                                          \
}                                                                              \

/*============================================================================*/

static const char          *WindowTitle       = "Mandelbrot";
static const unsigned int   WindowWidth       = 800;
static const unsigned int   WindowHeight      = 600;
static const float          WindowWidthFloat  = (float)WindowWidth;
static const float          WindowHeightFloat = (float)WindowHeight;
static const float          DefaultDX         = DefaultScale /
                                                WindowWidthFloat;
static const float          DefaultDY         = DefaultScale /
                                                WindowHeightFloat;
static const unsigned int   MaxIters          = 256;
static const float          PointOutRadiusSq  = 100.f;
static const float          DeltaTime         = 5.0;

/*============================================================================*/

static err_state_t          prog_ctor              (ctx_t          *ctx,
                                                    int             argc,
                                                    const char     *argv[]);

static err_state_t          prog_dtor              (ctx_t          *ctx);

static inline void          update_window          (ctx_t          *ctx);

static err_state_t          update_position        (ctx_t          *ctx);

static err_state_t          run_normal_mode        (ctx_t          *ctx);

static err_state_t          run_testing_mode       (ctx_t          *ctx,
                                                    size_t          point_iters,
                                                    double         *time);

static err_state_t          get_render_time        (ctx_t          *ctx);

static sf::Uint32    get_point_color        (unsigned int    iters);

static unsigned int  get_mandelbrot_iters   (float           x0,
                                                    float           y0);

/*============================================================================*/

int main(int argc, const char *argv[]) {
    /*------------------------------------------------------------------------*/
    /* Initializing context                                                   */
    ctx_t ctx = {};
    _EXIT_IF_ERROR(prog_ctor(&ctx, argc, argv));

    /*------------------------------------------------------------------------*/
    /* Running program depending on mode                                      */
    if(ctx.testing_mode) {
        _EXIT_IF_ERROR(get_render_time(&ctx));
    }
    else {
        _EXIT_IF_ERROR(run_normal_mode(&ctx));
    }
    /*------------------------------------------------------------------------*/
    /* Calling the destructor of context and exit                             */
    prog_dtor(&ctx);
    return EXIT_SUCCESS;
}

/*============================================================================*/

err_state_t get_render_time(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Avarage values to use least squares method                             */
    double avg_x  = 0;
    double avg_y  = 0;
    double avg_xx = 0;
    double avg_xy = 0;
    double avg_yy = 0;
    /*------------------------------------------------------------------------*/
    /* Current number of iterations is stored in variable iters               */
    size_t iters = ctx->iters_point_min;
    /*------------------------------------------------------------------------*/
    /* Opening file to write points                                           */
    FILE *output = fopen(ctx->output_file, "w");
    /*------------------------------------------------------------------------*/
    /* Checking for errors while opening file                                 */
    if(output == NULL) {
        err_msg("Error while opening file '%s'\n",
                ctx->output_file);
        return STATE_FILE_OPENING_ERROR;
    }
    /*------------------------------------------------------------------------*/
    /* Running rendering window with each point rendering iters time for each */
    /* value with requested step                                              */
    for(size_t step = 0; step < ctx->steps_num; step++) {
        /*--------------------------------------------------------------------*/
        /* Variable which will containg rendering time                        */
        double time = 0;
        /*--------------------------------------------------------------------*/
        /* Running rendering                                                  */
        _RETURN_IF_ERROR(run_testing_mode(ctx, iters, &time));
        /*--------------------------------------------------------------------*/
        log_msg("Test %lu completed\n", step + 1);
        /*--------------------------------------------------------------------*/
        /* Adding value of current point x, y, xx, xy and yy to sums          */
        avg_x  += (double)iters;
        avg_y  +=         time;
        avg_xx += (double)iters * (double)iters;
        avg_yy +=         time  *         time;
        avg_xy += (double)iters *         time;
        /*--------------------------------------------------------------------*/
        /* Writing point to file                                              */
        fprintf(output, "%lu, %f\n", iters, time);
        /*--------------------------------------------------------------------*/
        /* Updating iterations number                                         */
        iters += ctx->step_value;
    }
    /*------------------------------------------------------------------------*/
    /* Closing output file                                                    */
    fclose(output);
    /*------------------------------------------------------------------------*/
    /* Getting avarage values for least squares method                        */
    avg_x  /= (double)ctx->steps_num;
    avg_y  /= (double)ctx->steps_num;
    avg_xx /= (double)ctx->steps_num;
    avg_xy /= (double)ctx->steps_num;
    avg_yy /= (double)ctx->steps_num;
    /*------------------------------------------------------------------------*/
    /* Determining angle coefficient and its error using least squares method */
    /*  k = (<xy> - <x><y>) / (<xx> - <x><x>)                                 */
    /* dk = (1 / sqrt(n)) * sqrt((<yy> - <y><y>) / (<xx> - <x><x>) - k * k)   */
    double render_time_val = (avg_xy - avg_x * avg_y) /
                             (avg_xx - avg_x * avg_x);
    double render_time_err = sqrt(((avg_yy - avg_y * avg_y) /
                                   (avg_xx - avg_x * avg_x) -
                                   render_time_val * render_time_val) /
                                  (double)ctx->steps_num);
    /*------------------------------------------------------------------------*/
    /* Outputting value                                                       */
    color_printf(MAGENTA_TEXT, BOLD_TEXT, DEFAULT_BACKGROUND,
                 "╭─────────────────────────────────────────╮\n"
                 "│ Screen render time: %.6f ± %.6f │\n"
                 "╰─────────────────────────────────────────╯\n",
                 render_time_val,
                 render_time_err);
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/

err_state_t run_testing_mode(ctx_t *ctx, size_t point_iters, double *time) {
    /*------------------------------------------------------------------------*/
    /* Getting start time using clock_gettime() which is only for Linux but   */
    /* it has less error than clock() or time                                 */
    timespec start;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    /*------------------------------------------------------------------------*/
    /* Y-coordinate of first point on screen                                  */
    float y0 = -0.5f * DefaultScale - DefaultOffsetY;
    /*------------------------------------------------------------------------*/
    /* Running through points from top to bottom                              */
    for(unsigned int yi = 0; yi < WindowHeight; yi++, y0 += DefaultDY) {
        /*--------------------------------------------------------------------*/
        /* X-coordinate of first point for each line determined by Y0         */
        float x0 = -0.5f * DefaultScale - DefaultOffsetX;
        /*--------------------------------------------------------------------*/
        /* Running through points on line                                     */
        for(unsigned int xi = 0; xi < WindowWidth; xi++, x0 += DefaultDX) {
            /*----------------------------------------------------------------*/
            /* Variable to store itarations number                            */
            unsigned int iters = 0;
            /*----------------------------------------------------------------*/
            /* Getting point itarations number for point_iters times          */
            for(size_t i = 0; i < point_iters; i++) {
                /*------------------------------------------------------------*/
                /* Running function to get coordinates                        */
                iters = get_mandelbrot_iters(x0, y0);
                /*------------------------------------------------------------*/
                /* Writing color to array to make compiler save this function */
                ctx->image[xi + WindowWidth * yi] = get_point_color(iters);
                /*------------------------------------------------------------*/
            }
            /*----------------------------------------------------------------*/
        }
        /*--------------------------------------------------------------------*/
    }
    /*------------------------------------------------------------------------*/
    /* Getting end time of rendering                                          */
    timespec end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    /*------------------------------------------------------------------------*/
    /* Determining render time                                                */
    timespec render_time;
    /*------------------------------------------------------------------------*/
    /* If end nanosecond are less than start nanoseconds we add 10e9 to their */
    /* difference and subtructing 1 from seconds difference                   */
    if (end.tv_nsec < start.tv_nsec)
    {
        render_time.tv_sec  = end.tv_sec - start.tv_sec - 1;
        render_time.tv_nsec = 1000000000 + end.tv_nsec  - start.tv_nsec;
    }
    /*------------------------------------------------------------------------*/
    /* Else we just get difference of seconds and nanoseconds                 */
    else
    {
        render_time.tv_sec  = end.tv_sec  - start.tv_sec;
        render_time.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
    /*------------------------------------------------------------------------*/
    /* Writing answer as double in seconds                                    */
    *time = (double)render_time.tv_sec + (double)render_time.tv_nsec / 1e9;
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/

err_state_t run_normal_mode(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Running screen until it is closed                                      */
    while(true) {
        /*--------------------------------------------------------------------*/
        log_msg("\t-Frame start\n");
        /*--------------------------------------------------------------------*/
        /* Checking keys and window closed event                              */
        _RETURN_IF_ERROR(update_position(ctx));
        /*--------------------------------------------------------------------*/
        /* First point Y-coordinate                                           */
        float y0 = -0.5f * ctx->scale - ctx->offset_y;
        /*--------------------------------------------------------------------*/
        /* Addition to point Y-coordinate for each                            */
        float dy = ctx->scale / WindowHeightFloat;
        /*--------------------------------------------------------------------*/
        /* Running through all points                                         */
        for(unsigned int yi = 0; yi < WindowHeight; yi++, y0 += dy) {
            /*----------------------------------------------------------------*/
            /* X-coordinate of first point in line determined by y0           */
            float x0 = -0.5f * ctx->scale - ctx->offset_x;
            /*----------------------------------------------------------------*/
            /* Addition to point X-coordinate                                 */
            float dx = ctx->scale / WindowWidthFloat;
            /*----------------------------------------------------------------*/
            /* Running through line                                           */
            for(unsigned int xi = 0; xi < WindowWidth; xi++, x0 += dx) {
                /*------------------------------------------------------------*/
                /* Number of iterations to get out of circle with radius      */
                /* sqrt(PointOutRadiusSq)                                     */
                unsigned int iters = get_mandelbrot_iters(x0, y0);
                /*------------------------------------------------------------*/
                /* Printing point color depending on its iterations number    */
                ctx->image[xi + WindowWidth * yi] = get_point_color(iters);
                /*------------------------------------------------------------*/
            }
            /*----------------------------------------------------------------*/
        }
        /*--------------------------------------------------------------------*/
        log_msg("\t-Updated points\n");
        /*--------------------------------------------------------------------*/
        /* Updating window                                                    */
        update_window(ctx);
        /*--------------------------------------------------------------------*/
        log_msg("\t-Frame end\n");
        /*--------------------------------------------------------------------*/
    }
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/

sf::Uint32 get_point_color(unsigned int iters) {
    /*------------------------------------------------------------------------*/
    /* Red component of color                                                 */
    sf::Uint32 color_red   = 20 + (150 - 20) * (150 - 20) * (150 - 20) *
                                   iters     *  iters     *  iters /
                                  (MaxIters  *  MaxIters  *  MaxIters);
    /*------------------------------------------------------------------------*/
    /* Green component of color                                               */
    sf::Uint32 color_green = 10 + (100 - 10) * (100 - 10) *
                                   iters     *  iters /
                                  (MaxIters  *  MaxIters);
    /*------------------------------------------------------------------------*/
    /* Blue component of color                                                */
    sf::Uint32 color_blue  = 50 * (iters % 2) + (100 - 50) * (100 - 50) *
                                                 iters     *  iters /
                                                (MaxIters  *  MaxIters);
    /*------------------------------------------------------------------------*/
    /* Alpha is always 100%                                                   */
    sf::Uint32 alpha = 255;
    /*------------------------------------------------------------------------*/
    /* If point did not reach MaxIters we assume that it left forever         */
    if(iters != MaxIters) {
        return (color_red   << 0 ) +
               (color_green << 8 ) +
               (color_blue  << 16) +
               (alpha       << 24);
    }
    /*------------------------------------------------------------------------*/
    /* Setting color to black if it did not reached the circle                */
    return (sf::Uint32)(255 << 24);
    /*------------------------------------------------------------------------*/
}

/*============================================================================*/

unsigned int get_mandelbrot_iters(float x0, float y0) {
    /*------------------------------------------------------------------------*/
    /* Number of iterations                                                   */
    unsigned int n = 0;
    /*------------------------------------------------------------------------*/
    /* Current position                                                       */
    float x = x0;
    float y = y0;
    /*------------------------------------------------------------------------*/
    /* Running iterations                                                     */
    for(; n < MaxIters; n++) {
        /*--------------------------------------------------------------------*/
        /* Getting new point value                                            */
        float x_new = x * x - y * y + x0;
        y = 2 * x * y + y0;
        x = x_new;
        /*--------------------------------------------------------------------*/
        /* Ending if it reached circle with radius sqrt(PointOutRadiusSq)     */
        if(x * x + y * y > PointOutRadiusSq) {
            break;
        }
        /*--------------------------------------------------------------------*/
    }
    /*------------------------------------------------------------------------*/
    /* Returning number of iterations                                         */
    return n;
}

/*============================================================================*/

inline void update_window(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Updating texture with points array                                     */
    ctx->texture.update((sf::Uint8 *)ctx->image);
    /*------------------------------------------------------------------------*/
    /* Drawing sprite with this texture which fulls the screen                */
    ctx->window.draw(ctx->box);
    /*------------------------------------------------------------------------*/
    /* Updating window                                                        */
    ctx->window.display();
    /*------------------------------------------------------------------------*/
}

/*============================================================================*/

err_state_t prog_ctor(ctx_t *ctx, int argc, const char *argv[]) {
    /*------------------------------------------------------------------------*/
    log_msg("Starting initializing\n");
    /*------------------------------------------------------------------------*/
    /* Parsing command line flags                                             */
    _RETURN_IF_ERROR(parse_flags(ctx, argc, argv));
    /*------------------------------------------------------------------------*/
    log_msg("Successfully parsed flags\n");
    /*------------------------------------------------------------------------*/
    /* Graphics mode                                                          */
    if(!ctx->testing_mode) {
        /*--------------------------------------------------------------------*/
        /* Creating window                                                    */
        ctx->window.create(sf::VideoMode(WindowWidth, WindowHeight),
                           WindowTitle);
        /*--------------------------------------------------------------------*/
        /* Creating texture which will fill the window                        */
        ctx->texture.create(WindowWidth, WindowHeight);
        /*--------------------------------------------------------------------*/
        /* Setting sprite position and texture                                */
        ctx->box.setPosition(0, 0);
        ctx->box.setTexture(ctx->texture);
        ctx->box.setTextureRect(sf::IntRect(0, 0, WindowWidth, WindowHeight));
        /*--------------------------------------------------------------------*/
        log_msg("Successfully created window\n");
        /*--------------------------------------------------------------------*/
    }
    /*------------------------------------------------------------------------*/
    /* Allocating points array which is used both in normal and testing mode  */
    ctx->image = (sf::Uint32 *)calloc(WindowWidth * WindowHeight,
                                      sizeof(sf::Uint32));
    /*------------------------------------------------------------------------*/
    /* Checking for allocation errors                                         */
    if(ctx->image == NULL) {
        err_msg("Error while allocating memory for screen buffer.\n");
        return STATE_MEMORY_ERROR;
    }
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/

err_state_t prog_dtor(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Free of points array which is allocated in all modes                   */
    free(ctx->image);
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/

int log_msg(const char *format, ...) {
    /*------------------------------------------------------------------------*/
    /* Getting args list                                                      */
    va_list args;
    va_start(args, format);
    /*------------------------------------------------------------------------*/
    /* Printing log message                                                   */
    int printed_symbols = color_vprintf(YELLOW_TEXT,
                                        NORMAL_TEXT,
                                        DEFAULT_BACKGROUND,
                                        format,
                                        args);
    /*------------------------------------------------------------------------*/
    /* End of args                                                            */
    va_end(args);
    /*------------------------------------------------------------------------*/
    /* Printing stdout as X11 writes to stderr and its messages get colored   */
    fflush(stdout);
    /*------------------------------------------------------------------------*/
    return printed_symbols;
}

/*============================================================================*/

int err_msg(const char *format, ...) {
    /*------------------------------------------------------------------------*/
    /* Getting args list                                                      */
    va_list args;
    va_start(args, format);
    /*------------------------------------------------------------------------*/
    /* Printing error message                                                 */
    int printed_symbols = color_vprintf(RED_TEXT,
                                        BOLD_TEXT,
                                        DEFAULT_BACKGROUND,
                                        format,
                                        args);
    /*------------------------------------------------------------------------*/
    /* End of args                                                            */
    va_end(args);
    /*------------------------------------------------------------------------*/
    /* Printing stdout as X11 writes to stderr and its messages get colored   */
    fflush(stdout);
    /*------------------------------------------------------------------------*/
    return printed_symbols;
}

/*============================================================================*/

err_state_t update_position(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Checking if window is closed                                           */
    sf::Event event;
    while(ctx->window.pollEvent(event)) {
        if(event.type == sf::Event::Closed) {
            ctx->window.close();
            return STATE_EXIT_SUCCESS;
        }
    }
    /*------------------------------------------------------------------------*/
    /* Getting arrows pressed values                                          */
    /* Getting P and O pressed values (P to increase and O to decrease size)  */
    int  left    = sf::Keyboard::isKeyPressed(sf::Keyboard::Left );
    int  right   = sf::Keyboard::isKeyPressed(sf::Keyboard::Right);
    int  up      = sf::Keyboard::isKeyPressed(sf::Keyboard::Up   );
    int  down    = sf::Keyboard::isKeyPressed(sf::Keyboard::Down );
    bool larger  = sf::Keyboard::isKeyPressed(sf::Keyboard::P    );
    bool smaller = sf::Keyboard::isKeyPressed(sf::Keyboard::O    );
    /*------------------------------------------------------------------------*/
    /* Updating position                                                      */
    ctx->offset_x -= ctx->scale / DeltaTime * (float)(right - left);
    ctx->offset_y -= ctx->scale / DeltaTime * (float)(down  - up  );
    /*------------------------------------------------------------------------*/
    /* Updating scale                                                         */
    if(larger != smaller) {
        if(larger) {
            ctx->scale /= 1.5f;
        }
        else {
            ctx->scale *= 1.5f;
        }
    }
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
