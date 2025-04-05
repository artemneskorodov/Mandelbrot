/*============================================================================*/
/* This defines a number of points in one packed vector                       */
/* Supported constants are 1, 4 and 8. They use floats, XMM's and YMM's       */
/* respectively                                                               */

#define RENDER_VECTOR_4
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
#if defined(__x86_64__) or defined(_M_X64)
    #include <xmmintrin.h>
    #include <immintrin.h>
#elif defined(__ARM_NEON__)
    #include <arm_neon.h>
#else
    #error This program obly supports x86-64 or arm with arm neon
#endif
/*----------------------------------------------------------------------------*/
#include "colors.h"
#include "mandelbrot.h"
#include "parse_flags.h"
/*============================================================================*/

static const char          *WindowTitle         = "Mandelbrot";
static const unsigned int   WindowWidth         = 800;
static const unsigned int   WindowHeight        = 600;
static const float          WindowWidthFloat    = (float)WindowWidth;
static const float          WindowHeightFloat   = (float)WindowHeight;
static const unsigned int   MaxIters            = 256;
static const float          PointOutRadiusSq    = 100.f;
static const float          DeltaTime           = 100.f;
static const float          ScaleMult           = 1.2f;
static const size_t         FPSMeasurementIters = 1;
static const size_t         ProgressBarLength   = 70;

/*============================================================================*/

static err_state_t      prog_ctor              (ctx_t          *ctx,
                                                int             argc,
                                                const char     *argv[]);

static err_state_t      prog_dtor              (ctx_t          *ctx);

static inline void      update_window          (ctx_t          *ctx);

static err_state_t      update_position        (ctx_t          *ctx);

static err_state_t      run_normal_mode        (ctx_t          *ctx);

static err_state_t      run_testing_mode       (ctx_t          *ctx,
                                                size_t          render_iters,
                                                unsigned long  *ticks);

static err_state_t      get_render_time        (ctx_t          *ctx);

static sf::Uint32       get_point_color        (unsigned int    iters);

static err_state_t      render_mandelbrot      (ctx_t          *ctx,
                                                size_t          render_iters);

static err_state_t      set_colors             (ctx_t          *ctx);

double                  get_duration           (timespec       *start,
                                                timespec       *end);

static inline void      get_time_real          (timespec       *time);

static void             write_progress_bar     ();

static void             update_progress_bar    (size_t          max,
                                                size_t          current);

inline unsigned long    get_ticks              (void);

/*============================================================================*/
/* This macro is used in main to avoid writing of checking errors and calling */
/* context destructor a lot                                                   */

#define _EXIT_IF_ERROR(...) {                                                  \
    err_state_t _error_code = (__VA_ARGS__);                                   \
    if(_error_code != STATE_SUCCESS) {                                         \
        prog_dtor(&ctx);                                                       \
        return (int)_error_code;                                               \
    }                                                                          \
}                                                                              \

/*----------------------------------------------------------------------------*/
/* Main of program, creates context and read user input. Function runs modes  */
/* which are defined by flags (testing mode which gets render time and normal */
/* mode which draws Mandelbrot set on screen                                  */

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

/*----------------------------------------------------------------------------*/
#undef _EXIT_IF_ERROR
/*============================================================================*/
/* Function runs rendering with different number of iterations which is       */
/* defined by user input (look parse flags) and counts render time with least */
/* squares method. It also writes points (number of iterations and time) to   */
/* file which is also pushed from user with flag '--output'                   */

err_state_t get_render_time(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Avarage values to use least squares method                             */
    unsigned long lu_avg_x  = 0;
    unsigned long lu_avg_y  = 0;
    unsigned long lu_avg_xx = 0;
    unsigned long lu_avg_xy = 0;
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
    /* Initializing progress bar and total number of iterations.              */
    write_progress_bar();
    size_t max = (2 * iters + (ctx->steps_num - 1) * ctx->step_value) *
                 ctx->steps_num / 2;
    size_t current = 0;
    /*------------------------------------------------------------------------*/
    /* Running rendering window with each point rendering iters time for each */
    /* value with requested step                                              */
    for(size_t step = 0; step < ctx->steps_num; step++) {
        /*--------------------------------------------------------------------*/
        /* Variable which will containg rendering time                        */
        unsigned long ticks = 0;
        /*--------------------------------------------------------------------*/
        /* Running rendering                                                  */
        _RETURN_IF_ERROR(run_testing_mode(ctx, iters, &ticks));
        /*--------------------------------------------------------------------*/
        /* Adding value of current point x, y, xx, xy and yy to sums.         */
        lu_avg_x  += iters;
        lu_avg_y  += ticks;
        lu_avg_xx += iters * iters;
        lu_avg_xy += iters * ticks;
        /*--------------------------------------------------------------------*/
        /* Writing point to file.                                             */
        fprintf(output, "%lu, %lu\n", iters, ticks);
        /*--------------------------------------------------------------------*/
        /* Updating iterations number.                                        */
        iters += ctx->step_value;
        /*--------------------------------------------------------------------*/
        /* Updating progress bar.                                             */
        current += iters;
        update_progress_bar(max, current);
    }
    /*------------------------------------------------------------------------*/
    /* Writing new line character after progress bar finished its work.       */
    putchar('\n');
    /*------------------------------------------------------------------------*/
    /* Closing output file                                                    */
    fclose(output);
    /*------------------------------------------------------------------------*/
    /* Getting avarage values for least squares method                        */
    double avg_x  = (double)lu_avg_x  / (double)ctx->steps_num;
    double avg_y  = (double)lu_avg_y  / (double)ctx->steps_num;
    double avg_xx = (double)lu_avg_xx / (double)ctx->steps_num;
    double avg_xy = (double)lu_avg_xy / (double)ctx->steps_num;
    /*------------------------------------------------------------------------*/
    /* Determining angle coefficient and its error using least squares method */
    /*  k = (<xy> - <x><y>) / (<xx> - <x><x>)                                 */
    double render_time_val = (avg_xy - avg_x * avg_y) /
                             (avg_xx - avg_x * avg_x);
    /*------------------------------------------------------------------------*/
    /* Outputting value                                                       */
    color_printf(MAGENTA_TEXT, BOLD_TEXT, DEFAULT_BACKGROUND,
                 "Screen render ticks: %.0f\n",
                 render_time_val);
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
/* Function runs render of screen for number of times specified by parameter  */
/* 'render_iters'. It gets start time and end time using clock_gettime() with */
/* CLOCK_PROCESS_CPUTIME_ID parameter which allows to get procces time, so    */
/* the results does not depend on other proccesses runned on machine.         */
/* Funciton writes render time to 'time' in seconds.                          */

err_state_t run_testing_mode(ctx_t         *ctx,
                             size_t         render_iters,
                             unsigned long *ticks) {
    /*------------------------------------------------------------------------*/
    /* Getting start time using clock_gettime() which is only for Linux but   */
    /* it has less error than clock() or time                                 */
    unsigned long start_ticks = get_ticks();
    /*------------------------------------------------------------------------*/
    _RETURN_IF_ERROR(render_mandelbrot(ctx, render_iters));
    /*------------------------------------------------------------------------*/
    /* Getting end time of rendering                                          */
    unsigned long end_ticks = get_ticks();
    /*------------------------------------------------------------------------*/
    /* Writing render time as double in seconds                               */
    *ticks = end_ticks - start_ticks;
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
/* Function runs loop, which renders Mandelbrot set with current screen       */
/* position and scale. It also allows to move, by calling update_position()   */

err_state_t run_normal_mode(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Running screen until it is closed                                      */
    while(true) {
        timespec start;
        get_time_real(&start);
        for(size_t i = 0; i < FPSMeasurementIters; i++) {
            /*----------------------------------------------------------------*/
            /* Checking keys and events, updating view                        */
            _RETURN_IF_ERROR(update_position(ctx));
            /*----------------------------------------------------------------*/
            /* Counting number of iterations for each point to escape from set*/
            _RETURN_IF_ERROR(render_mandelbrot(ctx, 1));
            /*----------------------------------------------------------------*/
            _RETURN_IF_ERROR(set_colors(ctx));
            /*----------------------------------------------------------------*/
            /* Updating window                                                */
            update_window(ctx);
            /*----------------------------------------------------------------*/
        }
        timespec end;
        get_time_real(&end);
        ctx->fps = (double)FPSMeasurementIters / get_duration(&start, &end);
    }
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
#if defined(RENDER_VECTOR_1)
/* This is conditional compilation part for rendering Mandelbrot set without  */
/* any optimization.                                                          */

/*============================================================================*/
/* This function renders Mandelbrot set with current position for render_iters*/
/* times. It does not use packed float numbers so this is the slowest way of  */
/* rendering.                                                                 */

err_state_t render_mandelbrot(ctx_t *ctx, size_t render_iters) {
    /*------------------------------------------------------------------------*/
    /* Creating constant that used to move point to next                      */
    const float dx = ctx->scale / WindowWidthFloat;
    const float dy = ctx->scale / WindowHeightFloat;
    /*------------------------------------------------------------------------*/
    /* Rendering screen render_iters times                                    */
    for(size_t iteration = 0; iteration < render_iters; iteration++) {
        /*--------------------------------------------------------------------*/
        /* Setting y0 to first line                                           */
        float y0 = -0.5f * ctx->scale - ctx->offset_y;
        /*--------------------------------------------------------------------*/
        /* Creating pointer to current pixel                                  */
        sf::Uint32 *image = ctx->image;
        /*--------------------------------------------------------------------*/
        /* Running through lines                                              */
        for(unsigned int yi = 0; yi < WindowHeight; yi++, y0 += dy) {
            /*----------------------------------------------------------------*/
            /* Setting x0 to first pixel in line                              */
            float x0 = -0.5f * ctx->scale - ctx->offset_x;
            /*----------------------------------------------------------------*/
            /* Running through points in line                                 */
            for(unsigned int xi = 0; xi < WindowWidth; xi++, x0 += dx) {
                /*------------------------------------------------------------*/
                /* Variable to store number of iterations to escape           */
                unsigned int iters = 0;
                /*------------------------------------------------------------*/
                /* Current point position                                     */
                float x = x0;
                float y = y0;
                /*------------------------------------------------------------*/
                for(; iters < MaxIters; iters++) {
                    /*--------------------------------------------------------*/
                    /* Creating x^2, y^2 and 2xy                              */
                    float x_squared = x * x;
                    float y_squared = y * y;
                    float two_x_y = x * y * 2.f;
                    /*--------------------------------------------------------*/
                    /* Checking if escaped                                    */
                    if(x_squared + y_squared > PointOutRadiusSq) {
                        break;
                    }
                    /*--------------------------------------------------------*/
                    /* Updating point coordinates                             */
                    /* x_new = x^2 - y^2 + x0                                 */
                    /* y_new = 2xy + y0                                       */
                    x = x_squared - y_squared + x0;
                    y = two_x_y + y0;
                    /*--------------------------------------------------------*/
                }
                /*------------------------------------------------------------*/
                /* Setting pixel to number of iteration to escape which can   */
                /* be converted to color later                                */
                *image = iters;
                image++;
                /*------------------------------------------------------------*/
            }
            /*----------------------------------------------------------------*/
        }
        /*--------------------------------------------------------------------*/
    }
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
#elif defined(RENDER_VECTOR_4) or defined(RENDER_VECTOR_8)
/* This is conditional compilation part with rendering mandelbrot set using   */
/* packed float numbers                                                       */

/*============================================================================*/
    #if defined(RENDER_VECTOR_4) and (defined(__x86_64__) or defined(_M_X64))
    /* Macroses that are used in render function when rendering with 4 packed */
    /* float numbers is turned on. They are provided to avoid copying of code */
    /* for different packing sizes as rendering is used with similar functions*/

    /*========================================================================*/
        typedef __m128 vector_t;
        /*--------------------------------------------------------------------*/
        typedef __m128i vectorI_t;
        /*--------------------------------------------------------------------*/
        typedef __m128i *pixel_ptr_t;
        /*--------------------------------------------------------------------*/
        static vector_t mm_create_initializer(float dx);
        /*--------------------------------------------------------------------*/
        static const unsigned int PackedSize = 4;
        /*--------------------------------------------------------------------*/
        #define _MM_SET1_PS(_value)             _mm_set1_ps       ((_value ))
        /*--------------------------------------------------------------------*/
        #define _MM_SET1_EPI32(_value)          _mm_set1_epi32    ((_value ))
        /*--------------------------------------------------------------------*/
        #define _MM_ADD_PS(_value1, _value2)    _mm_add_ps        ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_MUL_PS(_value1, _value2)    _mm_mul_ps        ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_SUB_PS(_value1, _value2)    _mm_sub_ps        ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_ADD_EPI32(_value1, _value2) _mm_add_epi32     ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_CMPLE_PS(_value1, _value2)  _mm_cmple_ps      ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_CHECK_ZERO(_value)          _mm_testz_si128   ((_value),   \
                                                                   (_value))
        /*--------------------------------------------------------------------*/
        #define _MM_AND_SI(_value1, _value2)    _mm_and_si128     ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_STORE_SI(_dst, _value)      _mm_store_si128   ((_dst   ),  \
                                                                   (_value ))
        /*--------------------------------------------------------------------*/
        vector_t mm_create_initializer(float x0, float dx) {
            return _MM_ADD_PS(_MM_SET1_PS(x0),
                              __mm_set_ps(3.f * dx,
                                          2.f * dx,
                                          1.f * dx,
                                          0.f * dx));
        }
    /*========================================================================*/
    #elif defined(RENDER_VECTOR_8) and (defined(__x86_64__) or defined(_M_X64))
    /* Macroses that are used in render function when rendering with 8 packed */
    /* float numbers is turned on. They are provided to avoid copying of code */
    /* for different packing sizes as rendering is used with similar functions*/

    /*========================================================================*/
        typedef __m256 vector_t;
        /*--------------------------------------------------------------------*/
        typedef __m256i vectorI_t;
        /*--------------------------------------------------------------------*/
        typedef __m256i *pixel_ptr_t;
        /*--------------------------------------------------------------------*/
        static const unsigned int PackedSize = 8;
        /*--------------------------------------------------------------------*/
        static vector_t mm_create_initializer(float dx);
        /*--------------------------------------------------------------------*/
        #define _MM_SET1_PS(_value)             _mm256_set1_ps    ((_value ))
        /*--------------------------------------------------------------------*/
        #define _MM_SET1_EPI32(_value)          _mm256_set1_epi32 ((_value ))
        /*--------------------------------------------------------------------*/
        #define _MM_ADD_PS(_value1, _value2)    _mm256_add_ps     ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_MUL_PS(_value1, _value2)    _mm256_mul_ps     ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_SUB_PS(_value1, _value2)    _mm256_sub_ps     ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_ADD_EPI32(_value1, _value2) _mm256_add_epi32  ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_CMPLE_PS(_value1, _value2)  _mm256_cmp_ps     ((_value1),  \
                                                                   (_value2),  \
                                                                    18)
        /*--------------------------------------------------------------------*/
        #define _MM_CHECK_ZERO(_value)          _mm256_testz_si256((_value),   \
                                                                   (_value))
        /*--------------------------------------------------------------------*/
        #define _MM_AND_SI(_value1, _value2)    _mm256_and_si256  ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_STORE_SI(_dst, _value)      _mm256_store_si256((_dst   ),  \
                                                                   (_value ))
        /*--------------------------------------------------------------------*/
        vector_t mm_create_initializer(float x0, float dx) {
            return _MM_ADD_PS(_MM_SET1_PS   (x0),
                              __mm256_set_ps(7.f * dx,
                                             6.f * dx,
                                             5.f * dx,
                                             4.f * dx,
                                             3.f * dx,
                                             2.f * dx,
                                             1.f * dx,
                                             0.f * dx));
        }
    /*========================================================================*/
    #elif defined(RENDER_VECTOR_4) and defined(__ARM_NEON__)

    /*========================================================================*/
        typedef float32x4_t vector_t;
        /*--------------------------------------------------------------------*/
        typedef uint32x4_t vectorI_t;
        /*--------------------------------------------------------------------*/
        typedef uint32_t *pixel_ptr_t;
        /*--------------------------------------------------------------------*/
        static const unsigned int PackedSize = 4;
        /*--------------------------------------------------------------------*/
        static vector_t mm_create_x0_init(float x0, float dx);
        /*--------------------------------------------------------------------*/
        #define _MM_SET1_PS(_value)             vdupq_n_f32       ((_value ))
        /*--------------------------------------------------------------------*/
        #define _MM_SET1_EPI32(_value)          vdupq_n_u32       ((_value ))
        /*--------------------------------------------------------------------*/
        #define _MM_ADD_PS(_value1, _value2)    vaddq_f32         ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_MUL_PS(_value1, _value2)    vmulq_f32         ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_SUB_PS(_value1, _value2)    vsubq_f32         ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_ADD_EPI32(_value1, _value2) vaddq_u32         ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_CMPLE_PS(_value1, _value2)  vcleq_f32         ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_CHECK_ZERO(_cmp)          (!vgetq_lane_u64((_cmp), 0) &&   \
                                               !vgetq_lane_u64((_cmp), 1))
        /*--------------------------------------------------------------------*/
        #define _MM_AND_SI(_value1, _value2)    vandq_u32         ((_value1),  \
                                                                   (_value2))
        /*--------------------------------------------------------------------*/
        #define _MM_STORE_SI(_dst, _value)      vst1q_u32         ((_dst   ),  \
                                                                   (_value ))
        /*--------------------------------------------------------------------*/
        vector_t mm_create_x0_init(float x0, float dx) {
            float32_t dx_temp_arr[4] = {0.f * dx, 1.f * dx, 2.f * dx, 3.f * dx};
            vector_t dx_temp = vld1q_f32(dx_temp_arr);
            return _MM_ADD_PS(_MM_SET1_PS(x0), dx_temp);
        }
    /*========================================================================*/
    #endif
/*============================================================================*/
/* This function renders Mandelbrot set when rendering with inrinsics is      */
/* turned on. It uses defines that depend on number of floats in one packed   */
/* vector.                                                                    */

err_state_t render_mandelbrot(ctx_t *ctx, size_t render_iters) {
    /*------------------------------------------------------------------------*/
    /* dx_float variable is used in different places, so saving it once       */
    float     dx_float = ctx->scale / WindowWidthFloat;
    /*------------------------------------------------------------------------*/
    /* Delat for point for each loop iteration                                */
    vector_t  dy       = _MM_SET1_PS(ctx->scale / WindowHeightFloat);
    vector_t  dx       = _MM_SET1_PS((float)PackedSize * dx_float);
    /*------------------------------------------------------------------------*/
    /* Radius to compare with to check that point is escaped                  */
    vector_t  escape_r = _MM_SET1_PS(PointOutRadiusSq);
    /*------------------------------------------------------------------------*/
    /* Mask to use result of _mm_cmple_ps() as addition to iters variable     */
    vectorI_t mask     = _MM_SET1_EPI32(1);
    /*------------------------------------------------------------------------*/
    /* Creating initializer constant for y0 (with PackedSize points)          */
    vector_t  y0_init  = _MM_SET1_PS(-0.5f * ctx->scale - ctx->offset_y);
    /*------------------------------------------------------------------------*/
    /* Creating initializer constant for x0 (with PackedSize points)          */
    /* Temp variables use to avoid using of registers                         */
    vector_t  x0_init  = mm_create_x0_init(-0.5f * ctx->scale - ctx->offset_x,
                                           dx_float);
    /*------------------------------------------------------------------------*/
    /* Rendering the screen render_iters times                                */
    for(size_t iteration = 0; iteration < render_iters; iteration++) {
        /*--------------------------------------------------------------------*/
        /* Pointer to current position in image                               */
        uint32_t *image = ctx->image;
        /*--------------------------------------------------------------------*/
        /* Current y0 vector                                                  */
        vector_t y0 = y0_init;
        /*--------------------------------------------------------------------*/
        /* Running through lines                                              */
        for(unsigned int yi = 0; yi < WindowHeight; yi++) {
            /*----------------------------------------------------------------*/
            /* Current x0 vector                                              */
            vector_t x0 = x0_init;
            /*----------------------------------------------------------------*/
            /* Running through all points on line                             */
            for(unsigned int xi = 0; xi < WindowWidth; xi += PackedSize) {
                /*------------------------------------------------------------*/
                /* This is vector with number of iterations to escape         */
                vectorI_t iters = {0};
                /*------------------------------------------------------------*/
                /* Current points vector                                      */
                vector_t x = x0;
                vector_t y = y0;
                /*------------------------------------------------------------*/
                for(unsigned int n = 0; n < MaxIters; n++) {
                    /*--------------------------------------------------------*/
                    /* Creating x^2, y^2 and 2 * x * y in the start when we   */
                    /* are sure that x and y are in registers                 */
                    vector_t x_square = _MM_MUL_PS(x, x);
                    vector_t y_square = _MM_MUL_PS(y, y);
                    vector_t two_xy = _MM_MUL_PS(x, y);
                    two_xy = _MM_ADD_PS(two_xy, two_xy);
                    /*--------------------------------------------------------*/
                    /* Getting compare result for points and escape radius    */
                    vector_t r_square = _MM_ADD_PS(x_square, y_square);
                    vectorI_t cmp_result = (vectorI_t)_MM_CMPLE_PS(r_square,
                                                                   escape_r);
                    /*--------------------------------------------------------*/
                    /* Aplying bitmask to cmp_result that sets each non zero  */
                    /* element to 1 in low bit                                */
                    cmp_result = _MM_AND_SI(cmp_result, mask);
                    /*--------------------------------------------------------*/
                    /* Checking if all points escaped                         */
                    /*--------------------------------------------------------*/
                    if(_MM_CHECK_ZERO(cmp_result)) {
                        break;
                    }
                    /*--------------------------------------------------------*/
                    /* Applying bitmask that sets not escaped points elements */
                    /* to ones and adding them to iters vector                */
                    iters = _MM_ADD_EPI32(iters, cmp_result);
                    /*--------------------------------------------------------*/
                    /* Next point coordinates                                 */
                    /* x_new = x^2 - y^2 + x0                                 */
                    /* y_new = 2xy + y0                                       */
                    x = _MM_SUB_PS(x_square, y_square);
                    x = _MM_ADD_PS(x, x0);
                    y = _MM_ADD_PS(two_xy, y0);
                    /*--------------------------------------------------------*/
                }
                /*------------------------------------------------------------*/
                /* Setting screen element to escape iters that can be         */
                /* converted to color later and moving screen pointer         */
                _MM_STORE_SI((pixel_ptr_t)image, iters);
                image += PackedSize;
                /*------------------------------------------------------------*/
                /* Moving x0 to next set of PackedSize points                 */
                x0 = _MM_ADD_PS(x0, dx);
                /*------------------------------------------------------------*/
            }
            /*----------------------------------------------------------------*/
            /* Moving y0 to next line                                         */
            y0 = _MM_ADD_PS(y0, dy);
            /*----------------------------------------------------------------*/
        }
        /*--------------------------------------------------------------------*/
    }
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
    #undef _MM_CREATE_INITIALIZER
    #undef _MM_SET1_PS
    #undef _MM_SET1_EPI32
    #undef _MM_ADD_PS
    #undef _MM_MUL_PS
    #undef _MM_SUB_PS
    #undef _MM_ADD_EPI32
    #undef _MM_CMPLE_PS
    #undef _MM_TEST_SI
    #undef _MM_AND_SI
    #undef _MM_STORE_SI
/*============================================================================*/
#endif

/*============================================================================*/
/* Function converts context image to colors. It is expected that every item  */
/* is set to number of iterations for this specific point to escape from      */
/* Mandelbrot set.                                                            */

err_state_t set_colors(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Converting each pixel iterations number to escape to its color         */
    for(unsigned int elem = 0; elem < WindowHeight * WindowWidth; elem++) {
        ctx->image[elem] = get_point_color(ctx->image[elem]);
    }
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
/* Returns point color depending on number of iterations to leave from        */
/* Mandelbrot set. It is expected that parameter is less than MaxIters        */

sf::Uint32 get_point_color(unsigned int iters) {
    /*------------------------------------------------------------------------*/
    /* Red component of color                                                 */
    sf::Uint32 color_red   = 30 + (150 - 20) * (150 - 20)   *
                                   iters     *  iters       /
                                  (MaxIters  *  MaxIters);
    /*------------------------------------------------------------------------*/
    /* Green component of color                                               */
    sf::Uint32 color_green = 10 + (100 - 10) * (100 - 10) *
                                   iters     *  iters     /
                                  (MaxIters  *  MaxIters);
    /*------------------------------------------------------------------------*/
    /* Blue component of color                                                */
    sf::Uint32 color_blue  = 50 * (iters % 2) + (100 - 50) *
                                                 iters     /
                                                (MaxIters);
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
/* Drawing image on screen.                                                   */

inline void update_window(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Updating texture with points array                                     */
    ctx->texture.update((sf::Uint8 *)ctx->image);
    /*------------------------------------------------------------------------*/
    /* Drawing sprite with this texture which fulls the screen                */
    ctx->window.draw(ctx->box);
    /*------------------------------------------------------------------------*/
    /* Updating FPS value                                                     */
    snprintf(ctx->fps_buffer, FPSBufferSize, "%.3f FPS", ctx->fps);
    ctx->fps_text.setString(ctx->fps_buffer);
    ctx->window.draw(ctx->fps_text);
    /*------------------------------------------------------------------------*/
    /* Updating window                                                        */
    ctx->window.display();
    /*------------------------------------------------------------------------*/
}

/*============================================================================*/
/* Mandelbrot context constructor which gets user input, creates image and    */
/* window to draw on.                                                         */

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
        ctx->window.create             (sf::VideoMode(WindowWidth,
                                                      WindowHeight),
                                        WindowTitle);
        /*--------------------------------------------------------------------*/
        /* Creating texture which will fill the window                        */
        ctx->texture.create            (WindowWidth, WindowHeight);
        /*--------------------------------------------------------------------*/
        /* Setting sprite position and texture                                */
        ctx->box.setPosition           (0, 0);
        ctx->box.setTexture            (ctx->texture);
        ctx->box.setTextureRect        (sf::IntRect(0,
                                                    0,
                                                    WindowWidth,
                                                    WindowHeight));
        /*--------------------------------------------------------------------*/
        /* Setting FPS string info                                            */
        ctx->font.loadFromFile         ("font.ttf");
        ctx->fps_text.setFont          (ctx->font);
        ctx->fps_text.setCharacterSize (30);
        ctx->fps_text.setStyle         (sf::Text::Bold);
        ctx->fps_text.setFillColor     (sf::Color::Green);
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
/* Mandelbrot context destructor. It frees image which is allocated           */

err_state_t prog_dtor(ctx_t *ctx) {
    /*------------------------------------------------------------------------*/
    /* Free of points array which is allocated in all modes                   */
    free(ctx->image);
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
/* Writes log message. This function works for long time as it calls flush()  */

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
/* Writes error message. This function works for long time as it calls flush()*/

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
/* Checking if the window is closed and updating position depending on pressed*/
/* keys.                                                                      */

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
            ctx->scale /= ScaleMult;
        }
        else {
            ctx->scale *= ScaleMult;
        }
    }
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
/* This function writes current real time to 'time' structure. It is only     */
/* expected to use as time difference.                                        */
inline void get_time_real(timespec *time) {
    clock_gettime(CLOCK_REALTIME, time);
}

/*============================================================================*/
/* This function returns difference of times as seconds.                      */

double get_duration(timespec *start, timespec *end) {
    /*------------------------------------------------------------------------*/
    /* Determining render time                                                */
    timespec render_time;
    /*------------------------------------------------------------------------*/
    /* If end nanosecond are less than start nanoseconds we add 10e9 to their */
    /* difference and subtructing 1 from seconds difference                   */
    if (end->tv_nsec < start->tv_nsec)
    {
        render_time.tv_sec  = end->tv_sec - start->tv_sec - 1;
        render_time.tv_nsec = 1000000000 + end->tv_nsec  - start->tv_nsec;
    }
    /*------------------------------------------------------------------------*/
    /* Else we just get difference of seconds and nanoseconds                 */
    else
    {
        render_time.tv_sec  = end->tv_sec  - start->tv_sec;
        render_time.tv_nsec = end->tv_nsec - start->tv_nsec;
    }
    /*------------------------------------------------------------------------*/
    /* Returning duration in seconds                                          */
    return (double)render_time.tv_sec + (double)render_time.tv_nsec / 1e9;
    /*------------------------------------------------------------------------*/
}

/*============================================================================*/
/* Function writes progress bar for first time. It is expected that no other  */
/* text will be written to stdout or stderr untill progress reaches its end.  */
/* Use update_progress_bar(size_t, size_t) to update progress bar.            */

void write_progress_bar() {
    /*------------------------------------------------------------------------*/
    /* Printing progress bar.                                                 */
    putchar('[');
    for(size_t i = 0; i < ProgressBarLength; i++) {
        putchar('.');
    }
    putchar(']');
    /*------------------------------------------------------------------------*/
    /* Printing 0 percentage of progress.                                     */
    printf(" 0 %%");
    /*------------------------------------------------------------------------*/
    fflush(stdout);
}

/*============================================================================*/
/* Function updates progress bar. It is expected that no other function will  */
/* write to stdout ot stderr untill progress bar ended. Use function          */
/* write_progress_bar() first to initialize it. Also don't forget to print    */
/* new line character '\n' after progress bar reached its end.                */

void update_progress_bar(size_t max, size_t current) {
    /*------------------------------------------------------------------------*/
    /* Number of symbols in progress bar which are on.                        */
    size_t pixels_on = current * ProgressBarLength / max;
    /*------------------------------------------------------------------------*/
    /* If reached end setting to maximum.                                     */
    if(pixels_on > ProgressBarLength) {
        pixels_on = ProgressBarLength;
    }
    /*------------------------------------------------------------------------*/
    /* Deleting old progress bar from console.                                */
    putchar('\r');
    /*------------------------------------------------------------------------*/
    /* Writing new progress bar.                                              */
    putchar('[');
    for(size_t i = 0; i < pixels_on; i++) {
        putchar('*');
    }
    for(size_t i = 0; i < ProgressBarLength - pixels_on; i++) {
        putchar('.');
    }
    putchar(']');
    /*------------------------------------------------------------------------*/
    /* Writing percentage value of progress.                                  */
    size_t percentage = (size_t)(100.f * (float)current / (float)max) % 100;
    printf(" %lu %%", percentage);
    /*------------------------------------------------------------------------*/
    fflush(stdout);
}

/*============================================================================*/

inline unsigned long get_ticks(void) {
    #if defined(__ARM_NEON__)
        return __builtin_readcyclecounter();
    #elif defined(__x86_64__) or defined(_M_X64)
        return __rdtsc();
    #endif
}

/*============================================================================*/
