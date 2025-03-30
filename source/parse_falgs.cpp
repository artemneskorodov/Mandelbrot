/*============================================================================*/
#include <string.h>
/*----------------------------------------------------------------------------*/
#include "mandelbrot.h"
#include "parse_flags.h"
#include "colors.h"
/*============================================================================*/

err_state_t graphics_handler    (ctx_t         *ctx,
                                 const char    *argv[],
                                 int            argc,
                                 int           *current_flag);

err_state_t test_handler        (ctx_t         *ctx,
                                 const char    *argv[],
                                 int            argc,
                                 int           *current_flag);

err_state_t output_handler      (ctx_t         *ctx,
                                 const char    *argv[],
                                 int            argc,
                                 int           *current_flag);

/*============================================================================*/

struct flag_handler_t {
    err_state_t   (*handler)(ctx_t *, const char *[], int, int *);
    const char     *long_name;
    const char     *short_name;
};

/*============================================================================*/

static const flag_handler_t FlagHandlers[] = {
    {graphics_handler, "--graphics", "-g"},
    {    test_handler,     "--test", "-t"},
    {  output_handler,   "--output", "-o"}};

/*----------------------------------------------------------------------------*/

static const size_t FlagHandlersSize = sizeof(FlagHandlers) /
                                       sizeof(FlagHandlers[0]);

/*============================================================================*/

err_state_t parse_flags(ctx_t *ctx, int argc, const char *argv[]) {
    /*------------------------------------------------------------------------*/
    /* Current parsing flag. This index is moved by other functions           */
    int cur_flag = 1;
    /*------------------------------------------------------------------------*/
    /* Running graphics as default if no flags passed                         */
    if(argc == 1) {
        _RETURN_IF_ERROR(graphics_handler(ctx, argv, argc, &cur_flag));
    }
    /*------------------------------------------------------------------------*/
    /* Checking for strage errors                                             */
    else if(argc < 1) {
        err_msg("Unexpected amount of flags: %d\n", argc);
        return STATE_UNEXPECTED_PARAMETER;
    }
    /*------------------------------------------------------------------------*/
    /* Parsing flags                                                          */
    else {
        /*--------------------------------------------------------------------*/
        /* Running through flags                                              */
        while(cur_flag < argc) {
            /*----------------------------------------------------------------*/
            bool flag_found = false;
            /*----------------------------------------------------------------*/
            /* Running through known flags                                    */
            for(size_t i = 0; i < FlagHandlersSize; i++) {
                /*------------------------------------------------------------*/
                /* Checking if found particular flag                          */
                if(strcmp(FlagHandlers[i].long_name,  argv[cur_flag]) == 0 ||
                   strcmp(FlagHandlers[i].short_name, argv[cur_flag]) == 0) {
                    /*--------------------------------------------------------*/
                    /* Running handler                                        */
                    _RETURN_IF_ERROR(FlagHandlers[i].handler(ctx,
                                                             argv,
                                                             argc,
                                                             &cur_flag));
                    /*--------------------------------------------------------*/
                    /* Updating flags index                                   */
                    cur_flag++;
                    /*--------------------------------------------------------*/
                    /* Setting found flag on true and going out               */
                    flag_found = true;
                    break;
                    /*--------------------------------------------------------*/
                }
                /*------------------------------------------------------------*/
            }
            /*----------------------------------------------------------------*/
            /* Checking that flag was found                                   */
            if(!flag_found) {
                err_msg("Unknown flag: %s\n", argv[cur_flag]);
                return STATE_UNKNOWN_FLAG;
            }
            /*----------------------------------------------------------------*/
        }
        /*--------------------------------------------------------------------*/
    }
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/

err_state_t graphics_handler(ctx_t      *ctx,
                             const char */*argv*/[],
                             int         /*argc*/,
                             int        */*current_flag*/) {
    /*------------------------------------------------------------------------*/
    /* Setting context for graphics                                           */
    ctx->testing_mode = false;
    ctx->offset_x   = DefaultOffsetX;
    ctx->offset_y   = DefaultOffsetY;
    ctx->scale      = DefaultScale;
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/

err_state_t test_handler(ctx_t      *ctx,
                         const char *argv[],
                         int         argc,
                         int        *current_flag) {
    /*------------------------------------------------------------------------*/
    /* Checking if additional parameters appered                              */
    if(*current_flag + 3 >= argc) {
        err_msg("Flag '%s' expected to have 3 parameter:\n"
                "\t- minimal iterations per point\n"
                "\t- step\n"
                "\t- steps number\n",
                argv[*current_flag]);
        return STATE_UNEXPECTED_PARAMETER;
    }
    /*------------------------------------------------------------------------*/
    /* Reading additional parameters and setting context for testing          */
    ctx->iters_point_min = strtoul(argv[*current_flag + 1], NULL, 10);
    ctx->step_value      = strtoul(argv[*current_flag + 2], NULL, 10);
    ctx->steps_num       = strtoul(argv[*current_flag + 3], NULL, 10);
    *current_flag += 3;
    ctx->testing_mode = true;
    ctx->offset_x   = DefaultOffsetX;
    ctx->offset_y   = DefaultOffsetY;
    ctx->scale      = DefaultScale;
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/

err_state_t output_handler(ctx_t      *ctx,
                           const char *argv[],
                           int         argc,
                           int        *current_flag) {
    /*------------------------------------------------------------------------*/
    /* Checking that additional parameter appeared                            */
    if(*current_flag + 1 >= argc) {
        err_msg("Flag '%s' expected to have 1 parameter: output file\n",
                argv[*current_flag]);
    }
    /*------------------------------------------------------------------------*/
    /* Reading additional parameter                                           */
    ctx->output_file = argv[*current_flag + 1];
    *current_flag += 1;
    /*------------------------------------------------------------------------*/
    return STATE_SUCCESS;
}

/*============================================================================*/
