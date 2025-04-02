/*============================================================================*/
#ifndef MANDELBROT_H
#define MANDELBROT_H
/*============================================================================*/
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
/*============================================================================*/

enum err_state_t {
    STATE_SUCCESS               = 0,
    STATE_UNEXPECTED_PARAMETER  = 1,
    STATE_MEMORY_ERROR          = 2,
    STATE_UNKNOWN_FLAG          = 3,
    STATE_EXIT_SUCCESS          = 4,
    STATE_FILE_OPENING_ERROR    = 5,
};

/*============================================================================*/

#define _RETURN_IF_ERROR(...) {                     \
    err_state_t _error_code = (__VA_ARGS__);        \
    if(_error_code != STATE_SUCCESS) {              \
        return _error_code;                         \
    }                                               \
}                                                   \
/*============================================================================*/

static const float  DefaultScale    = 3.f;
static const float  DefaultOffsetX  = 0.5f;
static const float  DefaultOffsetY  = 0.f;
static const size_t DefaultIters    = 10;
static const size_t FPSBufferSize   = 32;

/*============================================================================*/

struct ctx_t {
    sf::RenderWindow            window;
    sf::Texture                 texture;
    sf::Sprite                  box;
    bool                        testing_mode;
    size_t                      iters_point_min;
    size_t                      step_value;
    size_t                      steps_num;
    const char                 *output_file;
    float                       offset_x;
    float                       offset_y;
    float                       scale;
    alignas(16) sf::Uint32     *image;
    sf::Text                    fps_text;
    char                        fps_buffer[FPSBufferSize];
    double                      fps;
    sf::Font                    font;
};

/*============================================================================*/

int log_msg(const char *format, ...);
int err_msg(const char *format, ...);

/*============================================================================*/
#endif
/*============================================================================*/
