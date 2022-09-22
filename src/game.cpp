#include "game.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <SDL2/SDL.h>
#include "common.h"

#include "beam.h"
#include "renderer.h"
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#include "stb_image.h"

struct Game
{
    bool is_running;

    double base_frequency;

    bool is_paused;

    bool    fullscreen;

    Patch    patch;
    Renderer renderer;
    Beam     beam;
};


void game_update_and_render(Game* game, double frame_sec)
{
    if(game->is_paused) return;

    BeamData beam_data;
    beam_simulate(&game->beam, &beam_data, &game->patch, game->base_frequency, frame_sec);
    gl_renderer_draw_beam_points(game->renderer, game->beam, beam_data);
}

Patch butterfly(void)
{
    Patch phi = pt_fmod(pt_time(), 1.0 / pt_frequency()) * pt_frequency();

    Patch a = 0.25 * pt_lift1(sin, cos, nullptr, (24.0 * M_PI) * phi);
    Patch p = a * (
        pt_lift1(exp, pt_lift1(cos, (24.0 * M_PI) * phi)) - 
        2.0 * pt_lift1(cos, (96.0 * M_PI) * phi) - 
        pt_pow(pt_lift1(sin, (2.0 * M_PI) * phi), 5)
    );

    return p;
}

static inline int sqr(int x)
{
    return x*x;
}

Patch crt(size_t num_scanlines, double refresh_rate)
{
    int w,h,n;
    unsigned char *image_data = stbi_load("../image.jpg", &w, &h, &n, 0);

    Patch signal;
    {
        size_t idx = 0;
        size_t n_points = 240 * 240;
        Patch* patches = (Patch*) malloc(sizeof(Patch) * (n_points - 1));
        //double last_intensity = 0.0;
        for(size_t i = 0; i < 240; ++i)
        {
            for(size_t j = 0; j < 240; ++j)
            {
                double intensity = image_data[240*(239-i)+j] / 255.0;
                // @todo: The are better approaches than using pt_line.
                patches[idx++] = pt_point(0, 0, intensity);//pt_line(0.0, 0.0, last_intensity, 0.0, 0.0, intensity);
                //last_intensity = intensity;
            }
        }
        signal = pt_seq(n_points - 1, patches);
    }

    stbi_image_free(image_data);

    auto phaser = 2.0 * pt_fmod(pt_time(), 1.0 / pt_frequency()) * pt_frequency() - 1.0;
    Patch p = 
    /* Horizontal deflector */ pt_point(1.0, 0.0, 0.0) * pt_fset(phaser, num_scanlines * refresh_rate) +
    /*  Vertical deflector  */ pt_point(0.0, 1.0, 0.0) * (pt_lift1(floor, 0.5 * num_scanlines * pt_fset(phaser, refresh_rate)) + 0.5) * (2.0 / num_scanlines) +
    ///*         Signal       */ pt_point(0.0, 0.0, 1.0);// * (pt_fset(phaser, (num_scanlines / 2) * refresh_rate) > 0.0);
    /*         Signal       */ pt_fset(signal, refresh_rate);// * (pt_fset(phaser, (num_scanlines / 2) * refresh_rate) > 0.0);
    // @todo: Perhaps use line_strip for generating the signal.
    return p;
}

void crt2_call(void* userdata, double t, double f, double* x, double* y, double* z)
{
    double* intensities = (double*) userdata;
    double m = 60.0 * 240.0 * 240.0 * t;
    double a = fmod(m, 240.0);
    int b = (int(m) / 240) % 240;
    *x = (a - 120.0) / 120.0;
    *y = (b - 120.0) / 120.0;
    *z = intensities[240 * b + int(a)];
}

Patch crt2(size_t num_scanlines, double refresh_rate)
{
    int w,h,n;
    unsigned char *image_data = stbi_load("../image.jpg", &w, &h, &n, 0);

    double* intensities = (double*) malloc(sizeof(double) * 240*240);
    {
        for(size_t i = 0; i < 240; ++i)
            for(size_t j = 0; j < 240; ++j)
                intensities[240 * i + j] = pow(image_data[240*(239-i)+j] / 255.0, 2.2);
    }

    stbi_image_free(image_data);

    Patch p;
    p.call = crt2_call;
    p.type = 2;
    p.userdata = (void*)intensities;
    return p;
}

Game* game_create(GameVars vars)
{
    Game* game = new Game;
    {
        game->is_running     = true;
        game->base_frequency = 60.0;
        game->is_paused      = false;
        game->fullscreen     = vars.fullscreen;

        // Beam setup
        {
            game->beam.num_edges  = 60000;
            game->beam.decay_time = 4e-1;
            game->beam.radius     = 1.0e-2;

            game->beam.sim_time = 0.0;

            game->beam.x = 0.0;
            game->beam.y = 0.0;
            game->beam.z = 0.0;

            game->beam.color[0]  = 1.0f;
            game->beam.color[1]  = 1.0f;
            game->beam.color[2]  = 1.0f;

            game->beam.intensity = 1500.0;
        }

        game->patch = crt2(240, 60.0);
    }

    // Create renderer struct
    if(!gl_renderer_init(&game->renderer, vars.window_width, vars.window_height)) return nullptr;
    gl_renderer_set_beam_parameters(game->renderer, game->beam);

    return game;
}

void game_destroy(Game* game)
{
    delete game;
}

void game_window_resize(Game* game, int width, int height)
{
    gl_renderer_resize(&game->renderer, width, height);

    // Keep beam at fixed size regardless of resolution.
    game->beam.radius = 1.0e-2 * 700.0 / fmin(width, height);
    gl_renderer_set_beam_parameters(game->renderer, game->beam);
}

bool game_is_running(const Game* game)
{
    return game->is_running;
}

void game_quit(Game* game)
{
    game->is_running = false;
}

bool game_fullscreen_get(const Game* game)
{
    return game->fullscreen;
}

void game_pause(Game* game, bool pause_or_unpause)
{
    game->is_paused = pause_or_unpause;
}
