#include "game.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <SDL2/SDL.h>

#include "beam.h"
#include "renderer.h"
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#include "stb_image.h"

struct Game
{
    bool     is_running;
    double   base_frequency;
    bool     is_paused;
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

void crt_call(void* userdata, double t, double f, double* x, double* y, double* z)
{
    double* intensities = (double*) userdata;
    double m = 60.0 * 240.0 * t;
    double a = m - int(m);
    int b = int(m) % 240;
    *x = 2.0 * (a - 0.5);
    *y = -(b - 120.0) / 120.0;
    *z = intensities[240 * b + int(a*240.0)];
}

Patch crt(size_t num_scanlines, double refresh_rate)
{
    int w,h,n;
    unsigned char *image_data = stbi_load("../image.jpg", &w, &h, &n, 0);

    double* intensities = (double*) malloc(sizeof(double) * 240*240);
    {
        for(size_t i = 0; i < 240; ++i)
            for(size_t j = 0; j < 240; ++j)
                intensities[240 * i + j] = pow(image_data[240*i+j] / 255.0, 2.2);
    }

    stbi_image_free(image_data);

    Patch p;
    p.call = crt_call;
    p.userdata = (void*)intensities;
    return p;
}

Game* game_create(int window_width, int window_height)
{
    Game* game = new Game;
    {
        game->is_running     = true;
        game->base_frequency = 60.0;
        game->is_paused      = false;

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

        game->patch = crt(240, 60.0);
    }

    // Create renderer struct
    if(!gl_renderer_init(&game->renderer, window_width, window_height)) return nullptr;
    gl_renderer_set_beam_parameters(game->renderer, game->beam);

    return game;
}

void game_destroy(Game* game)
{
    delete game;
}

void game_window_resize(Game* game, int drawable_width, int drawable_height, int width, int height)
{
    gl_renderer_resize(&game->renderer, drawable_width, drawable_height);

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

void game_pause(Game* game, bool pause_or_unpause)
{
    game->is_paused = pause_or_unpause;
}
