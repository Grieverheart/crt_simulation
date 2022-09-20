#include "game.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <SDL2/SDL.h>
#include "common.h"

#include "beam.h"
#include "renderer.h"

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

    Patch a = 0.25 * pt_lift1(sin, cos, (24.0 * M_PI) * phi);
    Patch p = a * (
        pt_lift1(exp, pt_lift1(cos, (24.0 * M_PI) * phi)) - 
        2.0 * pt_lift1(cos, (96.0 * M_PI) * phi) - 
        pt_pow(pt_lift1(sin, (2.0 * M_PI) * phi), 5)
    );

    return p;
}

Patch crt(size_t num_scanlines, double refresh_rate)
{
    auto phaser = 2.0 * pt_fmod(pt_time(), 1.0 / pt_frequency()) * pt_frequency() - 1.0;
    Patch p = 
    /* Horizontal deflector */ pt_point(1.0, 0.0) * pt_fset(phaser, num_scanlines * refresh_rate) +
    /*  Vertical deflector  */ pt_point(0.0, 1.0) * (pt_lift1(floor, 0.5 * num_scanlines * pt_fset(phaser, refresh_rate)) + 0.5) * (2.0 / num_scanlines);
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
            game->beam.num_edges  = 5000;
            game->beam.decay_time = 4e-2;
            game->beam.radius     = 1.0e-2;

            game->beam.sim_time = 0.0;

            game->beam.x = 0.0;
            game->beam.y = 0.0;

            game->beam.color[0]  = 0.05f;
            game->beam.color[1]  = 1.00f;
            game->beam.color[2]  = 0.05f;

            game->beam.intensity = 500.0;
        }

        game->patch = crt(240, 60.0);
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
