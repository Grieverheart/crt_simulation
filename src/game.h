#pragma once

#include <cstdint>

struct Game;

struct GameVars
{
    bool fullscreen;
    int window_width, window_height;
};

Game* game_create(GameVars vars);
void game_destroy(Game* game);
void game_window_resize(Game* game, int width, int height);
void game_update_and_render(Game* game, double frame_sec);
void game_pause(Game* game, bool pause_or_unpause);
void game_quit(Game* game);
bool game_is_running(const Game* game);
bool game_fullscreen_get(const Game* game);

