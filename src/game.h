#pragma once

#include <cstdint>

struct Game;

Game* game_create(int window_width, int window_height);
void game_destroy(Game* game);
void game_window_resize(Game* game, int drawable_width, int drawable_height, int width, int height);
void game_update_and_render(Game* game, double frame_sec);
void game_pause(Game* game, bool pause_or_unpause);
void game_quit(Game* game);
bool game_is_running(const Game* game);

