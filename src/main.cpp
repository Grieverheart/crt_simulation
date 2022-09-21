#include <SDL2/SDL.h>
#include <GL/glew.h>
#include <SDL2/SDL_opengl.h>

#include "common.h"
#include "game.h"

int main(int argc, char* argv[])
{
    // Create window and gl context
    int window_width = 700;
    int window_height = 700;

    GameVars game_vars
    {
        false,
        window_width, window_height
    };

    SDL_Window* window;
    SDL_GLContext gl_context;
    {
        if(SDL_Init(SDL_INIT_VIDEO) < 0)
        {
            fprintf(stderr, "Error initializing SDL. SDL_Error: %s\n", SDL_GetError());
            return -1;
        }


        // Use OpenGL 3.1 core
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
        SDL_GL_SetAttribute(SDL_GL_FRAMEBUFFER_SRGB_CAPABLE, 1);

        window = SDL_CreateWindow(
            "Vectron",
            SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
            window_width, window_height,
            SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL | SDL_WINDOW_ALLOW_HIGHDPI |
            (game_vars.fullscreen? SDL_WINDOW_FULLSCREEN_DESKTOP: SDL_WINDOW_RESIZABLE)
        );

        if(!window)
        {
            fprintf(stderr, "Error creating SDL window. SDL_Error: %s\n", SDL_GetError());
            SDL_Quit();
            return -1;
        }

        gl_context = SDL_GL_CreateContext(window);
        if(!gl_context)
        {
            fprintf(stderr, "Error creating SDL GL context. SDL_Error: %s\n", SDL_GetError());
            SDL_DestroyWindow(window);
            SDL_Quit();
            return -1;
        }

        int r, g, b, a, m, s;
        SDL_GL_GetAttribute(SDL_GL_RED_SIZE, &r);
        SDL_GL_GetAttribute(SDL_GL_GREEN_SIZE, &g);
        SDL_GL_GetAttribute(SDL_GL_BLUE_SIZE, &b);
        SDL_GL_GetAttribute(SDL_GL_ALPHA_SIZE, &a);
        SDL_GL_GetAttribute(SDL_GL_MULTISAMPLESAMPLES, &m);
        SDL_GL_GetAttribute(SDL_GL_FRAMEBUFFER_SRGB_CAPABLE, &s);

        SDL_GL_SetSwapInterval(1);
    }

    SDL_GL_GetDrawableSize(window, &game_vars.window_width, &game_vars.window_height);
    Game* game = game_create(game_vars);
    if(!game)
    {
        SDL_GL_DeleteContext(gl_context);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return -1;
    }


    uint64_t cpu_frequency = SDL_GetPerformanceFrequency();
    uint64_t current_clock = SDL_GetPerformanceCounter();

    bool paused = false;
    while(game_is_running(game))
    {
        // Process input
        SDL_Event sdl_event;
        while(SDL_PollEvent(&sdl_event) != 0)
        {
            if(sdl_event.type == SDL_KEYDOWN)
            {
                switch(sdl_event.key.keysym.sym){
                case SDLK_p:
                    paused = !paused;
                    game_pause(game, paused);
                    break;
                case SDLK_ESCAPE:
                    game_quit(game);
                    break;
                default:
                    break;
                }
            }
            else if(sdl_event.type == SDL_QUIT)
            {
                game_quit(game);
            }
            else if(sdl_event.type == SDL_WINDOWEVENT)
            {
                if(sdl_event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
                {
                    int drawable_width, drawable_height;
                    SDL_GL_GetDrawableSize(window, &drawable_width, &drawable_height);
                    game_window_resize(game, drawable_width, drawable_height);
                }
                else if(sdl_event.window.event == SDL_WINDOWEVENT_FOCUS_LOST)
                {
                    game_pause(game, true);
                }
                else if(sdl_event.window.event == SDL_WINDOWEVENT_FOCUS_GAINED)
                {
                    game_pause(game, false);
                }
            }
        }

        int drawable_width, drawable_height;
        SDL_GL_GetDrawableSize(window, &drawable_width, &drawable_height);

        uint64_t new_clock = SDL_GetPerformanceCounter();
        double frame_sec = double(new_clock - current_clock) / cpu_frequency;
        current_clock = new_clock;

        // Make sure we don't use a very small frame duration
        // after resuming.
        //frame_sec = fmax(1.0 / 200.0, frame_sec);
        game_update_and_render(game, frame_sec);

        SDL_GL_SwapWindow(window);
    }

    game_destroy(game);
    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
