#pragma once

#include "beam.h"
#include <GL/glew.h>

struct Renderer
{
    GLuint fbo;
    GLuint fbo_texture;

    GLuint fullscreen_triangle_vao;
    GLuint beam_point_vao;
    GLuint beam_point_vbo;
    GLuint beam_point_ibo;

    GLuint fullscreen_color_shader_id;
    GLuint fullscreen_tex_shader_id;
    GLuint beam_shader_id;

    GLint  color_location;
    GLint  tex_location;
    GLint  fade_location;
    GLint  num_edge_location;
    GLint  decay_time_location;
    GLint  beam_dt_location;
    GLint  beam_radius_location;
    GLint  beam_color_location;
    GLint  beam_intensity_location;
    GLint  aspect_ratio_location;

    GLint viewport_x, viewport_y;
    GLsizei viewport_size;
};

// @todo: Split into:
// void gl_init(void);
// gl_beam_renderer_init(BeamRenderer* renderer, size_t width, size_t height);
// ...
// gl_stage_editor_renderer_init(StageEditorRenderer* renderer, size_t width, size_t height)

bool gl_renderer_init(Renderer* renderer, size_t width, size_t height);
bool gl_renderer_resize(Renderer* renderer, size_t width, size_t height);
void gl_renderer_draw_beam_points(const Renderer& renderer, const Beam& beam, const BeamData& beam_data);
void gl_renderer_set_beam_parameters(const Renderer& renderer, const Beam& beam);
