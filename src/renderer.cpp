#include "renderer.h"

#include <cmath>

#include "gl_utils.h"
#include "fullscreen_color_shader.h"
#include "fullscreen_tex_shader.h"
#include "beam_shader.h"

namespace
{
    void renderer_upload_beam_points(const Renderer& renderer, double* beam_points, size_t num_edges)
    {
        // @todo: Try glMapBuffer.
        float* beam_gl_points = new float[(num_edges + 1) * 4 * 3];
        for(size_t i = 0; i < num_edges + 1; ++i)
            for(size_t j = 0; j < 12; ++j)
                beam_gl_points[12 * i + j] = (float)beam_points[3 * i + (j % 3)];

        glBindBuffer(GL_ARRAY_BUFFER, renderer.beam_point_vbo);
        glBufferData(GL_ARRAY_BUFFER, (num_edges + 1) * 4 * 3 * sizeof(*beam_gl_points), beam_gl_points, GL_STREAM_DRAW);
        delete[] beam_gl_points;

        GLuint* indices = new GLuint[6 * num_edges];
        for(GLuint i = 0; i < (GLuint)num_edges; ++i)
        {
            indices[6 * i + 0] = 4 * i + 0;
            indices[6 * i + 1] = 4 * i + 1;
            indices[6 * i + 2] = 4 * i + 2;
            indices[6 * i + 3] = 4 * i + 0;
            indices[6 * i + 4] = 4 * i + 2;
            indices[6 * i + 5] = 4 * i + 3;
        }

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, renderer.beam_point_ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6 * num_edges * sizeof(*indices), indices, GL_STREAM_DRAW);
        delete[] indices;
    }

} // namespace

bool gl_renderer_init(Renderer* renderer, size_t width, size_t height)
{
    GLenum err = glewInit();
    if(err != GLEW_OK)
    {
        fprintf(stderr, "Error initializing GLEW.\n");
        return false;
    }

    int glVersion[2] = {-1, 1};
    glGetIntegerv(GL_MAJOR_VERSION, &glVersion[0]);
    glGetIntegerv(GL_MINOR_VERSION, &glVersion[1]);

    gl_debug(__FILE__, __LINE__);

    printf("Using OpenGL: %d.%d\n", glVersion[0], glVersion[1]);
    printf("Renderer used: %s\n", glGetString(GL_RENDERER));
    printf("Shading Language: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

    size_t size = width;
    if(width > height)
    {
        size = height;
        renderer->viewport_x = (width - height) / 2;
        renderer->viewport_y = 0;
    }
    else
    {
        renderer->viewport_x = 0;
        renderer->viewport_y = (height - width) / 2;
    }
    renderer->viewport_size = size;

    glGenFramebuffers(1, &renderer->fbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, renderer->fbo);
    glGenTextures(1, &renderer->fbo_texture);
    glBindTexture(GL_TEXTURE_2D, renderer->fbo_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, size, size, 0, GL_RGBA, GL_FLOAT, NULL);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, renderer->fbo_texture, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if(status != GL_FRAMEBUFFER_COMPLETE)
    {
        printf("FB error, status 0x%04x", status);
        return false;
    }

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    // Create vao for generating fullscreen triangle
    glGenVertexArrays(1, &renderer->fullscreen_triangle_vao);


    GLuint fullscreen_color_shader_id = glCreateProgram();
    {
        //Create vertex shader
        GLuint shader_vp = glCreateShader(GL_VERTEX_SHADER);

        glShaderSource(shader_vp, 1, &fullscreen_color_shader::vs, 0);
        glCompileShader(shader_vp);
        validate_shader(shader_vp, fullscreen_color_shader::vs);
        glAttachShader(fullscreen_color_shader_id, shader_vp);

        glDeleteShader(shader_vp);

        //Create fragment shader
        GLuint shader_fp = glCreateShader(GL_FRAGMENT_SHADER);

        glShaderSource(shader_fp, 1, &fullscreen_color_shader::fs, 0);
        glCompileShader(shader_fp);
        validate_shader(shader_fp, fullscreen_color_shader::fs);
        glAttachShader(fullscreen_color_shader_id, shader_fp);

        glDeleteShader(shader_fp);

        glLinkProgram(fullscreen_color_shader_id);

        if(!validate_program(fullscreen_color_shader_id)){
            fprintf(stderr, "Error while validating shader.\n");
            //glDeleteVertexArrays(1, &renderer->fullscreen_triangle_vao);
            return false;
        }
    }
    renderer->fullscreen_color_shader_id = fullscreen_color_shader_id;

    GLuint fullscreen_tex_shader_id = glCreateProgram();
    {
        //Create vertex shader
        GLuint shader_vp = glCreateShader(GL_VERTEX_SHADER);

        glShaderSource(shader_vp, 1, &fullscreen_tex_shader::vs, 0);
        glCompileShader(shader_vp);
        validate_shader(shader_vp, fullscreen_tex_shader::vs);
        glAttachShader(fullscreen_tex_shader_id, shader_vp);

        glDeleteShader(shader_vp);

        //Create fragment shader
        GLuint shader_fp = glCreateShader(GL_FRAGMENT_SHADER);

        glShaderSource(shader_fp, 1, &fullscreen_tex_shader::fs, 0);
        glCompileShader(shader_fp);
        validate_shader(shader_fp, fullscreen_tex_shader::fs);
        glAttachShader(fullscreen_tex_shader_id, shader_fp);

        glDeleteShader(shader_fp);

        glLinkProgram(fullscreen_tex_shader_id);

        if(!validate_program(fullscreen_tex_shader_id)){
            fprintf(stderr, "Error while validating shader.\n");
            //glDeleteVertexArrays(1, &renderer->fullscreen_triangle_vao);
            return false;
        }
    }
    renderer->fullscreen_tex_shader_id = fullscreen_tex_shader_id;

    GLuint beam_shader_id = glCreateProgram();
    {
        //Create vertex shader
        GLuint shader_vp = glCreateShader(GL_VERTEX_SHADER);

        glShaderSource(shader_vp, 1, &beam_shader::vs, 0);
        glCompileShader(shader_vp);
        validate_shader(shader_vp, beam_shader::vs);
        glAttachShader(beam_shader_id, shader_vp);

        glDeleteShader(shader_vp);

        //Create fragment shader
        GLuint shader_fp = glCreateShader(GL_FRAGMENT_SHADER);

        glShaderSource(shader_fp, 1, &beam_shader::fs, 0);
        glCompileShader(shader_fp);
        validate_shader(shader_fp, beam_shader::fs);
        glAttachShader(beam_shader_id, shader_fp);

        glDeleteShader(shader_fp);

        glLinkProgram(beam_shader_id);

        if(!validate_program(beam_shader_id)){
            fprintf(stderr, "Error while validating shader.\n");
            //glDeleteVertexArrays(1, &renderer->fullscreen_triangle_vao);
            return false;
        }
    }
    renderer->beam_shader_id = beam_shader_id;

    GLuint beam_point_vao, beam_point_vbo, beam_point_ibo;
    glGenVertexArrays(1, &beam_point_vao);
    glBindVertexArray(beam_point_vao);
    {
        glGenBuffers(1, &beam_point_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, beam_point_vbo);

        glEnableVertexAttribArray((GLuint)0);
        glVertexAttribPointer((GLuint)0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray((GLuint)1);
        glVertexAttribPointer((GLuint)1, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid*)(12 * sizeof(float)));

        glGenBuffers(1, &beam_point_ibo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, beam_point_ibo);
    }
    glBindVertexArray(0);
    renderer->beam_point_vao = beam_point_vao;
    renderer->beam_point_vbo = beam_point_vbo;
    renderer->beam_point_ibo = beam_point_ibo;


    glUseProgram(fullscreen_color_shader_id);
    renderer->color_location = glGetUniformLocation(fullscreen_color_shader_id, "color");

    glUseProgram(fullscreen_tex_shader_id);
    renderer->tex_location = glGetUniformLocation(fullscreen_tex_shader_id, "tex");

    glUseProgram(beam_shader_id);
    renderer->num_edge_location       = glGetUniformLocation(beam_shader_id, "beam_num_edges");
    renderer->decay_time_location     = glGetUniformLocation(beam_shader_id, "beam_decay_time");
    renderer->beam_dt_location        = glGetUniformLocation(beam_shader_id, "beam_dt");
    renderer->beam_radius_location    = glGetUniformLocation(beam_shader_id, "beam_radius");
    renderer->beam_color_location     = glGetUniformLocation(beam_shader_id, "beam_color");
    renderer->beam_intensity_location = glGetUniformLocation(beam_shader_id, "beam_intensity");
    renderer->aspect_ratio_location   = glGetUniformLocation(beam_shader_id, "aspect_ratio_correction");

    float aspect_ratio_correction[] = {1, 1};
    //if(width > height) aspect_ratio_correction[0] = float(height) / width;
    //else aspect_ratio_correction[1] = float(width) / height;
    glUniform2fv(renderer->aspect_ratio_location, 1, aspect_ratio_correction);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_BLEND);

    glClear(GL_COLOR_BUFFER_BIT);

    return true;
}

bool gl_renderer_resize(Renderer* renderer, size_t width, size_t height)
{
    size_t size = width;
    if(width > height)
    {
        size = height;
        renderer->viewport_x = (width - height) / 2;
        renderer->viewport_y = 0;
    }
    else
    {
        renderer->viewport_x = 0;
        renderer->viewport_y = (height - width) / 2;
    }
    renderer->viewport_size = size;

    glBindTexture(GL_TEXTURE_2D, renderer->fbo_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, size, size, 0, GL_RGBA, GL_FLOAT, NULL);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, renderer->fbo);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glUseProgram(renderer->beam_shader_id);

    float aspect_ratio_correction[] = {1, 1};
    //if(width > height) aspect_ratio_correction[0] = float(height) / width;
    //else aspect_ratio_correction[1] = float(width) / height;
    //glUniform2fv(renderer->aspect_ratio_location, 1, aspect_ratio_correction);

    glUseProgram(0);

    return true;
}

void gl_renderer_draw_beam_points(const Renderer& renderer, const Beam& beam, const BeamData& beam_data)
{
    glViewport(0, 0, renderer.viewport_size, renderer.viewport_size);

    GLuint queries[3] = {};
    glGenQueries(3, queries);

    glEnable(GL_BLEND);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    renderer_upload_beam_points(renderer, beam_data.points, beam.num_edges);
    {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, renderer.fbo);

        glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
        glUseProgram(renderer.fullscreen_color_shader_id);

        double decay_factor = exp(-beam_data.dt * beam.num_edges / beam.decay_time);
        glUniform4f(renderer.color_location, 0.0f, 0.0f, 0.0f, decay_factor);
        glBindVertexArray(renderer.fullscreen_triangle_vao);

        // Multiply previous frame with the decay factor
        glBeginQuery(GL_TIME_ELAPSED, queries[0]);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glEndQuery(GL_TIME_ELAPSED);

        glBlendFunc(GL_ONE, GL_ONE);
        glUseProgram(renderer.beam_shader_id);
        glUniform1f(renderer.beam_dt_location, beam_data.dt);
        glBindVertexArray(renderer.beam_point_vao);

        // Draw the beam trail
        glBeginQuery(GL_TIME_ELAPSED, queries[1]);
        glDrawElements(GL_TRIANGLES, 6 * beam.num_edges, GL_UNSIGNED_INT, 0);
        glEndQuery(GL_TIME_ELAPSED);

        glDisable(GL_BLEND);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glViewport(renderer.viewport_x, renderer.viewport_y, renderer.viewport_size, renderer.viewport_size);
        glClear(GL_COLOR_BUFFER_BIT);
        glBindTexture(GL_TEXTURE_2D, renderer.fbo_texture);
        glUseProgram(renderer.fullscreen_tex_shader_id);

        glEnable(GL_FRAMEBUFFER_SRGB);
        glBindVertexArray(renderer.fullscreen_triangle_vao);

        // Draw the result into the framebuffer
        glBeginQuery(GL_TIME_ELAPSED, queries[2]);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glEndQuery(GL_TIME_ELAPSED);

        glDisable(GL_FRAMEBUFFER_SRGB);
    }

    glDeleteQueries(3, queries);
}

void gl_renderer_set_beam_parameters(const Renderer& renderer, const Beam& beam)
{
    glUseProgram(renderer.beam_shader_id);
    glUniform1ui(renderer.num_edge_location, beam.num_edges);
    glUniform1f(renderer.beam_radius_location, beam.radius);
    glUniform1f(renderer.beam_intensity_location, beam.intensity);
    glUniform3fv(renderer.beam_color_location, 1, beam.color);
    glUniform1f(renderer.decay_time_location, beam.decay_time);
}
