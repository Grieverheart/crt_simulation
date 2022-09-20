
namespace beam_shader
{

static const char* fs = R"(
#version 330

#define SQRT2 1.4142135623730951
#define SQRT2PI 2.506628274631001

float erf(float x) {
    float s = sign(x), a = abs(x);
    x = 1.0 + (0.278393 + (0.230389 + 0.078108 * (a * a)) * a) * a;
    x *= x;
    return s - s / (x * x);
}

uniform vec3 beam_color;
uniform float beam_radius;
uniform float beam_dt;
uniform float beam_intensity;
uniform float beam_decay_time;
uniform uint beam_num_edges;

out vec4 outColor;

flat in float path_length;
flat in int edge_id;
smooth in vec2 beam_coord;

void main(void)
{
    float sigma = beam_radius / 5.0;
    float total_factor = beam_intensity * beam_dt;
    if(path_length > 1e-5)
    {
        float f = beam_dt * sigma / (SQRT2 * path_length * beam_decay_time);
        total_factor *= erf(f + beam_coord.x / (SQRT2 * sigma)) - erf(f + (beam_coord.x - path_length) / (SQRT2 * sigma));
        total_factor *= exp(
            f * f - beam_coord.y * beam_coord.y / (2.0 * sigma * sigma) -
            beam_dt * (float(beam_num_edges) - float(edge_id) - beam_coord.x / path_length) / beam_decay_time
        ) / (2.0 * path_length);
    }
    else
    {
        total_factor *= exp(-dot(beam_coord, beam_coord) / (2.0 * sigma * sigma)) / (SQRT2PI * sigma);
    }
    outColor = vec4(total_factor * beam_color, 1.0);
}
)";

static const char* vs = R"(
#version 330

layout(location = 0) in vec2 start;
layout(location = 1) in vec2 end;

// 3--2
// | /|
// |/ |
// 0--1

// @todo: Convert to a struct or something.
uniform float beam_radius;
uniform vec2 aspect_ratio_correction;

flat out float path_length;
flat out int edge_id;
smooth out vec2 beam_coord;

void main(void)
{
    int vid = gl_VertexID % 4;
    path_length = length(end - start);

    vec2 dir = beam_radius * ((path_length > 1e-5)? (end - start) / path_length: vec2(1.0, 0.0));

    vec2 orth = vec2(-dir.y, dir.x);
    vec2 pos = ((vid < 2)? start - dir: end + dir) +
               ((vid == 0 || vid == 3)? orth: -orth);

    beam_coord.x = ((vid < 2)? -beam_radius: path_length + beam_radius);
    beam_coord.y = (vid == 0 || vid == 3)? beam_radius: -beam_radius;

    edge_id = gl_VertexID / 4;

    gl_Position = vec4(pos * aspect_ratio_correction, 0.0, 1.0);
}
)";

};

