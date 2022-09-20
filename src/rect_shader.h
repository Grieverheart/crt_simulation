namespace rect_shader
{

static const char* fs = R"(
#version 330
out vec4 outColor;
uniform vec3 inColor;
//uniform vec4 coords;
//smooth in vec2 local_position;

void main(void)
{
    //float a = min(1 - local_position.x, local_position.x) * coords.z / 0.1;
    //float b = min(1 - local_position.y, local_position.y) * coords.w / 0.1;
    //outColor = vec4(0.5 * (2 - smoothstep(0.03, 0.05, min(a, b))) * inColor, 1);
    outColor = vec4(inColor, 1.0);
}
)";

//
// 1_____3
// |     |
// |     |
// |_____|
// 0     2
//

static const char* vs = R"(
#version 330

uniform vec4 coords;
//smooth out vec2 local_position;

void main(void)
{
    vec2 TexCoord;
    TexCoord.x = (gl_VertexID < 2)? 0.0: 1.0;
    TexCoord.y = (gl_VertexID % 2 == 0)? 0.0: 1.0;

    vec2 position = 2.0 * (TexCoord * coords.zw + coords.xy) - 1.0;

    //local_position = TexCoord;
    gl_Position = vec4(position, 0.0, 1.0);
}
)";

};
