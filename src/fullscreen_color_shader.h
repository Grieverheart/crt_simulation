namespace fullscreen_color_shader
{

static const char* fs = R"(
#version 330

uniform vec4 color;
out vec4 outColor;

void main(void){
    outColor = color;
}
)";

static const char* vs = R"(
#version 330

void main(void){

    vec2 TexCoord;
    TexCoord.x = (gl_VertexID == 2)? 2.0: 0.0;
    TexCoord.y = (gl_VertexID == 1)? 2.0: 0.0;

    gl_Position = vec4(2.0 * TexCoord - 1.0, 0.0, 1.0);
}
)";

};
