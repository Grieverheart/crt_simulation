namespace fullscreen_tex_shader
{

static const char* fs = R"(
#version 330

uniform sampler2D tex;
noperspective in vec2 TexCoord;

out vec3 outColor;

void main(void){
    outColor = texture(tex, TexCoord).rgb;
}
)";

static const char* vs = R"(
#version 330

noperspective out vec2 TexCoord;

void main(void){

    TexCoord.x = (gl_VertexID == 2)? 2.0: 0.0;
    TexCoord.y = (gl_VertexID == 1)? 2.0: 0.0;

    gl_Position = vec4(2.0 * TexCoord - 1.0, 0.0, 1.0);
}
)";

};
