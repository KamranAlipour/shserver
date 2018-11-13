#version 330 core
uniform sampler2D Texture0;
in      vec4      oColor;
in      vec2      oTexCoord;
out     vec4      FragColor;
void main()
{
   FragColor = oColor * texture2D(Texture0, oTexCoord);
}