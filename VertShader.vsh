#version 330 core
uniform mat4 matWVP;
uniform float LightScale; // Scale the diffuse color

layout (location = 0) in vec4 Position;
layout (location = 1) in vec4 Normal;
layout (location = 2) in vec4 Color;
layout (location = 3) in vec2 TexCoord;
out     vec2 oTexCoord;
out     vec4 oColor;

// coefficients from spherical harmonics
uniform vec3 L00 ;
uniform vec3 L1m1;
uniform vec3 L10 ;
uniform vec3 L11 ;
uniform vec3 L2m2;
uniform vec3 L2m1;
uniform vec3 L20 ;
uniform vec3 L21 ;
uniform vec3 L22 ;

const float C1 = 0.429043;
const float C2 = 0.511664;
const float C3 = 0.743125;
const float C4 = 0.886227;
const float C5 = 0.247708;

void main()
{
   vec4 tnorm    = Normal;
   vec3 DiffuseColor =  C1 * L22 * (tnorm.x * tnorm.x - tnorm.y * tnorm.y) +
                   C3 * L20 * tnorm.z * tnorm.z +
                   C4 * L00 -
                   C5 * L20 +
                   2.0 * C1 * L2m2 * tnorm.x * tnorm.y +
                   2.0 * C1 * L21  * tnorm.x * tnorm.z +
                   2.0 * C1 * L2m1 * tnorm.y * tnorm.z +
                   2.0 * C2 * L11  * tnorm.x +
                   2.0 * C2 * L1m1 * tnorm.y +   
                   2.0 * C2 * L10  * tnorm.z;
    
   DiffuseColor *= LightScale;
	
   oTexCoord   = TexCoord;
   
   //oColor.rgb  = pow(Color.rgb, vec3(2.2)); 
   //oColor.a    = Color.a;
   oColor.rgb = DiffuseColor;
   oColor.a = 1.0;
   
   gl_Position = (matWVP * Position);
}