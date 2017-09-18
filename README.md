## GLLA ##


GLLA (GL Linear Algebra) is a small C11 static library that provides types and functions for a 3-vector, a 3x3-matrix, and a 4x4 affine matrix. It is written to be exceedingly clear in both implementation and usage.

#### vec3 ####
This type is a 3-component float vector, using Clang's vector extensions.

Some of the usual vec3 linear algebra functions which Clang does not provide have been implemented.

#### mat3 ####
This type is a struct wrapper around an array of 3 vec3 vectors, to allow for easy assignment.

Some of the usual mat3 linear algebra functions have been implemented.

#### amat4 ####
This type is a struct wrapper around a mat3 (upper-left 3x3 matrix) and a vec3 (first three components of the last column).

Since it represents an affine 4x4 matrix, all provided amat4 functions behave as if the last row is <0, 0, 0, 1>, without actually carrying that data around. Where possible, operations which would multiply by zero have been exploited to reduce computation.

Access the upper-left inner mat3 with .a

Access the first three elements of the final column as a vec3 with .t

### Example Usage ###

Using Clang's vector extensions, components can be accessed similar to GLSL:

	vec3 v1 = {1.0, 2.0, 3.0};
	vec3 v2 = {v1.x * 4, v1.y, v1.z};

GLLA does not use any operator overloading, except as provided by Clang's vector extensions. Algebra is performed by calling functions on the defined types, passed by value. Functions return by value as well.

	vec3 v3 = vec3_lerp(v1, v2, 0.5);
	vec3 normal = vec3_cross(v3 - v1, v2 - v1);
	normal = vec3_normalize(normal);
    
For computer graphics, you can quickly and efficiently compose affine frames by assigning to the rotation and translation components independently.

	vec3 eye_target = {0.0, 0.0, 0.0};
	amat4 eye_frame = {.t = (vec3){0.0, 1.0, 5.0}};
	eye_frame.a = mat3_lookat(eye_frame.t, eye_target, (vec3){0, 1, 0});
	
### Compiling and Installing ###
The library is just one header file and one source file. You can simply copy it into your project if you like. Alternatively, you can use the provided Makefile to compile and install the library into your /usr/local/lib/ and /usr/local/include/ folders with this command:

	make install
	
Then you can use it in any future project by prepending this to your C source files:

	#include "glla.h"

	
and compiling with the flags:

	-L/usr/local/lib/ -lglla
	
### Motivation ###

I wrote a prior linear algebra library, glalgebra, in the course of making a 3D game engine in C (still a work in progress). I did not like the syntax used for other linear algebra libraries, so I made my own, which seems to be a rite of passage among people who write their own game engines.

After using it for some time, I became aware of Clang's vector extensions, which I decided to use for my successor to the original library. This revised library removes complexity I have not found useful with the previous library. GLLA stands for "GL Linear Algebra" or "Gavin's Library for Linear Algebra".
