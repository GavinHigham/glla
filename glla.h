#ifndef GLLA_H
#define GLLA_H
#include <inttypes.h>

typedef float vec2 __attribute__((ext_vector_type(2)));
typedef float vec3 __attribute__((ext_vector_type(3)));
typedef float vec4 __attribute__((ext_vector_type(4)));
typedef double dvec3 __attribute__((ext_vector_type(3)));
typedef int64_t qvec3 __attribute__((ext_vector_type(3)));
typedef int32_t ivec3 __attribute__((ext_vector_type(3)));
typedef int16_t svec3 __attribute__((ext_vector_type(3)));

typedef struct matrix3 {
	vec3 rows[3];
} mat3;

typedef struct affine_matrix4 {
	mat3 a;
	vec3 t;
} amat4;

#define MAT3_IDENT {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}
#define AMAT4_IDENT {MAT3_IDENT, {0.0, 0.0, 0.0}}

//Return the dot product of a and b.
float vec2_dot(vec2 a, vec2 b);
//Return the sum of the components of a.
float vec2_sum(vec2 a);

//Unpack a vec3, which might be larger than 3 floats because of alignment, into an array of 3 floats.
void vec3_unpack(float out[3], vec3 in);
//Returns a new vector pointing in the same direction which has been normalized (magnitude set to 1.0)
vec3 vec3_normalize(vec3 a) __attribute__ ((const));
//Returns a new vector that represents the cross product of a and b.
vec3 vec3_cross(vec3 a, vec3 b) __attribute__ ((const));
//Returns a new vector that is linearly interpolated between a and by by parameter alpha.
//(a*alpha + b*(1 - alpha)
vec3 vec3_lerp(vec3 a, vec3 b, float alpha) __attribute__ ((const));
//Return the dot product of a and b.
float vec3_dot(vec3 a, vec3 b) __attribute__ ((const));
//Return the magnitude of a.
float vec3_mag(vec3 a) __attribute__ ((const));
//Return the distance between a and b.
float vec3_dist(vec3 a, vec3 b) __attribute__ ((const));
//Return the sum of the components of a.
float vec3_sum(vec3 a) __attribute__ ((const));
//Prints a vec3 like so: "{x, y, z}" (no newline).
void vec3_print(vec3 a);
//Prints a vec3 like so: "{x, y, z}" (with newline).
void vec3_println(vec3 a);
//Prints a vec3 like so: "{x, y, z}" (no newline). Takes a printf format for printing each float.
void vec3_printf(char *fmt, vec3 a);

//Unpack a dvec3, which might be larger than 3 doubles because of alignment, into an array of 3 doubles.
void dvec3_unpack(double out[3], dvec3 in);
//Returns a new vector pointing in the same direction which has been normalized (magnitude set to 1.0)
dvec3 dvec3_normalize(dvec3 a) __attribute__ ((const));
//Returns a new vector that represents the cross product of a and b.
dvec3 dvec3_cross(dvec3 a, dvec3 b) __attribute__ ((const));
//Returns a new vector that is linearly interpolated between a and by by parameter alpha.
//(a*alpha + b*(1 - alpha)
dvec3 dvec3_lerp(dvec3 a, dvec3 b, double alpha) __attribute__ ((const));
//Return the dot product of a and b.
float dvec3_dot(dvec3 a, dvec3 b) __attribute__ ((const));
//Return the magnitude of a.
float dvec3_mag(dvec3 a) __attribute__ ((const));
//Return the distance between a and b.
float dvec3_dist(dvec3 a, dvec3 b) __attribute__ ((const));
//Return the sum of the components of a.
float dvec3_sum(dvec3 a) __attribute__ ((const));
//Return the sum of the components of a.
double dvec3_sumd(dvec3 a) __attribute__ ((const));
//Prints a dvec3 like so: "{x, y, z}" (no newline).
void dvec3_print(dvec3 a);
//Prints a dvec3 like so: "{x, y, z}" (no newline). Takes a printf format for printing each float.
void dvec3_printf(char *fmt, dvec3 a);

//Return the magnitude of a as a float.
float ivec3_magf(ivec3 a) __attribute__ ((const));
//Return the distance between a and b.
float ivec3_distf(ivec3 a, ivec3 b) __attribute__ ((const));
//Return the sum of the components of a.
int64_t ivec3_sumd(ivec3 a) __attribute__ ((const));

//Return the magnitude of a as a float.
float qvec3_magf(qvec3 a) __attribute__ ((const));
//Return the magnitude of a as a double.
double qvec3_magd(qvec3 a) __attribute__ ((const));
//Return the distance between a and b.
float qvec3_distf(qvec3 a, qvec3 b) __attribute__ ((const));
//Return the distance between a and b.
double qvec3_distd(qvec3 a, qvec3 b) __attribute__ ((const));
//Return the sum of the components of a. Could overflow!
int64_t qvec3_sum(qvec3 a) __attribute__ ((const));
//Return the sum of the components of a. Could lose precision!
double qvec3_sumd(qvec3 a) __attribute__ ((const));
//Prints a qvec3 origin like so: "[x, y, z]" (no newline).
void qvec3_print(qvec3 b);
//Prints a qvec3 origin like so: "[x, y, z]".
void qvec3_println(qvec3 b);
//Prints a qvec3 origin like so: "[x, y, z]" (no newline). Takes a printf format for printing each int64.
void qvec3_printf(char *fmt, qvec3 b);

//Unpack a svec3, which might be larger than 3 int16_t because of alignment, into an array of 3 int16_t.
void svec3_unpack(int16_t out[3], svec3 in);

//Return the dot product of a and b.
float vec4_dot(vec4 a, vec4 b);
//Return the sum of the components of a.
float vec4_sum(vec4 a);

//Create a new mat3 from an array of floats. Row-major order.
mat3 mat3_from_array(float *array) __attribute__ ((const));
//3x3 Identity matrix.
mat3 mat3_ident() __attribute__ ((const));
//Multiply a by b
mat3 mat3_mult(mat3 a, mat3 b) __attribute__ ((const));
//Multiply a by the column vector b.
vec3 mat3_multvec(mat3 a, vec3 b) __attribute__ ((const));
//Rotate a about <ux, uy, uz> by some angle, provided by s and c.
//s and c should be the sine and cosine of the angle, respectively.
mat3 mat3_rot(mat3 a, float ux, float uy, float uz, float s, float c) __attribute__ ((const));
//Produces the matrix for a rotation about <ux, uy, uz> by some angle, provided by s and c.
//s and c should be the sine and cosine of the angle, respectively.
mat3 mat3_rotmat(float ux, float uy, float uz, float s, float c) __attribute__ ((const));
//Produces a rotation matrix about the three basis vectors by some angle, provided by s and c.
//s and c should be the sine and cosine of the angle, respectively.
mat3 mat3_rotmatx(float s, float c) __attribute__ ((const));
mat3 mat3_rotmaty(float s, float c) __attribute__ ((const));
mat3 mat3_rotmatz(float s, float c) __attribute__ ((const));
//Scale a by x, y, z.
mat3 mat3_scale(mat3 a, float x, float y, float z) __attribute__ ((const));
//Produce a matrix that will scale by x, y, z.
mat3 mat3_scalemat(float x, float y, float z) __attribute__ ((const));
//Produce the transpose matrix of a.
mat3 mat3_transp(mat3 a) __attribute__ ((const));
//Produce a rotation matrix that will look from p to q, with u up.
mat3 mat3_lookat(vec3 p, vec3 q, vec3 u) __attribute__ ((const));
//Copy a into a buffer representing a true 3x3 row-major matrix.
//The buffer should be large enough to store 9 floats.
void mat3_to_array(mat3 a, float *buf);
//Copy a into a buffer representing a true 3x3 column-major matrix.
//The buffer should be large enough to store 9 floats.
void mat3_to_array_cm(mat3 a, float *buf);
//Takes a mat3 and a vec3, and copies them into a buffer representing a true, row-major 4x4 matrix.
//a becomes the rotation portion, and b becomes the translation.
void mat3_vec3_to_array(mat3 a, vec3 b, float *buf) __attribute__ ((const));
//Prints a mat3, row by row.
void mat3_print(mat3 a);
//Prints an array of 9 floats, row by row.
void mat3_buf_print(float *a);


//Multiply a by b and return the result as a new affine matrix.
amat4 amat4_mult(amat4 a, amat4 b) __attribute__ ((const));
//Multiply a by b and return the result as a new affine matrix.
//This version has branches to try to reduce multiplications. (Faster on ARM?)
amat4 amat4_mult_b(amat4 a, amat4 b) __attribute__ ((const));
//Multiply a by b as a column vector and return a new vector.
//b is implied to be a 4-vec with the form <x, y, z, 1>
vec3 amat4_multpoint(amat4 a, vec3 b) __attribute__ ((const));
//Multiply a by b as a column vector and return a new vector.
//b is implied to be a 4-vec with the form <x, y, z, 0>
//It's probably faster to do mat3_multvec(a.a, b), since that copies less.
vec3 amat4_multvec(amat4 a, vec3 b) __attribute__ ((const));
//Produce a new matrix that is the result of rotating a about the axis <ux, uy, uz> by angle, in radians.
//<ux, uy, uz> should be normalized beforehand.
//s and c should be the sine and cosine of the angle, respectively.
amat4 amat4_rot(amat4 a, float ux, float uy, float uz, float s, float c) __attribute__ ((const));
//Copy a into a buffer representing a true 4x4 row-major matrix.
//The buffer should be large enough to store 16 floats.
//The last row will be <0, 0, 0, 1>.
void amat4_to_array(amat4 a, float *buf);
//Produce a rotation matrix which represents a rotation about <ux, uy, uz> by some angle.
//s and c should be the sine and cosine of the angle, respectively.
amat4 amat4_rotmat(float ux, float uy, float uz, float s, float c) __attribute__ ((const));
//Produce a rotation matrix which represents a rotation about <ux, uy, uz> by angle, in radians.
//This version tries to reduce multiplications at the cost of more interdependant local variables.
amat4 amat4_rotmat_lomult(float ux, float uy, float uz, float s, float c) __attribute__ ((const));
//Produce a lookat matrix that points from point p to point q, with u as "up".
amat4 amat4_lookat(vec3 p, vec3 q, vec3 u) __attribute__ ((const));
//Produce the inverse matrix to a, provided that a represents a rotation and a translation.
//This is not a true inverse, which in practice would be quite slow, and should be avoided.
//Thus, this does not work if a represents a scale or skew.
amat4 amat4_inverse(amat4 a) __attribute__ ((const));
//Do a true 4x4 matrix multiply between a and b, and put the result into out.
//a, b and out should represent 4x4 row-major matrices, as arrays of floats.
//a, b and out must all be different 16-float buffers, or the behaviour of this function is undefined.
void amat4_buf_mult(float * restrict a, float * restrict b, float * restrict out);
//Do a true 4x4 matrix multiply between a and b, and put the result into out.
//b and out should represent 4x4 row-major matrices, as arrays of floats.
//b and out must all be different 16-float buffers, or the behaviour of this function is undefined.
void amat4_amat_buf_mult(amat4 a, float * restrict b, float * restrict out);
//Do a true 4x4 matrix multiply between a and b, and put the result into out.
//a and out should represent 4x4 row-major matrices, as arrays of floats.
//a and out must all be different 16-float buffers, or the behaviour of this function is undefined.
void amat4_buf_amat_mult(float * restrict a, amat4 b, float * restrict out);
//Multiply the 4x4 row-major matrix a by the column vector b and return a new vector.
//b is implied to be a 4-vec with the form <x, y, z, 1>
//a should be an array of 16 floats.
vec3 amat4_buf_multpoint(float *a, float *b, float *out);
//Prints an amat4, row by row.
void amat4_print(amat4 a);
//Prints an array of 16 floats, row by row.
void amat4_buf_print(float *a);


#endif