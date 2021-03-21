#include <stdio.h>
#include <math.h> //sqrt, sin, cos
#include <string.h> //memcpy
#include "glla.h"

float vec2_sum(vec2 a)
{
	return a.x + a.y;
}

float vec2_dot(vec2 u, vec2 v)
{
	return vec2_sum(u * v);
}

vec3 vec3_normalize(vec3 a)
{
	return a/vec3_mag(a);
}

vec3 vec3_cross(vec3 u, vec3 v)
{
	return (vec3){u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x};
}

vec3 vec3_lerp(vec3 a, vec3 b, float alpha)
{
	return ((1 - alpha) * a) + (alpha * b);
}

float vec3_dot(vec3 u, vec3 v)
{
	return vec3_sum(u * v);
}

float vec3_mag(vec3 a)
{
	return sqrt(vec3_sum(a * a));
}

float vec3_dist(vec3 a, vec3 b)
{
	return vec3_mag(a - b);
}

float vec3_sum(vec3 a)
{
	return a.x + a.y + a.z;
}

void vec3_unpack(float out[3], vec3 in)
{
	out[0] = in.x;
	out[1] = in.y;
	out[2] = in.z;
}

void vec3_print(vec3 a)
{
	printf("{%f, %f, %f}", a.x, a.y, a.z);
}

void vec3_println(vec3 a)
{
	printf("{%f, %f, %f}\n", a.x, a.y, a.z);
}

void vec3_printf(char *fmt, vec3 a)
{
	char format[strlen(fmt) * 3 + 7]; //2 braces + 2 commas + 2 spaces + 1 newline = 7
	sprintf(format, "{%s, %s, %s}", fmt, fmt, fmt);
	printf(format, a.x, a.y, a.z);
}

dvec3 dvec3_normalize(dvec3 a)
{
	return a/dvec3_mag(a);
}

dvec3 dvec3_cross(dvec3 u, dvec3 v)
{
	return (dvec3){u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x};
}

dvec3 dvec3_lerp(dvec3 a, dvec3 b, double alpha)
{
	return ((1 - alpha) * a) + (alpha * b);
}

float dvec3_dot(dvec3 u, dvec3 v)
{
	return dvec3_sum(u * v);
}

float dvec3_mag(dvec3 a)
{
	return sqrt(dvec3_sum(a * a));
}

float dvec3_dist(dvec3 a, dvec3 b)
{
	return dvec3_mag(a - b);
}

float dvec3_sum(dvec3 a)
{
	return a.x + a.y + a.z;
}

double dvec3_sumd(dvec3 a)
{
	return a.x + a.y + a.z;
}

void dvec3_unpack(double out[3], dvec3 in)
{
	out[0] = in.x;
	out[1] = in.y;
	out[2] = in.z;
}

void dvec3_print(dvec3 a)
{
	printf("{%lf, %lf, %lf}", a.x, a.y, a.z);
}

void dvec3_printf(char *fmt, dvec3 a)
{
	char format[strlen(fmt) * 3 + 7]; //2 braces + 2 commas + 2 spaces + 1 newline = 7
	sprintf(format, "{%s, %s, %s}", fmt, fmt, fmt);
	printf(format, a.x, a.y, a.z);
}

int64_t ivec3_sumd(ivec3 a)
{
	return a.x + a.y + a.z;
}

float ivec3_magf(ivec3 a)
{
	return sqrt(ivec3_sumd(a * a));
}

float ivec3_distf(ivec3 a, ivec3 b)
{
	return ivec3_magf(a - b);
}

float qvec3_magf(qvec3 a)
{
	return sqrt(qvec3_sum(a * a));
}

double qvec3_magd(qvec3 a)
{
	dvec3 tmp = {a.x, a.y, a.z};
	return sqrt(dvec3_sumd(tmp * tmp));
}

float qvec3_distf(qvec3 a, qvec3 b)
{
	return qvec3_magf(a - b);
}

double qvec3_distd(qvec3 a, qvec3 b)
{
	return qvec3_magd(a - b);
}

int64_t qvec3_sum(qvec3 a)
{
	return a.x + a.y + a.z;
}

double qvec3_sumd(qvec3 a)
{
	return (double)a.x + (double)a.y + (double)a.z;
}

void qvec3_print(qvec3 b)
{
	printf("[%lli, %lli, %lli]", b.x, b.y, b.z);
}

void qvec3_println(qvec3 b)
{
	printf("[%lli, %lli, %lli]\n", b.x, b.y, b.z);
}

void qvec3_printf(char *fmt, qvec3 b)
{
	char format[strlen(fmt) * 3 + 7]; //2 braces + 2 commas + 2 spaces + 1 newline = 7
	sprintf(format, "{%s, %s, %s}", fmt, fmt, fmt);
	printf(format, b.x, b.y, b.z);
}

void svec3_unpack(int16_t out[3], svec3 in)
{
	out[0] = in.x;
	out[1] = in.y;
	out[2] = in.z;
}

float vec4_sum(vec4 a)
{
	return a.x + a.y + a.z + a.w;
}

float vec4_dot(vec4 u, vec4 v)
{
	return vec4_sum(u * v);
}

mat3 mat3_mult(mat3 a, mat3 b)
{
	return (mat3)
	{
		{
			{ //First row
				vec3_dot(a.rows[0], (vec3){b.rows[0].x, b.rows[1].x, b.rows[2].x}),
				vec3_dot(a.rows[0], (vec3){b.rows[0].y, b.rows[1].y, b.rows[2].y}),
				vec3_dot(a.rows[0], (vec3){b.rows[0].z, b.rows[1].z, b.rows[2].z})
			},
			{ //Second row
				vec3_dot(a.rows[1], (vec3){b.rows[0].x, b.rows[1].x, b.rows[2].x}),
				vec3_dot(a.rows[1], (vec3){b.rows[0].y, b.rows[1].y, b.rows[2].y}),
				vec3_dot(a.rows[1], (vec3){b.rows[0].z, b.rows[1].z, b.rows[2].z})
			},
			{//Third row
				vec3_dot(a.rows[2], (vec3){b.rows[0].x, b.rows[1].x, b.rows[2].x}),
				vec3_dot(a.rows[2], (vec3){b.rows[0].y, b.rows[1].y, b.rows[2].y}),
				vec3_dot(a.rows[2], (vec3){b.rows[0].z, b.rows[1].z, b.rows[2].z})
			}
		}
	};
}

vec3 mat3_multvec(mat3 a, vec3 b)
{
	return (vec3){
		vec3_dot(a.rows[0], b), 
		vec3_dot(a.rows[1], b),
		vec3_dot(a.rows[2], b)
	};
}

mat3 mat3_rot(mat3 a, float ux, float uy, float uz, float s, float c)
{
	return mat3_mult(a, mat3_rotmat(ux, uy, uz, s, c));
}

mat3 mat3_rotmat(float ux, float uy, float uz, float s, float c)
{
	float c1 = 1-c;
	return (mat3){
		{
			{c + ux*ux*c1, ux*uy*c1 - uz*s, ux*uz*c1 + uy*s},
			{uy*ux*c1 + uz*s, c + uy*uy*c1, uy*uz*c1 - ux*s},
			{ux*uz*c1 - uy*s, uy*uz*c1 + ux*s, c + uz*uz*c1}
		}
	};
}

mat3 mat3_rotmatx(float s, float c)
{
	return (mat3){
		{
			{1, 0, 0},
			{0, c, -s},
			{0, s, c}
		}
	};
}
mat3 mat3_rotmaty(float s, float c)
{
	return (mat3){
		{
			{c, 0, s},
			{0, 1, 0},
			{-s, 0, c}
		}
	};
}
mat3 mat3_rotmatz(float s, float c)
{
	return (mat3){
		{
			{c, -s, 0},
			{s, c, 0},
			{0, 0, 1}
		}
	};
}

mat3 mat3_scale(mat3 a, float x, float y, float z)
{
	return (mat3){
		{
			x * a.rows[0],
			y * a.rows[1],
			z * a.rows[2]
		}
	};
};

mat3 mat3_scalemat(float x, float y, float z)
{
	return (mat3){
		{
			{x, 0, 0},
			{0, y, 0},
			{0, 0, z}
		}
	};
}

mat3 mat3_transp(mat3 a)
{
	return (mat3){
		{
			{a.rows[0].x, a.rows[1].x, a.rows[2].x},
			{a.rows[0].y, a.rows[1].y, a.rows[2].y},
			{a.rows[0].z, a.rows[1].z, a.rows[2].z}
		}
	};
}

mat3 mat3_lookat(vec3 p, vec3 q, vec3 u)
{
	vec3 z = vec3_normalize(p - q);
	vec3 x = vec3_normalize(vec3_cross(u, z));
	vec3 y = vec3_cross(z, x);
	return (mat3){{
		{x.x, y.x, z.x},
		{x.y, y.y, z.y},
		{x.z, y.z, z.z}
	}};
}

void mat3_to_array(mat3 a, float *buf)
{
	float tmp[] = {
		a.rows[0].x, a.rows[0].y, a.rows[0].z,
		a.rows[1].x, a.rows[1].y, a.rows[1].z,
		a.rows[2].x, a.rows[2].y, a.rows[2].z
	};
	memcpy(buf, tmp, sizeof(tmp));
}

void mat3_to_array_cm(mat3 a, float *buf)
{
	float tmp[] = {
		a.rows[0].x, a.rows[1].x, a.rows[2].x,
		a.rows[0].y, a.rows[1].y, a.rows[2].y,
		a.rows[0].z, a.rows[1].z, a.rows[2].z
	};
	memcpy(buf, tmp, sizeof(tmp));
}

void mat3_vec3_to_array(mat3 a, vec3 b, float *buf)
{
	float tmp[] = {
		a.rows[0].x, a.rows[0].y, a.rows[0].z, b.x,
		a.rows[1].x, a.rows[1].y, a.rows[1].z, b.y,
		a.rows[2].x, a.rows[2].y, a.rows[2].z, b.z,
		0,                     0,           0,   1
	};
	memcpy(buf, tmp, sizeof(tmp));
}

void mat3_print(mat3 a)
{
#define mrow(v) v[0], v[1], v[2]
	vec3 *r = a.rows;
	printf("%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n",
		mrow(r[0]),
		mrow(r[1]),
		mrow(r[2]));
#undef mrow
}

void mat3_buf_print(float *a)
{
#define mrow(a, r) a[3*r], a[3*r+1], a[3*r+2]
	printf("%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n",
		mrow(a, 0),
		mrow(a, 1),
		mrow(a, 2));
#undef mrow
}

amat4 amat4_mult(amat4 a, amat4 b)
{
	return (amat4){
		mat3_mult(a.a, b.a),
		amat4_multpoint(a, b.t)
	};
}

vec3 amat4_multpoint(amat4 a, vec3 b)
{
	return mat3_multvec(a.a, b) + a.t;
}

vec3 amat4_multvec(amat4 a, vec3 b)
{
	return mat3_multvec(a.a, b);
}

amat4 amat4_rot(amat4 a, float ux, float uy, float uz, float s, float c)
{
	return (amat4){
		mat3_rot(a.a, ux, uy, uz, s, c),
		a.t
	};
}

void amat4_to_array(amat4 a, float *buf)
{
	float tmp[] = {
		a.a.rows[0].x, a.a.rows[0].y, a.a.rows[0].z, a.t.x,
		a.a.rows[1].x, a.a.rows[1].y, a.a.rows[1].z, a.t.y,
		a.a.rows[2].x, a.a.rows[2].y, a.a.rows[2].z, a.t.z,
		            0,             0,             0,     1
	};
	memcpy(buf, tmp, sizeof(tmp));
};

amat4 amat4_rotmat(float ux, float uy, float uz, float s, float c)
{
	return (amat4){mat3_rotmat(ux, uy, uz, s, c), {0, 0, 0}};
}

amat4 amat4_rotmat_lomult(float ux, float uy, float uz, float s, float c)
{
	float c1 = 1-c;
	//This strategy may exhibit poorer instruction-level parallelism.
	float uxc1 = ux*c1;    //Costs one multiply, saves five
	float uyc1 = uy*c1;    //Costs one multiply, saves three
	float uyzc1 = uz*uyc1; //Costs one multiply, saves two
	float uxs = ux * s;    //Costs one multiply, saves two
	float uys = uy * s;    //Costs one multiply, saves two
	float uzs = uz * s;    //Costs one multiply, saves two
	return (amat4){
		{
			{
				{c + ux*uxc1, uy*uxc1 - uzs, uz*uxc1 + uys},
				{uy*uxc1 + uzs, c + uy*uyc1, uyzc1 - uxs},
				{uz*uxc1 - uys, uyzc1 + uxs, c + uz*uz*c1}
			}
		},
		{0, 0, 0}
	};
}

//"Fake" inverse, only handles translation and rotation (but cheap!)
amat4 amat4_inverse(amat4 a)
{
	return (amat4){
		mat3_transp(a.a),
		{
			vec3_dot((vec3){a.a.rows[0].x, a.a.rows[1].x, a.a.rows[2].x}, -a.t),
			vec3_dot((vec3){a.a.rows[0].y, a.a.rows[1].y, a.a.rows[2].y}, -a.t),
			vec3_dot((vec3){a.a.rows[0].z, a.a.rows[1].z, a.a.rows[2].z}, -a.t)
		}
	};
}

void amat4_buf_mult(float * restrict a, float * restrict b, float * restrict out)
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			out[i * 4 + j] =
			a[i * 4 + 0] * b[ 0 + j]+
			a[i * 4 + 1] * b[ 4 + j]+
			a[i * 4 + 2] * b[ 8 + j]+
			a[i * 4 + 3] * b[12 + j];
		}
	}
}

void amat4_mat_buf_mult(amat4 a, float * restrict b, float * restrict out)
{
	float tmp[16];
	amat4_to_array(a, tmp);
	amat4_buf_mult(tmp, b, out);
}

void amat4_buf_mat_mult(float * restrict a, amat4 b, float * restrict out)
{
	float tmp[16];
	amat4_to_array(b, tmp);
	amat4_buf_mult(a, tmp, out);
}

vec3 amat4_buf_multpoint(float *a, float *b, float *out)
{
	float tmp[4];
	for (int i = 0; i < 4; i++)
		tmp[i] = a[4*i+0]*b[0] + a[4*i+1]*b[1] + a[4*i+2]*b[2] + a[4*i+3]*b[3];
	if (out)
		memcpy(out, tmp, 4*sizeof(float));
	return (vec3){tmp[0], tmp[1], tmp[2]};
}

void amat4_print(amat4 a)
{
#define mrow(v, f) v[0], v[1], v[2], f
	vec3 *r = a.a.rows;
	printf("%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n",
		mrow(r[0], a.t[0]),
		mrow(r[1], a.t[1]),
		mrow(r[2], a.t[2]),
		0.0, 0.0, 0.0, 1.0);
#undef mrow
}

void amat4_buf_print(float *a)
{
#define mrow(a, r) a[4*r], a[4*r+1], a[4*r+2], a[4*r+3]
	printf("%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n",
		mrow(a, 0),
		mrow(a, 1),
		mrow(a, 2),
		mrow(a, 3));
#undef mrow
}
