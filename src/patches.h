#pragma once
// @note: Although we can get rid of the frequency argument for patches
// by passing instead x = t * f, this causes problems when we have patches
// that we want to animate with a speed that is independent of the
// frequency set in the patch sequence.
//
// Note that although it might not be productive to move the frequency
// argument from the patch instance to the patch sequencer, it might sound
// like a good idea to move it into the patch itself. The problem, then, is
// that each patch needs to explicitely take a frequency argument. Instead,
// we can use patches such as fmul or fset to only set the frequency for
// those patches that need to be different.

#include <cstdint>
#include <cstddef>

// @note: Random thought. We could make Patch work on 3d data instead,
// and do the perspective divide in the sequencer.

struct Patch
{
    uint8_t type;
    void* userdata;
    void (*call)(void* userdata, double t, double f, double* x, double* y, double* z);
};

typedef double (*DoubleFunc1)(double);
typedef double (*DoubleFunc2)(double, double);

Patch patch_clone(Patch p);
void patch_data_free(Patch p);

/* --- Core patches --- */

struct AudioData
{
    uint8_t* buffer;
    uint32_t num_samples;
    uint32_t sample_size;
};

Patch pt_lift1(DoubleFunc1 func_x, DoubleFunc1 func_y, DoubleFunc1 func_z, Patch patch);
Patch pt_lift2(DoubleFunc2 func_x, DoubleFunc2 func_y, DoubleFunc2 func_z, Patch patch_a, Patch patch_b);

Patch pt_time(void);
Patch pt_frequency(void);

Patch pt_point(double x, double y, double z);

Patch pt_random(uint32_t seed);
Patch pt_seq(size_t n_patches, Patch* patches);
Patch pt_flip(Patch x);
Patch pt_audio(const AudioData* audio_Data);

Patch pt_fmul(Patch patch, double factor);
Patch pt_fset(Patch patch, double frequency);
// @todo: offset?

/* --- Derivative patches --- */
Patch pt_lift1(DoubleFunc1 func, Patch patch);
Patch pt_lift2(DoubleFunc2 func_x, DoubleFunc2 func_y, DoubleFunc2 func_z, Patch patch, double value);
Patch pt_lift2(DoubleFunc2 func_x, DoubleFunc2 func_y, DoubleFunc2 func_z, double value, Patch patch);
Patch pt_lift2(DoubleFunc2 func, Patch patch, double value);
Patch pt_lift2(DoubleFunc2 func, double value, Patch patch);
Patch pt_lift2(DoubleFunc2 func, Patch a, Patch b);

Patch operator*(Patch a,  Patch b);
Patch operator*(Patch a, double b);
Patch operator*(double a, Patch b);

Patch operator/(Patch a,  Patch b);
Patch operator/(Patch a, double b);
Patch operator/(double a, Patch b);

Patch operator+(Patch a,  Patch b);
Patch operator+(Patch a, double b);
Patch operator+(double a, Patch b);

Patch operator-(Patch a,  Patch b);
Patch operator-(Patch a, double b);
Patch operator-(double a, Patch b);
                                  
Patch operator<(Patch a,  Patch b);
Patch operator<(Patch a, double b);
Patch operator<(double a, Patch b);

Patch operator>(Patch a,  Patch b);
Patch operator>(Patch a, double b);
Patch operator>(double a, Patch b);

Patch pt_point(double a);

Patch pt_line(double xa, double ya, double za, double xb, double yb, double zb);
Patch pt_line(double xa, double ya, double xb, double yb);
Patch pt_line(Patch a, Patch b);
Patch pt_line_strip(size_t n_points, const double* points);
Patch pt_line_loop(size_t n_points, const double* points);

Patch pt_fmod(Patch t, double value);
Patch pt_fmod(Patch a, Patch b);
Patch pt_fabs(Patch t);
Patch pt_pow(Patch x, double power);
Patch pt_circle(void);
Patch pt_square_wave(void);

enum class TextAlignment {LEFT, RIGHT, CENTER};
Patch pt_text(const char* text, TextAlignment alignment);

