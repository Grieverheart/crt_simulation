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

#include <limits>
#include <new>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>

#include "patches.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// @note: Random thought. We could make Patch work on 3d data instead,
// and do the perspective divide in the sequencer.

namespace
{
    double osc_phaser(double t, double frequency)
    {
        return fmod(t, 1.0 / frequency) * frequency;
    }

    double osc_sawtooth_wave(double t)
    {
        return 2.0 * (t - uint64_t(t)) - 1.0;
    }

    double osc_triangle_wave(double t)
    {
        return 2.0 * fabs(osc_sawtooth_wave(2.0 * t)) - 1.0;
    }

    double osc_square_wave(double t)
    {
        return 2.0 * (osc_sawtooth_wave(t) > 0.0) - 1.0;
    }

    enum PatchType: uint8_t
    {
        PT_LIFT1,
        PT_LIFT2,
        PT_FREQUENCY,
        PT_TIME,
        PT_POINT,
        PT_POINT_OBS,
        PT_RANDOM,
        PT_SEQ,
        PT_FLIP,
        PT_AUDIO,
        PT_FMUL,
        PT_FSET,
    };

    struct ScalarPatchData
    {
        Patch patch;
        double value;
    };

    struct Lift1Data
    {
        DoubleFunc1 func_x, func_y, func_z;
        Patch patch;
    };

    struct Lift2PartialData
    {
        DoubleFunc2 func_x, func_y, func_z;
        double value;
        Patch patch;
    };

    struct Lift2Data
    {
        DoubleFunc2 func_x, func_y, func_z;
        Patch patch_a, patch_b;
    };

    struct PatchArray
    {
        size_t n_patches;
        Patch* patches;
    };

    void lift1_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        const Lift1Data* data = (const Lift1Data*) userdata;
        data->patch.call(data->patch.userdata, t, frequency, x, y, z);
        if(data->func_x)
            *x = data->func_x(*x);
        if(data->func_y)
            *y = data->func_y(*y);
        if(data->func_z)
            *z = data->func_z(*z);
    }

    void lift2l_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        const Lift2PartialData* data = (const Lift2PartialData*) userdata;
        data->patch.call(data->patch.userdata, t, frequency, x, y, z);
        if(data->func_x)
            *x = data->func_x(*x, data->value);
        if(data->func_y)
            *y = data->func_y(*y, data->value);
        if(data->func_z)
            *z = data->func_z(*z, data->value);
    }

    void lift2r_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        const Lift2PartialData* data = (const Lift2PartialData*) userdata;
        data->patch.call(data->patch.userdata, t, frequency, x, y, z);
        *x = data->func_x(data->value, *x);
        *y = data->func_y(data->value, *y);
        *z = data->func_z(data->value, *z);
    }

    void lift2_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        const Lift2Data* data = (const Lift2Data*) userdata;
        double xa, ya, xb, yb, za, zb;
        data->patch_a.call(data->patch_a.userdata, t, frequency, &xa, &ya, &za);
        data->patch_b.call(data->patch_b.userdata, t, frequency, &xb, &yb, &zb);
        if(data->func_x)
            *x = data->func_x(xa, xb);
        if(data->func_y)
            *y = data->func_y(ya, yb);
        if(data->func_z)
            *z = data->func_z(za, zb);
    }

    void point_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        const double* coords = (double*)userdata;
        *x = coords[0];
        *y = coords[1];
        *z = coords[2];
    }

    void time_call(void* userdata, double t, double f, double* x, double* y, double* z)
    {
        *x = t;
        *y = t;
        *z = t;
    }

    void seq_call(void* userdata, double t, double f, double* x, double* y, double* z)
    {
        const PatchArray* seq_data = (const PatchArray*)userdata;

        double phi = osc_phaser(t, f);
        size_t pid = seq_data->n_patches * phi;

        if(pid < seq_data->n_patches)
        {
            const Patch* p = &seq_data->patches[pid];
            p->call(p->userdata, t, f * seq_data->n_patches, x, y, z);
        }
    }

    void flip_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        Patch p = *(const Patch*) userdata;
        p.call(p.userdata, t, frequency, y, x, z);
    }

    /* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
    uint32_t xorshift32(uint32_t* rng)
    {
        uint32_t x = *rng;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        *rng = x;
        return x;
    }

    double random(uint32_t* rng)
    {
        return (double)xorshift32(rng) / std::numeric_limits<uint32_t>::max();
    }

    void random_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        uint32_t* rng = (uint32_t*) userdata;
        *x = random(rng);
        *y = random(rng);
        *z = random(rng);
    }

    void wav_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        *z = 1.0;
        const AudioData* audio = (const AudioData*)userdata;

        double phi = osc_phaser(t, frequency);

        uint32_t sid      = uint32_t(phi * audio->num_samples);
        uint32_t sid_next = (sid != audio->num_samples - 1)? sid + 1: 0;

        double f = phi * audio->num_samples - uint32_t(phi * audio->num_samples);

        switch(audio->sample_size)
        {
            case 1:
            {
                const int8_t* buffer = (int8_t*)audio->buffer;
                *x = (1.0 - f) * (buffer[2 * sid + 0] / 128.0) + f * (buffer[2 * sid_next + 0] / 128.0);
                *y = (1.0 - f) * (buffer[2 * sid + 1] / 128.0) + f * (buffer[2 * sid_next + 1] / 128.0);
            } break;

            case 2:
            {
                const int16_t* buffer = (int16_t*)audio->buffer;
                *x = (1.0 - f) * (buffer[2 * sid + 0] / 32768.0) + f * (buffer[2 * sid_next + 0] / 32768.0);
                *y = (1.0 - f) * (buffer[2 * sid + 1] / 32768.0) + f * (buffer[2 * sid_next + 1] / 32768.0);
            } break;

            case 4:
            {
                const int32_t* buffer = (int32_t*)audio->buffer;
                *x = (1.0 - f) * (buffer[2 * sid + 0] / 2147483648.0) + f * (buffer[2 * sid_next + 0] / 2147483648.0);
                *y = (1.0 - f) * (buffer[2 * sid + 1] / 2147483648.0) + f * (buffer[2 * sid_next + 1] / 2147483648.0);
            } break;

            default: assert(false); break;
        }
    }

    void fmul_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        const ScalarPatchData* data = (const ScalarPatchData*)userdata;
        data->patch.call(data->patch.userdata, t, data->value * frequency, x, y, z);
    }

    void fset_call(void* userdata, double t, double frequency, double* x, double* y, double* z)
    {
        const ScalarPatchData* data = (const ScalarPatchData*)userdata;
        data->patch.call(data->patch.userdata, t, data->value, x, y, z);
    }

} // namespace


void patch_data_free(Patch p)
{
    switch(p.type)
    {
        case PT_LIFT1:
        {
            Lift1Data* data = (Lift1Data*) p.userdata;
            patch_data_free(data->patch);
            free(data);
        } break;

        case PT_LIFT2:
        {
            Lift2Data* data = (Lift2Data*) p.userdata;
            patch_data_free(data->patch_a);
            patch_data_free(data->patch_b);
            free(p.userdata);
        } break;

        case PT_FREQUENCY:
        case PT_TIME:
        {
        } break;

        case PT_POINT:
        {
            free(p.userdata);
        } break;

        case PT_RANDOM:
        {
            free(p.userdata);
        } break;

        case PT_SEQ:
        {
            PatchArray* data = (PatchArray*) p.userdata;
            for(size_t i = 0; i < data->n_patches; ++i)
                patch_data_free(data->patches[i]);
            free(data->patches);
            free(data);
        } break;

        case PT_FLIP:
        {
            Patch* data = (Patch*) p.userdata;
            patch_data_free(*data);
            free(data);
        } break;

        case PT_AUDIO: break;

        case PT_FSET:
        case PT_FMUL:
        {
            ScalarPatchData* data = (ScalarPatchData*) p.userdata;
            patch_data_free(data->patch);
            free(data);
        } break;

        default: assert(false); break;
    }
}

Patch patch_clone(Patch p)
{
    Patch r;
    r.type = p.type;
    r.call = p.call;

    switch(p.type)
    {
        case PT_LIFT1:
        {
            r.userdata = malloc(sizeof(Lift1Data));
            Lift1Data* data = (Lift1Data*) r.userdata;

            const Lift1Data* other_data = (const Lift1Data*) p.userdata;
            data->func_x = other_data->func_x;
            data->func_y = other_data->func_y;
            data->func_z = other_data->func_z;
            data->patch = patch_clone(other_data->patch);
        } break;

        case PT_LIFT2:
        {
            r.userdata = malloc(sizeof(Lift2Data));
            Lift2Data* data = (Lift2Data*) r.userdata;

            const Lift2Data* other_data = (const Lift2Data*) p.userdata;
            data->func_x = other_data->func_x;
            data->func_y = other_data->func_y;
            data->func_z = other_data->func_z;
            data->patch_a = patch_clone(other_data->patch_a);
            data->patch_b = patch_clone(other_data->patch_b);
        } break;

        case PT_FREQUENCY:
        case PT_TIME:
        {
            r.userdata = nullptr;
        } break;

        case PT_POINT:
        {
            r.userdata = malloc(3 * sizeof(double));
            memcpy(r.userdata, p.userdata, 3 * sizeof(double));
        } break;

        case PT_RANDOM:
        {
            r.userdata = malloc(sizeof(uint32_t));
            memcpy(r.userdata, p.userdata, sizeof(uint32_t));
        } break;

        case PT_SEQ:
        {
            r.userdata = malloc(sizeof(PatchArray));
            PatchArray* data = (PatchArray*) r.userdata;
            const PatchArray* other_data = (const PatchArray*) p.userdata;

            data->n_patches = other_data->n_patches;
            data->patches = (Patch*) malloc(other_data->n_patches * sizeof(Patch));

            for(size_t i = 0; i < other_data->n_patches; ++i)
                data->patches[i] = patch_clone(other_data->patches[i]);
        } break;

        case PT_FLIP:
        {
            r.userdata = malloc(sizeof(Patch));
            Patch* data = (Patch*) r.userdata;
            *data = patch_clone(*(const Patch*)p.userdata);
        } break;

        case PT_AUDIO:
        {
            r.userdata = p.userdata;
        } break;

        case PT_FSET:
        case PT_FMUL:
        {
            r.userdata = malloc(sizeof(ScalarPatchData));
            ScalarPatchData* data = (ScalarPatchData*) r.userdata;
            const ScalarPatchData* other_data = (const ScalarPatchData*) p.userdata;

            data->value = other_data->value;
            data->patch = patch_clone(other_data->patch);
        } break;

        default: assert(false); break;
    }

    return r;
}


Patch pt_lift1(DoubleFunc1 func_x, DoubleFunc1 func_y, DoubleFunc1 func_z, Patch patch)
{
    Patch p;
    p.userdata = malloc(sizeof(Lift1Data));
    p.type = PT_LIFT1;
    new (p.userdata) Lift1Data{func_x, func_y, func_z, patch};
    p.call = lift1_call;
    return p;
}

Patch pt_lift1(DoubleFunc1 func, Patch patch)
{
    return pt_lift1(func, func, func, patch);
}

Patch pt_lift2(DoubleFunc2 func_x, DoubleFunc2 func_y, DoubleFunc2 func_z, Patch patch_a, Patch patch_b)
{
    Patch p;
    p.type = PT_LIFT2;
    p.call = lift2_call;
    p.userdata = malloc(sizeof(Lift2Data));
    new (p.userdata) Lift2Data{func_x, func_y, func_z, patch_a, patch_b};
    return p;
}

Patch pt_lift2(DoubleFunc2 func_x, DoubleFunc2 func_y, DoubleFunc2 func_z, Patch patch, double value)
{
    return pt_lift2(func_x, func_y, func_z, patch, pt_point(value));
}

Patch pt_lift2(DoubleFunc2 func_x, DoubleFunc2 func_y, DoubleFunc2 func_z, double value, Patch patch)
{
    return pt_lift2(func_x, func_y, func_z, pt_point(value), patch);
}

Patch pt_lift2(DoubleFunc2 func, Patch patch, double value)
{
    return pt_lift2(func, func, func, patch, pt_point(value));
}

Patch pt_lift2(DoubleFunc2 func, double value, Patch patch)
{
    return pt_lift2(func, func, func, pt_point(value), patch);
}

Patch pt_lift2(DoubleFunc2 func, Patch a, Patch b)
{
    return pt_lift2(func, func, func, a, b);
}

Patch pt_fmod(Patch t, double value)
{
    return pt_lift2(fmod, t, value);
}

Patch pt_fmod(Patch a, Patch b)
{
    return pt_lift2(fmod, a, b);
}

Patch pt_fabs(Patch t)
{
    return pt_lift1(fabs, t);
}

double mul_double(double x, double y)
{
    return x * y;
}

double div_double(double x, double y)
{
    return x / y;
}

double add_double(double x, double y)
{
    return x + y;
}

double sub_double(double x, double y)
{
    return x - y;
}

double lesser_double(double x, double y)
{
    return x < y;
}

double greater_double(double x, double y)
{
    return x > y;
}

#define DEF_BINARY_OPERATOR(op_name, op, Ta, Tb)\
Patch operator op(Ta a, Tb b)\
{\
    return pt_lift2(op_name##_double, a, b);\
}

DEF_BINARY_OPERATOR(mul, *, Patch,  Patch);
DEF_BINARY_OPERATOR(mul, *, Patch, double);
DEF_BINARY_OPERATOR(mul, *, double, Patch);

DEF_BINARY_OPERATOR(div, /, Patch,  Patch);
DEF_BINARY_OPERATOR(div, /, Patch, double);
DEF_BINARY_OPERATOR(div, /, double, Patch);

DEF_BINARY_OPERATOR(add, +, Patch,  Patch);
DEF_BINARY_OPERATOR(add, +, Patch, double);
DEF_BINARY_OPERATOR(add, +, double, Patch);

DEF_BINARY_OPERATOR(sub, -, Patch,  Patch);
DEF_BINARY_OPERATOR(sub, -, Patch, double);
DEF_BINARY_OPERATOR(sub, -, double, Patch);

DEF_BINARY_OPERATOR(lesser, <, Patch,  Patch);
DEF_BINARY_OPERATOR(lesser, <, Patch, double);
DEF_BINARY_OPERATOR(lesser, <, double, Patch);

DEF_BINARY_OPERATOR(greater, >, Patch,  Patch);
DEF_BINARY_OPERATOR(greater, >, Patch, double);
DEF_BINARY_OPERATOR(greater, >, double, Patch);

#undef DEF_BINARY_OPERATOR

Patch pt_point(double x, double y, double z)
{
    double* data = (double*) malloc(3 * sizeof(double));
    data[0] = x;
    data[1] = y;
    data[2] = z;

    Patch p;
    p.userdata = (void*) data;
    p.type = PT_POINT;
    p.call = point_call;
    return p;
}

Patch pt_point(double a)
{
    return pt_point(a, a, a);
}

Patch pt_time(void)
{
    return Patch{PT_TIME, nullptr, time_call};
}

void frequency_call(void* userdata, double t, double f, double* x, double* y, double* z)
{
    *x = f;
    *y = f;
    *z = f;
}

Patch pt_frequency(void)
{
    return Patch{PT_FREQUENCY, nullptr, frequency_call};
}

Patch pt_pow(Patch x, double power)
{
    return pt_lift2(pow, x, power);
}

Patch pt_line(double xa, double ya, double za, double xb, double yb, double zb)
{
    return 0.5 * (
        pt_point(xb + xa, yb + ya, zb + za) -
        pt_point(xb - xa, yb - ya, zb - za) * pt_lift1(cos, M_PI * pt_fmod(pt_time(), 1.0 / pt_frequency()) * pt_frequency())
    );
}

Patch pt_line(double xa, double ya, double xb, double yb)
{
    return 0.5 * (
        pt_point(xb + xa, yb + ya, 2.0) -
        pt_point(xb - xa, yb - ya, 0.0) * pt_lift1(cos, M_PI * pt_fmod(pt_time(), 1.0 / pt_frequency()) * pt_frequency())
    );
}

Patch pt_line(Patch a, Patch b)
{
    return 0.5 * (
         b + a -
        (patch_clone(b) - patch_clone(a)) * pt_lift1(cos, M_PI * pt_fmod(pt_time(), 1.0 / pt_frequency()) * pt_frequency())
    );
}

Patch pt_seq(size_t n_patches, Patch* patches)
{
    PatchArray* seq_data = (PatchArray*) malloc(sizeof(PatchArray));
    seq_data->n_patches = n_patches;
    seq_data->patches = patches;

    Patch p;
    p.userdata = (void*) seq_data;
    p.type = PT_SEQ;
    p.call = seq_call;
    return p;
}

Patch pt_flip(Patch x)
{
    Patch* data = (Patch*) malloc(sizeof(Patch));
    *data = x;

    Patch p;
    p.userdata = (void*) data;
    p.type = PT_FLIP;
    p.call = flip_call;
    return p;
}

Patch pt_line_strip(size_t n_points, const double* points)
{
    Patch* patches = (Patch*) malloc(sizeof(Patch) * (n_points - 1));
    for(size_t i = 0; i < n_points - 1; ++i)
        patches[i] = pt_line(points[2 * i], points[2 * i + 1], points[2 * i + 2], points[2 * i + 3]);
    return pt_seq(n_points - 1, patches);
}

Patch pt_line_loop(size_t n_points, const double* points)
{
    Patch* patches = (Patch*) malloc(sizeof(Patch) * n_points);
    for(size_t i = 0; i < n_points - 1; ++i)
        patches[i] = pt_line(points[2 * i], points[2 * i + 1], points[2 * i + 2], points[2 * i + 3]);
    patches[n_points - 1] = pt_line(points[2 * n_points - 2], points[2 * n_points - 1], points[0], points[1]);
    return pt_seq(n_points, patches);
}

Patch pt_circle(void)
{
    return pt_lift1(sin, cos, nullptr, (2.0 * M_PI) * pt_time() * pt_frequency());
}

Patch pt_char_a(void)
{
    Patch* patches = (Patch*) malloc(5 * sizeof(Patch));
    new (patches) Patch[5]
    {
        pt_line( 1.0, -0.5, -1.0, -0.5),
        pt_line(-1.0, -1.0, -1.0,  0.5),
        pt_line(-1.0,  0.5,  0.0,  1.0),
        pt_line( 0.0,  1.0,  1.0,  0.5),
        pt_line( 1.0,  0.5,  1.0, -1.0)
    };

    return pt_seq(5, patches);
}

Patch pt_char_b(void)
{
    Patch* patches = (Patch*) malloc(8 * sizeof(Patch));
    new (patches) Patch[8]
    {
        pt_line(-1.0, -1.0, -1.0,  1.0),
        pt_line(-1.0,  1.0,  0.5,  1.0),
        pt_line( 0.5,  1.0,  1.0,  0.5),
        pt_line( 1.0,  0.5,  0.5,  0.0),
        pt_line( 0.5,  0.0, -1.0,  0.0),
        pt_line( 0.5,  0.0,  1.0, -0.5),
        pt_line( 1.0, -0.5,  0.5, -1.0),
        pt_line( 0.5, -1.0, -1.0, -1.0)
    };

    return pt_seq(8, patches);
}

Patch pt_char_c(void)
{
    double points[] =
    {
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
         1.0,  1.0,
    };

    return pt_line_strip(4, points);
}

Patch pt_char_d(void)
{
    double points[] =
    {
        -1.0, -1.0,
        -1.0,  1.0,
         0.5,  1.0,
         1.0,  0.0,
         0.5, -1.0,
    };

    return pt_line_loop(5, points);
}

Patch pt_char_e(void)
{
    Patch* patches = (Patch*) malloc(5 * sizeof(Patch));
    new (patches) Patch[5]
    {
        pt_line( 1.0, -1.0, -1.0, -1.0),
        pt_line(-1.0, -1.0, -1.0,  0.0),
        pt_line(-1.0,  0.0,  1.0,  0.0),
        pt_line(-1.0,  0.0, -1.0,  1.0),
        pt_line(-1.0,  1.0,  1.0,  1.0),
    };

    return pt_seq(5, patches);
}

Patch pt_char_f(void)
{
    Patch* patches = (Patch*) malloc(4 * sizeof(Patch));
    new (patches) Patch[4]
    {
        pt_line(-1.0, -1.0, -1.0,  0.0),
        pt_line(-1.0,  0.0,  1.0,  0.0),
        pt_line(-1.0,  0.0, -1.0,  1.0),
        pt_line(-1.0,  1.0,  1.0,  1.0),
    };

    return pt_seq(4, patches);
}

Patch pt_char_g(void)
{
    double points[] =
    {
         0.0, -0.0,
         1.0, -0.0,
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
         1.0,  1.0,
         1.0,  0.5,
    };

    return pt_line_strip(7, points);
}

Patch pt_char_h(void)
{
    Patch* patches = (Patch*) malloc(3 * sizeof(Patch));
    new (patches) Patch[3]
    {
        pt_line(-1.0, -1.0, -1.0,  1.0),
        pt_line(-1.0,  0.0,  1.0,  0.0),
        pt_line( 1.0, -1.0,  1.0,  1.0),
    };

    return pt_seq(3, patches);
}

Patch pt_char_i(void)
{
    Patch* patches = (Patch*) malloc(3 * sizeof(Patch));
    new (patches) Patch[3]
    {
        pt_line(-1.0, -1.0,  1.0, -1.0),
        pt_line( 0.0, -1.0,  0.0,  1.0),
        pt_line(-1.0,  1.0,  1.0,  1.0),
    };

    return pt_seq(3, patches);
}

Patch pt_char_j(void)
{
    double points[] =
    {
        -1.0, -0.5,
         0.0, -1.0,
         1.0, -1.0,
         1.0,  1.0,
    };

    return pt_line_strip(4, points);
}

Patch pt_char_k(void)
{
    Patch* patches = (Patch*) malloc(3 * sizeof(Patch));
    new (patches) Patch[3]
    {
        pt_line( 1.0, -1.0, -1.0,  0.0),
        pt_line(-1.0, -1.0, -1.0,  1.0),
        pt_line(-1.0,  0.0,  1.0,  1.0),
    };

    return pt_seq(3, patches);
}

Patch pt_char_l(void)
{
    double points[] =
    {
        -1.0,  1.0,
        -1.0, -1.0,
         1.0, -1.0,
    };

    return pt_line_strip(3, points);
}

Patch pt_char_m(void)
{
    double points[] =
    {
        -1.0, -1.0,
        -1.0,  1.0,
         0.0,  0.0,
         1.0,  1.0,
         1.0, -1.0,
    };

    return pt_line_strip(5, points);
}

Patch pt_char_n(void)
{
    double points[] =
    {
        -1.0, -1.0,
        -1.0,  1.0,
         1.0, -1.0,
         1.0,  1.0,
    };

    return pt_line_strip(4, points);
}

Patch pt_char_o(void)
{
    double points[] =
    {
         1.0,  1.0,
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
    };

    return pt_line_loop(4, points);
}

Patch pt_char_p(void)
{
    double points[] =
    {
        -1.0, -1.0,
        -1.0,  1.0,
         1.0,  1.0,
         1.0,  0.0,
        -1.0,  0.0,
    };

    return pt_line_strip(5, points);
}

Patch pt_char_q(void)
{
    Patch* patches = (Patch*) malloc(7 * sizeof(Patch));
    new (patches) Patch[7]
    {
        pt_line( 0.75, -0.75,  0.50, -1.00),
        pt_line( 0.50, -1.00, -1.00, -1.00),
        pt_line(-1.00, -1.00, -1.00,  1.00),
        pt_line(-1.00,  1.00,  1.00,  1.00),
        pt_line( 1.00,  1.00,  1.00, -0.50),
        pt_line( 1.00, -0.50,  0.75, -0.75),
        pt_line( 0.50, -0.50,  1.00, -1.00),
    };

    return pt_seq(7, patches);
}

Patch pt_char_r(void)
{
    double points[] =
    {
        -1.0, -1.0,
        -1.0,  1.0,
         1.0,  1.0,
         1.0,  0.0,
        -1.0,  0.0,
         1.0, -1.0,
    };

    return pt_line_strip(6, points);
}

Patch pt_char_s(void)
{
    double points[] =
    {
        -1.0, -1.0,
         1.0, -1.0,
         1.0,  0.0,
        -1.0,  0.0,
        -1.0,  1.0,
         1.0,  1.0,
    };

    return pt_line_strip(6, points);
}

Patch pt_char_t(void)
{
    Patch* patches = (Patch*) malloc(2 * sizeof(Patch));
    new (patches) Patch[2]
    {
        pt_line( 0.0, -1.0,  0.0,  1.0),
        pt_line(-1.0,  1.0,  1.0,  1.0),
    };

    return pt_seq(2, patches);
}

Patch pt_char_u(void)
{
    double points[] =
    {
        -1.0,  1.0,
        -1.0, -1.0,
         1.0, -1.0,
         1.0,  1.0,
    };

    return pt_line_strip(4, points);
}

Patch pt_char_v(void)
{
    double points[] =
    {
        -1.0,  1.0,
         0.0, -1.0,
         1.0,  1.0,
    };

    return pt_line_strip(3, points);
}

Patch pt_char_w(void)
{
    double points[] =
    {
        -1.0,  1.0,
        -1.0, -1.0,
         0.0,  0.0,
         1.0, -1.0,
         1.0,  1.0,
    };

    return pt_line_strip(5, points);
}

Patch pt_char_x(void)
{
    Patch* patches = (Patch*) malloc(4 * sizeof(Patch));
    new (patches) Patch[4]
    {
        pt_line(-1.0, -1.0,  0.0,  0.0),
        pt_line( 0.0,  0.0, -1.0,  1.0),

        pt_line( 1.0, -1.0,  0.0,  0.0),
        pt_line( 0.0,  0.0,  1.0,  1.0),
    };

    return pt_seq(4, patches);
}

Patch pt_char_y(void)
{
    Patch* patches = (Patch*) malloc(3 * sizeof(Patch));
    new (patches) Patch[3]
    {
        pt_line( 0.0, -1.0,  0.0,  0.0),
        pt_line(-1.0,  1.0,  0.0,  0.0),
        pt_line( 0.0,  0.0,  1.0,  1.0),
    };

    return pt_seq(3, patches);
}

Patch pt_char_z(void)
{
    double points[] =
    {
         1.0, -1.0,
        -1.0, -1.0,
         1.0,  1.0,
        -1.0,  1.0,
    };

    return pt_line_strip(4, points);
}

Patch pt_char_0(void)
{
    double points[] =
    {
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
         1.0,  1.0,
    };

    return pt_line_loop(4, points);
}

Patch pt_char_1(void)
{
    return pt_line(0.0, 1.0, 0.0, -1.0);
}

Patch pt_char_2(void)
{
    double points[] =
    {
        -1.0,  1.0,
         1.0,  1.0,
         1.0,  0.0,
        -1.0,  0.0,
        -1.0, -1.0,
         1.0, -1.0,
    };

    return pt_line_strip(6, points);
}

Patch pt_char_3(void)
{
    Patch* patches = (Patch*) malloc(5 * sizeof(Patch));
    new (patches) Patch[5]
    {
        pt_line(-1.0, -1.0,  1.0, -1.0),
        pt_line( 1.0, -1.0,  1.0,  0.0),
        pt_line( 1.0,  0.0, -1.0,  0.0),
        pt_line( 1.0,  0.0,  1.0,  1.0),
        pt_line( 1.0,  1.0, -1.0,  1.0),
    };

    return pt_seq(5, patches);
}

Patch pt_char_4(void)
{
    Patch* patches = (Patch*) malloc(3 * sizeof(Patch));
    new (patches) Patch[3]
    {
        pt_line(-1.0,  1.0, -1.0,  0.0),
        pt_line(-1.0,  0.0,  1.0,  0.0),
        pt_line( 1.0,  1.0,  1.0, -1.0),
    };

    return pt_seq(3, patches);
}

Patch pt_char_5(void)
{
    double points[] =
    {
        -1.0, -1.0,
         1.0, -1.0,
         1.0,  0.0,
        -1.0,  0.0,
        -1.0,  1.0,
         1.0,  1.0,
    };

    return pt_line_strip(6, points);
}

Patch pt_char_6(void)
{
    double points[] =
    {
        -1.0,  1.0,
        -1.0, -1.0,
         1.0, -1.0,
         1.0,  0.0,
        -1.0,  0.0,
    };

    return pt_line_strip(5, points);
}

Patch pt_char_7(void)
{
    double points[] =
    {
        -1.0,  1.0,
         1.0,  1.0,
         1.0, -1.0,
    };

    return pt_line_strip(3, points);
}

Patch pt_char_8(void)
{
    double points[] =
    {
        -1.0,  0.0,
         1.0,  0.0,
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
         1.0,  1.0,
         1.0,  0.0,
    };

    return pt_line_strip(7, points);
}

Patch pt_char_9(void)
{
    double points[] =
    {
         1.0,  0.0,
        -1.0,  0.0,
        -1.0,  1.0,
         1.0,  1.0,
         1.0, -1.0,
    };

    return pt_line_strip(5, points);
}

Patch pt_char_colon(void)
{
    Patch* patches = (Patch*) malloc(2 * sizeof(Patch));
    new (patches) Patch[2]
    {
        pt_point(0.0,  0.5, 1.0),
        pt_point(0.0, -0.5, 1.0),
    };
    return pt_seq(2, patches);
}

Patch pt_char_semicolon(void)
{
    Patch* patches = (Patch*) malloc(2 * sizeof(Patch));
    new (patches) Patch[2]
    {
        pt_point(0.0,  0.5, 1.0),
        pt_line(0.0, -0.5, 0.0, -1.0),
    };
    return pt_seq(2, patches);
}

Patch pt_char_less(void)
{
    double points[] =
    {
         1.0,  1.0,
        -1.0,  0.0,
         1.0, -1.0,
    };

    return pt_line_strip(3, points);
}

Patch pt_char_equal(void)
{
    Patch* patches = (Patch*) malloc(2 * sizeof(Patch));
    new (patches) Patch[2]
    {
        pt_line(-1.0,  0.5,  1.0,  0.5),
        pt_line(-1.0, -0.5,  1.0, -0.5),
    };

    return pt_seq(2, patches);
}

Patch pt_char_greater(void)
{
    double points[] =
    {
        -1.0,  1.0,
         1.0,  0.0,
        -1.0, -1.0,
    };

    return pt_line_strip(3, points);
}

Patch pt_char_period(void)
{
    return pt_point(0.0, -1.0, 1.0);
}

Patch pt_char_slash(void)
{
    return pt_line(-1.0, -1.0, 1.0, 1.0);
}

Patch pt_char_none(void)
{
    Patch* patches = (Patch*) malloc(6 * sizeof(Patch));
    new (patches) Patch[6]
    {
        pt_line( 1.0, -1.0, -1.0, -1.0),
        pt_line(-1.0, -1.0, -1.0,  1.0),
        pt_line(-1.0,  1.0,  1.0,  1.0),
        pt_line( 1.0,  1.0,  1.0, -1.0),
        pt_line( 1.0, -1.0, -1.0,  1.0),
        pt_line(-1.0, -1.0,  1.0,  1.0)
    };

    return pt_seq(6, patches);
}

static Patch (*characters[])(void) =
{
    pt_char_period,
    pt_char_slash,
    pt_char_0,
    pt_char_1,
    pt_char_2,
    pt_char_3,
    pt_char_4,
    pt_char_5,
    pt_char_6,
    pt_char_7,
    pt_char_8,
    pt_char_9,
    pt_char_colon,
    pt_char_semicolon,
    pt_char_less,
    pt_char_equal,
    pt_char_greater,
    pt_char_none,
    pt_char_none,
    pt_char_a,
    pt_char_b,
    pt_char_c,
    pt_char_d,
    pt_char_e,
    pt_char_f,
    pt_char_g,
    pt_char_h,
    pt_char_i,
    pt_char_j,
    pt_char_k,
    pt_char_l,
    pt_char_m,
    pt_char_n,
    pt_char_o,
    pt_char_p,
    pt_char_q,
    pt_char_r,
    pt_char_s,
    pt_char_t,
    pt_char_u,
    pt_char_v,
    pt_char_w,
    pt_char_x,
    pt_char_y,
    pt_char_z,
};

Patch pt_text(const char* text, TextAlignment alignment)
{
    const char* start_char = text;
    const char* current_char = text;

    size_t size_patches = 32;
    Patch* text_patches = (Patch*)malloc(size_patches * sizeof(Patch));
    size_t* line_widths = (size_t*)malloc(size_patches * sizeof(size_t));
    size_t* line_sizes  = (size_t*)malloc(size_patches * sizeof(size_t));

    size_t max_width = 0;
    size_t pt_char_count = 0;
    size_t line_count = 0;
    while(true)
    {
        while(*current_char != '\0' && *current_char != '\n') ++current_char;

        size_t line_width = current_char - start_char;

        line_sizes[line_count] = 0;

        if(line_width > 0)
        {
            line_widths[line_count] = line_width;

            if(line_width > max_width) max_width = line_width;

            for(size_t i = 0; i < line_width; ++i)
            {
                if(pt_char_count >= size_patches)
                {
                    size_patches *= 2;
                    text_patches = (Patch*)realloc(text_patches, size_patches * sizeof(Patch));
                    line_widths  = (size_t*)realloc(line_widths, size_patches * sizeof(size_t));
                    line_sizes   = (size_t*)realloc(line_sizes, size_patches * sizeof(size_t));
                }

                if(start_char[i] >= 46 && start_char[i] <= 90)
                {
                    uint8_t val = start_char[i] - '.';
                    text_patches[pt_char_count] = characters[val]();
                }
                else if(start_char[i] == 32)
                {
                    line_sizes[line_count] = line_count > 0?
                        pt_char_count - line_sizes[line_count - 1]:
                        pt_char_count;

                    continue;
                }
                else
                {
                    text_patches[pt_char_count] = pt_char_none();
                }

                text_patches[pt_char_count] = pt_point(2.5 * i, -2.5 * line_count, 1.0) + text_patches[pt_char_count];

                ++pt_char_count;
            }
        }

        ++line_count;

        if(*current_char == '\0') break;

        start_char = ++current_char;
    }

    if(!pt_char_count)
    {
        free(text_patches);
        free(line_widths);
        free(line_sizes);
        return pt_char_none();
    }

    text_patches = (Patch*)realloc(text_patches, pt_char_count * sizeof(Patch));

    if(alignment == TextAlignment::CENTER || alignment == TextAlignment::RIGHT)
    {
        for(size_t wi = 0, i = 0; wi < line_count; ++wi)
        {
            //@todo: Center alignement seems to be a bit off for e.g.
            //pt_text("120.24\nSECONDS", TextAlignment::CENTER)
            double offset = (alignment == TextAlignment::CENTER)?
                0.5 * max_width - 0.5 * line_widths[wi]:
                (double)max_width - (double)line_widths[wi];

            for(size_t li = 0; li < line_sizes[wi]; ++li)
            {
                text_patches[i] = pt_point(2.5 * offset, 0.0, 1.0) + text_patches[i];
                ++i;
            }
        }
    }

    free(line_widths);
    free(line_sizes);

    Patch p = pt_seq(pt_char_count, text_patches);
    //@note: The center of each letter is actually 0,0.
    p = p + pt_point(
        -0.5 * 2.5 * (max_width - 1.0),
         0.5 * 2.5 * (line_count - 1.0),
         1.0
    );

    double a = 0.08 * 2.0 / 3.0;
    double b = 0.08 * 1.0;

    return pt_point(a, b, 1.0) * p;
}

Patch pt_cardioid(void)
{
    Patch patch =
        pt_point(1.0, 1.0, 0.0) * (
            2.0 * (1.0 - pt_lift1(cos, (2.00 * M_PI) * pt_time() * pt_frequency())) *
            pt_lift1(sin, cos, nullptr, (2.00 * M_PI) * pt_time() * pt_frequency())
        ) +
        pt_point(0.0, 2.0, 1.0);

    return patch;
}

Patch pt_random(uint32_t seed)
{
    uint32_t* data = (uint32_t*) malloc(sizeof(uint32_t));
    *data = seed;

    Patch p;
    p.userdata = (void*) data;
    p.type = PT_RANDOM;
    p.call = random_call;
    return p;
}

Patch pt_square_wave(void)
{
    return 2.0 * ((pt_frequency() * pt_fmod(pt_time(), 1.0 / pt_frequency())) > 0.5) - 1.0;
}

Patch pt_audio(const AudioData* audio_data)
{
    Patch p;
    p.userdata = (void*) audio_data;
    p.type = PT_AUDIO;
    p.call = wav_call;

    return p;
}

Patch pt_fmul(Patch patch, double factor)
{
    Patch p;
    p.userdata = malloc(sizeof(ScalarPatchData));
    p.type = PT_FMUL;
    new (p.userdata) ScalarPatchData{patch, factor};
    p.call = fmul_call;
    return p;
}

Patch pt_fset(Patch patch, double frequency)
{
    Patch p;
    p.userdata = malloc(sizeof(ScalarPatchData));
    p.type = PT_FSET;
    new (p.userdata) ScalarPatchData{patch, frequency};
    p.call = fset_call;
    return p;
}


