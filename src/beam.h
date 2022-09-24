#pragma once

#include <cstddef>

struct Patch
{
    void* userdata;
    void (*call)(void* userdata, double t, double f, double* x, double* y, double* z);
};

struct BeamData
{
    size_t num_points;
    double* points;
    double dt;
};

struct Beam
{
    size_t num_edges;
    double decay_time;
    double radius;

    float intensity;
    float color[3];

    double sim_time;

    double x, y, z;
};

void beam_simulate(Beam* beam, BeamData* beam_data, const Patch* p, double base_frequency, double frame_sec);

