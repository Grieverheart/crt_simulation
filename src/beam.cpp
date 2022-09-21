#include "beam.h"
#include <cstdio>

void beam_simulate(Beam* beam, BeamData* beam_data, const Patch* p, double base_frequency, double frame_sec)
{
    beam_data->dt = frame_sec / beam->num_edges;
    // @note: we actually report one less point, but we use the fact that the last
    // point is repeated into the next batch so that we can more easily do the
    // lerp in the audio thread.
    beam_data->num_points = beam->num_edges;
    beam_data->points = new double[3 * (beam->num_edges + 1)];

    for(size_t n = 0; n < beam->num_edges; ++n)
    {
        beam_data->points[3 * n + 0] = beam->x;
        beam_data->points[3 * n + 1] = beam->y;
        beam_data->points[3 * n + 2] = beam->z;

        beam->x = 0.0;
        beam->y = 0.0;
        beam->z = 0.0;

        p->call(p->userdata, beam->sim_time + n * beam_data->dt, base_frequency, &beam->x, &beam->y, &beam->z);
    }

    beam_data->points[3 * beam->num_edges + 0] = beam->x;
    beam_data->points[3 * beam->num_edges + 1] = beam->y;
    beam_data->points[3 * beam->num_edges + 2] = beam->z;

    beam->sim_time += frame_sec;
}
