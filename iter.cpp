#include <iostream>
#include <fstream>
#include <random>
#include "CImg.h"

using namespace cimg_library;
using namespace std;

struct PARTICLE
{
    long double mass;
    long double position_x;
    long double position_y;
    long double velocity_x;
    long double velocity_y;
};

int main()
{

    /* Units */
    double G = 6.674E-11;
    double c = 299792458;
    double s_in_year = 31557600;
    double lightyear = c * s_in_year;
    double kg_in_Ms = 1.989E30;

    /* Simulations basic vars*/
    int end_step = 10000;
    double time_scale = 1E5 * s_in_year; // 100k years per step
    int n_particles = 1000;                        // Messier 67 500-1000
    long double sim_radius = 30000 * lightyear;    // 10 * Messier 67 2600-2900

    /* Vars */
    vector<PARTICLE> particles;
    long double Fx;
    long double Fy;
    long double delta_x;
    long double delta_y;
    long double dist;
    long double sq_dist;

    /* Output*/
    int img_size = 1000;
    int pixel_x;
    int pixel_y;
    const unsigned char color_white[] = {255, 255, 255};
    CImg<unsigned char> img(img_size, img_size, 1, 1, 0);
    CImgDisplay dsp(img_size, img_size, "Iter", 0);
    dsp.display(img);

    /* RNG */
    long unsigned int seed = 1;
    mt19937 rng{seed};
    normal_distribution<> nd{0, 1};

    printf("Seed: %d\n", seed);

    /* Generate particles */

    for (size_t i = 0; i < n_particles; i++)
    {
        PARTICLE p;
        p.mass = abs(nd(rng)) * kg_in_Ms;
        p.position_x = nd(rng) * sim_radius / 10;
        p.position_y = nd(rng) * sim_radius / 10; 
        p.velocity_x = nd(rng) * 5E3;
        p.velocity_y = nd(rng) * 5E3;

        /* Cheat to make spiral galaxy */
        if (p.position_x > 0 and p.position_y > 0)
        {
            p.velocity_y = -1 * abs(nd(rng) * 5E5 * p.position_x / sim_radius);
        }     
        if (p.position_x > 0 and p.position_y < 0)
        {
            p.velocity_x = -1 * abs(nd(rng) * 5E5 * p.position_y / sim_radius);
        }
        if (p.position_x < 0 and p.position_y < 0)
        {
            p.velocity_y = +1 * abs(nd(rng) * 5E5 * p.position_x / sim_radius);
        }
        if (p.position_x < 0 and p.position_y > 0)
        {
            p.velocity_x = +1 * abs(nd(rng) * 5E5 * p.position_y / sim_radius);
        }
        /* End cheat to make spiral galaxy */

        particles.push_back(p);
    }

    /* Simulate*/
    for (size_t i = 0; i < end_step; i++)
    {

        /* Interate for every particle*/
        for (auto &&p0 : particles)
        {

            /* Interate it's interactions with other particles*/
            Fx = 0;
            Fy = 0;
            for (auto &&pn : particles)
            {
                delta_x = pn.position_x - p0.position_x;
                delta_y = pn.position_y - p0.position_y;
                sq_dist = delta_x * delta_x + delta_y * delta_y;
                Fx += (G * p0.mass * pn.mass * delta_x) / sq_dist;
                Fy += (G * p0.mass * pn.mass * delta_y) / sq_dist;
                /* Fix division */
                if (isnan(Fx))
                {
                    Fx = 0;
                }
                if (isnan(Fy))
                {
                    Fy = 0;
                }
            }

            p0.velocity_x += (Fx / p0.mass);
            p0.velocity_y += (Fy / p0.mass);
            p0.position_x += p0.velocity_x * time_scale;
            p0.position_y += p0.velocity_y * time_scale;
        }

        /* Out put current step */
        cout << "Step #" << i + 1 << "/" << end_step << "\n"; /* Print to console*/
        img.fill(0);                                          /* Clean image*/
        for (auto &&p : particles)                            /* Plot each star */
        {
            /* pos   =    img_size     *  ( ( 1 + ( position / radius ) ) / 2 )
             *         scaling to img      push -1 to 0, 0 to 0.5, and 1 to 1                      
             */
            pixel_x = img_size * ((1 + (p.position_x / sim_radius)) / 2);
            pixel_y = img_size * ((1 + (p.position_y / sim_radius)) / 2);
            img.draw_point(pixel_x, pixel_y, 0, color_white);
        }

        /* Display full image*/
        dsp.display(img);
    }

    return 0;
}
