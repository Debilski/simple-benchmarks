#include <iostream>
#include <fstream>

#include <boost/generator_iterator.hpp>
#include <boost/program_options.hpp>
#include <boost/random.hpp>

#include <boost/timer/timer.hpp>
#include "progress.hpp"


typedef boost::mt19937 RNGType;
RNGType rng;

boost::uniform_01<> rand01;
boost::variate_generator<RNGType, boost::uniform_01<>> dice(rng, rand01);

double lambda = 1;
double kappa = -0.05;
double tau = 10;
double sigma = 1;

void r_step(double dt, int i, int j, double **uu, double **vv)
{
    double u = uu[i][j];
    double v = vv[i][j];

    double du = (lambda * u - u * u * u + kappa - sigma * v) * dt;
    double dv = (1. / tau * (u - v)) * dt;

    uu[i][j] = u + du;
    vv[i][j] = v + dv;
}

inline int positive_modulo(int i, int n)
{
    return (i % n + n) % n;
}

void diff_step(double dt, int size_x, int size_y, double **uu, double **vv, double **diff_u, double **diff_v)
{
    for (int i = 0; i < size_x; ++i)
    {
        for (int j = 0; j < size_y; ++j)
        {
            diff_u[i][j] = (-4 * uu[i][j] + uu[positive_modulo(i - 1, size_x)][j] + uu[positive_modulo(i + 1, size_x)][j] + uu[i][positive_modulo(j - 1, size_y)] + uu[i][positive_modulo(j + 1, size_y)]) * dt;
            diff_v[i][j] = (-4 * vv[i][j] + vv[positive_modulo(i - 1, size_x)][j] + vv[positive_modulo(i + 1, size_x)][j] + vv[i][positive_modulo(j - 1, size_y)] + vv[i][positive_modulo(j + 1, size_y)]) * dt;
        }
    }

    for (int i = 0; i < size_x; ++i)
    {
        for (int j = 0; j < size_y; ++j)
        {
            uu[i][j] = uu[i][j] + 0.0000028 * diff_u[i][j];
            vv[i][j] = vv[i][j] + 0.005 * diff_v[i][j];
        }
    }
}

int main(int argc, char *argv[])
{

    using namespace boost::program_options;

int size_x, size_y, num_steps;
double dt;

      try
  {
    options_description desc{"Reaction–Diffusion system"};
    desc.add_options()
      ("help,h", "Help screen")
      ("size_x", value<int>())
      ("size_y", value<int>())
      ("steps", value<int>())
      ("dt", value<double>());

    variables_map result;
    store(parse_command_line(argc, argv, desc), result);

    if (result.count("help"))
      std::cout << desc << '\n';

    size_x = result["size_x"].as<int>();
    size_y = result["size_y"].as<int>();

    num_steps = result["steps"].as<int>();
    dt = result["dt"].as<double>();

  }
  catch (const error &ex)
  {
    std::cerr << ex.what() << '\n';
  }

    double **uu;
    double **vv;

    double **diff_u;
    double **diff_v;
    uu = new double *[size_x];
    vv = new double *[size_x];

    diff_u = new double *[size_x];
    diff_v = new double *[size_x];

    double fix_u = - pow(abs(kappa), 1./3.);
    double fix_v = - pow(abs(kappa), 1./3.);

    for (int i = 0; i < size_x; ++i)
    {
        uu[i] = new double[size_y];
        vv[i] = new double[size_y];

        diff_u[i] = new double[size_y];
        diff_v[i] = new double[size_y];
        // each i-th pointer is now pointing to dynamic array (size 10) of actual int values
        for (int j = 0; j < size_y; ++j)
        { // for each column
            uu[i][j] = fix_u;
            vv[i][j] = fix_v;

            diff_u[i][j] = 1;
            diff_v[i][j] = 1;
        }
    }

    std::ofstream mydata;
    mydata.open("data2.txt");
    uu[0][0] = -0.1;

    boost::progress_display show_progress(num_steps);

    {
        boost::timer::auto_cpu_timer t;

        for (int step = 0; step < num_steps; ++step)
        {
            ++show_progress;

            mydata << uu[0][0] << " ";
            vv[0][0] += 0.01 * (dice() - 0.5);

            for (int i = 0; i < size_x; ++i)
            {
                for (int j = 0; j < size_y; ++j)
                { // for each column
                    if (step < 1000)
                        vv[i][j] += 0.01 * (dice() - 0.5);
                    r_step(dt, i, j, uu, vv);
                }
            }

            diff_step(dt, size_x, size_y, uu, vv, diff_u, diff_v);
        }
    }

    std::ofstream uufile, vvfile;
    uufile.open("uu.txt");
    vvfile.open("vv.txt");

    for (int i = 0; i < size_x; ++i)
    {
        for (int j = 0; j < size_y; ++j)
        { // for each column
            uufile << uu[i][j] << " ";
            vvfile << vv[i][j] << " ";
        }
        uufile << std::endl;
        vvfile << std::endl;
    }
}
