#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <cmath>
#include <tuple>
#include <vector>

#include "EquazioneDifferenzialeBase.hpp"
#include "DoublePendulum.hpp"
#include "RungeKutta.hpp"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMultiGraph.h"
#include "TColor.h"

// TODO: add phase space
// TODO: make the lenghts different and the masses different, e con i marker di dimensione diversa
// TODO: aggiungere l'attrito

using namespace std;

std::tuple<float, float, float> HSVtoRGB(float h, float s, float v);

int main(int argc, char **argv)
{
    if (argc > 6){
        cout << "Usage: " << argv[0]
             << "<number of pendulums, default: 1> "
             << "<show trace (1 true/0 false, default: 1)> "
             << "<show pendulums (options: [0: don't draw, 1: draw just the point, 2: draw all], default: 2)> "
             << "<trace limit (default: 100)>"
             << "<configuration space or phase space (options: [0: configuration space, 1: phase space], default: 0)>"
             << ", note if you want to set trace limit, you also have to set all other 3 parameters.\n"
             << "If you want to keep the trace use a big number."
             << endl;
        return -1;
    }

    int npend = argc >= 2 ? atoi(argv[1]) : 1;
    int show_pendulums = argc >= 3 ? atoi(argv[2]) : 2;
    bool show_trace = argc >= 4 ? atoi(argv[3]) : true;
    int trace_limit = argc >= 5 ? atoi(argv[4]) : 100;

    if (!show_pendulums and !show_trace){
        cout << "WHY WOULD YOU WANT TO PLOT AN EMPTY GRAPH?" << endl;
        return -104;
    }

    if (npend < 1){
        cout << "Number of pendulums must be at least 1." << endl;
        return -1;
    }

    if (show_pendulums > 2 or show_pendulums < 0){
        cout << "Show pendulums should be 0, 1 or 2." << endl;
        return -1;
    }

    if (trace_limit <= 0){
        show_trace = false;
    }

    // simulation parameters
    double spread = 0.01;
    double h = 1E-2;


    TApplication myApp("myApp", &argc, argv);
    TCanvas *c1 = new TCanvas("c1", "Doppio pendolo", 800, 800);

    RungeKutta myRK;
    DoublePendulum pend;

    double tmax = 100; // [s]
    double t = 0.;

    vector<vector<double>> positions(npend, vector<double>(4));
    for (int i = 0; i < npend; i++){
        positions[i] = {(1 - 2 * spread) * M_PI, (1 - spread /2.) * M_PI + spread * M_PI * i / npend, 0., 0.};
    }


    vector<TGraph> pendulums(npend);
    vector<TGraph> traces(npend);

    TMultiGraph myGraph;

    // generate the graphs
    float k;
    int colorIndex;

    for (int i = 0; i < npend; i++){
        // color the pendulum
        k = double(360.0 / npend) * i;
        auto [r, g, b] = HSVtoRGB(k, 1.0, 1.0);
        colorIndex = 5000 + i;
        new TColor(colorIndex, r, g, b);

        // graph pendulum
        if (show_pendulums){
            pendulums[i].SetLineColor(colorIndex);
            pendulums[i].SetMarkerColor(colorIndex);
            pendulums[i].SetLineWidth(4);
            pendulums[i].SetMarkerStyle(8);
            pendulums[i].SetMarkerSize(1);

            // add to multigraph
            myGraph.Add(&pendulums[i], "LP");
        }

        // graph pendulum trace
        if (show_trace){
            traces[i].SetLineColor(colorIndex);
            traces[i].SetLineWidth(2);

            // add to multigraph
            myGraph.Add(&traces[i], "L");
        }
    }

    // wait time (to record the screen)
    //std::this_thread::sleep_for(std::chrono::milliseconds(5000));


    int nstep = tmax / h + 0.5;

    // simulation loop
    for (int step = 0; step < nstep; step++)
    {
        c1->Clear();
        c1->DrawFrame(-2.2, -2.2, 2.2, 2.2, "Double Pendulum Animation;X-axis;Y-axis");

        for (int j = 0; j < npend; j++){
            double x1 = sin(positions[j][0]);
            double y1 = -cos(positions[j][0]);
            double x2 = x1 + sin(positions[j][1]);
            double y2 = y1 - cos(positions[j][1]);

            if (show_pendulums){
                pendulums[j].SetPoint(0, 0, 0);
                pendulums[j].SetPoint(1, x1, y1);
                pendulums[j].SetPoint(2, x2, y2);

                pendulums[j].Draw("LP SAME");
            }

            if (show_trace){
                // Add the current point to the trace
                if (step < trace_limit) {
                    traces[j].SetPoint(step, x2, y2);
                } else {
                    // Shift the points: remove the first point and add the new point
                    for (int k = 1; k < trace_limit; k++) {
                        double px, py;
                        traces[j].GetPoint(k, px, py);
                        traces[j].SetPoint(k - 1, px, py);
                    }
                    traces[j].SetPoint(trace_limit - 1, x2, y2);
                }

                traces[j].Draw("L SAME");
            }

            positions[j] = myRK.Passo(t, positions[j], h, pend);
        }

        c1->Update();

        t += h;

        gSystem->ProcessEvents();

        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    myApp.Run();

    return 0;
}

// HSV to RGB conversion function
std::tuple<float, float, float> HSVtoRGB(float h, float s, float v) {
    float c = v * s;
    float x = c * (1 - fabs(fmod(h / 60.0, 2) - 1));
    float m = v - c;
    float r, g, b;
    if (h >= 0 && h < 60) { r = c; g = x; b = 0; }
    else if (h >= 60 && h < 120) { r = x; g = c; b = 0; }
    else if (h >= 120 && h < 180) { r = 0; g = c; b = x; }
    else if (h >= 180 && h < 240) { r = 0; g = x; b = c; }
    else if (h >= 240 && h < 300) { r = x; g = 0; b = c; }
    else { r = c; g = 0; b = x; }
    return std::make_tuple((r + m), (g + m), (b + m));
}