#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <cmath>
#include <tuple>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <unordered_map>
#include <functional>

#include "EquazioneDifferenzialeBase.hpp"
#include "DoublePendulum.hpp"
#include "RungeKutta.hpp"
#include "EulerMethod.hpp"

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

tuple<float, float, float> HSVtoRGB(float h, float s, float v);

int main(int argc, char **argv)
{
    // ======= Settings =======
    YAML::Node config = YAML::LoadFile("SETTINGS.yaml");

    cout << "Settings recap: " << endl;

    // Accessing 'visualization' parameters
    auto visualization       = config["visualization"];
    bool showConfigSpace     = visualization["show_config_space"].as<bool>();
    bool showPhaseSpace      = visualization["show_phase_space"].as<bool>();
    bool showGrid            = visualization["show_grid"].as<bool>();
    bool showPendulums       = visualization["show_pendulums"].as<bool>();
    bool showTrace           = visualization["show_trace"].as<bool>();
    unsigned int traceLength = visualization["trace_length"].as<unsigned int>();

    cout << "  Visualization Parameters:\n";
    cout << "    Show Config Space: " << showConfigSpace << "\n";
    cout << "    Show Phase Space: "  << showPhaseSpace  << "\n";
    cout << "    Show Grid: "         << showGrid        << "\n";
    cout << "    Trace Length: "      << traceLength     << "\n";

    // Accessing 'simulation_parameters'
    auto simulation         = config["simulation_parameters"];
    unsigned int nPendulums = simulation["number_of_pendulums"].as<unsigned int>();
    double timeStep         = simulation["time_step"].as<double>();
    double totalTime        = simulation["total_time"].as<double>();

    cout << "\n  Simulation Parameters:\n";
    cout << "    Number of Pendulums: " << nPendulums << "\n";
    cout << "    Time Step: "           << timeStep   << "\n";
    cout << "    Total Time: "          << totalTime  << "\n";

    // Accessing 'advanced' settings
    auto advanced                  = config["advanced"];
    string integrator              = advanced["integration_method"].as<string>();
    bool energy_log                = advanced["energy_logging"].as<bool>();
    bool multi_threading           = advanced["multi_threading"].as<bool>();
    unsigned int number_of_threads = advanced["thread_count"].as<unsigned int>();

    cout << "\n  Advanced Parameters:\n";
    cout << "    Integration Method: " << integrator        << "\n";
    cout << "    Energy Log: "         << energy_log        << "\n";
    cout << "    Multi Threading: "    << multi_threading   << "\n";
    cout << "    Number of Threads: "  << number_of_threads << "\n";

    // Accessing 'pendulums' settings
    auto pendulums_conf = config["pendulums"];
    double gravity      = pendulums_conf["gravity"].as<double>();
    double mass1        = pendulums_conf["mass1"].as<double>();
    double mass2        = pendulums_conf["mass2"].as<double>();
    double length1      = pendulums_conf["length1"].as<double>();
    double length2      = pendulums_conf["length2"].as<double>();
    double start_theta1 = pendulums_conf["theta1"].as<double>();
    double start_theta2 = pendulums_conf["theta2"].as<double>();
    double damping1     = pendulums_conf["damping1"].as<double>();
    double damping2     = pendulums_conf["damping2"].as<double>();
    double spread       = pendulums_conf["spread"].as<double>();

    cout << "\n  Pendulums Settings:\n";
    cout << "    Gravity: "      << gravity      << "\n";
    cout << "    Mass1: "        << mass1        << "\n";
    cout << "    Mass2: "        << mass2        << "\n";
    cout << "    Length1: "      << length1      << "\n";
    cout << "    Length2: "      << length2      << "\n";
    cout << "    Start_theta1: " << start_theta1 << "\n";
    cout << "    Start_theta2: " << start_theta2 << "\n";
    cout << "    Damping1: "     << damping1     << "\n";
    cout << "    Damping2: "     << damping2     << "\n";
    cout << "    Spread: "       << spread       << "\n";

    // check if settings are correct
    if (!(showConfigSpace || showPhaseSpace))                              { cout << "Config space and Phase space not showing why do you want to even run the simulation?!?\n"; exit(1); }
    if (!(showPendulums || showTrace))                                     { cout << "Both pendulums and trace are deactivated!\n"; exit(1); }
    if (traceLength < 2 && showTrace)                                      { cout << "Trace length must be greater than 2\nIf you don't want it deactivate it."; exit(1); }
    if (!nPendulums)                                                       { cout << "Number of pendulums cannot be zero!"; exit(1); }
    if (timeStep > totalTime)                                              { cout << "Time step is greater than total time!"; exit(1); }
    if (timeStep <= 0)                                                     { cout << "Time step cannot be less than or equal to zero!"; exit(1); }
    if (!(integrator == "eu2" || integrator == "rk4" || integrator == "v")){ cout << "Integration method not recognized!\n"; exit(1); }
    if (number_of_threads < 2)                                             { cout << "Number of Threads must be greater than 2\n"; exit(1); }

    // create the pendulums and positions
    vector<DoublePendulumDamped> pendulums;
    vector<vector<double>> positions;
    for (unsigned int i = 0; i < nPendulums; ++i) {
        pendulums.push_back(DoublePendulumDamped(gravity, mass1, mass2, length1, length2, damping1, damping2));
        positions.push_back(vector<double> {start_theta1, start_theta2 * (1 - spread / 2. + spread * i/nPendulums), 0, 0});
    }

    // Map each setting to a function
    unordered_map<string, function<void(DoublePendulumDamped&, double)>> pendulumActions;
    pendulumActions["gravity"]  = [](DoublePendulumDamped &pend, double value) { pend.SetG(value); };
    pendulumActions["mass1"]    = [](DoublePendulumDamped &pend, double value) { pend.SetM1(value); };
    pendulumActions["mass2"]    = [](DoublePendulumDamped &pend, double value) { pend.SetM2(value); };
    pendulumActions["length1"]  = [](DoublePendulumDamped &pend, double value) { pend.SetL1(value); };
    pendulumActions["length2"]  = [](DoublePendulumDamped &pend, double value) { pend.SetL2(value); };
    pendulumActions["damping1"] = [](DoublePendulumDamped &pend, double value) { pend.SetGamma1(value); };
    pendulumActions["damping2"] = [](DoublePendulumDamped &pend, double value) { pend.SetGamma2(value); };
    unordered_map<string, function<void(vector<double>&, double)>> positionsActions;
    positionsActions["theta1"] = [](vector<double> &pos, double value) { pos[0] = value; };
    positionsActions["theta2"] = [](vector<double> &pos, double value) { pos[1] = value; };

   	// overwrite single pendulums settings using the maps provided earlier
    for (const auto& node : config) {
        string section = node.first.as<string>();

        if (section.find("pendulum_") != 0)
            continue;

        int pendulumIndex = stoi(section.substr(9)); // Extract the index

        // check if index is valid
        if (pendulumIndex > nPendulums || pendulumIndex < 1) {
            cout << "Pendulum Index " << pendulumIndex << " is out of range\n";
            exit(3);
        }

        cout << "\n  Override for Pendulum " << pendulumIndex << ":\n";

        // Loop through settings for the specific pendulum
        for (const auto& setting : node.second) {
            string key = setting.first.as<string>();
            double value = setting.second.as<double>();
            cout << "    " << key << ": " << value << "\n";

            // Call the corresponding function based on the key
    		if (pendulumActions.find(key) != pendulumActions.end())
        		pendulumActions[key](pendulums[pendulumIndex -1], value);
            else if (positionsActions.find(key) != positionsActions.end())
                positionsActions[key](positions[pendulumIndex -1], value);
            else {
        		cout << "Unknown pendulum setting: " << key << "\n";
        		exit(2);
    		}
        }
    }

    // simulation
    TApplication myApp("myApp", &argc, argv);
    TCanvas *c1 = new TCanvas("c1", "Doppio pendolo", 800, 800);

    EquazioneDifferenzialeBase *myInteg = nullptr;
    if (integrator == "rk4")
        myInteg = new RungeKutta();
    if (integrator == "eu2")
        myInteg = new Eulero();
    if (integrator == "v"){
        cout << "Verlet not yet implemented\n";
        exit(4);
    }

    double t = 0.;

    vector<TGraph> Gpends(nPendulums);
    vector<TGraph> traces(nPendulums);

    TMultiGraph myGraph;

    // generate the graphs
    float k;
    int colorIndex;

    for (int i = 0; i < nPendulums; i++){
        // color the pendulum
        k = double(360.0 / nPendulums) * i;
        auto [r, g, b] = HSVtoRGB(k, 1.0, 1.0);
        colorIndex = 5000 + i;
        new TColor(colorIndex, r, g, b);

        // graph pendulum
        if (showPendulums){
            Gpends[i].SetLineColor(colorIndex);
            Gpends[i].SetMarkerColor(colorIndex);
            Gpends[i].SetLineWidth(4);
            Gpends[i].SetMarkerStyle(8);
            Gpends[i].SetMarkerSize(1);

            // add to multigraph
            myGraph.Add(&Gpends[i], "LP");
        }

        // graph pendulum trace
        if (showTrace){
            traces[i].SetLineColor(colorIndex);
            traces[i].SetLineWidth(2);

            // add to multigraph
            myGraph.Add(&traces[i], "L");
        }
    }

    // wait time (to record the screen)
    //this_thread::sleep_for(chrono::milliseconds(1000));


    int nstep = totalTime / timeStep + 0.5;
    double max_pend_size = 0;

    // simulation loop
    for (int j = 0; j < nPendulums; j++){
        if (showPendulums)
       		Gpends[j].SetPoint(0, 0, 0);
        max_pend_size = max(max_pend_size, pendulums[j].GetL1() + pendulums[j].GetL2());
    }
    max_pend_size *= 1.01; // add some padding



    for (int step = 0; step < nstep; step++)
    {
        c1->Clear();
        c1->DrawFrame(-max_pend_size, -max_pend_size, max_pend_size, max_pend_size, "Double Pendulum Animation;X-axis;Y-axis");

        for (int j = 0; j < nPendulums; j++){
            double x1 =      pendulums[j].GetL1() * sin(positions[j][0]);
            double y1 =    - pendulums[j].GetL1() * cos(positions[j][0]);
            double x2 = x1 + pendulums[j].GetL2() * sin(positions[j][1]);
            double y2 = y1 - pendulums[j].GetL2() * cos(positions[j][1]);

            if (showPendulums){
                Gpends[j].SetPoint(1, x1, y1);
                Gpends[j].SetPoint(2, x2, y2);

                Gpends[j].Draw("LP SAME");
            }

            if (showTrace){
                // Add the current point to the trace
                if (step < traceLength) {
                    traces[j].SetPoint(step, x2, y2);
                } else {
                    // Shift the points: remove the first point and add the new point
                    for (int k = 1; k < traceLength; k++) {
                        double px, py;
                        traces[j].GetPoint(k, px, py);
                        traces[j].SetPoint(k - 1, px, py);
                    }
                    traces[j].SetPoint(traceLength - 1, x2, y2);
                }

                traces[j].Draw("L SAME");
            }

            positions[j] = myInteg->Passo(t, positions[j], timeStep, pendulums[j]);
        }

        c1->Update();

        t += timeStep;

        gSystem->ProcessEvents();

        //this_thread::sleep_for(chrono::milliseconds(50));
    }
    myApp.Run();

    return 0;
}

// HSV to RGB conversion function
tuple<float, float, float> HSVtoRGB(float h, float s, float v) {
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
    return make_tuple((r + m), (g + m), (b + m));
}