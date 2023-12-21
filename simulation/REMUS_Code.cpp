// C++ Coursework Motion of a Rocket Simulation
// The purpose of this code is to simulate and output the basic
// flight characteristics for a model rocket that can have multiple booster stages
//User Note: In atom the default setting creates trailing whitespace, a new line, this code works
//when there are only the number of lines entered, no additional blank lines when saved
// i.e parameter file for two stage booster had 3 lines only, not 4, Atom sometimes has a default to add another line

//Standard header files
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
using namespace std;

//Environmental Class to store Environental Parameters
class EnvironmentVar {
public:
    double rho, g, Cd, v0, h0;
};

//Rocket Class to store Rocket Parameters
//Includes functions to store the total final mass and the stage burn time
class RocketVar {
public:
    double Rocket();
    vector<double> massfunc();
    vector<double> m0, mR, A, mf, ue, t;
    int nstage;
    double m_tot;

    //Vector function creates and stores the masses vector for subsequent rocket stages
    vector<double> mass() {
        vector<double> massvect;
        for (int j = 0; j < nstage; j++) {
            m_tot = m0[j] + m0[j + 1];
            massvect.push_back(m_tot);
        }
        massvect.push_back(mR[nstage - 1]);
        return massvect;
    }

    //Vector function creates and stores the stage burn time vector so that the model can loop through each stage for the duration of the burn
    vector<double> burn() {
        vector<double> burntime;
        double stagetime; //temporary variable used to create the vector for burn times
        for (int i = 0; i < nstage; i++) {
            stagetime = (m0[i] - mR[i]) / mf[i];
            burntime.push_back(stagetime);
        }
        return burntime;
    }
};
//--------------------------------------------------------------------------
int main() {
    //Calling and allocating variables to access class objects
    EnvironmentVar envp;
    RocketVar rockp;
    //Vector stores the number of lines of the input file
    //This is used to identify the number of rocket stages
    vector<string> vectparam;
    ifstream vParamfile("parameters.txt");
    string line; //Temporary variable created to read number of lines in input file

    if (vParamfile.good()) {
        while (!vParamfile.eof()) {
            getline(vParamfile, line);
            vectparam.push_back(line);
        }
    }
    else {
        cout << "File inaccessible, please make sure there is no trailing whitespace in your document file, try again! " << endl; //Console ouput in case the file cannot be opened
    }
    rockp.nstage = vectparam.size() - 1; //Identifies the number of rocket stages

    //Environmental Parameters are read from the input file and stored
    stringstream Temp(vectparam[0]);
    Temp >> envp.rho; Temp >> envp.g; Temp >> envp.Cd;
    Temp >> envp.v0; Temp >> envp.h0;
    vParamfile.close();

    //Defining Temporary variables to store rocket variables read from input file
    string m_init, m_rocket, m_area, massf, ex_vel, timestep;
    string temp;

    //Opens the file to allocate the rocket parameters
    //Each entry into each parameter vector corresponds to a new stage
    ifstream param("parameters.txt");
    getline(param, temp);
    if (param.good()) {
        while (!param.eof()) {
            getline(param, m_init, ' '); rockp.m0.push_back(stod(m_init)); //Each variable is stored in a corresponding class vector object
            getline(param, m_rocket, ' '); rockp.mR.push_back(stod(m_rocket));
            getline(param, m_area, ' '); rockp.A.push_back(stod(m_area));
            getline(param, massf, ' '); rockp.mf.push_back(stod(massf));
            getline(param, ex_vel, ' '); rockp.ue.push_back(stod(ex_vel));
            getline(param, timestep, '\n'); rockp.t.push_back(stod(timestep));
        }
        param.close();
    }
    else {
        cout << "File inaccessible, please make sure there is no trailing whitespace in your document file, try again! " << endl; //Console ouput in case the file cannot be opened
    }
    //Calling class member functions and reassinging output to variables
    vector<double> stageburnT;
    stageburnT = rockp.burn();
    vector<double> massT;
    massT = rockp.mass();

    //Opens the output file to create columns to fill with data
    ofstream outfile("output.txt", ios::out | ios::trunc);
    double m;
    int i = 0; //used for while loop to iterate through
    double h = envp.h0; //Assigning height variable
    double v = envp.v0; //Assigning velocity variable
    double dt = rockp.t[0]; ///Assigning timestep
    double maxh = h; //Assigning maximum altitude storage variable
    double duration = 0; //Assigining duration of flight variable
    outfile << setw(3) << "T" << setw(8) << "H" << setw(10) << "V" << setw(8) << "M" << endl;

    //Defining Runge Kutta Variables
    long double dragCf, mk1, vk1, vk2, vk3, vk4, hk1, hk2, hk3, hk4, freefall, term_vel;

    //Rocket Equations of motion evaluated for each stage
    for (int q = 0; q < rockp.nstage; q++) {
        // Rocket motion for each stage, mass vector called derived from class
        m = massT[q];
        for (int i = 0; dt * i <= stageburnT[q]; i++) {
            duration = duration + dt;
            mk1 = dt * -rockp.mf[q];
            //Drag coefficient Calcualtion (excluding mass)
            dragCf = 0.5 * envp.rho * envp.Cd * rockp.A[q];
            //Velocity RK4 k Values Evaluated
            vk1 = dt * (-envp.g - (dragCf * v * abs(v)) / m + (v * rockp.mf[q] * rockp.ue[q]) / (abs(v) * m));
            vk2 = dt * (-envp.g - (dragCf * (v + 0.5 * vk1) * abs(v + 0.5 * vk1)) / m + ((v + 0.5 * vk1) * rockp.mf[q] * rockp.ue[q]) / (abs(v + 0.5 * vk1) * m));
            vk3 = dt * (-envp.g - (dragCf * (v + 0.5 * vk2) * abs(v + 0.5 * vk2)) / m + ((v + 0.5 * vk2) * rockp.mf[q] * rockp.ue[q]) / (abs(v + 0.5 * vk2) * m));
            vk4 = dt * (-envp.g - (dragCf * (v + vk3) * abs(v + vk3)) / m + ((v + vk3) * rockp.mf[q] * rockp.ue[q]) / (abs(v + vk3) * m));
            //Height RK4 Values Evaluated
            hk1 = dt * v; hk2 = dt * (v + 0.5 * hk1); hk3 = dt * (v + 0.5 * hk2); hk4 = dt * (v + hk3);
            //Updating equations of motion with RK4 for every timestep
            m = m + mk1; //total mass of the rocket at any one point in time (with the stages)
            v = v + (vk1 + 2.0 * (vk2 + vk3) + vk4) / 6.0;
            h = h + (hk1 + 2.0 * (hk2 + hk3) + hk4) / 6.0;
            if ((h > maxh)) {
              maxh = h;
            }
            if ((h <0)) {
              break;
            }
            //Ouputs all data to text file with appropriate precision
            outfile << setprecision(3) << fixed << i * dt + dt << " " << setprecision(5) << fixed << h << " " << setprecision(5) << fixed << v << " " << setprecision(3) << fixed << m << " " << endl;
        }
    }
    //Rocket Equations of Motion for freefall, equations simplified
    //mdot = 0; vdot = -g - drag with m equal to rocket mass; hdot = v
    while(h>0) {
      //Drag coefficient (excluding mass)
        dragCf = 0.5 * envp.rho * envp.Cd * rockp.A.back();
        //Freefall Velocity RK4 k Values Evaluated
        vk1 = dt * (-envp.g - (dragCf * v * abs(v)) / m);
        vk2 = dt * (-envp.g - (dragCf * (v + 0.5 * vk1) * abs(v + 0.5 * vk1)) / m);
        vk3 = dt * (-envp.g - (dragCf * (v + 0.5 * vk2) * abs(v + 0.5 * vk2)) / m);
        vk4 = dt * (-envp.g - (dragCf * (v + vk3) * abs(v + vk3)) / m);
        //Freefall Height RK4 k Values Evaluated
        hk1 = dt * v; hk2 = dt * (v + 0.5 * hk1); hk3 = dt * (v + 0.5 * hk2); hk4 = dt * (v + hk3);
        //Updating equations of motion with RK4 forevery timestep
        h = h + (hk1 + 2.0 * (hk2 + hk3) + hk4) / 6.0;
        v = v + (vk1 + 2.0 * (vk2 + vk3) + vk4) / 6.0;
        m = massT.back(); // mass of the rocket only
        duration = duration + dt;
        term_vel = v;
        //Ouputs all data to text file with appropriate precision
        outfile << setprecision(3) << fixed << duration << " " << setprecision(5) << fixed << h << " " << setprecision(5) << fixed << v << " " << setprecision(3) << fixed << m << " " << endl;
        i++;
        if ((h > maxh)) {
          maxh = h;
        }
    }
    //Console Ouput providing useful data of flight profile
    cout << "Flight Profile Characteristics: " << endl;
    cout << "Number of rocket stages: " << rockp.nstage << endl;
    cout << "Simulation step size: " << dt << endl;
    cout << "Maximum Altitude reached: " << maxh << " m" << endl;
    cout << "Terminal Velocity: " << term_vel << " m/s" << endl;
    cout << "Time of flight " << duration << " s" <<endl;
    outfile.close(); //Closes the output file
}