//
//  World.cpp
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 08/04/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

//Header files
#include "Animal.h"
#include "Sensor.h"
#include "World.h"



std::vector<Sensor*>  World::trappinggrid(int NoSensors, double cam_interval, std::ofstream &Sensors, std::vector<double> SensorWidth, std::vector<double> SensorRadius, int NoAnimals){
    
    int no_cam_locations_y = sqrt(NoSensors); int no_cam_locations_x = sqrt(NoSensors);

    std::vector<Sensor*> AllSensors(pow(no_cam_locations_y,2));
    
    int sensorcount=0;
    int LengthSW = SensorWidth.size();
    int LengthSR = SensorRadius.size();
    
    
    for(int grid_number_x = 0; grid_number_x < no_cam_locations_x; grid_number_x ++){
        double x_location = grid_number_x*cam_interval;
        
        for(int grid_number_y = 0; grid_number_y < no_cam_locations_y; grid_number_y ++){
            double y_location = grid_number_y*cam_interval;
            
            for(int sensor1 =0; sensor1<LengthSW; sensor1 ++){
                for(int sensor=0; sensor<LengthSR; sensor++){
                    //std::cout<< "sensor: " << sensorcount << "/"<< NoSensors <<std::endl;
                    AllSensors[sensorcount] =new Sensor(NoAnimals, x_location, y_location, SensorRadius[sensor],0, SensorWidth[sensor1],sensorcount);
                    //Saves the locations and the angle of the Sensor
                    Sensors << sensorcount << //1st column
                        "," << x_location << //2nd column
                        "," << y_location << //...
                        "," << 0 << //...
                        "," << SensorWidth[sensor1] << //4th column
                        "," << SensorRadius[sensor] << //4th column
                        "," << "Grid" <<
                        "," << grid_number_x <<
                        "," << grid_number_y <<
                        "\n";
                    sensorcount +=1;
                };
            };
        };
    };
    return AllSensors;
};


void World::AnimalMovement(int id, double randomstart,
                           int Iteration, int SaveMovement,
                           double CentreHome_r,double MaximumDistance,
                           double no_of_move_states,std::vector<double> probability,
                           std::vector<std::vector<double>> mean_vector,std::vector<std::vector<double>> variance_vector,
                           std::vector<std::vector<double>> transitions,
                           int NoSteps, std::ofstream &Movement ,  std::vector<Sensor*> AllSensors , std::ofstream &Captures,
                           double cam_interval,int  NoSensors, int Sex,
                           double buffer){
    srand(randomstart);
    double seed1, seed2, seed3,seed4,seed5, seed6, temp;
    for(int j=0; j<251; j++){
        if(j==50){seed1=double(rand());}
        else if(j==100){seed2=double(rand());}
        else if(j==150){seed3=double(rand());}
        else if(j==200){seed4=double(rand());}
        else if(j==250){seed5=double(rand());}
        else {temp=double(rand());};
    };
    
    //sets random location
    double xlocation = StartLocation(NoSensors, cam_interval,seed1, CentreHome_r);
    double ylocation = StartLocation(NoSensors, cam_interval,seed2, CentreHome_r);
    
    // To choose a start angle, sets up a random number class
    // Uses a radom number from stream RandomNumberStreamAnimalAngle for a random seed
    srand(seed3);double CurrentAngleTemp = ((double) rand()/RAND_MAX)*2*M_PI;
    
    // New animal given start locations - at the centre of the home range
    //  Inputs are: ID & Starts location (x,y) &  Initial angle
    Animal Animal1(id,Iteration, SaveMovement,
                   CentreHome_r,MaximumDistance,
                   no_of_move_states,probability,mean_vector,variance_vector,
                   transitions,
                   xlocation, ylocation, CurrentAngleTemp, seed4,
                   Movement,  AllSensors , Captures, Sex);
    
    //Sets seed for a random number
    //Random number stream for the movemnet of the animal
    srand(seed5);
    for(int j=0; j<NoSteps; j++){
        seed6 =double(rand());
        for(int extra=0; extra<999; extra++){temp=double(rand());};
        Animal1.UpdateLocation(seed6, Movement,  AllSensors, Captures);
    }; //End of j loop for Steps

}

double World::StartLocation(int NoSensors, double cam_interval, double seed, double CentreHome_r){
    
    int no_locations = sqrt(NoSensors);
    double width = (no_locations-1)*cam_interval;
    
    srand(seed); double output= ((double) rand()/RAND_MAX)*(width + CentreHome_r*2)-CentreHome_r;
    
    return(output);
};




std::string World::make_filename( const std::string& directory ,const std::string& basename){
    std::ostringstream result;
    result << directory <<basename;
    return result.str();
};



void World::oneiteration(double cam_interval,
                int NoSteps,
                int Iteration, int SaveMovement,
                std::vector<double> CentreHome_r,double MaximumDistance,
                std::vector<double>  no_of_move_states,std::vector<std::vector<double>> probability,
                std::vector<std::vector<std::vector<double>>> mean_vector,std::vector<std::vector<std::vector<double>>> variance_vector,
                 std::vector<std::vector<std::vector<double>>> transitions,
                double seed,
                std::ofstream &Captures, std::ofstream &Movement,std::vector<Sensor*> AllSensors, int NoAnimal, double propMale,int NoSensors,
                         std::ofstream &Settings){
    
    double seed1, temp;
    int sex;
    int nomales = propMale*NoAnimal;
    int lengthsensors = AllSensors.size();
    double buffer ;
   // if(lengthsensors==2){ buffer = fmax(CentreHome_r[0],CentreHome_r[1]); } else if(lengthsensors==1){buffer = CentreHome_r[0];} else{std::cout<<"buffer problem"; EXIT_FAILURE;};
    buffer = fmax(CentreHome_r[0],CentreHome_r[1]);
    std::cout<<"nomales " <<nomales <<std::endl;
    double area = pow(buffer*2 + cam_interval*(sqrt(NoSensors)-1),2);
    double density = NoAnimal/area;
    Settings << NoAnimal << "," << NoSteps <<"," <<Iteration << "," << nomales <<"," << area <<","<< density << "\n";
    //run animals
    for(int i=0; i<NoAnimal; i++){
        
        srand(seed);
        for(int k=0; k<NoAnimal;k++){temp=double(rand());}
        seed1=double(rand());
        for(int k=0; k<NoAnimal;k++){temp=double(rand());}
        seed= double(rand());

        if(i<nomales){sex = 1;} else{sex =0;};
        
        std::cout<< "Animal " <<i+1 <<" / "<<NoAnimal<< " r " << CentreHome_r[sex] <<std::endl;
        AnimalMovement(i, seed1, Iteration, SaveMovement,
                        CentreHome_r[sex],MaximumDistance,no_of_move_states[sex],probability[sex],mean_vector[sex],variance_vector[sex], transitions[sex],
                       NoSteps, Movement ,  AllSensors , Captures,cam_interval, lengthsensors, sex, buffer);
    }

};

void World::MultipleIterations(int NoSensors, int NoInterations, std::string savevalue,
                                 double cam_interval, std::vector<double> SensorWidth, std::vector<double> SensorRadius,
                                 int NoSteps,
                                 int SaveMovement,
                               std::vector<double> CentreHome_r,double MaximumDistance,
                               std::vector<double>  no_of_move_states,std::vector<std::vector<double>> probability,
                               std::vector<std::vector<std::vector<double>>> mean_vector,std::vector<std::vector<std::vector<double>>> variance_vector,
                               std::vector<std::vector<std::vector<double>>> transitions,
                               int NoAnimals,
                                double propMale
                                 ){
    
    //Creates file for Sensors (CSV file) and writes in the header
    std::ofstream Sensors;
    Sensors.open(make_filename(savevalue, ",Sensors.csv" ).c_str());
    Sensors << "sensor.ID" <<
    "," << "x.location" <<
    "," << "y.location" <<
    "," << "centre.angle" <<
    "," << "half.width.angle" <<
    "," << "radius" <<
    "," << "placement" <<
    "," << "grid.row" <<
    "," << "grid.column" <<
    "\n";
    
    //Creates file for Captures (CSV file) and writes in the header
    std::ofstream CapturesNotRef;
    std::ofstream &Captures = CapturesNotRef;
    Captures.open(make_filename(savevalue, ",Captures.csv" ).c_str());
    Captures << "animal.number" <<
    "," << "step.number" <<
    "," << "sensor.ID" <<
    "," << "iteration.number" <<
    "," << "x.location" <<
    "," << "y.location" <<
    "," << "time" <<
    "," << "angle.from.camera.to.animal" <<
    "," << "angle.from.animal.to.camera" <<
    "," << "distance.to.camera" <<
    "," << "in.range"
    "," << "just.entered"
    "\n";
    
    //Creates file for Movement (CSV file) and writes in the header
    std::ofstream MovementNotRef;
    std::ofstream &Movement = MovementNotRef;
    Movement.open(make_filename(savevalue, ",Movement.csv" ).c_str());
    Movement << "animal.number" <<
    "," << "step.number" <<
    "," << "x.location" <<
    "," << "y.location" <<
    "," << "movement.angle" <<
    "," << "total.distance" <<
    "," << "current.state" <<
    "," << "edge.interaction" <<
    "," << "distance.moved" <<
    "," << "iteration.number" <<
    "," << "sex" <<
    "\n";
    
    //Creates file for Movement (CSV file) and writes in the header
    std::ofstream SettingsNotRef;
    std::ofstream &Settings = SettingsNotRef;
    Settings.open(make_filename(savevalue, ",Settings.csv" ).c_str());
    Settings << "number.of.animals" <<
    "," << "number.of.steps" <<
    "," <<"iteration.number" <<
    "," << "number.of.males" <<
    "," << "area" <<
    ","<< "density" <<
    "\n";
    
    
    //set up grid
    std::vector<Sensor*>  AllSensors = trappinggrid(NoSensors, cam_interval, Sensors,  SensorWidth, SensorRadius, NoAnimals);
    
 
    double seed, temp;
    for(int iteration =0; iteration<NoInterations;iteration++){
        
        srand(iteration);
        for(int j=0; j<NoInterations; j++){
            seed=double(rand());
        };
        
        std::cout << "iteration "  << iteration << "/" << NoInterations << std::endl;
        
        oneiteration(cam_interval,
                     NoSteps,
                     iteration, SaveMovement,
                     CentreHome_r,MaximumDistance,
                     no_of_move_states, probability,
                     mean_vector, variance_vector,transitions,
                     seed,
                     Captures, Movement, AllSensors, NoAnimals, propMale, NoSensors, Settings);
      
    };
    
    // Closes files that are open in all iterations
    Captures.close();
    Movement.close();
    CapturesNotRef.close();
    MovementNotRef.close();
    for(int i=0; i<NoSensors; i++) {delete AllSensors[i];};
};


double World::string_to_double( const std::string& s )
{
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        return 0;
    return x;
};

std::vector<std::vector <std::vector <double>>> World::readparamters(std::string nameoffile){
    
    std::vector<std::vector<double>> mean;
    std::vector<std::vector<double>> variance;
    std::vector<double> probability;
    std::vector<std::vector<double>> transistions;
    
    std::vector<double> meantemp(2);
    std::vector<double> vartemp(4);
    std::vector<double> probtemp(10);
    
    std::vector<double>  noofclusters(1);
    std::vector<double>  maxarea(1);
    
    std::string line;
    std::ifstream input;
    input.open(nameoffile);
    
    int linecounter =0;
    if (input.is_open()) {
        
        while (getline(input, line)) {
            //std::cout <<"linecounter: " <<linecounter << std::endl;
            std::stringstream ss(line);
            std::string l;
            
            if(linecounter ==0){ // No clusters
                getline(ss, l, ' ');
                noofclusters[0] = string_to_double(l);

            }
            else if(linecounter < (noofclusters[0]+1)){
                int count=0;
                while (count<2) {
                    getline(ss, l, ' ');
                    //std::cout <<"mean l: " <<l <<" ,count: "<< count << std::endl;
                    meantemp[count] = string_to_double(l);
                    count+=1;
                }
                mean.push_back (meantemp);
            }
            else if(linecounter < (2*noofclusters[0]+1)){
                int count=0;
                while (count<4) {
                    getline(ss, l, ' ');
                    //std::cout <<"var l: " <<l <<" ,count: "<< count << std::endl;
                    vartemp[count] = string_to_double(l);
                    count+=1;
                }
                variance.push_back (vartemp);
            }
            else if(linecounter < (2*noofclusters[0]+2)){
                int count=0;
                while (count<noofclusters[0]) {
                    getline(ss, l, ' ');
                    //std::cout <<"prob l: " <<l <<" ,count: "<< count << std::endl;
                    probtemp[count] = string_to_double(l);
                    count+=1;
                }
            }
            else if(linecounter < (2*noofclusters[0]+3)){
                int count=0;
                while (count<1) {
                    getline(ss, l, ' ');
                    //std::cout <<"prob l: " <<l <<" ,count: "<< count << std::endl;
                    maxarea[0] = string_to_double(l);
                    count+=1;
                }
            }
            else {
                int count=0;
                std::vector<double>  transtemp(noofclusters[0]);
                while (count<noofclusters[0]) {
                    getline(ss, l, ' ');
                    //std::cout <<"trans l: " <<l <<" ,count: "<< count << std::endl;
                    transtemp[count] = string_to_double(l);
                    count+=1;
                }
                //std::cout <<"transtemp: " <<transtemp[0] <<" "<< transtemp[1] << std::endl;
                transistions.push_back(transtemp);
            };
            
            linecounter+=1;
        }
        
        input.close();
    } else{std::cout<<"input probs" <<std::endl;}// END OF CAMERA CHECK
    
    
    std::vector<std::vector <std::vector <double>>> parameters(6);
    std::vector <std::vector <double>> clusters = {noofclusters};
    std::vector <std::vector <double>> maximumarea = {maxarea};
    std::vector <std::vector <double>> states = {probtemp};
    
    parameters[0] = clusters;
    parameters[1] = mean;
    parameters[2] = variance;
    parameters[3] = states;
    parameters[4] = transistions;
    parameters[5] = maximumarea;
    
    return parameters;
    
}

void World::runsimulation(std::vector<std::string> inputs, int NoSensors, int NoInterations, std::string savevalue,
                          double cam_interval, std::vector<double> SensorWidth, std::vector<double> SensorRadius,
                          int NoSteps,
                          int SaveMovement,
                          double MaximumDistance,
                          int NoAnimals, double propMale){
    int nostrings= inputs.size();
    
    std::vector<double>  no_of_move_states;
    std::vector<std::vector<double>> probability;
    std::vector<std::vector<std::vector<double>>> mean_vector;
    std::vector<std::vector<std::vector<double>>> variance_vector;
    std::vector<std::vector<std::vector<double>>> transitions;
    std::vector<double> CentreHome_r;
    
    for(int f=0; f<nostrings;f++){
        
        std::vector<std::vector <std::vector <double>>> files = readparamters(inputs[f]);
        //std::cout<<"f "<< f<< " " << files[0][0][0]<<std::endl;
        no_of_move_states.push_back(files[0][0][0]);
        probability.push_back(files[3][0]);
        mean_vector.push_back(files[1]);
        variance_vector.push_back(files[2]);
        transitions.push_back(files[4]);
        double temp = sqrt(files[5][0][0]/M_PI);
        CentreHome_r.push_back(temp);
    };
    
    //std::cout<< " no states "<<no_of_move_states[0]<<std::endl;
    
    MultipleIterations(NoSensors, NoInterations, savevalue,
                       cam_interval, SensorWidth, SensorRadius,
                       NoSteps,
                       SaveMovement,
                        CentreHome_r, MaximumDistance,
                       no_of_move_states, probability,
                        mean_vector, variance_vector,
                       transitions,
                      NoAnimals,
                        propMale
                       );
};

