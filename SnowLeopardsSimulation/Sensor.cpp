//
//  Sensor.cpp
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 27/03/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#include "Sensor.h"
#include <random>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>

Sensor::Sensor(){};
Sensor::Sensor(double NoAnimal, double x, double y, double r, double direct, double halfwidth, int id){
    //creates vectors
    std::vector<int> temp(NoAnimal,0);
    animals_in_sensor_range =temp;
    animals_just_enter_range =temp;
    // sets values
    set_location_x(x);
    set_location_y(y);
    set_radius(r);
    set_angle_direction(direct);
    set_angle_halfwidth(halfwidth);
    SensorID = id;
    
};


void Sensor::calculate_sensoredges(){
    lhs_angle = RangeAngle(angle_direction - angle_halfwidth);
    rhs_angle = RangeAngle(angle_direction + angle_halfwidth);
    
    if(fabs(lhs_angle-M_PI)<0.0001 || fabs(lhs_angle)<0.0001){
        lhs_vert=1;lhs_horz=0;
        lhs_gradient = NAN; lhs_intercept = NAN;
    } else if(fabs(lhs_angle-M_PI_2)<0.0001 || fabs(lhs_angle-3*M_PI_2)<0.0001){
        lhs_vert=0;lhs_horz=1;
        lhs_gradient = NAN; lhs_intercept = NAN;
    } else {
        lhs_gradient = GradientFromAngle(lhs_angle);
        lhs_intercept = location_y - location_x*lhs_gradient;
        lhs_vert = NAN; lhs_horz = NAN;
    }
    
    if(fabs(rhs_angle-M_PI)<0.0001 || fabs(rhs_angle)<0.0001){
        rhs_vert=1;rhs_horz=0;
        rhs_gradient = NAN; rhs_intercept = NAN;
    } else if(fabs(rhs_angle-M_PI_2)<0.0001 || fabs(rhs_angle-3*M_PI_2)<0.0001){
        rhs_vert=0;rhs_horz=1;
        rhs_gradient = NAN; rhs_intercept = NAN;
    } else {
        rhs_gradient = GradientFromAngle(rhs_angle);
        rhs_intercept = location_y - location_x*rhs_gradient;
        rhs_vert = NAN; rhs_horz = NAN;
        //std::cout<< "rhs_vert should be NAN " << rhs_vert<<std::endl;
    }
    //std::cout<< "rhs_vert" << rhs_vert<<std::endl;
 };
double Sensor::GradientFromAngle(double angle){
    // The Gradient is calculated as = delta(y)/delta(x)
    //  tan(angle) = Opposite/Adjacent = delta (x)/delta(y)
    // therefore gradient = 1/tan(angle)
    double Gradient = 1/tan(angle);
    return(Gradient);
};
std::vector <double> Sensor::AngleAndAngleInteraction(double m1_Angle, double c1_Angle, double m2_Angle, double c2_Angle){
    
    double Coordinate1 = (c2_Angle - c1_Angle)/(m1_Angle-m2_Angle);
    double Coordinate2 = m1_Angle*Coordinate1 +c1_Angle;
    
    std::vector <double> Coord(2);
    Coord[0] = Coordinate1;
    Coord[1] = Coordinate2;
    
    return(Coord);
};
std::vector <double> Sensor::HorzAndCircInteraction(double Horz){

    std::vector <double> Coord(4);
    double temp = pow(radius,2) - pow((Horz - location_y),2);
    if(temp>-0.0001){
        if(temp<0){temp=0;};
        double Coordinate1 = sqrt(temp) +location_x;
        double Coordinate2 = - sqrt(temp) +location_x;
        Coord[0] = Coordinate1;
        Coord[1] = Horz;
        Coord[2] = Coordinate2;
        Coord[3] = Horz;
        
    };

    return(Coord);
};
std::vector <double> Sensor::VertAndCircInteraction(double Vert){
    std::vector <double> Coord(4);
    double temp = pow(radius,2) - pow((Vert - location_x),2);
    if(temp>-0.0001){
        if(temp<0){temp=0;};
        double Coordinate1 = sqrt(temp)+ location_y;
        double Coordinate2 = -sqrt(temp)+ location_y;
        Coord[0] = Vert;
        Coord[1] = Coordinate1;
        Coord[2] = Vert;
        Coord[3] = Coordinate2;
    }
    return(Coord);
};
double Sensor::VertAndAngleInteraction(double Vert, double m_Angle, double c_Angle){
    double YCoordinate = m_Angle*Vert +c_Angle;
    return(YCoordinate);
};
double Sensor::HorzAndAngleInteraction(double Horz, double m_Angle, double c_Angle){
    double XCoordinate = (Horz - c_Angle)/m_Angle;
    return(XCoordinate);
};
std::vector <double> Sensor::AngleAndCircInteraction(double m_Angle, double c_Angle){
    
    std::vector <double> Coord(4);
    double circ_term1 = pow(m_Angle,2)+1;
    double temp_A =  location_y - c_Angle;
    double circ_term2 = -(m_Angle*temp_A +location_x)*2;
    double circ_term3 = pow(temp_A,2)+pow(location_x,2) - pow(radius,2);
    double temp = pow(circ_term2,2)-(4*circ_term1*circ_term3);
    
    if(temp>-0.00001){
        if(temp<0){temp=0;} // This needs to be included for tangents
        
        double circ_solsqrt = sqrt(temp);
        double circ_xsol1 = (-circ_term2 - circ_solsqrt)/(2*circ_term1);
        double circ_xsol2 = (-circ_term2 + circ_solsqrt)/(2*circ_term1);
        double circ_ysol1 = m_Angle*circ_xsol1 + c_Angle;
        double circ_ysol2 = m_Angle*circ_xsol2 + c_Angle;
        // Vector to return
        Coord[0] = circ_xsol1;
        Coord[1] = circ_ysol1;
        Coord[2] = circ_xsol2;
        Coord[3] = circ_ysol2;
    }
    return(Coord);
};



void Sensor::SensorEdgeAndMovement(double location_x_animal, double location_y_animal,
                               double previous_x_animal, double previous_y_animal,
                               int Individual_ID,
                               double move_angle,
                               int itnumber,
                               double m_animal,double c_animal,
                               double m_detector, double c_detector,double g_detector, double vert, double horz,
                               double disttotal,std::ofstream &Captures, int stepnumber){
    
    double currentlocx =location_x_animal;
    double currentlocy =location_y_animal;
    double previouslocx =previous_x_animal;
    double previouslocy =previous_y_animal;
    
    double X1, Y1;
    std::vector<double> renameXandY(2);std::vector<double> renameTandA(2);
    double mindect, maxdect;
    if(currentlocx == previouslocx){//Animal has vertical movement
        if(horz==1){ //Horzontal detector
            X1 = previous_x_animal; Y1 = location_y;
            if(fabs(g_detector)<0.0001){mindect = location_y; maxdect = location_y+radius;} else{mindect = location_y-radius; maxdect = location_y;};
            if((Y1-mindect)<-0.0001 || (Y1-maxdect)>0.0001){X1 = NAN; Y1 = NAN;}
        } else if(vert==1 && location_x==previouslocx){ //Vertical detector
            X1 = previous_x_animal;
            //could occur anywhere on line so: min possible location and max possible location
            if(fabs(g_detector-M_PI)<0.0001){mindect = location_x-radius; maxdect = location_x;} else{mindect=location_x; maxdect=location_x+radius;};
            if(fabs(move_angle-M_PI)<0.0001){Y1 = std::min(previouslocy,maxdect);} else{Y1 = std::max(previouslocy,mindect);};
            if(fabs(g_detector)<0.0001){mindect = location_y; maxdect = location_y+radius;} else{mindect = location_y-radius; maxdect = location_y;};
            if((Y1-mindect)<-0.0001 || (Y1-maxdect)>0.0001){X1 = NAN; Y1 = NAN;}
        } else{ //Detector at angle
            X1 = currentlocx; Y1 = VertAndAngleInteraction(currentlocx, m_detector, c_detector);
            if(fabs(g_detector)<0.0001){mindect = location_y; maxdect = location_y+radius;} else{mindect = location_y-radius; maxdect = location_y;};
            if((Y1-mindect)<-0.0001 || (Y1-maxdect)>0.0001){X1 = NAN; Y1 = NAN;}
        };
    } else if(currentlocy == previouslocy){ //Animal has horizontal movement
        if(horz==1 && location_y== previouslocy){ //Horzontal detector (share the same y-coord)
            Y1 = previouslocy;
            // I need to idenfy the max. and min. x-value of the detector
            if(fabs(g_detector-M_PI/2)<0.0001){mindect=location_y; maxdect=location_y+radius;} else{mindect = location_y-radius; maxdect = location_y;};
            if(fabs(move_angle-M_PI/2)<0.0001){X1 = std::max(previouslocx,mindect);} else{X1 = std::min(previouslocx,maxdect);};
        } else if(vert==1){// Vertical detector
            X1 = location_x; Y1= previouslocy;
            if(fabs(g_detector)<0.0001){mindect = location_y; maxdect = location_y+radius;} else{mindect = location_y-radius; maxdect = location_y;};
                   // std::cout<< " mindect " << mindect << " maxdect "<< maxdect<< " Y1 " << Y1 <<" X1 "<< X1 << std::endl;
            if((Y1-mindect)<-0.0001 || (Y1-maxdect)>0.0001){X1 = NAN; Y1 = NAN;}
        } else {// Detector at angle
            X1 = HorzAndAngleInteraction(currentlocy, m_detector, c_detector); Y1 = currentlocy;
        };
        //std::cout<< " radius " <<radius<<" mindect " << mindect << " maxdect "<< maxdect<<" m_animal " <<m_animal<< " c_animal " <<c_animal<< " Y1 " << Y1 <<" X1 "<< X1 << std::endl;
    } else if(fabs(m_detector-m_animal)<0.0001 && fabs(c_detector-c_animal)<0.0001){// same intercept and angle
        if(g_detector<M_PI){mindect=location_x; maxdect=location_x+radius*sin(g_detector);} else{mindect = location_x+radius*sin(g_detector); maxdect = location_x;};
        // => to the left => is the right
        if(move_angle<M_PI){X1 = std::max(previouslocx,mindect);} else{X1 = std::min(previouslocx,maxdect);};
        
        // Calculation of y location
        if(g_detector<M_PI/2 ||g_detector>3*M_PI/2 ){
            mindect=location_y; maxdect=location_y+radius*cos(g_detector);
        } else{
            mindect = location_y+radius*cos(g_detector); maxdect = location_y;
        };
        // => up => down
        if(move_angle<M_PI/2 ||move_angle>3*M_PI/2 ){ Y1 = std::max(previouslocy,mindect);} else{Y1 = std::min(previouslocy,maxdect);};
    } else if(fabs(m_detector-m_animal)<0.0001 && c_detector!=c_animal){ // same angle BUT different intercept
        X1 = NAN;
        Y1 = NAN;
    }  else if(m_detector != m_animal){
        if(horz==1){ //Horzontal detector
            X1 = HorzAndAngleInteraction(location_y, m_animal, c_animal); Y1 = location_y;
            if(fabs(g_detector-M_PI/2)<0.0001){mindect=location_y; maxdect=location_y+radius;} else{mindect = location_y-radius; maxdect = location_y;};
            if((X1-mindect)<-0.0001 || (X1-maxdect)>0.0001){X1 = NAN; Y1 = NAN;}

        }  else if(vert==1){ //Vertical detector
            X1 = location_x; Y1 = VertAndAngleInteraction(location_x, m_animal, c_animal);
            
            if(fabs(g_detector-M_PI)<0.0001){mindect = location_x-radius; maxdect = location_x;} else{mindect=location_x; maxdect=location_x+radius;};
            if((Y1-mindect)<-0.0001 || (Y1-maxdect)>0.0001){X1 = NAN; Y1 = NAN;}
            //std::cout<<" mindect " << mindect << " maxdect "<< maxdect<<" m_animal " <<m_animal<< " c_animal " <<c_animal<< " Y1 " << Y1 <<"X1 "<< X1 << std::endl;
        }else{
            renameXandY =AngleAndAngleInteraction(m_detector, c_detector, m_animal, c_animal);
            X1 = renameXandY[0];Y1 = renameXandY[1];
        };
    } else{ //ALL other cases exit failure
        std::cout<<"Something is very wrong with Sensor::SensorAndMovement " << "gradient of detector: " <<m_detector <<"gradient of animal: "<< m_animal << std::endl;
        exit (EXIT_FAILURE);
    };
    
    //Calculate time and the angle of the capture
     //std::cout<<"X1 " <<X1<< " Y1 "  << Y1 <<" disttotal " << disttotal<<std::endl;
    renameTandA = TimeAndAngleCal(Y1, X1, previouslocy, previouslocx, disttotal);
    double time = renameTandA[0]; double angle = renameTandA[1];
    double distance = sqrt(pow(Y1-previouslocy,2)+pow(X1-previouslocx,2));
    
    if((time<=1 && fabs(angle-move_angle)<0.0001)||(time==0 && angle==0)){
        capture_current = 1;
        capture_x = X1; capture_y = Y1;
        capture_t = time; capture_a = angle;
        capture_d = distance;
        CapturesInsideCameraAngle(Individual_ID, move_angle, itnumber, Captures, stepnumber);
    } else{
        capture_current = 0;
        capture_x = NAN; capture_y = NAN;
        capture_t = NAN; capture_a = NAN;
    };
}; // END OF FUNCTION
void Sensor::SensorCircAndMovement(double location_x_animal, double location_y_animal, 
                                   double previous_x_animal, double previous_y_animal,
                                   int Individual_ID,
                                   double move_angle,
                                   int itnumber,
                                   double m_animal, double c_animal,
                                   double disttotal,
                                   std::ofstream &Captures, int stepnumber){
    
    std::vector<double> TandA(2);
    std::vector<double> XandY(4);
    double X,Y;
    
    
    if(location_x_animal == previous_x_animal){ //Animal has vertical movement
        XandY = VertAndCircInteraction(location_x_animal);
    } else if(location_y_animal == previous_y_animal){ //Animal has horizontal movement
        XandY = HorzAndCircInteraction(location_y_animal);
    } else{
        XandY = AngleAndCircInteraction(m_animal, c_animal);
    }//END OF ANGLE ANIMAL MOVEMENT

    for(int v=0; v<2; v++){
        X=XandY[v*2];Y=XandY[(v*2)+1];
        
        TandA = TimeAndAngleCal(Y, X, previous_y_animal, previous_x_animal, disttotal);//time and angle of the interscept
        double time = TandA[0]; double angle = TandA[1];
       //std::cout<<"Circ"<< " location_x_animal "<< previous_x_animal<< " X " <<X <<" Y " <<Y <<" time " <<time <<" angle " <<angle <<std::endl;
        if(((time<=1 || fabs(time-1)<0.0001) && fabs(angle-move_angle)<0.0001)||(time==0 && angle==0)){
            capture_current = 1;
            capture_x=X;
            capture_y=Y;
            capture_t=time;
            capture_a=angle;
            CapturesInsideCameraAngle(Individual_ID, move_angle, itnumber, Captures, stepnumber);
        }; // END OF IF LOOP
    }; // END OF FOR v loop
};

void Sensor::CapturesInsideCameraAngle(int Individual_ID,
                                double move_angle,
                                int itnumber,
                                std::ofstream &Captures,
                                int stepnumber
                                ){
    // If on the exact same spot as the Sensor assume it will be captured
    if(fabs(capture_x-location_x)<0.0001 && fabs(capture_y-location_y)<0.0001){
        capture_a_c2a =0; capture_a_a2c =0;
        UpdateCaptures(Individual_ID,itnumber, Captures, stepnumber);
    } else{
        
        double AngleFromSensor = AngleTwoPoints(location_x, capture_x, location_y, capture_y);
        capture_a_c2a = RangeAngle(AngleFromSensor- angle_direction); // the minus angle so from centre
        double movementcorrection;
        if(move_angle>M_PI){movementcorrection=move_angle-2*M_PI;}else{movementcorrection=move_angle;};
        capture_a_a2c = RangeAngle(AngleFromSensor+M_PI - movementcorrection);
        
        if(angle_halfwidth == M_PI){// If the Sensor is a circle
            UpdateCaptures(Individual_ID,itnumber, Captures, stepnumber);
        } else if((lhs_angle>rhs_angle && (AngleFromSensor <= rhs_angle || AngleFromSensor >= lhs_angle) )|| //if lhs is <0 (and recal as #>0)
                fabs(lhs_angle-AngleFromSensor) <0.0001  || fabs(rhs_angle-AngleFromSensor) <0.0001  || //approx equal to
                (AngleFromSensor >= lhs_angle && AngleFromSensor <= rhs_angle)){ //between both edges
            UpdateCaptures(Individual_ID,itnumber, Captures, stepnumber);
        } else{
            capture_a_c2a =NAN; capture_a_a2c =NAN;
        };
        //End of IF - "in width of detector"
    }; //END of else - not directly on same spot as Sensor
    
}; // End of function

///////////////////////
///////////////////////

void Sensor::CapturesIntersection(double location_x_animal, double location_y_animal,   // Current locations
                                  double previous_x_animal, double previous_y_animal,   // Previous locations
                                  int Individual_ID,                                    // Animal ID
                                  double move_angle,                                    // Call direction & width
                                  int itnumber,                                         // Iteration number
                                  std::ofstream &Captures,
                                  int stepnumber){

    
    double disttotal = DistTwoPoints(previous_x_animal,location_x_animal,previous_y_animal,location_y_animal);
    
    //If the animal movement was a line on a graph with a gradient and a intercept, Y=mX+c, then:
    //  - gadient would be, m=(change x/change y)
    //  - intercpet would be: y-mx=c (where y and x are known)
    double m_animal  = GradientFromAngle(move_angle);
    double c_animal  = location_y_animal-location_x_animal*m_animal;
    //std::cout<< "location_y_animal " << location_y_animal <<" location_x_animal " << location_x_animal <<" location_x " << location_x<<" location_y " << location_y<<std::endl;
    // Checks the beginning of the step
    capture_x=previous_x_animal;
    capture_y=previous_y_animal;
    capture_t=0;
    capture_a=0;
    double disttocam =DistTwoPoints(previous_x_animal,location_x,previous_y_animal,location_y);
    if(disttocam<=radius){CapturesInsideCameraAngle(Individual_ID, move_angle, itnumber, Captures,stepnumber);};
    
    
    if(angle_halfwidth < M_PI){ // If the Sensor width is 360˚ then the straight line edges are not important
        // Checks for crossing the boundaries for detector 1
        SensorEdgeAndMovement(location_x_animal, location_y_animal, previous_x_animal, previous_y_animal,
                          Individual_ID, move_angle,itnumber, m_animal,c_animal,
                          lhs_gradient, lhs_intercept, lhs_angle, lhs_vert, lhs_horz,
                          disttotal,Captures,stepnumber);
        // Checks for crossing the boundaries for detector 2
        SensorEdgeAndMovement(location_x_animal, location_y_animal, previous_x_animal, previous_y_animal,
                          Individual_ID, move_angle,itnumber, m_animal,c_animal,
                          rhs_gradient, rhs_intercept, rhs_angle, rhs_vert, rhs_horz,
                          disttotal,Captures,stepnumber);
    };
    // Checks for crossing the boundaries for circular edge of detector
    SensorCircAndMovement(location_x_animal,location_y_animal, previous_x_animal, previous_y_animal,
                          Individual_ID,move_angle,itnumber,m_animal, c_animal, disttotal,Captures,stepnumber);
};



double Sensor::DistTwoPoints(double X1, double X2, double Y1, double Y2){
    double distance = sqrt(pow(X1-X2,2) + pow(Y1-Y2,2));
    return distance;
};
double Sensor::AngleTwoPoints(double X1, double X2, double Y1, double Y2){
    // The angle from X1,Y1 to X2,Y2 is atan(d(y)/d(x))
    // Don't want always want theta, I want the baring from north
    double diffx = (X2-X1);
    double diffy = (Y2-Y1);
    // Gives the angle from north
    double theta = atan(diffx/diffy);// - M_PI_2;
    
    if(diffy<0){theta+= M_PI; }
    //Gets rid of negative values
    theta = RangeAngle(theta);
    
    return theta;
};
double Sensor::RangeAngle(double angle){
    while(angle<0){angle += 2*M_PI;}
    while(angle>=2*M_PI){angle -= 2*M_PI;};
    return angle;
};
std::vector <double> Sensor::TimeAndAngleCal(double Y, double X, double previous_y_animal, double previous_x_animal, double disttotal){
    
    std::vector <double> returnvalues(2);
    
    double distedge = DistTwoPoints(X,previous_x_animal,Y,previous_y_animal);
    double time = distedge/disttotal;
    //std::cout<< "distedge: " << distedge<< " disttotal: " << disttotal << " time: " << time <<std::endl;
    //Find the angle
    
    double AngleAnimalCap = AngleTwoPoints(previous_x_animal,X,previous_y_animal,Y);
    
    
    if((fabs(Y-previous_y_animal)<0.0001 && fabs(X-previous_x_animal)<0.0001)||fabs(AngleAnimalCap-2*M_PI)<0.0001 ){AngleAnimalCap=0;};
    if(time<0 && time>-0.0001){time=0;};
    
    returnvalues[0] = time;
    returnvalues[1] = AngleAnimalCap;
    return(returnvalues);
};
 
void Sensor::UpdateCaptures(int Animal_ID, int itnumber, std::ofstream &Captures, int stepnumber){
    capture_count +=1;
    //std::cout<<"Animal_ID " << Animal_ID <<" stepnumber " << stepnumber <<std::endl;
    animal_in_out_range(Animal_ID, capture_t, stepnumber); //runs algorithm to find if moving into/out of camera zone
    Captures << Animal_ID <<","
            << stepnumber <<","
            << SensorID <<","
            << itnumber << ","
            << capture_x <<","
            << capture_y <<","
            << capture_t <<","
            << capture_a_c2a <<","
            << capture_a_a2c <<","
            << capture_d <<","
            << animals_in_sensor_range[Animal_ID] <<","
            << animals_just_enter_range[Animal_ID] << "\n";
    
};
void Sensor::animal_in_out_range(int animal_id, double time, int stepnumber){
    if(time==0 && stepnumber==0){ // if starts sim in camera shot
        animals_in_sensor_range[animal_id]=1;
        animals_just_enter_range[animal_id]=1;
    } else if(time != 0 && time!=1) { // moves into or out of camera
        if(animals_in_sensor_range[animal_id]==1){ // already been captured => moving out
            animals_in_sensor_range[animal_id]=0;
            animals_just_enter_range[animal_id]=0;
        }else{ // => moving in
            animals_in_sensor_range[animal_id]=1;
            animals_just_enter_range[animal_id]=1;
        }
    } else{ // in range at start or end of movement => remaining in zone
        animals_in_sensor_range[animal_id]=1;
        animals_just_enter_range[animal_id]=0;
    };
};
