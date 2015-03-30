//
//  testAnimal.cpp
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 25/03/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#include "catch.hpp"

#include "Animal.h"


#include <stdio.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <iostream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <random>



TEST_CASE("Testing Animal functions") {
    Animal sc;
    
    SECTION("pick state - when probability is 1") {
        INFO("Using TestRN") // Only appears on a FAIL
        
        CHECK(sc.PickState(1, 2, {0,1}) == 1);
        CHECK(sc.PickState(1, 2, {1,0}) == 0);
        
        // Displays this variable on a FAIL
        int N=(1000000);
        std::vector<double> temp(N);
        for(int i = 0; i < N; i ++) {
            srand(i);
            double t_seed = (double) rand();
            temp[i] = sc.PickState(t_seed, 2, {0.5,0.5});
        };
        double sum = std::accumulate(temp.begin(),temp.end(),0);
        double mean =sum/N;
        
        CAPTURE(sum);
        CHECK(mean == Approx(0.5));
        
    }
    
    SECTION("matrix multiplication") {
        INFO("Using TestRN") // Only appears on a FAIL
        std::vector<double> returnedvector = sc.MatrixMultip({1,2},{1,3,2,4});
        REQUIRE(returnedvector[0] == 7);
        REQUIRE(returnedvector[1] == 10);
        
    }
    
    SECTION("Range angle "){
        INFO("Range angle")
        CHECK(sc.RangeAngle(3*M_PI)   == M_PI*1);
        CHECK(sc.RangeAngle(3.5*M_PI) == M_PI*1.5);
        CHECK(sc.RangeAngle(-.5*M_PI) == M_PI*1.5);
        CHECK(sc.RangeAngle(2*M_PI) == 0);
        CHECK(sc.RangeAngle(-1.5*M_PI) == M_PI*0.5);
    }
    
    SECTION("CalNext_x "){
        INFO("Calculating next x location")
        sc.set_Current_x(0);
        CHECK(sc.CalNext_X(1,0)==0);
        CHECK(sc.CalNext_X(sqrt(2),M_PI/4)==1);
        CHECK(sc.CalNext_X(1,M_PI/2)==1);
    }
    
    SECTION("CalNext_y "){
        INFO("Calculating next y location")
        sc.set_Current_y(0);
        CHECK(sc.CalNext_Y(1,0)==1);
        CHECK(sc.CalNext_Y(sqrt(2),M_PI/4)==Approx(1));
        CHECK(sc.CalNext_Y(1,M_PI/2)==Approx(0));
    }
    
    SECTION("Normal"){
        INFO("normal")
        // Displays this variable on a FAIL
        int N=(10000000);
        std::vector<double> temp(N);
        std::vector<double> temp_1sd(N);
        std::vector<double> temp_2sd(N);
        double t_seed;
        for(int i = 0; i < N; i ++) {
            t_seed =i;
            for(int j = 0; j < 10; j ++) {srand(t_seed); t_seed = (double) rand();}
            temp[i] = sc.Normal(t_seed, 0, 1);
            if(temp[i]>-1 && temp[i]<1){temp_1sd[i]=1;}else{temp_1sd[i]=0;};
            if(temp[i]>-2 && temp[i]<2){temp_2sd[i]=1;}else{temp_2sd[i]=0;};

        };
        double sum = 0; double sum_1sd =0; double sum_2sd =0;
        for(int i=0; i<N; i++){ sum += temp[i]; sum_1sd += temp_1sd[i];sum_2sd += temp_2sd[i];};
        double mean = sum/N; double mean_1sd = sum_1sd/N; double mean_2sd = sum_2sd/N;
        
       // double sum =  std::accumulate(temp.begin(),temp.end(),0); double mean = sum/N;
       // double sum_1sd =  std::accumulate(temp_1sd.begin(),temp_1sd.end(),0); double mean_1sd = sum_1sd/N;
       // double sum_2sd =  std::accumulate(temp_2sd.begin(),temp_2sd.end(),0); double mean_2sd = sum_2sd/N;
        CHECK(mean == Approx(0));
        CHECK(mean_1sd ==Approx(0.68268));
        CHECK(mean_2sd ==Approx(0.95449));
    }
    
    SECTION("Multivaraite Normal"){
        INFO("Multinormal")
        // Displays this variable on a FAIL
        sc.set_MaximumDistance(exp(20));
        int N=(10000000);
        std::vector<double> temp(2);
        std::vector<double> temp_x(N);
        std::vector<double> temp_y(N);
        double t_seed;
        for(int i = 0; i < N; i ++) {
            t_seed =i;
            for(int j = 0; j < 10; j ++) {srand(t_seed); t_seed = (double) rand();}
            temp = sc.MultivariteNorm (t_seed, {10, M_PI},{1,0,0,0.1});
            //std::cout<< i <<" temp: "<< temp[0] <<", " <<temp[1] <<std::endl;
            temp_x[i] = temp[0];
            temp_y[i] = temp[1];
            
        };
        
        double sum_x = 0; double sum_y = 0;
        for(int i=0; i<N; i++){ sum_x += temp_x[i]; sum_y += temp_y[i];};
        double mean_x = sum_x/N; double mean_y = sum_y/N;

        CHECK(mean_x == Approx(10));
        CHECK(mean_y == Approx(M_PI));

    }
    
    SECTION("New location from multivariate normal"){
        INFO("NEWLOCATION")
        double t_seed;
        int N=(100000);
        
        sc.set_MaximumDistance(exp(20));
        sc.set_Current_x(0);
        sc.set_Current_y(0);
        
        std::vector<double> temp(2);
        std::vector<double> temp_x(N);
        std::vector<double> temp_y(N);
        
        for(int i = 0; i < N; i ++) {
            t_seed =i;
            for(int j = 0; j < 10; j ++) {srand(t_seed); t_seed = (double) rand();}
            temp=sc.NewLocationFromMult(t_seed, {3,0}, {0.001,0,0,0.0001});
            //std::cout<< i<<"/" << N<<" temp: "<< temp[0] <<", " <<temp[1] <<std::endl;
            temp_x[i] = temp[0];
            temp_y[i] = temp[1];
        }
        
        double sum_x = 0; double sum_y = 0;
        for(int i=0; i<N; i++){ sum_x += temp_x[i]; sum_y += temp_y[i];};
        double mean_x = sum_x/N; double mean_y = sum_y/N;
       
        CHECK(mean_x == Approx(0));
        CHECK(mean_y == Approx(exp(3)));
    }

    SECTION("New location from multivariate normal"){
        INFO("NEWLOCATION_1")
        double t_seed;
        int N=(100000);
        
        sc.set_MaximumDistance(exp(20));
        sc.set_Current_x(10);
        sc.set_Current_y(10);
        
        std::vector<double> temp(2);
        std::vector<double> temp_x(N);
        std::vector<double> temp_y(N);
        
        for(int i = 0; i < N; i ++) {
            t_seed =i;
            for(int j = 0; j < 10; j ++) {srand(t_seed); t_seed = (double) rand();}
            temp=sc.NewLocationFromMult(t_seed, {3,M_PI}, {0.001,0,0,0.0001});
            //std::cout<< i<<"/" << N<<" temp: "<< temp[0] <<", " <<temp[1] <<std::endl;
            temp_x[i] = temp[0];
            temp_y[i] = temp[1];
        }
        
        double sum_x = 0; double sum_y = 0;
        for(int i=0; i<N; i++){ sum_x += temp_x[i]; sum_y += temp_y[i];};
        double mean_x = sum_x/N; double mean_y = sum_y/N;
        
        CHECK(mean_x == Approx(10));
        CHECK(mean_y == Approx(+10-exp(3)));
    }
    
    SECTION("New location from multivariate normal"){
        INFO("NEWLOCATION_2")
        double t_seed;
        int N=(100000);
        
        sc.set_MaximumDistance(exp(20));
        sc.set_Current_x(0);
        sc.set_Current_y(0);
        
        std::vector<double> temp(2);
        std::vector<double> temp_x(N);
        std::vector<double> temp_y(N);
        
        for(int i = 0; i < N; i ++) {
            t_seed =i;
            for(int j = 0; j < 10; j ++) {srand(t_seed); t_seed = (double) rand();}
            temp=sc.NewLocationFromMult(t_seed, {3,M_PI/2}, {0.001,0,0,0.0001});
            //std::cout<< i<<"/" << N<<" temp: "<< temp[0] <<", " <<temp[1] <<std::endl;
            temp_x[i] = temp[0];
            temp_y[i] = temp[1];
        }
        
        double sum_x = 0; double sum_y = 0;
        for(int i=0; i<N; i++){ sum_x += temp_x[i]; sum_y += temp_y[i];};
        double mean_x = sum_x/N; double mean_y = sum_y/N;
        
        CHECK(mean_x == Approx(exp(3)));
        CHECK(mean_y == Approx(0));
    }
    
    SECTION("New location from multivariate normal"){
        INFO("NEWLOCATION_3")
        double t_seed;
        int N=(100000);
        
        sc.set_MaximumDistance(exp(20));
        sc.set_Current_x(0);
        sc.set_Current_y(0);
        
        
        std::vector<double> temp(2);
        std::vector<double> temp_x(N);
        std::vector<double> temp_y(N);
        
        for(int i = 0; i < N; i ++) {
            t_seed =i;
            for(int j = 0; j < 10; j ++) {srand(t_seed); t_seed = (double) rand();}
            temp=sc.NewLocationFromMult(t_seed, {3,M_PI/4}, {0.001,0,0,0.0001});
            //std::cout<< i<<"/" << N<<" temp: "<< temp[0] <<", " <<temp[1] <<std::endl;
            temp_x[i] = temp[0];
            temp_y[i] = temp[1];
        }
        
        double sum_x = 0; double sum_y = 0;
        for(int i=0; i<N; i++){ sum_x += temp_x[i]; sum_y += temp_y[i];};
        double mean_x = sum_x/N; double mean_y = sum_y/N;
        
        CHECK(mean_x == Approx(sqrt(pow(exp(3),2)/2)));
        CHECK(mean_y == Approx(sqrt(pow(exp(3),2)/2)));
    }
    
    SECTION("New location"){
        INFO("NEWLOCATION_cal")
        double t_seed;
        int N=(200000);
        

        sc.set_CentreHome_x(0);sc.set_CentreHome_y(0);sc.set_CentreHome_r(10000000);
        sc.set_Current_x(0);sc.set_Current_y(0);sc.set_Current_state(0);
        sc.set_MaximumDistance(15000);
        sc.set_no_of_move_states(2);
        sc.set_probability({0.5,0.5});
        sc.set_mean_vector({{3,M_PI},{3,M_PI}});
        sc.set_variance_vector({{0.1,0,0,0.1},{0.1,0,0,0.1}});

        
        std::vector<double> temp_s(N);
        std::vector<double> temp_x(N);
        std::vector<double> temp_y(N);
        std::vector<double> temp_d(N);
        std::vector<double> temp_m(N);
        
        for(int i = 0; i < N; i ++) {
            t_seed =i;
            for(int j = 0; j < 10; j ++) {srand(t_seed); t_seed = (double) rand();}

            sc.NewLocation (t_seed);
            temp_s[i] = sc.get_Current_state();
            temp_x[i] = sc.get_Current_x();
            temp_y[i] = sc.get_Current_y();
            //std::cout<< i<<"/" << N<<" temp: "<< temp_x[i] <<", " <<temp_y[i] <<std::endl;
            if(i>2){temp_m[i] = sqrt(pow(temp_x[i]-temp_x[i-1],2)+pow(temp_y[i]-temp_y[i-1],2));}
            temp_d[i] =sqrt(pow(temp_x[i],2)+pow(temp_y[i],2));

        }
        
        double sum_s = 0; double sum_m = 0;
        for(int i=0; i<N; i++){ sum_s += temp_s[i]; sum_m += temp_m[i]; };
        double mean_s = sum_s/N;  double mean_m = sum_m/(N-1);
        double radius_check =0;
        for(int i=0; i<N; i++){ if(temp_d[i]>10000000){radius_check +=1;}};
        
        //CHECK(sc.get_Total_distance() == Approx(exp(3)*N));
        CHECK(mean_s == Approx(0.5));
        CHECK(radius_check == 0);
    }

    
}