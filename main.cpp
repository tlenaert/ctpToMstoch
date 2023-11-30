//
//  main.cpp
//  signalling
//
//  Created by Tom Lenaerts on 23/01/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <chrono>
#include <gsl/gsl_matrix.h>

#include "kernel.hpp"
#include "rangen.h"
#include "ctpgame.hpp"
#include "strategy.hpp"
#include "ctpdata.hpp"

using namespace std;
using namespace std::chrono;
#define FIXED_FLOAT(x) std::fixed <<setprecision(9)<<(x)


bool runNetwork(unsigned approach, unsigned psize, double cost, unsigned levels, unsigned maxlevel, double betas, double mut, double eps, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, CtpData& data, ofstream& nodef, ofstream& edgef, ofstream& belf, ofstream& levf){
    
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;

    high_resolution_clock::time_point start = high_resolution_clock::now();
    Kernel run(psize, strategies.size(), game);
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost, approach);
    run.execute(&strategies, betas, mut, scaling);
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    duration<double> duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " s" << endl;
    run.showEvoRobustStrategies(&strategies,1.0);
    //belief transitions
    for(unsigned first=0; first <= game.length(); first++){
        for(unsigned second=0; second <= game.length(); second++){
            if(first != second){
                double trans = run.calcBeliefTransition(&strategies,5, first, second);
                belf << first << "\t" << second << "\t" << trans <<endl;
            }
        }
    }
    //level transitions
    for(unsigned first=0; first <= levels; first++){
        for(unsigned second=0; second <= levels; second++){
            if(first != second){
                double trans = run.calcLevelTransition(&strategies,5, first, second);
                levf << first << "\t" << second << "\t" << trans <<endl;
            }
        }
    }
    //nodes output
    run.printStationary(&strategies, 0.0, nodef); // include all SD
    run.printLevelStationary(&strategies, nodef); // include all SD
    run.printBeliefStationary(&strategies,nodef); // include all SD
    //edges output (fixation probabilities)
    run.printFixation(&strategies,1.0/double(psize),edgef); // scale links with 1/Z, i.e. neutral selection.
    return true;
}

void generateMixedStrategies(unsigned approach, unsigned psize, unsigned levels, unsigned maxlevel, double eps, unsigned repeats, CtpGame& game, RanGen& ran, ofstream& of){
    
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
//    strategies.createEquilibriumStrategies(game); // only (1,0) and (5,3) for L=4

    cout << strategies << endl;

    high_resolution_clock::time_point start = high_resolution_clock::now();
    Kernel run(psize, strategies.size(), game);
    run.generateMixedStrategies(&strategies, levels, eps, repeats, &ran, approach);
    
    for(unsigned i=0; i < strategies.size(); i++){
        Strategy* s = (strategies)[i];
        vector<double> mixed = s->getMixed();
        
        of << *s << "-" << "D" << "\t";
        for(unsigned j = 0; j < mixed.size(); j++){
            of << mixed[j];
            if (j < (mixed.size()-1))
                of << "\t";
        }
        of<< endl;
    }
    
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    duration<double> duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " s" << endl;
}

void generateSingleMixedStrategy(unsigned approach, unsigned psize, unsigned levels, unsigned maxlevel, double eps, unsigned repeats, CtpGame& game, RanGen& ran, ofstream& of){
    
    Strategy test(2,2,2);
    unsigned role=1;
    cout << test << endl;

    high_resolution_clock::time_point start = high_resolution_clock::now();
    Kernel run(psize, 1, game);
    run.generateSingleMixedStrategy(&test,1, levels, eps, repeats, &ran, approach);
    
    vector<double> mixed = test.getMixed();
    
    of << test << "\t";
    for(unsigned j = 0; j < mixed.size(); j++){
        of << mixed[j];
        if (j < (mixed.size()-1))
            of << "\t";
    }
    of<< endl;
    
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    duration<double> duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " s" << endl;
}


void runAnalytical(unsigned approach, unsigned psize, double cost, unsigned levels, unsigned maxlevel, double betas, double mut, double eps, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, CtpData& data, ofstream& of, ofstream& lf, ofstream& sf, ofstream& ef, ofstream& df, ofstream& levf, ofstream& mbf){
    
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;

    CtpEntry* elm1 = data["MKP4avg"]; //change to MKP6avg for L=6 ICG

    high_resolution_clock::time_point start = high_resolution_clock::now();
    Kernel run(psize, strategies.size(), game);
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost, approach);
    
    run.execute(&strategies, betas, mut, scaling);
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    duration<double> duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " s" << endl;
    
    run.levelKdistribution(&strategies, betas, levels+1, lf);
    run.stepsdistribution(&strategies, betas, eps, levels+1, repeats, approach, &ran, sf);
    run.decperkDistribution(&strategies, betas, levels+1, eps, repeats, approach, &ran, levf);
    run.misbeliefperkDistribution(&strategies, betas, levels+1, eps, repeats, approach, &ran, mbf);

    run.beliefDistribution(&strategies, betas, of);
    run.decisionDistribution(&strategies, betas, eps, repeats, approach, &ran, df);
    //compare to actual data
    double MSE1 = run.calcMSE(&strategies,elm1, eps, repeats, &ran, approach);
    ef << FIXED_FLOAT(betas) << "\t" << FIXED_FLOAT(eps) << "\t" << FIXED_FLOAT(MSE1) << "\t" << FIXED_FLOAT(sqrt(MSE1));
    ef << endl;
    
}


void runBeta(unsigned approach, unsigned psize, double cost, unsigned levels, unsigned maxlevel, double mut, double eps, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, CtpData& data, ofstream& of, ofstream& lf, ofstream& sf, ofstream& ef, ofstream& df, ofstream& levf){
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;
    double tmp=0.0001;
    double basebeta = tmp;
    double beta=basebeta;
    unsigned step=0;
    unsigned interval = 1;
    CtpEntry* elm = data["MKP6avg"];

    Kernel run(psize, strategies.size(), game);
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost, approach);

    while(beta<=10.0){
        cout << interval << " - " << step << " - " << FIXED_FLOAT(beta) << endl;
        high_resolution_clock::time_point start = high_resolution_clock::now();
        run.execute(&strategies, beta, mut, scaling);

        //don't forget to uncomment the correct inference method in the output functions
        run.levelKdistribution(&strategies, beta, levels+1, lf);
        run.stepsdistribution(&strategies, beta, eps, levels+1, repeats, approach, &ran, sf);
        run.decisionDistribution(&strategies, beta, eps, repeats, approach, &ran, df);
        run.beliefDistribution(&strategies, beta, of);
        run.showEvoRobustStrategies(&strategies, 1.0);
        //compare to actual data
        double MSE = run.calcMSE(&strategies,elm, eps, repeats, &ran, approach);
//        double avgp = run.averagePrediction(&strategies);


        ef << FIXED_FLOAT(beta) << "\t" << FIXED_FLOAT(eps) << "\t" << FIXED_FLOAT(MSE) << "\t" << FIXED_FLOAT(sqrt(MSE)) << endl; // << "\t" << avgp << endl;
        high_resolution_clock::time_point stop = high_resolution_clock::now();
        duration<double> duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by run function: "
             << fixed << duration.count() << " seconds" << endl;

        
        interval +=1;
        beta = (basebeta*interval);
        if(interval==10){
            interval=1;
            step +=1;
            basebeta = (tmp * pow(10, step));
        }
    }
}

void runEpsilon(unsigned approach, unsigned psize, double cost, unsigned levels, unsigned maxlevel, double mut, double beta, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, ofstream& of, ofstream& lf, ofstream& sf){
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;
    
    for(double eps=0.0 ; eps < 0.5  ; eps+=0.02){
        Kernel run(psize, strategies.size(), game);
        cout << FIXED_FLOAT(eps) << " - " << FIXED_FLOAT(beta) << endl;
        high_resolution_clock::time_point start = high_resolution_clock::now();
        run.calcPayoffs(&strategies, eps, repeats, &ran, cost, approach);
        run.execute(&strategies, beta, mut, scaling);
        high_resolution_clock::time_point stop = high_resolution_clock::now();
        duration<double> duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by run function: "
             << fixed << duration.count() << " seconds" << endl;
        run.levelKdistribution(&strategies, eps, levels+1, lf);
        run.beliefDistribution(&strategies, eps, of);
        run.stepsdistribution(&strategies, beta, eps, levels+1, repeats, approach, &ran, sf);
    }
}


void runBoth(unsigned approach, unsigned psize, double cost, unsigned levels, unsigned maxlevel, double mut, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, CtpData& data, ofstream& ef){
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;
    CtpEntry* elm1 = data["MKP4avg"];

    for(double eps=0.15; eps <= 0.23; eps+=0.01){
        Kernel run(psize, strategies.size(), game);
        run.calcPayoffs(&strategies, eps, repeats, &ran,cost, approach); //use the same payoff matrix for all beta

        for(double beta=0.9; beta <= 1.2; beta+=0.01){
            cout << "BETA=" << beta << "EPS=" << eps << endl;
            high_resolution_clock::time_point start = high_resolution_clock::now();
            run.execute(&strategies, beta, mut,scaling);
            high_resolution_clock::time_point stop = high_resolution_clock::now();
            duration<double> duration = duration_cast<microseconds>(stop - start);
            cout << "Time taken by function: "
                 << duration.count() << " s" << endl;

            //compare to actual data
            double MSE1 = run.calcMSE(&strategies,elm1, eps, repeats, &ran, approach);
            double avgb = run.averageBelief(&strategies);
            double avgl = run.averageLevel(&strategies);
            double avgf = run.averageFitness(&strategies);
            double avgp = run.averagePrediction(&strategies);
            double avgd = run.averageDecision(&strategies, repeats, eps, &ran,approach);

            ef << FIXED_FLOAT(beta) << "\t" << FIXED_FLOAT(eps) << "\t" << FIXED_FLOAT(MSE1) << "\t" << FIXED_FLOAT(sqrt(MSE1));
            ef << "\t" << avgb << "\t" << avgl << "\t" << avgd << "\t" << (avgb-avgd) << "\t" << avgf << "\t" << avgp << endl;
        }
    }
}


void runGradient(unsigned approach, unsigned psize, double beta, double eps, double mut, unsigned repeats, double cost, CtpGame& game, RanGen& ran, ofstream& gf){
    //redo implementation of run.gradient as it now assumes that each action is present
    StrategySpace strategies;
    strategies.createEquilibriumStrategies(game);
    Strategy neq(3,4,4);
    Kernel run(psize, strategies.size(), game);

    high_resolution_clock::time_point start = high_resolution_clock::now();
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost, approach);
    for(unsigned i = 1; i < psize; i+=1){
        double grad = run.gradientSimple(&strategies,&neq, i, beta, mut);
        gf << i << "\t" << grad << endl;
    }
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    duration<double> duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " s" << endl;
}



void runGradientEpsilon(unsigned approach, unsigned psize, double beta, double mut, unsigned repeats, double cost, CtpGame& game, RanGen& ran, ofstream& gf){
    //redo implementation of run.gradient as it now assumes that each action is present
    cout << "Beta =" << beta << endl;
    StrategySpace strategies;
    strategies.createEquilibriumStrategies(game);
//    strategies.createSymmetricUTStrategies(game, 4, 4);
    cout << strategies << endl;

    Strategy neq(3,4,4);
    vector<double> epsilons={0.05,0.06,0.07,0.08};//{0.0,0.07,0.14,0.21,0.28};
    unsigned num_lines = 4;
    gsl_matrix *results =  gsl_matrix_alloc(psize+1, num_lines);
    gsl_matrix_set_zero(results);
    unsigned count=0;
    for(double seps = 0; seps < num_lines; seps++){ // only look at the eps=0 and eps=0,18 case for now
        high_resolution_clock::time_point start = high_resolution_clock::now();
        Kernel run(psize, strategies.size(), game);
        run.calcPayoffs(&strategies, epsilons[seps], repeats, &ran,cost,approach);
        for(unsigned i = 1; i < psize; i+=1){
            double grad = run.gradient(&strategies,&neq, i, beta, mut); //C=true
            gsl_matrix_set(results, i, count, grad);
        }
        run.execute(&strategies, beta, mut, 1.0);
        cout << epsilons[seps] << "\t" << run << endl;

        high_resolution_clock::time_point stop = high_resolution_clock::now();
        duration<double> duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by function: "
             << duration.count() << " s" << endl;
        count +=1;
    }
    for(unsigned i=0; i <= psize; i++){
        gf << i << "\t";
        for(unsigned j=0; j < num_lines; j++){
            gf << gsl_matrix_get(results, i, j);
            if (j<(num_lines-1))
                gf << "\t";
        }
        gf << endl;
    }
}


int main(int argc, char * argv[]) {
    // relevant variables to set
    unsigned length = 4; // or 6
    double first = 0.4;
    double second = 0.1;
    double factor= 2.0;
    unsigned levels = 4;// or 6
    double epsilon=0.19; //0.25 for L=6
    unsigned repeats=50000;
    unsigned psize = 500;
    double betas=0.31; //0.49 for L=6
    double mut = 0.0;
    unsigned maxlevel=4; // or 6
    double scaling = 1.0;
    double cost = 0.0;
    unsigned approach=0; //0 = conditional with intertia, 1= only conditional, 2= unconditional

    //Random number generator
    RanGen ran;
    
    // outputfiles, activat based on called function
    string lfname("ctplevelsLall.txt");
    string dfname("ctpbeliefsLall.txt");
    string sfname("ctpstepsLall.txt");
    string efname("ctperrorLall.txt");
    string decfname("ctpdecisionLall.txt");
    string levfname("ctpldecperKLall.txt");
    string mbfname ("ctpmbperKLall.txt");
    string belfname("ctpltransitbeliefs.txt");
    string nodefname("ctpnodes.txt");
    string edgefname("ctpedges.txt");
    string gradfname("ctpgradient.txt");
    string sdfname("ctpsd.txt");
    string msfname("mixed.txt");

    //create output files
    ofstream df(dfname);
    ofstream lf(lfname);
    ofstream sf(sfname);
    ofstream errf(efname);

    ofstream decf(decfname);
    ofstream mbf(mbfname);
    ofstream levf(levfname);

    ofstream belf(belfname);

    ofstream nodef(nodefname);
    ofstream edgef(edgefname);

    ofstream gf(gradfname);
    ofstream sd(sdfname);
    
    ofstream ms(msfname);

    //create incremental centipede game
    CtpGame game(first, second, factor, length);
    
    //data from McKelvey, R. D., & Palfrey, T. R. (1992). An experimental study of the centipede game. Econometrica: Journal of the Econometric Society, 803-836.
    CtpData data;
    CtpEntry MKP4avg(4, {0.071, 0.356, 0.370, 0.153, 0.049});
    data.add("MKP4avg", &MKP4avg);
    CtpEntry MKP6avg(6, {0.007, 0.064, 0.199, 0.384, 0.253, 0.078, 0.014});
    data.add("MKP6avg", &MKP6avg);


    //Possible function (see above)
    
    //Producing the results for 1 specific beta and epsilon setting; see results in Figure 2 in the manuscript
    runAnalytical(approach, psize, cost, levels, maxlevel, betas, mut, epsilon, repeats, scaling, game, ran, data, df, lf, sf, errf, decf, levf,mbf);
    
    //Producing the results for a range of beta values with given epsilon (or vice versa), See Figure 1B and 1C
//    runBeta(approach, psize, cost, levels, maxlevel, mut, epsilon, repeats, scaling, game, ran, data, df, lf, sf, errf, decf,levf);
//    runEpsilon(approach, psize, cost, levels, maxlevel, mut, betas, repeats, scaling, game, ran, df, lf, sf);

    //Producing contour data for varying values of both epsilon and beta (see e.g. Figure 1A).
//    runBoth(approach, psize, cost, levels, maxlevel, mut, repeats, scaling, game, ran, data, errf);

    //Producing the markov chain for the population states, see Extended Data Figure 5
//    runNetwork(approach, psize, cost, levels, maxlevel, betas, mut, epsilon, repeats, scaling, game, ran, data, nodef, edgef, belf, levf);

    //Functions producing the gradient of selection for one or multiple epsilon values, see Figure 3.
//    runGradient(approach, psize, betas, epsilon, mut, repeats, cost, game, ran, gf);
//    runGradientEpsilon(approach, psize, betas, mut, repeats, cost, game, ran, gf);

    //other functions
//    generateSingleMixedStrategy(approach, psize, levels, maxlevel, epsilon, 100*repeats, game, ran, ms);
    errf.close();
    gf.close();
    lf.close();
    sf.close();
    df.close();
    decf.close();
    mbf.close();
    levf.close();
    belf.close();
    nodef.close();
    edgef.close();
    sd.close();
    ms.close();
    return 0;
}

