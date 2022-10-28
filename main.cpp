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


bool runNetwork(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double betas, double mut, double eps, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, CtpData& data, ofstream& nodef, ofstream& edgef, ofstream& belf, ofstream& levf){
    
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;

    high_resolution_clock::time_point start = high_resolution_clock::now();
    Kernel run(psize, strategies.size(), game);
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost);
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


void runAnalytical(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double betas, double mut, double eps, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, CtpData& data, ofstream& of, ofstream& lf, ofstream& sf, ofstream& ef, ofstream& df, ofstream& levf, ofstream& mbf){
    
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;

    CtpEntry* elm1 = data["MKP6avg"]; //change tpo MKP4avg for L=4 ICG

    high_resolution_clock::time_point start = high_resolution_clock::now();
    Kernel run(psize, strategies.size(), game);
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost);
    
    run.execute(&strategies, betas, mut, scaling);
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    duration<double> duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " s" << endl;
    
    run.levelKdistribution(&strategies, betas, levels+1, lf);
    run.stepsdistribution(&strategies, betas, eps, levels+1, repeats, &ran, sf);
    run.decperkDistribution(&strategies, betas, levels+1, eps, repeats, &ran, levf);
    run.misbeliefperkDistribution(&strategies, betas, levels+1, eps, repeats, &ran, mbf);

    run.beliefDistribution(&strategies, betas, of);
    run.decisionDistribution(&strategies, betas, eps, repeats, &ran, df);
    //compare to actual data
    double MSE1 = run.calcMSE(&strategies,elm1, eps, repeats, &ran);
    ef << FIXED_FLOAT(betas) << "\t" << FIXED_FLOAT(eps) << "\t" << FIXED_FLOAT(MSE1) << "\t" << FIXED_FLOAT(sqrt(MSE1));
    ef << endl;
    
}


void runBeta(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double mut, double eps, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, CtpData& data, ofstream& of, ofstream& lf, ofstream& sf, ofstream& ef, ofstream& df, ofstream& levf){
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;
    double tmp=0.0001;
    double basebeta = tmp;
    double beta=basebeta;
    unsigned step=0;
    unsigned interval = 1;
    CtpEntry* elm = data["MKP4avg"];

    Kernel run(psize, strategies.size(), game);
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost);

    while(beta<=10.0){
        cout << interval << " - " << step << " - " << FIXED_FLOAT(beta) << endl;
        high_resolution_clock::time_point start = high_resolution_clock::now();
        run.execute(&strategies, beta, mut, scaling);

        run.levelKdistribution(&strategies, beta, levels+1, lf);
        run.stepsdistribution(&strategies, beta, eps, levels+1, repeats, &ran, sf);
        run.decisionDistribution(&strategies, beta, eps, repeats, &ran, df);
        run.beliefDistribution(&strategies, beta, of);
        run.showEvoRobustStrategies(&strategies, 1.0);
        //compare to actual data
        double MSE = run.calcMSE(&strategies,elm, eps, repeats, &ran);

        ef << FIXED_FLOAT(beta) << "\t" << FIXED_FLOAT(eps) << "\t" << FIXED_FLOAT(MSE) << "\t" << FIXED_FLOAT(sqrt(MSE)) << endl;
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

void runEpsilon(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double mut, double beta, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, ofstream& of, ofstream& lf, ofstream& sf){
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;
    
    for(double eps=0.0 ; eps < 0.5  ; eps+=0.02){
        Kernel run(psize, strategies.size(), game);
        cout << FIXED_FLOAT(eps) << " - " << FIXED_FLOAT(beta) << endl;
        high_resolution_clock::time_point start = high_resolution_clock::now();
        run.calcPayoffs(&strategies, eps, repeats, &ran, cost);
        run.execute(&strategies, beta, mut, scaling);
        high_resolution_clock::time_point stop = high_resolution_clock::now();
        duration<double> duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by run function: "
             << fixed << duration.count() << " seconds" << endl;
        run.levelKdistribution(&strategies, eps, levels+1, lf);
        run.beliefDistribution(&strategies, eps, of);
        run.stepsdistribution(&strategies, beta, eps, levels+1, repeats, &ran, sf);
    }
}

void runBoth(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double mut, unsigned repeats, double scaling, CtpGame& game, RanGen& ran, CtpData& data, ofstream& ef){
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;
    CtpEntry* elm1 = data["MKP4avg"];

    for(double eps=0.0; eps <= 0.31; eps+=0.02){
        Kernel run(psize, strategies.size(), game);
        run.calcPayoffs(&strategies, eps, repeats, &ran,cost); //use the same payoff matrix for all beta

        for(double beta=0.0; beta <= 1.01; beta+=0.02){
            cout << "BETA=" << beta << "EPS=" << eps << endl;
            high_resolution_clock::time_point start = high_resolution_clock::now();
            run.execute(&strategies, beta, mut,scaling);
            high_resolution_clock::time_point stop = high_resolution_clock::now();
            duration<double> duration = duration_cast<microseconds>(stop - start);
            cout << "Time taken by function: "
                 << duration.count() << " s" << endl;

            //compare to actual data
            double MSE1 = run.calcMSE(&strategies,elm1, eps, repeats, &ran);
            double avgb = run.averageBelief(&strategies);
            double avgl = run.averageLevel(&strategies);
            double avgf = run.averageFitness(&strategies);
            double avgd = run.averageDecision(&strategies, repeats, eps, &ran);

            ef << FIXED_FLOAT(beta) << "\t" << FIXED_FLOAT(eps) << "\t" << FIXED_FLOAT(MSE1) << "\t" << FIXED_FLOAT(sqrt(MSE1));
            ef << "\t" << avgb << "\t" << avgl << "\t" << avgd << "\t" << (avgb-avgd) << "\t" << avgf << endl;
        }
    }
}


void runGradient(unsigned psize, double beta, double eps, double mut, unsigned repeats, double cost, CtpGame& game, RanGen& ran, ofstream& gf){
    StrategySpace strategies;
    strategies.createEquilibriumStrategies(game);
    Strategy neq(3,4,4);
    Kernel run(psize, strategies.size(), game);

    high_resolution_clock::time_point start = high_resolution_clock::now();
    run.calcPayoffs(&strategies, eps, repeats, &ran,cost); 
    for(unsigned i = 1; i < psize; i+=1){
        double grad = run.gradientSimple(&strategies,&neq, i, beta, mut);
        gf << i << "\t" << grad << endl;
    }
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    duration<double> duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " s" << endl;
}

void runGradientEpsilon(unsigned psize, double beta, double mut, unsigned repeats, double cost, CtpGame& game, RanGen& ran, ofstream& gf){
    cout << "Beta =" << beta << endl;
    StrategySpace strategies;
    strategies.createEquilibriumStrategies(game);
    Strategy neq(3,4,4);
    vector<double> epsilons={0.09, 0.10, 0.11, 0.12, 0.13};
    unsigned num_lines = 5;
    gsl_matrix *results =  gsl_matrix_alloc(psize+1, num_lines);
    gsl_matrix_set_zero(results);
    unsigned count=0;
    for(double seps = 0; seps < num_lines; seps++){ // only look at the eps=0 and eps=0,18 case for now
        high_resolution_clock::time_point start = high_resolution_clock::now();
        Kernel run(psize, strategies.size(), game);
        run.calcPayoffs(&strategies, epsilons[seps], repeats, &ran,cost);
        for(unsigned i = 1; i < psize; i+=1){
            double grad = run.gradient(&strategies,&neq, i, beta, mut); //C=true
            gsl_matrix_set(results, i, count, grad);
        }
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
    double epsilon=0.18;
    unsigned repeats=50000;
    unsigned psize = 500;
    double betas=0.3;
    double mut = 0.0;
    unsigned maxlevel=4;
    double scaling = 1.0;
    double cost = 0.0;

    //Random number generator
    RanGen ran;
    
    // outputfiles, activat based on called function
    string lfname("./ctplevelsLall.txt");
    string dfname("./ctpbeliefsLall.txt");
    string sfname("./ctpstepsLall.txt");
    string efname("./ctperrorLall.txt");
    string decfname("./ctpdecisionLall.txt");
    string levfname("./ctpldecperKLall.txt");
    string mbfname ("./ctpmbperKLall.txt");
    string belfname("./ctpltransitbeliefs.txt");
    string nodefname("./ctpnodes.txt");
    string edgefname("./ctpedges.txt");
    string gradfname("./ctpgradient.txt");
    
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
    runAnalytical(psize, cost, levels, maxlevel, betas, mut, epsilon, repeats, scaling, game, ran, data, df, lf, sf, errf, decf, levf,mbf);
    
    //Producing the results for a range of beta values with given epsilon (or vice versa), See Figure 1B and 1C
//    runBeta(psize, cost, levels, maxlevel, mut, epsilon, repeats, scaling, game, ran, data, df, lf, sf, errf, decf,levf);
//    runEpsilon(psize, cost, levels, maxlevel, mut, betas, repeats, scaling, game, ran, df, lf, sf);

    //Producing contour data for varying values of both epsilon and beta (see e.g. Figure 1A).
//    runBoth(psize, cost, levels, maxlevel, mut, repeats, scaling, game, ran, data, errf);

    //Producing the markov chain for the population states, see Extended Data Figure 5
//    runNetwork(psize, cost, levels, maxlevel, betas, mut, epsilon, repeats, scaling, game, ran, data, nodef, edgef, belf, levf);

    //Functions producing the gradient of selection for one or multiple epsilon values, see Figure 3.
//    runGradient(psize, betas, epsilon, mut, repeats, cost, game, ran, gf);
//    runGradientEpsilon(psize, betas, mut, repeats, cost, game, ran, gf);


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
    return 0;
}

