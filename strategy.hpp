//
//  strategy.hpp
//  signal
//
//  Created by Tom Lenaerts on 03/02/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#ifndef strategy_hpp
#define strategy_hpp

#include <cmath>
#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include "ctpgame.hpp"
#include "rangen.h"

using namespace std;

class Strategy {
public:
    Strategy():_level(0),_belief1(0), _belief2(0),_decision1(0), _decision2(0){};
    Strategy(unsigned level, unsigned belf1, unsigned belf2):_level(level),_belief1(belf1), _belief2(belf2),_decision1(belf1), _decision2(belf2){};
    Strategy(const Strategy& other);
    Strategy& operator=(const Strategy& other);

    unsigned level() const {return _level;}
    bool setLevel(unsigned value) {_level=value; return (_level == value);}
    unsigned beliefR1() const {return _belief1;}
    bool setBeliefR1(unsigned value) {_belief1=value; return (_belief1 == value);}
    unsigned beliefR2() const {return _belief2;}
    bool setBeliefR2(unsigned value) {_belief2=value; return (_belief2 == value);}
    unsigned decisionR1() const {return _decision1;}
    bool setDecisionR1(unsigned value) {_decision1=value; return (_decision1 == value);}
    unsigned decisionR2() const {return _decision2;}
    bool setDecisionR2(unsigned value) {_decision2=value; return (_decision2 == value);}

    unsigned belief(unsigned role) const;
    bool setBelief(unsigned role, unsigned value);
    unsigned decision(unsigned role) const;
    bool setDecision(unsigned role, unsigned value);

    bool inferDecision(CtpGame* game);
    bool stochasticInferDecision(CtpGame* game, double epsilon, RanGen* ran);

    bool operator==(const Strategy& other) const;
    bool operator!=(const Strategy& other) const;

    friend ostream & operator<<(ostream &o, Strategy& s)  {return s.display(o);}

protected:
    virtual ostream& display(ostream& os) const ;
    bool inferR1Decision(CtpGame* game, unsigned& R1, unsigned& R2);
    bool inferR2Decision(CtpGame* game, unsigned& R1, unsigned& R2);
    int addNoise(RanGen* ran, CtpGame* game, int base, double epsilon);
    bool firstOpportunity(unsigned& dec1, unsigned& dec2, unsigned length);

    //strategy is a baseline belief + a level of reasoning.
    unsigned _level;
    unsigned _belief1;
    unsigned _belief2;
    //after reasoning, whent take the money
    unsigned _decision1;
    unsigned _decision2;
};



class StrategySpace {
public:
    StrategySpace(){
        _cleared = true;
    }
    ~StrategySpace(){
        while(_space.size() > 0){
            Strategy* s = _space.back();
            _space.pop_back();
            delete s;
        }
    }
    bool createRandNowakStrategies(CtpGame& game,unsigned levels);
    bool createSubRandNowakStrategies(CtpGame& game,unsigned levels);
    bool createAltruistStrategies(CtpGame& game,unsigned levels);
    bool createAllSymmetricStrategies(CtpGame& game,unsigned levels);
    bool createNecessaryStrategies(CtpGame& game,unsigned levels);
    bool createRestrictedStrategies(CtpGame& game,unsigned levels);
    bool createSymmetricStrategies(CtpGame& game, unsigned beliefs, unsigned level);
    bool createSymmetricUTStrategies(CtpGame& game, unsigned beliefs, unsigned level);
    bool createSymmetricReducedStrategies(CtpGame& game, unsigned beliefs, unsigned level);
    bool createEquilibriumStrategies(CtpGame& game);

    bool clearStrategies(){
        while(_space.size() > 0){
            Strategy* s = _space.back();
            _space.pop_back();
            delete s;
        }
        _cleared=true;
        return _cleared;
    }

    unsigned int size() const {return (unsigned)_space.size();}
    Strategy* operator[](unsigned pos) const;
    int find(Strategy& s);
    int find(unsigned level, unsigned t1, unsigned t2);

    std::vector<Strategy*>::iterator begin() {return _space.begin();}
    std::vector<Strategy*>::iterator end() {return _space.end();}
    std::vector<Strategy*>::const_iterator begin() const{return _space.begin();}
    std::vector<Strategy*>::const_iterator end() const {return _space.end();}

    friend std::ostream & operator<<(std::ostream &o, StrategySpace& s){return s.display(o);}

protected:
    virtual std::ostream& display(std::ostream& os) const ;

    vector<Strategy*> _space;
    bool _cleared;
};


#endif /* strategy_hpp */
