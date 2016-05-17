/*
 * Splitter.h
 *
 *  Created on: Feb 6, 2015
 *      Author: zhonghua
 */

#ifndef SPLITTER_H_
#define SPLITTER_H_

#include "model.h"
#include "config.h"

namespace masc {

class Splitter {
public:
    virtual void measure(model* m);
    virtual ~Splitter(){};
    virtual vector<float> assignWeights(model* m, const Config& config)=0;
    static Splitter* createSplitter(CutHeuristic heruistic);
protected:
    Splitter();
    float m_min_edge_length;
    float m_max_edge_length;
};

class FlatTreeSpliiter : public Splitter {
public:
    FlatTreeSpliiter(){};
    virtual ~FlatTreeSpliiter(){};
    virtual vector<float> assignWeights(model* m, const Config& config) override;
};

class MinimumPerimeterSpliiter : public Splitter {
public:
    MinimumPerimeterSpliiter(){}
    virtual ~MinimumPerimeterSpliiter(){}
    virtual vector<float> assignWeights(model* m, const Config& config) override;
};

class SteepestEdgeSplitter : public Splitter {
public:
    SteepestEdgeSplitter(){}
    virtual ~SteepestEdgeSplitter() {}
    virtual vector<float> assignWeights(model *m, const Config& config) override;
};

// assign random weights on edges
class RandomSplitter : public Splitter {
public:
    RandomSplitter(){}
    virtual ~RandomSplitter() {}
    virtual vector<float> assignWeights(model *m, const Config& config) override;
};

// try all possible random tree
class BruteForceSplitter : public Splitter {
public:
    BruteForceSplitter():m_inited(false) {}
    virtual ~BruteForceSplitter() {}
    virtual vector<float> assignWeights(model *m, const Config& config) override;
protected:
    void init(int edges);

    vector<float> m_weights;
    bool m_inited;
};


// 
// TODO #1: This is your first splitter
//
class MySplitter01 : public Splitter {
public:
	MySplitter01(){}
	virtual ~MySplitter01() {}
	virtual vector<float> assignWeights(model *m, const Config& config) override;
};


// 
// TODO #2: This is your second splitter
//
class MySplitter02 : public Splitter {
public:
	MySplitter02(){}
	virtual ~MySplitter02() {}
	virtual vector<float> assignWeights(model *m, const Config& config) override;
};


} /* namespace masc */

#endif /* SPLITTER_H_ */
