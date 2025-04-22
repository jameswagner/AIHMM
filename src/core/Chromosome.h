/*
 * Chromosome.h
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <vector>
#include "ExpressionInfo.h"

class Chromosome {
public:
    Chromosome();
    ~Chromosome();

    std::vector<ExpressionInfo*> getFullExpressions(bool useLow);
    void addExpression(ExpressionInfo* ei);
    void setExpression(const std::vector<ExpressionInfo*>& ei);

private:
    std::vector<ExpressionInfo*> expressions;
    std::vector<ExpressionInfo*> heterozygotes;
};

#endif /* CHROMOSOME_H_ */
