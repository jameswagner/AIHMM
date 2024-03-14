#include "Chromosome.h"

Chromosome::Chromosome() {}

void Chromosome::setExpression(const std::vector<ExpressionInfo*>& ei) {
    expressions.clear();
    heterozygotes.clear();
    expressions = ei;
}

Chromosome::~Chromosome() {
    heterozygotes.clear();
    expressions.clear();
}

std::vector<ExpressionInfo*> Chromosome::getFullExpressions(bool useLow) {
    if (useLow) {
        return expressions;
    }
    else {
        // Uncomment if needed
        // return getHighOnlyExpressions();
        return expressions; // Return all expressions for now
    }
}

void Chromosome::addExpression(ExpressionInfo* ei) {
    expressions.push_back(ei);
}
