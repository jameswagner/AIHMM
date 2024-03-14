#include "ExpressionInfo.h"

ExpressionInfo::ExpressionInfo() : snp(nullptr), ratioOfRatios(0.0f), expressionRatio(0.0f), heterozygote(false) {}

ExpressionInfo::ExpressionInfo(SNP* new_snp, bool new_heterozygote, float new_ratioOfRatios, float new_expressionRatio)
    : snp(new_snp), ratioOfRatios(new_ratioOfRatios), expressionRatio(new_expressionRatio), heterozygote(new_heterozygote) {

    // Ensure ratioOfRatios is within bounds
    if (ratioOfRatios < -5.0f) {
        ratioOfRatios = -5.0f;
    }
    else if (ratioOfRatios > 5.0f) {
        ratioOfRatios = 5.0f;
    }
}


ExpressionInfo::~ExpressionInfo() {
    delete snp;
}

bool ExpressionInfo::isHeterozygote() {
    return heterozygote;
}

void ExpressionInfo::setHeterozygote(bool isHeterozygote) {
    heterozygote = isHeterozygote;
}

SNP* ExpressionInfo::getSNP() {
    return snp;
}

float ExpressionInfo::getExpressionRatio() {
    return expressionRatio;
}

void ExpressionInfo::setExpressionRatio(float new_ratio) {
    expressionRatio = new_ratio;
}

float ExpressionInfo::getExpressionLogRatio() {
    if (expressionRatio <= -10.0f) {
        expressionRatio = -10.0f;
    }
    else if (expressionRatio >= 4.0f) {
        expressionRatio = 1.0f;
    }
    return expressionRatio;
}

float ExpressionInfo::getRatioOfRatios() {
    return ratioOfRatios;
}
