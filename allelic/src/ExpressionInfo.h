/*
 * ExpressionInfo.h
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#ifndef EXPRESSIONINFO_H_
#define EXPRESSIONINFO_H_
#include "SNP.h"
class ExpressionInfo {
public:
	ExpressionInfo();
	
	ExpressionInfo(SNP* new_snp, bool new_heterozygote, float new_ratioOfRatios, float new_expressionRatio);
	

	virtual ~ExpressionInfo();

	SNP* getSNP();
	bool isHeterozygote();
	void setHeterozygote(bool isHeterozygote);
	void setExpressionRatio(float new_ratio);
	
	float getExpressionRatio();
	float getRatioOfRatios();
	void setRatioOfRatios(float ratio);
	float getExpressionLogRatio();
private:
	SNP *snp;
	float ratioOfRatios;
	float expressionRatio;
	bool heterozygote;
};

#endif /* EXPRESSIONINFO_H_ */
