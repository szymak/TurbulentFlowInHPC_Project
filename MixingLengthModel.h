#ifndef _MIXING_LENGTH_MODEL_H_
#define _MIXING_LENGTH_MODEL_H_

#include "Definitions.h"
#include "Parameters.h"
#include "TurbulentFlowField.h"

class MixingLengthModel {
	
	protected:
		const Parameters & _parameters;
	public:
		MixingLengthModel(const Parameters & parameters);
		virtual ~MixingLengthModel();
		// Return the value of lm at position (i, j, k)
		virtual FLOAT at(TurbulentFlowField & flowField, int i, int j) = 0;
		virtual FLOAT at(TurbulentFlowField & flowField, int i, int j, int k) = 0;
};

class LmKh : public MixingLengthModel {
	
	public:
		LmKh(const Parameters & parameters);
		virtual FLOAT at(TurbulentFlowField & flowField, int i, int j);
		virtual FLOAT at(TurbulentFlowField & flowField, int i, int j, int k);

};

class LmLaminarFlatPlate : public MixingLengthModel {
	
	public:
		LmLaminarFlatPlate(const Parameters & parameters);
		virtual FLOAT at(TurbulentFlowField & flowField, int i, int j);
		virtual FLOAT at(TurbulentFlowField & flowField, int i, int j, int k);
	
};

class LmTurbulentFlatPlate : public MixingLengthModel {
	
	public:
		LmTurbulentFlatPlate(const Parameters & parameters);
		virtual FLOAT at(TurbulentFlowField & flowField, int i, int j);
		virtual FLOAT at(TurbulentFlowField & flowField, int i, int j, int k);
	
};

#endif
