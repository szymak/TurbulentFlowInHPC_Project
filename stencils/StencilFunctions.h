#ifndef _DERIVATIVES_H_
#define _DERIVATIVES_H_

#include <math.h>
#include "../Definitions.h"
#include "../Parameters.h"
#include "../FlowField.h"
#include "../TurbulentFlowField.h"


// Load the local velocity cube with relevant velocities of the 2D plane
inline void loadLocalVelocity2D(FlowField & flowField, FLOAT * const localVelocity, int i, int j){
    for (int row = -1; row <= 1; row++ ){
        for ( int column = -1; column <= 1; column ++ ){
            const FLOAT * const point = flowField.getVelocity().getVector(i + column, j + row);
            localVelocity[39 + 9*row + 3*column]     = point[0]; // x-component
            localVelocity[39 + 9*row + 3*column + 1] = point[1]; // y-component
        }
    }
}

// Load the local velocity cube with surrounding velocities
inline void loadLocalVelocity3D(FlowField & flowField, FLOAT * const localVelocity, int i, int j, int k){
    for ( int layer = -1; layer <= 1; layer ++ ){
        for ( int row = -1; row <= 1; row++ ){
            for ( int column = -1; column <= 1; column ++ ){
                const FLOAT * const point = flowField.getVelocity().getVector(i + column, j + row, k + layer);
                localVelocity[39 + 27*layer + 9*row + 3*column    ] = point[0]; // x-component
                localVelocity[39 + 27*layer + 9*row + 3*column + 1] = point[1]; // y-component
                localVelocity[39 + 27*layer + 9*row + 3*column + 2] = point[2]; // z-component
            }
        }
    }
}


// load local meshsize for 2D -> same as loadLocalVelocity2D, but invoking call to meshsize-ptr
inline void loadLocalMeshsize2D(const Parameters& parameters, FLOAT * const localMeshsize, int i, int j){
    for (int row = -1; row <= 1; row++ ){
        for ( int column = -1; column <= 1; column ++ ){
            localMeshsize[39 + 9*row + 3*column]     = parameters.meshsize->getDx(i+column,j+row);
            localMeshsize[39 + 9*row + 3*column + 1] = parameters.meshsize->getDy(i+column,j+row);
        }
    }
}

// load local meshsize for 3D
inline void loadLocalMeshsize3D(const Parameters& parameters, FLOAT * const localMeshsize, int i, int j, int k){
    for ( int layer = -1; layer <= 1; layer ++ ){
        for ( int row = -1; row <= 1; row++ ){
            for ( int column = -1; column <= 1; column ++ ){
                localMeshsize[39 + 27*layer + 9*row + 3*column    ] = parameters.meshsize->getDx(i+column,j+row,k+layer);
                localMeshsize[39 + 27*layer + 9*row + 3*column + 1] = parameters.meshsize->getDy(i+column,j+row,k+layer);
                localMeshsize[39 + 27*layer + 9*row + 3*column + 2] = parameters.meshsize->getDz(i+column,j+row,k+layer);
            }
        }
    }
}

// load local Turbulent Viscosity for 2D
inline void loadLocalTurbulentViscosity2D(TurbulentFlowField & flowField, FLOAT * const localTurbulentViscosity, int i, int j){
    for (int row = -1; row <= 1; row++ ){
        for ( int column = -1; column <= 1; column ++ ){
            localTurbulentViscosity[39 + 9*row + 3*column] = flowField.getTurbulentViscosity().getScalar(i + column, j + row);
        }
    }
}

// load local Turbulent Viscosity for 3D
inline void loadLocalTurbulentViscosity3D(TurbulentFlowField & flowField, FLOAT * const localTurbulentViscosity, int i, int j, int k){
    for ( int layer = -1; layer <= 1; layer ++ ){
        for ( int row = -1; row <= 1; row++ ){
            for ( int column = -1; column <= 1; column ++ ){
                localTurbulentViscosity[39 + 27*layer + 9*row + 3*column] = flowField.getTurbulentViscosity().getScalar(i + column, j + row, k + layer);
            }
        }
    }
}


// Maps an index and a component to the corresponding value in the cube.
inline int mapd (int i, int j, int k, int component){
   return 39 + 27*k + 9*j + 3*i + component;
}

// Derivative functions. They are applied to a cube of 3x3x3 cells. lv stands for the local velocity, lm represents the local mesh sizes
// dudx <-> first derivative of u-component of velocity field w.r.t. x-direction
inline FLOAT dudx ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv [mapd(0,0,0,0)] - lv [mapd(-1,0,0,0)] ) / GeometricParameters::dx;

    // evaluate dudx in the cell center by a central difference
    const int index0 = mapd(0,0,0,0);
    const int index1 = mapd(-1,0,0,0);
    return  ( lv [index0] - lv [index1] ) / lm[index0];
    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "dudx");}

    return tmp2;*/
}

inline FLOAT dudy ( const FLOAT * const lv, const FLOAT * const lm ) {

    // evaluate dudy in the cell center by interpolation
    
    /*
    	0 1 dyP
    	2 3 dy
    	4 5 dyM
    */
    
	//cells vel up
    FLOAT u0 = lv[mapd(-1,1,0,0)];
    FLOAT u1 = lv[mapd(0,1,0,0)];
    
	//cells center
    FLOAT u2 = lv[mapd(-1,0,0,0)];
    FLOAT u3 = lv[mapd(0,0,0,0)];
    
	//cells vel down
    FLOAT u4 = lv[mapd(-1,-1,0,0)];
    FLOAT u5 = lv[mapd(0,-1,0,0)];

	// cells length y
    FLOAT dy = lm[mapd(0,0,0,1)];
    FLOAT dyM = lm[mapd(0,-1,0,1)];
    FLOAT dyP = lm[mapd(0,1,0,1)];
    
	FLOAT veldown = ((u0+u1)*dy + (u2+u3)*dyP) / (2*(dy+dyP));
	FLOAT velup = ((u4+u5)*dy + (u2+u3)*dyM) / (2*(dy+dyM));

	return (velup-veldown) / dy;
	
}

inline FLOAT dudz ( const FLOAT * const lv, const FLOAT * const lm ) {

	//cells vel front
    FLOAT u0 = lv[mapd(-1,0,1,0)];
    FLOAT u1 = lv[mapd(0,0,1,0)];
    
	//cells center
    FLOAT u2 = lv[mapd(-1,0,0,0)];
    FLOAT u3 = lv[mapd(0,0,0,0)];
    
	//cells vel back
    FLOAT u4 = lv[mapd(-1,0,-1,0)];
    FLOAT u5 = lv[mapd(0,0,-1,0)];

	// cells length z
    FLOAT dz = lm[mapd(0,0,0,2)];
    FLOAT dzM = lm[mapd(0,0,-1,2)];
    FLOAT dzP = lm[mapd(0,0,1,2)];
    
	FLOAT velback = ((u0+u1)*dz + (u2+u3)*dzP) / (2*(dz+dzP));
	FLOAT velfront = ((u4+u5)*dz + (u2+u3)*dzM) / (2*(dz+dzM));

	return (velfront-velback) / dz;

}


inline FLOAT dvdy ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv [mapd(0,0,0,1)] - lv [mapd(0,-1,0,1)] ) / GeometricParameters::dy;
    const int index0 = mapd(0, 0,0,1);
    const int index1 = mapd(0,-1,0,1);
    return ( lv [index0] - lv [index1] ) / lm[index0];

    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "dvdy");}

    return tmp2;*/
}

inline FLOAT dvdx ( const FLOAT * const lv, const FLOAT * const lm ) {

	//cells vel right
    FLOAT v0 = lv[mapd(1,-1,0,1)];
    FLOAT v1 = lv[mapd(1,0,0,1)];
    
	//cells center
    FLOAT v2 = lv[mapd(0,-1,0,1)];
    FLOAT v3 = lv[mapd(0,0,0,1)];
    
	//cells vel left
    FLOAT v4 = lv[mapd(-1,-1,0,1)];
    FLOAT v5 = lv[mapd(-1,0,0,1)];

	// cells length x
    FLOAT dx = lm[mapd(0,0,0,0)];
    FLOAT dxM = lm[mapd(-1,0,0,0)];
    FLOAT dxP = lm[mapd(1,0,0,0)];
    
	FLOAT velleft = ((v0+v1)*dx + (v2+v3)*dxP) / (2*(dx+dxP));
	FLOAT velright = ((v4+v5)*dx + (v2+v3)*dxM) / (2*(dx+dxM));

	return (velright-velleft) / dx;
	
}

inline FLOAT dvdz ( const FLOAT * const lv, const FLOAT * const lm ) {

	//cells vel front
    FLOAT v0 = lv[mapd(0,-1,1,1)];
    FLOAT v1 = lv[mapd(0,0,1,1)];
    
	//cells center
    FLOAT v2 = lv[mapd(0,-1,0,1)];
    FLOAT v3 = lv[mapd(0,0,0,1)];
    
	//cells vel back
    FLOAT v4 = lv[mapd(0,-1,-1,1)];
    FLOAT v5 = lv[mapd(0,0,-1,1)];

	// cells length z
    FLOAT dx = lm[mapd(0,0,0,0)];
    FLOAT dxM = lm[mapd(-1,0,0,0)];
    FLOAT dxP = lm[mapd(1,0,0,0)];
    
	FLOAT velback = ((v0+v1)*dx + (v2+v3)*dxP) / (2*(dx+dxP));
	FLOAT velfront = ((v4+v5)*dx + (v2+v3)*dxM) / (2*(dx+dxM));

	return (velfront-velback) / dx;
	
}

inline FLOAT dwdz ( const FLOAT * const lv, const FLOAT * const lm ) {
    const int index0 = mapd(0,0, 0,2);
    const int index1 = mapd(0,0,-1,2);
    return ( lv [index0] - lv [index1] ) / lm[index0];
}

inline FLOAT dwdx ( const FLOAT * const lv, const FLOAT * const lm ) {

	//cells vel right
    FLOAT w0 = lv[mapd(1,0,-1,2)];
    FLOAT w1 = lv[mapd(1,0,0,2)];
    
	//cells center
    FLOAT w2 = lv[mapd(0,0,-1,2)];
    FLOAT w3 = lv[mapd(0,0,0,2)];
    
	//cells vel left
    FLOAT w4 = lv[mapd(-1,0,-1,2)];
    FLOAT w5 = lv[mapd(-1,0,0,2)];

	// cells length x
    FLOAT dx = lm[mapd(0,0,0,0)];
    FLOAT dxM = lm[mapd(-1,0,0,0)];
    FLOAT dxP = lm[mapd(1,0,0,0)];
    
	FLOAT velleft = ((w0+w1)*dx + (w2+w3)*dxP) / (2*(dx+dxP));
	FLOAT velright = ((w4+w5)*dx + (w2+w3)*dxM) / (2*(dx+dxM));

	return (velright-velleft) / dx;
	
}

inline FLOAT dwdy ( const FLOAT * const lv, const FLOAT * const lm ) {

	//cells vel up
    FLOAT w0 = lv[mapd(0,1,-1,2)];
    FLOAT w1 = lv[mapd(0,1,0,2)];
    
	//cells center
    FLOAT w2 = lv[mapd(0,0,-1,2)];
    FLOAT w3 = lv[mapd(0,0,0,2)];
    
	//cells vel down
    FLOAT w4 = lv[mapd(0,-1,-1,2)];
    FLOAT w5 = lv[mapd(0,-1,0,2)];

	// cells length y
    FLOAT dy = lm[mapd(0,0,0,1)];
    FLOAT dyM = lm[mapd(0,-1,0,1)];
    FLOAT dyP = lm[mapd(0,1,0,1)];
    
	FLOAT veldown = ((w0+w1)*dy + (w2+w3)*dyP) / (2*(dy+dyP));
	FLOAT velup = ((w4+w5)*dy + (w2+w3)*dyM) / (2*(dy+dyM));

	return (velup-veldown) / dy;
	
}
// second derivative of u-component w.r.t. x-direction, evaluated at the location of the u-component
inline FLOAT d2udx2 ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv[mapd(1,0,0,0)] - 2*lv[mapd(0,0,0,0)] + lv[mapd(-1,0,0,0)] )
    //    / ( GeometricParameters::dx * GeometricParameters::dx );

    // evaluate the second derivative at the location of the u-component of the velocity field;
    // we therefore use the two neighbouring u-components and assume arbitrary mesh sizes in both
    // directions -> the formula arises from a straight-forward taylor expansion.
    // -> for equal meshsizes, we obtain the usual [1 -2 1]-like stencil
    const int index_M1    = mapd(-1,0,0,0);
    const int index_0     = mapd(0,0,0,0);
    const int index_P1    = mapd(1,0,0,0);

    const FLOAT dx0   = lm[index_0];
    const FLOAT dx1   = lm[index_P1];
    const FLOAT dxSum = dx0+dx1;
    return 2.0*(lv[index_P1]/(dx1*dxSum) - lv[index_0]/(dx1*dx0) + lv[index_M1]/(dx0*dxSum) );

    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2udx2");}

    return tmp2;*/
}

inline FLOAT d2udy2 ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1=( lv[mapd(0,1,0,0)] - 2*lv[mapd(0,0,0,0)] + lv[mapd(0,-1,0,0)] )
    //    / ( GeometricParameters::dy * GeometricParameters::dy );
    // average mesh sizes, since the component u is located in the middle of the cell's face
    const FLOAT dy_M1 = lm[mapd(0,-1,0,1)];
    const FLOAT dy_0  = lm[mapd(0, 0,0,1)];
    const FLOAT dy_P1 = lm[mapd(0, 1,0,1)];
    const FLOAT dy0 = 0.5*(dy_0+dy_M1);
    const FLOAT dy1 = 0.5*(dy_0+dy_P1);
    const FLOAT dySum = dy0+dy1;
    return 2.0*(lv[mapd(0,1,0,0)]/(dy1*dySum) - lv[mapd(0,0,0,0)]/(dy1*dy0) + lv[mapd(0,-1,0,0)]/(dy0*dySum) );

    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2udy2");}

    return tmp2;*/
}

inline FLOAT d2udz2 ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv[mapd(0,0,1,0)] - 2*lv[mapd(0,0,0,0)] + lv[mapd(0,0,-1,0)] )
    //    / ( GeometricParameters::dz * GeometricParameters::dz );
    const FLOAT dz_M1 = lm[mapd(0, 0,-1,2)];
    const FLOAT dz_0  = lm[mapd(0, 0, 0,2)];
    const FLOAT dz_P1 = lm[mapd(0, 0 ,1,2)];
    const FLOAT dz0 = 0.5*(dz_0+dz_M1);
    const FLOAT dz1 = 0.5*(dz_0+dz_P1);
    const FLOAT dzSum = dz0+dz1;
    return 2.0*(lv[mapd(0,0,1,0)]/(dz1*dzSum) - lv[mapd(0,0,0,0)]/(dz1*dz0) + lv[mapd(0,0,-1,0)]/(dz0*dzSum) );
    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2udz2");}

    return tmp2;*/
}


/** second derivative of the v-component, evaluated at the location of the v-component */
inline FLOAT d2vdx2 ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv[mapd(1,0,0,1)] - 2*lv[mapd(0,0,0,1)] + lv[mapd(-1,0,0,1)] )
    //    / ( GeometricParameters::dx * GeometricParameters::dx );
    const FLOAT dx_M1 = lm[mapd(-1,0,0,0)];
    const FLOAT dx_0  = lm[mapd( 0,0,0,0)];
    const FLOAT dx_P1 = lm[mapd( 1,0,0,0)];
    const FLOAT dx0 = 0.5*(dx_0+dx_M1);
    const FLOAT dx1 = 0.5*(dx_0+dx_P1);
    const FLOAT dxSum = dx0+dx1;
    return 2.0*(lv[mapd(1,0,0,1)]/(dx1*dxSum) - lv[mapd(0,0,0,1)]/(dx1*dx0) + lv[mapd(-1,0,0,1)]/(dx0*dxSum) );

    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2vdx2");}

    return tmp2;*/
}

inline FLOAT d2vdy2 ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv[mapd(0,1,0,1)] - 2*lv[mapd(0,0,0,1)] + lv[mapd(0,-1,0,1)] )
    //    / ( GeometricParameters::dy * GeometricParameters::dy );
    const int index_M1    = mapd(0,-1,0,1);
    const int index_0     = mapd(0, 0,0,1);
    const int index_P1    = mapd(0, 1,0,1);

    const FLOAT dy0   = lm[index_0];
    const FLOAT dy1   = lm[index_P1];
    const FLOAT dySum = dy0+dy1;
    return 2.0*(lv[index_P1]/(dy1*dySum) - lv[index_0]/(dy1*dy0) + lv[index_M1]/(dy0*dySum) );

    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2vdy2");}

    return tmp2;*/
}

inline FLOAT d2vdz2 ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv[mapd(0,0,1,1)] - 2*lv[mapd(0,0,0,1)] + lv[mapd(0,0,-1,1)] )
    //    / ( GeometricParameters::dz * GeometricParameters::dz );
    const FLOAT dz_M1 = lm[mapd(0,0,-1,2)];
    const FLOAT dz_0  = lm[mapd(0,0, 0,2)];
    const FLOAT dz_P1 = lm[mapd(0,0, 1,2)];
    const FLOAT dz0 = 0.5*(dz_0+dz_M1);
    const FLOAT dz1 = 0.5*(dz_0+dz_P1);
    const FLOAT dzSum = dz0+dz1;
    return 2.0*(lv[mapd(0,0,1,1)]/(dz1*dzSum) - lv[mapd(0,0,0,1)]/(dz1*dz0) + lv[mapd(0,0,-1,1)]/(dz0*dzSum) );

    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2vdz2");}

    return tmp2;*/
}


/** second derivative of the w-component, evaluated at the location of the w-component */
inline FLOAT d2wdx2 ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv[mapd(1,0,0,2)] - 2*lv[mapd(0,0,0,2)] + lv[mapd(-1,0,0,2)] )
    //    / ( GeometricParameters::dx * GeometricParameters::dx );
    const FLOAT dx_M1 = lm[mapd(-1,0,0,0)];
    const FLOAT dx_0  = lm[mapd( 0,0,0,0)];
    const FLOAT dx_P1 = lm[mapd( 1,0,0,0)];
    const FLOAT dx0 = 0.5*(dx_0+dx_M1);
    const FLOAT dx1 = 0.5*(dx_0+dx_P1);
    const FLOAT dxSum = dx0+dx1;
    return 2.0*(lv[mapd(1,0,0,2)]/(dx1*dxSum) - lv[mapd(0,0,0,2)]/(dx1*dx0) + lv[mapd(-1,0,0,2)]/(dx0*dxSum) );

    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2wdx2");}

    return tmp2;*/
}

inline FLOAT d2wdy2 ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv[mapd(0,1,0,2)] - 2*lv[mapd(0,0,0,2)] + lv[mapd(0,-1,0,2)] )
    //    / ( GeometricParameters::dy * GeometricParameters::dy );
    const FLOAT dy_M1 = lm[mapd(0,-1,0,1)];
    const FLOAT dy_0  = lm[mapd(0, 0,0,1)];
    const FLOAT dy_P1 = lm[mapd(0, 1,0,1)];
    const FLOAT dy0 = 0.5*(dy_0+dy_M1);
    const FLOAT dy1 = 0.5*(dy_0+dy_P1);
    const FLOAT dySum = dy0+dy1;
    return 2.0*(lv[mapd(0,1,0,2)]/(dy1*dySum) - lv[mapd(0,0,0,2)]/(dy1*dy0) + lv[mapd(0,-1,0,2)]/(dy0*dySum) );

    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2wdy2");}

    return tmp2;*/
}

inline FLOAT d2wdz2 ( const FLOAT * const lv, const FLOAT * const lm ) {
    //double tmp1= ( lv[mapd(0,0,1,2)] - 2*lv[mapd(0,0,0,2)] + lv[mapd(0,0,-1,2)] )
    //    / ( GeometricParameters::dz * GeometricParameters::dz );
    const int index_M1    = mapd(0,0,-1,2);
    const int index_0     = mapd(0,0, 0,2);
    const int index_P1    = mapd(0,0, 1,2);

    const FLOAT dz0   = lm[index_0];
    const FLOAT dz1   = lm[index_P1];
    const FLOAT dzSum = dz0+dz1;
    return 2.0*(lv[index_P1]/(dz1*dzSum) - lv[index_0]/(dz1*dz0) + lv[index_M1]/(dz0*dzSum) );

    /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2wdz2");}

    return tmp2;*/
}


/** first-derivative of product (u*v), evaluated at the location of the v-component */
inline FLOAT duvdx ( const FLOAT * const lv, const Parameters & parameters, const FLOAT * const lm ) {
/*
    const FLOAT tmp1= 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,0)] + lv [mapd(0,1,0,0)] ) *
                         ( lv [mapd(0,0,0,1)] + lv [mapd(1,0,0,1)] ) ) -
                       ( ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,1,0,0)] ) *
                         ( lv [mapd(-1,0,0,1)] + lv [mapd(0,0,0,1)] ) ) )
      + parameters.solver.gamma *( ( fabs ( lv [mapd(0,0,0,0)] + lv [mapd(0,1,0,0)] ) *
                              ( lv [mapd(0,0,0,1)] - lv [mapd(1,0,0,1)] ) ) -
                       ( fabs ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,1,0,0)] ) *
                              ( lv [mapd(-1,0,0,1)] - lv [mapd(0,0,0,1)] ) ) )
                       ) / lm[mapd(0,0,0,0)];
*/

    const FLOAT hxShort = 0.5*lm[mapd( 0,0,0,0)];                       // distance of corner points in x-direction from center v-value
    const FLOAT hxLong0 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd(-1,0,0,0)]); // distance between center and west v-value
    const FLOAT hxLong1 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd( 1,0,0,0)]); // distance between center and east v-value
    const FLOAT hyShort = 0.5*lm[mapd(0,0,0,1)];                        // distance of center u-value from upper edge of cell
    const FLOAT hyLong  = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd(0,1,0,1)]);  // distance of north and center u-value

    const FLOAT u00  = lv[mapd( 0, 0, 0, 0)];
    const FLOAT u01  = lv[mapd( 0, 1, 0, 0)];
    const FLOAT v00  = lv[mapd( 0, 0, 0, 1)];
    const FLOAT v10  = lv[mapd( 1, 0, 0, 1)];

    const FLOAT uM10 = lv[mapd(-1, 0, 0, 0)];
    const FLOAT uM11 = lv[mapd(-1, 1, 0, 0)];
    const FLOAT vM10 = lv[mapd(-1, 0, 0, 1)];

    // this a central difference expression for the first-derivative. We therefore linearly interpolate u*v onto the surface of the
    // current cell (in 2D: upper left and upper right corner) and then take the central difference
    const FLOAT secondOrder = (  ((hyLong-hyShort)/hyLong*u00 +hyShort/hyLong*u01) * ((hxLong1-hxShort)/hxLong1*v00+hxShort/hxLong1*v10)
                               - ((hyLong-hyShort)/hyLong*uM10+hyShort/hyLong*uM11) *((hxLong0-hxShort)/hxLong0*v00+hxShort/hxLong0*vM10) )/ (2.0*hxShort);


    // this is a forward-difference in donor-cell style. We apply donor cell and again interpolate the velocity values (u-comp.)
    // onto the surface of the cell. We then apply the standard donor cell scheme. This will, however, result in non-equal
    // mesh spacing evaluations (in case of stretched meshes)
    const FLOAT kr = (hyLong-hyShort)/hyLong*u00 +hyShort/hyLong*u01;
    const FLOAT kl = (hyLong-hyShort)/hyLong*uM10+hyShort/hyLong*uM11;

    const FLOAT firstOrder  = 1.0/(4.0*hxShort)* (
                                kr*(v00+v10) - kl*(vM10+v00) + fabs(kr)*(v00 - v10) - fabs(kl)*(vM10 - v00)
                              );

    // return linear combination of central and donor-cell difference
    const FLOAT tmp2 = (1.0-parameters.solver.gamma)*secondOrder + parameters.solver.gamma*firstOrder;

//    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "Error duv_dx"); }
    return tmp2;
}


/** evaluates first derivative w.r.t. y for u*v at location of u-component. For details on implementation, see duvdx */
inline FLOAT duvdy ( const FLOAT * const lv, const Parameters & parameters, const FLOAT * const lm ) {
/*    const FLOAT tmp1 = 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,1)] + lv [mapd(1,0,0,1)] ) *
                         ( lv [mapd(0,0,0,0)] + lv [mapd(0,1,0,0)] ) ) -
                       ( ( lv [mapd(0,-1,0,1)] + lv [mapd(1,-1,0,1)] ) *
                         ( lv [mapd(0,-1,0,0)] + lv [mapd(0,0,0,0)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,1)] + lv [mapd(1,0,0,1)] ) *
                              ( lv [mapd(0,0,0,0)] - lv [mapd(0,1,0,0)] ) ) -
                       ( fabs ( lv [mapd(0,-1,0,1)] + lv [mapd(1,-1,0,1)] ) *
                              ( lv [mapd(0,-1,0,0)] - lv [mapd(0,0,0,0)] ) ) ) ) /
                       lm[mapd(0,0,0,1)];
*/
    const FLOAT hyShort = 0.5*lm[mapd( 0,0,0,1)];                       // distance of corner points in x-direction from center v-value
    const FLOAT hyLong0 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd( 0,-1,0,1)]); // distance between center and west v-value
    const FLOAT hyLong1 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd( 0,1,0,1)]); // distance between center and east v-value
    const FLOAT hxShort = 0.5*lm[mapd(0,0,0,0)];                        // distance of center u-value from upper edge of cell
    const FLOAT hxLong  = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd(1,0,0,0)]);  // distance of north and center u-value

    const FLOAT v00  = lv[mapd( 0, 0, 0, 1)];
    const FLOAT v10  = lv[mapd( 1, 0, 0, 1)];
    const FLOAT u00  = lv[mapd( 0, 0, 0, 0)];
    const FLOAT u01  = lv[mapd( 0, 1, 0, 0)];

    const FLOAT v0M1 = lv[mapd( 0,-1, 0, 1)];
    const FLOAT v1M1 = lv[mapd( 1,-1, 0, 1)];
    const FLOAT u0M1 = lv[mapd( 0,-1, 0, 0)];

    const FLOAT secondOrder = (  ((hxLong-hxShort)/hxLong*v00 +hxShort/hxLong*v10) * ((hyLong1-hyShort)/hyLong1*u00+hyShort/hyLong1*u01)
                               - ((hxLong-hxShort)/hxLong*v0M1+hxShort/hxLong*v1M1) *((hyLong0-hyShort)/hyLong0*u00+hyShort/hyLong0*u0M1) )/ (2.0*hyShort);


    const FLOAT kr = (hxLong-hxShort)/hxLong*v00 +hxShort/hxLong*v10;
    const FLOAT kl = (hxLong-hxShort)/hxLong*v0M1+hxShort/hxLong*v1M1;

    const FLOAT firstOrder  = 1.0/(4.0*hyShort)* (
                                kr*(u00+u01) - kl*(u0M1+u00) + fabs(kr)*(u00 - u01) - fabs(kl)*(u0M1 - u00)
                              );
    const FLOAT tmp2 = (1.0-parameters.solver.gamma)*secondOrder + parameters.solver.gamma*firstOrder;
//    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"Error duvdy"); }
    return tmp2;
}


/** evaluates first derivative w.r.t. x for u*w at location of w-component. For details on implementation, see duvdx */
inline FLOAT duwdx ( const FLOAT * const lv, const Parameters & parameters, const FLOAT * const lm ) {
/*    const FLOAT tmp1 = 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,0)] + lv [mapd(0,0,1,0)] ) *
                         ( lv [mapd(0,0,0,2)] + lv [mapd(1,0,0,2)] ) ) -
                       ( ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,0,1,0)] ) *
                         ( lv [mapd(-1,0,0,2)] + lv [mapd(0,0,0,2)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,0)] + lv [mapd(0,0,1,0)] ) *
                              ( lv [mapd(0,0,0,2)] - lv [mapd(1,0,0,2)] ) ) -
                       ( fabs ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,0,1,0)] ) *
                              ( lv [mapd(-1,0,0,2)] - lv [mapd(0,0,0,2)] ) ) ) ) /
                       lm[mapd(0,0,0,0)];
*/
    const FLOAT hxShort = 0.5*lm[mapd( 0,0,0,0)];                       // distance of corner points in x-direction from center v-value
    const FLOAT hxLong0 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd(-1,0,0,0)]); // distance between center and west v-value
    const FLOAT hxLong1 = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd( 1,0,0,0)]); // distance between center and east v-value
    const FLOAT hzShort = 0.5*lm[mapd(0,0,0,2)];                        // distance of center u-value from upper edge of cell
    const FLOAT hzLong  = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd(0,0,1,2)]);  // distance of north and center u-value

    const FLOAT u00  = lv[mapd( 0, 0, 0, 0)];
    const FLOAT u01  = lv[mapd( 0, 0, 1, 0)];
    const FLOAT w00  = lv[mapd( 0, 0, 0, 2)];
    const FLOAT w10  = lv[mapd( 1, 0, 0, 2)];

    const FLOAT uM10 = lv[mapd(-1, 0, 0, 0)];
    const FLOAT uM11 = lv[mapd(-1, 0, 1, 0)];
    const FLOAT wM10 = lv[mapd(-1, 0, 0, 2)];

    const FLOAT secondOrder = (  ((hzLong-hzShort)/hzLong*u00 +hzShort/hzLong*u01) * ((hxLong1-hxShort)/hxLong1*w00+hxShort/hxLong1*w10)
                               - ((hzLong-hzShort)/hzLong*uM10+hzShort/hzLong*uM11) *((hxLong0-hxShort)/hxLong0*w00+hxShort/hxLong0*wM10) )/ (2.0*hxShort);


    const FLOAT kr = (hzLong-hzShort)/hzLong*u00 +hzShort/hzLong*u01;
    const FLOAT kl = (hzLong-hzShort)/hzLong*uM10+hzShort/hzLong*uM11;

    const FLOAT firstOrder  = 1.0/(4.0*hxShort)* (
                                kr*(w00+w10) - kl*(wM10+w00) + fabs(kr)*(w00 - w10) - fabs(kl)*(wM10 - w00)
                              );
    const FLOAT tmp2 = (1.0-parameters.solver.gamma)*secondOrder + parameters.solver.gamma*firstOrder;
//    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"Error duwdx");}
    return tmp2;
}


/** evaluates first derivative w.r.t. z for u*w at location of u-component. For details on implementation, see duvdx */
inline FLOAT duwdz ( const FLOAT * const lv, const Parameters & parameters, const FLOAT * const lm ) {
/*    const FLOAT tmp1= 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,2)] + lv [mapd(1,0,0,2)] ) *
                         ( lv [mapd(0,0,0,0)] + lv [mapd(0,0,1,0)] ) ) -
                       ( ( lv [mapd(0,0,-1,2)] + lv [mapd(1,0,-1,2)] ) *
                         ( lv [mapd(0,0,-1,0)] + lv [mapd(0,0,0,0)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,2)] + lv [mapd(1,0,0,2)] ) *
                              ( lv [mapd(0,0,0,0)] - lv [mapd(0,0,1,0)] ) ) -
                       ( fabs ( lv [mapd(0,0,-1,2)] + lv [mapd(1,0,-1,2)] ) *
                              ( lv [mapd(0,0,-1,0)] - lv [mapd(0,0,0,0)] ) ) ) ) /
                       lm[mapd(0,0,0,2)];
*/
    const FLOAT hzShort = 0.5*lm[mapd( 0,0,0,2)];                       // distance of corner points in x-direction from center v-value
    const FLOAT hzLong0 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd( 0,0,-1,2)]); // distance between center and west v-value
    const FLOAT hzLong1 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd( 0,0, 1,2)]); // distance between center and east v-value
    const FLOAT hxShort = 0.5*lm[mapd(0,0,0,0)];                        // distance of center u-value from upper edge of cell
    const FLOAT hxLong  = 0.5*(lm[mapd(0,0,0,0)] + lm[mapd(1,0,0,0)]);  // distance of north and center u-value

    const FLOAT w00  = lv[mapd( 0, 0, 0, 2)];
    const FLOAT w10  = lv[mapd( 1, 0, 0, 2)];
    const FLOAT u00  = lv[mapd( 0, 0, 0, 0)];
    const FLOAT u01  = lv[mapd( 0, 0, 1, 0)];

    const FLOAT w0M1 = lv[mapd( 0, 0,-1, 2)];
    const FLOAT w1M1 = lv[mapd( 1, 0,-1, 2)];
    const FLOAT u0M1 = lv[mapd( 0, 0,-1, 0)];

    const FLOAT secondOrder = (  ((hxLong-hxShort)/hxLong*w00 +hxShort/hxLong*w10) * ((hzLong1-hzShort)/hzLong1*u00+hzShort/hzLong1*u01)
                               - ((hxLong-hxShort)/hxLong*w0M1+hxShort/hxLong*w1M1) *((hzLong0-hzShort)/hzLong0*u00+hzShort/hzLong0*u0M1) )/ (2.0*hzShort);


    const FLOAT kr = (hxLong-hxShort)/hxLong*w00 +hxShort/hxLong*w10;
    const FLOAT kl = (hxLong-hxShort)/hxLong*w0M1+hxShort/hxLong*w1M1;

    const FLOAT firstOrder  = 1.0/(4.0*hzShort)* (
                                kr*(u00+u01) - kl*(u0M1+u00) + fabs(kr)*(u00 - u01) - fabs(kl)*(u0M1 - u00)
                              );
    const FLOAT tmp2 = (1.0-parameters.solver.gamma)*secondOrder + parameters.solver.gamma*firstOrder;

//    if (fabs(tmp1-tmp2)> 1.0e-12){handleError(1,"duwdz");}
    return tmp2;
}


/** evaluates first derivative w.r.t. y for v*w at location of w-component. For details on implementation, see duvdx */
inline FLOAT dvwdy ( const FLOAT * const lv, const Parameters & parameters,const FLOAT * const lm ) {
/*    const FLOAT tmp1 =  1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,1)] + lv [mapd(0,0,1,1)] ) *
                         ( lv [mapd(0,0,0,2)] + lv [mapd(0,1,0,2)] ) ) -
                       ( ( lv [mapd(0,-1,0,1)] + lv [mapd(0,-1,1,1)] ) *
                         ( lv [mapd(0,-1,0,2)] + lv [mapd(0,0,0,2)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,1)] + lv [mapd(0,0,1,1)] ) *
                              ( lv [mapd(0,0,0,2)] - lv [mapd(0,1,0,2)] ) ) -
                       ( fabs ( lv [mapd(0,-1,0,1)] + lv [mapd(0,-1,1,1)] ) *
                              ( lv [mapd(0,-1,0,2)] - lv [mapd(0,0,0,2)] ) ) ) ) /
                       lm[mapd(0,0,0,1)];
*/
    const FLOAT hyShort = 0.5*lm[mapd( 0,0,0,1)];                       // distance of corner points in x-direction from center v-value
    const FLOAT hyLong0 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd(0,-1,0,1)]); // distance between center and west v-value
    const FLOAT hyLong1 = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd( 0,1,0,1)]); // distance between center and east v-value
    const FLOAT hzShort = 0.5*lm[mapd(0,0,0,2)];                        // distance of center u-value from upper edge of cell
    const FLOAT hzLong  = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd(0,0,1,2)]);  // distance of north and center u-value

    const FLOAT v00  = lv[mapd( 0, 0, 0, 1)];
    const FLOAT v01  = lv[mapd( 0, 0, 1, 1)];
    const FLOAT w00  = lv[mapd( 0, 0, 0, 2)];
    const FLOAT w10  = lv[mapd( 0, 1, 0, 2)];

    const FLOAT vM10 = lv[mapd( 0,-1, 0, 1)];
    const FLOAT vM11 = lv[mapd( 0,-1, 1, 1)];
    const FLOAT wM10 = lv[mapd( 0,-1, 0, 2)];

    const FLOAT secondOrder = (  ((hzLong-hzShort)/hzLong*v00 +hzShort/hzLong*v01) * ((hyLong1-hyShort)/hyLong1*w00+hyShort/hyLong1*w10)
                               - ((hzLong-hzShort)/hzLong*vM10+hzShort/hzLong*vM11) *((hyLong0-hyShort)/hyLong0*w00+hyShort/hyLong0*wM10) )/ (2.0*hyShort);


    const FLOAT kr = (hzLong-hzShort)/hzLong*v00 +hzShort/hzLong*v01;
    const FLOAT kl = (hzLong-hzShort)/hzLong*vM10+hzShort/hzLong*vM11;

    const FLOAT firstOrder  = 1.0/(4.0*hyShort)* (
                                kr*(w00+w10) - kl*(wM10+w00) + fabs(kr)*(w00 - w10) - fabs(kl)*(wM10 - w00)
                              );
    const FLOAT tmp2 = (1.0-parameters.solver.gamma)*secondOrder + parameters.solver.gamma*firstOrder;
//    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dvwdy");}
    return tmp2;
}


/** evaluates first derivative w.r.t. z for v*w at location of v-component. For details on implementation, see duvdx */
inline FLOAT dvwdz ( const FLOAT * const lv, const Parameters & parameters, const FLOAT * const lm ) {
/*    const FLOAT tmp1 = 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,2)] + lv [mapd(0,1,0,2)] ) *
                         ( lv [mapd(0,0,0,1)] + lv [mapd(0,0,1,1)] ) ) -
                       ( ( lv [mapd(0,0,-1,2)] + lv [mapd(0,1,-1,2)] ) *
                         ( lv [mapd(0,0,-1,1)] + lv [mapd(0,0,0,1)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,2)] + lv [mapd(0,1,0,2)] ) *
                              ( lv [mapd(0,0,0,1)] - lv [mapd(0,0,1,1)] ) ) -
                       ( fabs ( lv [mapd(0,0,-1,2)] + lv [mapd(0,1,-1,2)] ) *
                              ( lv [mapd(0,0,-1,1)] - lv [mapd(0,0,0,1)] ) ) ) ) /
                       lm[mapd(0,0,0,2)];
*/
    const FLOAT hzShort = 0.5*lm[mapd( 0,0,0,2)];                       // distance of corner points in x-direction from center v-value
    const FLOAT hzLong0 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd( 0,0,-1,2)]); // distance between center and west v-value
    const FLOAT hzLong1 = 0.5*(lm[mapd(0,0,0,2)] + lm[mapd( 0,0, 1,2)]); // distance between center and east v-value
    const FLOAT hyShort = 0.5*lm[mapd(0,0,0,1)];                        // distance of center u-value from upper edge of cell
    const FLOAT hyLong  = 0.5*(lm[mapd(0,0,0,1)] + lm[mapd(0,1,0,1)]);  // distance of north and center u-value

    const FLOAT w00  = lv[mapd( 0, 0, 0, 2)];
    const FLOAT w10  = lv[mapd( 0, 1, 0, 2)];
    const FLOAT v00  = lv[mapd( 0, 0, 0, 1)];
    const FLOAT v01  = lv[mapd( 0, 0, 1, 1)];

    const FLOAT w0M1 = lv[mapd( 0, 0,-1, 2)];
    const FLOAT w1M1 = lv[mapd( 0, 1,-1, 2)];
    const FLOAT v0M1 = lv[mapd( 0, 0,-1, 1)];

    const FLOAT secondOrder = (  ((hyLong-hyShort)/hyLong*w00 +hyShort/hyLong*w10) * ((hzLong1-hzShort)/hzLong1*v00+hzShort/hzLong1*v01)
                               - ((hyLong-hyShort)/hyLong*w0M1+hyShort/hyLong*w1M1) *((hzLong0-hzShort)/hzLong0*v00+hzShort/hzLong0*v0M1) )/ (2.0*hzShort);


    const FLOAT kr = (hyLong-hyShort)/hyLong*w00 +hyShort/hyLong*w10;
    const FLOAT kl = (hyLong-hyShort)/hyLong*w0M1+hyShort/hyLong*w1M1;

    const FLOAT firstOrder  = 1.0/(4.0*hzShort)* (
                                kr*(v00+v01) - kl*(v0M1+v00) + fabs(kr)*(v00 - v01) - fabs(kl)*(v0M1 - v00)
                              );
    const FLOAT tmp2 = (1.0-parameters.solver.gamma)*secondOrder + parameters.solver.gamma*firstOrder;
//    if (fabs(tmp1-tmp2) > 1.0e-12){std::cout << tmp1 << ", " << tmp2 << std::endl;handleError(1,"dvwdz");}
    return tmp2;
}


/** first derivative of u*u w.r.t. x, evaluated at location of u-component. */
inline FLOAT du2dx ( const FLOAT * const lv, const Parameters & parameters, const FLOAT * const lm ) {
/*    const FLOAT tmp1 = 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,0)] + lv [mapd(1,0,0,0)] ) *
                         ( lv [mapd(0,0,0,0)] + lv [mapd(1,0,0,0)] ) ) -
                       ( ( lv [mapd(-1,0,0,0)] + lv [mapd(0,0,0,0)] ) *
                         ( lv [mapd(-1,0,0,0)] + lv [mapd(0,0,0,0)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,0)] + lv [mapd(1,0,0,0)] ) *
                              ( lv [mapd(0,0,0,0)] - lv [mapd(1,0,0,0)] ) ) -
                       ( fabs ( lv [mapd(-1,0,0,0)] + lv [mapd(0,0,0,0)] ) *
                              ( lv [mapd(-1,0,0,0)] - lv [mapd(0,0,0,0)] ) ) ) ) /
                       lm[mapd(0,0,0,0)];
*/
    const FLOAT dxShort = 0.5*lm[mapd(0,0,0,0)];
    const FLOAT dxLong0 = 0.5*(lm[mapd(-1,0,0,0)] + lm[mapd(0,0,0,0)]);
    const FLOAT dxLong1 = 0.5*(lm[mapd( 0,0,0,0)] + lm[mapd(1,0,0,0)]);

    const FLOAT u0 = lv[mapd(0,0,0,0)];
    const FLOAT uM1= lv[mapd(-1,0,0,0)];
    const FLOAT u1 = lv[mapd(1,0,0,0)];

    const FLOAT kr = (dxLong1-dxShort)/dxLong1*u0 + dxShort/dxLong1*u1;
    const FLOAT kl = (dxLong0-dxShort)/dxLong0*u0 + dxShort/dxLong0*uM1;

    // central difference expression which is second-order accurate for uniform meshes. We interpolate u half-way between
    // neighboured u-component values and afterwards build the central difference for u*u
    const FLOAT secondOrder = (  ((dxLong1-dxShort)/dxLong1*u0 + dxShort/dxLong1*u1 )*((dxLong1-dxShort)/dxLong1*u0 + dxShort/dxLong1*u1 )
                               - ((dxLong0-dxShort)/dxLong0*u0 + dxShort/dxLong0*uM1)*((dxLong0-dxShort)/dxLong0*u0 + dxShort/dxLong0*uM1)
                              )/(2.0*dxShort);

    // donor-cell like derivative expression. We evaluate u half-way between neighboured u-components and use this as a prediction
    // of the transport direction
    const FLOAT firstOrder = 1.0/(4.0*dxShort)* (
                               kr*(u0+u1) - kl*(uM1+u0) + fabs(kr)*(u0 - u1) - fabs(kl)*(uM1 - u0)
                             );

    // return linear combination of central- and upwind difference
    const FLOAT tmp2 = (1.0-parameters.solver.gamma)*secondOrder + parameters.solver.gamma*firstOrder;
//    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"du2dx");}
    return tmp2;
}


/** first derivative of v*v w.r.t. y, evaluated at location of v-component; for details, see du2dx */
inline FLOAT dv2dy ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm ) {
/*    const FLOAT tmp1= 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,1)] + lv [mapd(0,1,0,1)] ) *
                         ( lv [mapd(0,0,0,1)] + lv [mapd(0,1,0,1)] ) ) -
                       ( ( lv [mapd(0,-1,0,1)] + lv [mapd(0,0,0,1)] ) *
                         ( lv [mapd(0,-1,0,1)] + lv [mapd(0,0,0,1)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,1)] + lv [mapd(0,1,0,1)] ) *
                              ( lv [mapd(0,0,0,1)] - lv [mapd(0,1,0,1)] ) ) -
                       ( fabs ( lv [mapd(0,-1,0,1)] + lv [mapd(0,0,0,1)] ) *
                              ( lv [mapd(0,-1,0,1)] - lv [mapd(0,0,0,1)] ) ) ) ) /
                       lm[mapd(0,0,0,1)];
*/
    const FLOAT dyShort = 0.5*lm[mapd(0,0,0,1)];
    const FLOAT dyLong0 = 0.5*(lm[mapd(0,-1,0,1)] + lm[mapd(0,0,0,1)]);
    const FLOAT dyLong1 = 0.5*(lm[mapd( 0,0,0,1)] + lm[mapd(0,1,0,1)]);

    const FLOAT v0 = lv[mapd(0,0,0,1)];
    const FLOAT vM1= lv[mapd(0,-1,0,1)];
    const FLOAT v1 = lv[mapd(0,1,0,1)];

    const FLOAT kr = (dyLong1-dyShort)/dyLong1*v0 + dyShort/dyLong1*v1;
    const FLOAT kl = (dyLong0-dyShort)/dyLong0*v0 + dyShort/dyLong0*vM1;

    const FLOAT secondOrder = (  ((dyLong1-dyShort)/dyLong1*v0 + dyShort/dyLong1*v1 )*((dyLong1-dyShort)/dyLong1*v0 + dyShort/dyLong1*v1 )
                               - ((dyLong0-dyShort)/dyLong0*v0 + dyShort/dyLong0*vM1)*((dyLong0-dyShort)/dyLong0*v0 + dyShort/dyLong0*vM1)
                              )/(2.0*dyShort);
    const FLOAT firstOrder = 1.0/(4.0*dyShort)* (
                               kr*(v0+v1) - kl*(vM1+v0) + fabs(kr)*(v0 - v1) - fabs(kl)*(vM1 - v0)
                             );
    const FLOAT tmp2 = (1.0-parameters.solver.gamma)*secondOrder + parameters.solver.gamma*firstOrder;
//    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dv2dy");}
    return tmp2;
}


/** first derivative of w*w w.r.t. z, evaluated at location of w-component; for details, see du2dx */
inline FLOAT dw2dz ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm ) {
/*    const FLOAT tmp1= 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,2)] + lv [mapd(0,0,1,2)] ) *
                         ( lv [mapd(0,0,0,2)] + lv [mapd(0,0,1,2)] ) ) -
                       ( ( lv [mapd(0,0,-1,2)] + lv [mapd(0,0,0,2)] ) *
                         ( lv [mapd(0,0,-1,2)] + lv [mapd(0,0,0,2)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,2)] + lv [mapd(0,0,1,2)] ) *
                              ( lv [mapd(0,0,0,2)] - lv [mapd(0,0,1,2)] ) ) -
                       ( fabs ( lv [mapd(0,0,-1,2)] + lv [mapd(0,0,0,2)] ) *
                              ( lv [mapd(0,0,-1,2)] - lv [mapd(0,0,0,2)] ) ) ) ) /
                       lm[mapd(0,0,0,2)];
*/
    const FLOAT dzShort = 0.5*lm[mapd(0,0,0,2)];
    const FLOAT dzLong0 = 0.5*(lm[mapd(0,0,-1,2)] + lm[mapd(0,0,0,2)]);
    const FLOAT dzLong1 = 0.5*(lm[mapd( 0,0,0,2)] + lm[mapd(0,0,1,2)]);

    const FLOAT w0 = lv[mapd(0,0,0,2)];
    const FLOAT wM1= lv[mapd(0,0,-1,2)];
    const FLOAT w1 = lv[mapd(0,0,1,2)];

    const FLOAT kr = (dzLong1-dzShort)/dzLong1*w0 + dzShort/dzLong1*w1;
    const FLOAT kl = (dzLong0-dzShort)/dzLong0*w0 + dzShort/dzLong0*wM1;

    const FLOAT secondOrder = (  ((dzLong1-dzShort)/dzLong1*w0 + dzShort/dzLong1*w1 )*((dzLong1-dzShort)/dzLong1*w0 + dzShort/dzLong1*w1 )
                               - ((dzLong0-dzShort)/dzLong0*w0 + dzShort/dzLong0*wM1)*((dzLong0-dzShort)/dzLong0*w0 + dzShort/dzLong0*wM1)
                              )/(2.0*dzShort);
    const FLOAT firstOrder = 1.0/(4.0*dzShort)* (
                               kr*(w0+w1) - kl*(wM1+w0) + fabs(kr)*(w0 - w1) - fabs(kl)*(wM1 - w0)
                             );
    const FLOAT tmp2 = (1.0-parameters.solver.gamma)*secondOrder + parameters.solver.gamma*firstOrder;
//    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dw2dz");}
    return tmp2;
}


inline FLOAT dxdudx ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm, const FLOAT* lt ) {

    const FLOAT dx = lm[mapd(0,0,0,0)];
    const FLOAT dxP1= lm[mapd(1,0,0,0)];

    const FLOAT u0 = lv[mapd(0,0,0,0)];
    const FLOAT uM1= lv[mapd(-1,0,0,0)];
    const FLOAT uP1 = lv[mapd(1,0,0,0)];

    // Kinematic viscosity
    const FLOAT kni = 1.0 / parameters.flow.Re;
    // Turbulent viscosity
    const FLOAT tni0 = lt[mapd(0,0,0,0)];
    const FLOAT tni1 = lt[mapd(1,0,0,0)];
    // Total viscosity
    const FLOAT ni0 = kni + tni0;
    const FLOAT ni1 = kni + tni1;

    return 2.0 /(dx+dxP1)*(ni1*(uP1-u0)/dxP1-ni0*(u0-uM1)/dx);

}

inline FLOAT dydvdy ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm, const FLOAT* lt  ) {

    const FLOAT dy = lm[mapd(0,0,0,1)];
    const FLOAT dyP1= lm[mapd(0,1,0,1)];

    const FLOAT v0 = lv[mapd(0,0,0,1)];
    const FLOAT vM1= lv[mapd(0,-1,0,1)];
    const FLOAT vP1 = lv[mapd(0,1,0,1)];

    // Kinematic viscosity
    const FLOAT kni = 1.0 / parameters.flow.Re;
    // Turbulent viscosity
    const FLOAT tni0 = lt[mapd(0,0,0,0)];
    const FLOAT tni1 = lt[mapd(0,1,0,0)];
    // Total viscosity
    const FLOAT ni0 = kni + tni0;
    const FLOAT ni1 = kni + tni1;

return 2.0 /(dy+dyP1)*(ni1*(vP1-v0)/dyP1-ni0*(v0-vM1)/dy);

}

inline FLOAT dzdwdz ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm, const FLOAT* lt  ) {

    const FLOAT dz = lm[mapd(0,0,0,2)];
    const FLOAT dzP1= lm[mapd(0,0,1,2)];

    const FLOAT w0 = lv[mapd(0,0,0,2)];
    const FLOAT wM1= lv[mapd(0,0,-1,2)];
    const FLOAT wP1 = lv[mapd(0,0,1,2)];

    // Kinematic viscosity
    const FLOAT kni = 1.0 / parameters.flow.Re;
    // Turbulent viscosity
    const FLOAT tni0 = lt[mapd(0,0,0,0)];
    const FLOAT tni1 = lt[mapd(0,0,1,0)];
    // Total viscosity
    const FLOAT ni0 = kni + tni0;
    const FLOAT ni1 = kni + tni1;

return 2.0 /(dz+dzP1)*(ni1*(wP1-w0)/dzP1-ni0*(w0-wM1)/dz);

}

inline FLOAT diffusive_term(
	const FLOAT r1, const FLOAT r2, const FLOAT r3, 
	const FLOAT p1, const FLOAT p2, const FLOAT p3, const FLOAT p4, 
	const FLOAT nu0, const FLOAT nu1, const FLOAT nu2, const FLOAT nu3, const FLOAT nu4, const FLOAT nu5, const FLOAT nu6,
	const FLOAT d1, const FLOAT d2,
	const FLOAT h1, const FLOAT h2, const FLOAT h3)
{
	
	const FLOAT drlow = (r2 - r1) / (h1 + h2) * 2;
	const FLOAT dplow = (p2 - p1) / (d1 + d2) * 2;
	const FLOAT nulow = ((nu1*h2 + nu3*h1)*d2 + (nu2*h2 + nu4*h1)*d1) / ((h1 + h2) * (d1 + d2));
	
	const FLOAT drhigh = (r3 - r2) / (h2 + h3) * 2;
	const FLOAT dphigh = (p4 - p3) / (d1 + d2) * 2;
	const FLOAT nuhigh = ((nu3*h3 + nu5*h2)*d2 + (nu4*h3 + nu6*h2)*d1) / ((h2 + h3) * (d1 + d2));
	
	return ((nuhigh+nu0)*(drhigh + dphigh) - (nulow+nu0)*(drlow + dplow)) / h2;
	
}

inline FLOAT dy_dudy_dvdx ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm, const FLOAT* lt  ) {
	
	// u velocities j=(-1,0,1)
	const FLOAT r1 = lv[mapd( 0,-1, 0, 0)];
	const FLOAT r2 = lv[mapd( 0, 0, 0, 0)];
	const FLOAT r3 = lv[mapd( 0, 1, 0, 0)];
	
	// v velocities i=(0,1) x j=(-1,0)
	const FLOAT p1 = lv[mapd( 0,-1, 0, 1)];
	const FLOAT p2 = lv[mapd( 1,-1, 0, 1)];
	const FLOAT p3 = lv[mapd( 0, 0, 0, 1)];
	const FLOAT p4 = lv[mapd( 1, 0, 0, 1)];
	
	// dx i=(0,1)
	const FLOAT d1 = lm[mapd( 0, 0, 0, 0)];
	const FLOAT d2 = lm[mapd( 1, 0, 0, 0)];
	
	// dy j=(-1,0,1)
	const FLOAT h1 = lm[mapd( 0,-1, 0, 1)];
	const FLOAT h2 = lm[mapd( 0, 0, 0, 1)];
	const FLOAT h3 = lm[mapd( 0, 1, 0, 1)];
	
	// viscosity i=(0,1) x j=(-1,0,1)
	const FLOAT nu0 = 1.0 / parameters.flow.Re;
	const FLOAT nu1 = lt[mapd( 0,-1, 0, 0)];
	const FLOAT nu2 = lt[mapd( 1,-1, 0, 0)];
	const FLOAT nu3 = lt[mapd( 0, 0, 0, 0)];
	const FLOAT nu4 = lt[mapd( 1, 0, 0, 0)];
	const FLOAT nu5 = lt[mapd( 0, 1, 0, 0)];
	const FLOAT nu6 = lt[mapd( 1, 1, 0, 0)];
	
	return diffusive_term(r1, r2, r3, p1, p2, p3, p4, nu0, nu1, nu2, nu3, nu4, nu5, nu6, d1, d2, h1, h2, h3);

}

inline FLOAT dz_dudz_dwdx ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm, const FLOAT* lt  ) {

	// u velocities k=(-1,0,1)
	const FLOAT r1 = lv[mapd( 0, 0,-1, 0)];
	const FLOAT r2 = lv[mapd( 0, 0, 0, 0)];
	const FLOAT r3 = lv[mapd( 0, 0, 1, 0)];
	
	// w velocities i=(0,1) x k=(-1,0)
	const FLOAT p1 = lv[mapd( 0, 0,-1, 2)];
	const FLOAT p2 = lv[mapd( 1, 0,-1, 2)];
	const FLOAT p3 = lv[mapd( 0, 0, 0, 2)];
	const FLOAT p4 = lv[mapd( 1, 0, 0, 2)];
	
	// dx i=(0,1)
	const FLOAT d1 = lm[mapd( 0, 0, 0, 0)];
	const FLOAT d2 = lm[mapd( 1, 0, 0, 0)];
	
	// dz k=(-1,0,1)
	const FLOAT h1 = lm[mapd( 0, 0,-1, 2)];
	const FLOAT h2 = lm[mapd( 0, 0, 0, 2)];
	const FLOAT h3 = lm[mapd( 0, 0, 1, 2)];
	
	// viscosity i=(0,1) x k=(-1,0,1)
	const FLOAT nu0 = 1.0 / parameters.flow.Re;
	const FLOAT nu1 = lt[mapd( 0, 0,-1, 0)];
	const FLOAT nu2 = lt[mapd( 1, 0,-1, 0)];
	const FLOAT nu3 = lt[mapd( 0, 0, 0, 0)];
	const FLOAT nu4 = lt[mapd( 1, 0, 0, 0)];
	const FLOAT nu5 = lt[mapd( 0, 0, 1, 0)];
	const FLOAT nu6 = lt[mapd( 1, 0, 1, 0)];
	
	return diffusive_term(r1, r2, r3, p1, p2, p3, p4, nu0, nu1, nu2, nu3, nu4, nu5, nu6, d1, d2, h1, h2, h3);

}

inline FLOAT dx_dudy_dvdx ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm, const FLOAT* lt  ) {

	// v velocities i=(-1,0,1)
	const FLOAT r1 = lv[mapd(-1, 0, 0, 1)];
	const FLOAT r2 = lv[mapd( 0, 0, 0, 1)];
	const FLOAT r3 = lv[mapd( 1, 0, 0, 1)];
	
	// u velocities i=(-1,0) x j=(0,1)
	const FLOAT p1 = lv[mapd(-1, 0, 0, 0)];
	const FLOAT p2 = lv[mapd(-1, 1, 0, 0)];
	const FLOAT p3 = lv[mapd( 0, 0, 0, 0)];
	const FLOAT p4 = lv[mapd( 0, 1, 0, 0)];
	
	// dy j=(0,1)
	const FLOAT d1 = lm[mapd( 0, 0, 0, 1)];
	const FLOAT d2 = lm[mapd( 0, 1, 0, 1)];
	
	// dx i=(-1,0,1)
	const FLOAT h1 = lm[mapd(-1, 0, 0, 0)];
	const FLOAT h2 = lm[mapd( 0, 0, 0, 0)];
	const FLOAT h3 = lm[mapd( 1, 0, 0, 0)];
	
	// viscosity j=(0,1) x i=(-1,0,1)
	const FLOAT nu0 = 1.0 / parameters.flow.Re;
	const FLOAT nu1 = lt[mapd(-1, 0, 0, 0)];
	const FLOAT nu2 = lt[mapd(-1, 1, 0, 0)];
	const FLOAT nu3 = lt[mapd( 0, 0, 0, 0)];
	const FLOAT nu4 = lt[mapd( 0, 1, 0, 0)];
	const FLOAT nu5 = lt[mapd( 1, 0, 0, 0)];
	const FLOAT nu6 = lt[mapd( 1, 1, 0, 0)];
	
	return diffusive_term(r1, r2, r3, p1, p2, p3, p4, nu0, nu1, nu2, nu3, nu4, nu5, nu6, d1, d2, h1, h2, h3);

}

inline FLOAT dz_dvdz_dwdy ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm, const FLOAT* lt  ) {

	// v velocities k=(-1,0,1)
	const FLOAT r1 = lv[mapd( 0, 0,-1, 1)];
	const FLOAT r2 = lv[mapd( 0, 0, 0, 1)];
	const FLOAT r3 = lv[mapd( 0, 0, 1, 1)];
	
	// w velocities k=(-1,0) x j=(0,1)
	const FLOAT p1 = lv[mapd( 0, 0,-1, 2)];
	const FLOAT p2 = lv[mapd( 0, 1,-1, 2)];
	const FLOAT p3 = lv[mapd( 0, 0, 0, 2)];
	const FLOAT p4 = lv[mapd( 0, 1, 0, 2)];
	
	// dy j=(0,1)
	const FLOAT d1 = lm[mapd( 0, 0, 0, 1)];
	const FLOAT d2 = lm[mapd( 0, 1, 0, 1)];
	
	// dz k=(-1,0,1)
	const FLOAT h1 = lm[mapd( 0, 0,-1, 2)];
	const FLOAT h2 = lm[mapd( 0, 0, 0, 2)];
	const FLOAT h3 = lm[mapd( 0, 0, 1, 2)];
	
	// viscosity j=(0,1) x k=(-1,0,1)
	const FLOAT nu0 = 1.0 / parameters.flow.Re;
	const FLOAT nu1 = lt[mapd( 0, 0,-1, 0)];
	const FLOAT nu2 = lt[mapd( 0, 1,-1, 0)];
	const FLOAT nu3 = lt[mapd( 0, 0, 0, 0)];
	const FLOAT nu4 = lt[mapd( 0, 1, 0, 0)];
	const FLOAT nu5 = lt[mapd( 0, 0, 1, 0)];
	const FLOAT nu6 = lt[mapd( 0, 1, 1, 0)];
	
	return diffusive_term(r1, r2, r3, p1, p2, p3, p4, nu0, nu1, nu2, nu3, nu4, nu5, nu6, d1, d2, h1, h2, h3);

}

inline FLOAT dx_dudz_dwdx ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm, const FLOAT* lt  ) {

	// w velocities i=(-1,0,1)
	const FLOAT r1 = lv[mapd(-1, 0, 0, 2)];
	const FLOAT r2 = lv[mapd( 0, 0, 0, 2)];
	const FLOAT r3 = lv[mapd( 1, 0, 0, 2)];
	
	// u velocities i=(-1,0) x k=(0,1)
	const FLOAT p1 = lv[mapd(-1, 0, 0, 0)];
	const FLOAT p2 = lv[mapd(-1, 0, 1, 0)];
	const FLOAT p3 = lv[mapd( 0, 0, 0, 0)];
	const FLOAT p4 = lv[mapd( 0, 0, 1, 0)];
	
	// dz k=(0,1)
	const FLOAT d1 = lm[mapd( 0, 0, 0, 2)];
	const FLOAT d2 = lm[mapd( 0, 0, 1, 2)];
	
	// dx i=(-1,0,1)
	const FLOAT h1 = lm[mapd(-1, 0, 0, 0)];
	const FLOAT h2 = lm[mapd( 0, 0, 0, 0)];
	const FLOAT h3 = lm[mapd( 1, 0, 0, 0)];
	
	// viscosity k=(0,1) x i=(-1,0,1)
	const FLOAT nu0 = 1.0 / parameters.flow.Re;
	const FLOAT nu1 = lt[mapd(-1, 0, 0, 0)];
	const FLOAT nu2 = lt[mapd(-1, 0, 1, 0)];
	const FLOAT nu3 = lt[mapd( 0, 0, 0, 0)];
	const FLOAT nu4 = lt[mapd( 0, 0, 1, 0)];
	const FLOAT nu5 = lt[mapd( 1, 0, 0, 0)];
	const FLOAT nu6 = lt[mapd( 1, 0, 1, 0)];
	
	return diffusive_term(r1, r2, r3, p1, p2, p3, p4, nu0, nu1, nu2, nu3, nu4, nu5, nu6, d1, d2, h1, h2, h3);

}

inline FLOAT dy_dvdz_dwdy ( const FLOAT * const lv, const Parameters & parameters, const FLOAT* const lm, const FLOAT* lt  ) {

	// w velocities j=(-1,0,1)
	const FLOAT r1 = lv[mapd( 0,-1, 0, 2)];
	const FLOAT r2 = lv[mapd( 0, 0, 0, 2)];
	const FLOAT r3 = lv[mapd( 0, 1, 0, 2)];
	
	// v velocities j=(-1,0) x k=(0,1)
	const FLOAT p1 = lv[mapd( 0,-1, 0, 0)];
	const FLOAT p2 = lv[mapd( 0,-1, 1, 0)];
	const FLOAT p3 = lv[mapd( 0, 0, 0, 0)];
	const FLOAT p4 = lv[mapd( 0, 0, 1, 0)];
	
	// dz k=(0,1)
	const FLOAT d1 = lm[mapd( 0, 0, 0, 2)];
	const FLOAT d2 = lm[mapd( 0, 0, 1, 2)];
	
	// dy j=(-1,0,1)
	const FLOAT h1 = lm[mapd( 0,-1, 0, 1)];
	const FLOAT h2 = lm[mapd( 0, 0, 0, 1)];
	const FLOAT h3 = lm[mapd( 0, 1, 0, 1)];
	
	// viscosity k=(0,1) x j=(-1,0,1)
	const FLOAT nu0 = 1.0 / parameters.flow.Re;
	const FLOAT nu1 = lt[mapd( 0,-1, 0, 0)];
	const FLOAT nu2 = lt[mapd( 0,-1, 1, 0)];
	const FLOAT nu3 = lt[mapd( 0, 0, 0, 0)];
	const FLOAT nu4 = lt[mapd( 0, 0, 1, 0)];
	const FLOAT nu5 = lt[mapd( 0, 1, 0, 0)];
	const FLOAT nu6 = lt[mapd( 0, 1, 1, 0)];
	
	return diffusive_term(r1, r2, r3, p1, p2, p3, p4, nu0, nu1, nu2, nu3, nu4, nu5, nu6, d1, d2, h1, h2, h3);


}

inline FLOAT computeF2D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,0)]
        + dt * ( 1 / parameters.flow.Re * ( d2udx2 ( localVelocity, localMeshsize )
                    + d2udy2(localVelocity, localMeshsize)) - du2dx (localVelocity, parameters, localMeshsize)
                    - duvdy (localVelocity, parameters, localMeshsize) + parameters.environment.gx);
}

inline FLOAT computeG2D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,1)]
        + dt * ( 1 / parameters.flow.Re * ( d2vdx2 ( localVelocity, localMeshsize )
                    + d2vdy2(localVelocity, localMeshsize)) - duvdx (localVelocity, parameters, localMeshsize)
                    - dv2dy (localVelocity, parameters, localMeshsize) + parameters.environment.gy);
}


inline FLOAT computeF3D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,0)]
                +  dt * ( 1 / parameters.flow.Re * ( d2udx2 ( localVelocity, localMeshsize )
                + d2udy2 ( localVelocity, localMeshsize ) + d2udz2 ( localVelocity, localMeshsize ) )
                - du2dx ( localVelocity, parameters, localMeshsize ) - duvdy ( localVelocity, parameters, localMeshsize )
                - duwdz ( localVelocity, parameters, localMeshsize ) + parameters.environment.gx );
}


inline FLOAT computeG3D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,1)]
                +  dt * ( 1 / parameters.flow.Re * ( d2vdx2 ( localVelocity, localMeshsize )
                + d2vdy2 ( localVelocity, localMeshsize ) + d2vdz2 ( localVelocity, localMeshsize ) )
                - dv2dy ( localVelocity, parameters, localMeshsize ) - duvdx ( localVelocity, parameters, localMeshsize )
                - dvwdz ( localVelocity, parameters, localMeshsize ) + parameters.environment.gy );
}


inline FLOAT computeH3D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,2)]
                +  dt * ( 1 / parameters.flow.Re * ( d2wdx2 ( localVelocity, localMeshsize )
                + d2wdy2 ( localVelocity, localMeshsize ) + d2wdz2 ( localVelocity, localMeshsize ) )
                - dw2dz ( localVelocity, parameters, localMeshsize ) - duwdx ( localVelocity, parameters, localMeshsize )
                - dvwdy ( localVelocity, parameters, localMeshsize ) + parameters.environment.gz );
}

inline FLOAT computeF2DTurbulence(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const FLOAT * const localTurbulentViscosity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,0)]
        + dt * (     2* dxdudx       (localVelocity, parameters, localMeshsize, localTurbulentViscosity) + dy_dudy_dvdx (localVelocity, parameters, localMeshsize, localTurbulentViscosity)
                    - du2dx        (localVelocity, parameters, localMeshsize)                          - duvdy        (localVelocity, parameters, localMeshsize)
                     + parameters.environment.gx);
}

inline FLOAT computeG2DTurbulence(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const FLOAT * const localTurbulentViscosity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,1)]
        + dt * (      dx_dudy_dvdx (localVelocity, parameters, localMeshsize, localTurbulentViscosity) +2* dydvdy       (localVelocity, parameters, localMeshsize, localTurbulentViscosity)
                    - duvdx        (localVelocity, parameters, localMeshsize)                          - dv2dy        (localVelocity, parameters, localMeshsize)
                    + parameters.environment.gy);
}

inline FLOAT computeF3DTurbulence(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const FLOAT * const localTurbulentViscosity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,0)]
        + dt * (     2* dxdudx       (localVelocity, parameters, localMeshsize, localTurbulentViscosity) + dy_dudy_dvdx (localVelocity, parameters, localMeshsize, localTurbulentViscosity)
                    + dz_dudz_dwdx (localVelocity, parameters, localMeshsize, localTurbulentViscosity)
                    - du2dx        (localVelocity, parameters, localMeshsize)                          - duvdy        (localVelocity, parameters, localMeshsize)
                    - duwdz        (localVelocity, parameters, localMeshsize)
                    + parameters.environment.gx );
}

inline FLOAT computeG3DTurbulence(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const FLOAT * const localTurbulentViscosity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,1)]
        + dt * (      dx_dudy_dvdx (localVelocity, parameters, localMeshsize, localTurbulentViscosity) + 2*dydvdy       (localVelocity, parameters, localMeshsize, localTurbulentViscosity)
                    + dz_dvdz_dwdy (localVelocity, parameters, localMeshsize, localTurbulentViscosity)
                    - dv2dy        (localVelocity, parameters, localMeshsize)                          - duvdx        (localVelocity, parameters, localMeshsize)
                    - dvwdz        (localVelocity, parameters, localMeshsize)
                    + parameters.environment.gy );
}

inline FLOAT computeH3DTurbulence(const FLOAT * const localVelocity, const FLOAT * const localMeshsize, const FLOAT * const localTurbulentViscosity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,2)]
        +  dt * (     dx_dudz_dwdx (localVelocity, parameters, localMeshsize, localTurbulentViscosity) + dy_dvdz_dwdy (localVelocity, parameters, localMeshsize, localTurbulentViscosity)
                    +2* dzdwdz       (localVelocity, parameters, localMeshsize, localTurbulentViscosity)
                    - dw2dz        (localVelocity, parameters, localMeshsize)                          - duwdx        (localVelocity, parameters, localMeshsize)
                    - dvwdy        (localVelocity, parameters, localMeshsize)
                    + parameters.environment.gz );
}

inline FLOAT computeSdotS3D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize){
return 
	(
		+dudx(localVelocity, localMeshsize)*dudx(localVelocity, localMeshsize)
		+dvdy(localVelocity, localMeshsize)*dvdy(localVelocity, localMeshsize)
		+dwdz(localVelocity, localMeshsize)*dwdz(localVelocity, localMeshsize)
	)
	+0.5*(
		+(dudy(localVelocity, localMeshsize)+dvdx(localVelocity, localMeshsize))*(dudy(localVelocity, localMeshsize)+dvdx(localVelocity, localMeshsize))
		+(dudz(localVelocity, localMeshsize)+dwdx(localVelocity, localMeshsize))*(dudz(localVelocity, localMeshsize)+dwdx(localVelocity, localMeshsize))
		+(dvdz(localVelocity, localMeshsize)+dwdy(localVelocity, localMeshsize))*(dvdz(localVelocity, localMeshsize)+dwdy(localVelocity, localMeshsize))
	);
}

inline FLOAT computeSdotS2D(const FLOAT * const localVelocity, const FLOAT * const localMeshsize){
return 
	(
		+dudx(localVelocity, localMeshsize)*dudx(localVelocity, localMeshsize)
		+dvdy(localVelocity, localMeshsize)*dvdy(localVelocity, localMeshsize))
	+0.5*(
		(dudy(localVelocity, localMeshsize)+dvdx(localVelocity, localMeshsize))*(dudy(localVelocity, localMeshsize)+dvdx(localVelocity, localMeshsize)));
}



#endif
