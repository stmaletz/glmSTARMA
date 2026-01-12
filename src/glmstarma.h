/* 
-----------------------------------------------------------------------------
    File: glmstarma.h
    Purpose: Main header file for glmstarma package
    Author: Steffen Maletz
    Last modified: 2025-12-06
-----------------------------------------------------------------------------
*/
#ifndef GLMSTARMA_H
#define GLMSTARMA_H

#include <cmath>
#include <iostream>
#include <iterator>
#include <list>
#include <RcppArmadillo.h>

class Family;
class Covariate;
class Orders;
class Neighborhood;
class CovariateList;

#include "family.h"
#include "neighbors.h"
#include "covariates.h"
#include "orders.h"
#include "design.h"
#include "parameter.h"
#include "model.h"
#include "fitting.h"

#endif
