#include <QDebug>
#include "generalglm.h"

QString GeneralGLM::getConditionString()
{
    QString string="";
    for (int jString=0; jString<_conditionList.size(); jString++)
    {
        string.append(_conditionList.at(jString));
        string.append(" ");
    }
    return string;
}

double GeneralGLM::getWeight(int iTime)
{
    if ( _nTime == 0 || iTime >= _nTime )
        return 1.;
    else
        return _weight_t[iTime];
};


void GeneralGLM::calculateCovarianceMatrix()
{ // call this after conditions have been defined
    for (int jCondition=0; jCondition<_conditionList.size(); jCondition++)
    {
        // Define the F convariance matrix for each condition.
        _fCovar[jCondition].resize(_nContrasts[jCondition]);
        _fCovarInv[jCondition].resize(_nContrasts[jCondition]);
        for (int jc1=0; jc1<_nContrasts[jCondition]; jc1++)
        {
            _fCovar[jCondition][jc1].fill(0.,_nContrasts[jCondition]);
            _fCovarInv[jCondition][jc1].fill(0.,_nContrasts[jCondition]);
            for (int jc2=0; jc2<_nContrasts[jCondition]; jc2++)
            {
                _fCovar[jCondition][jc1][jc2] = 0.;
                for (int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
                {
                    for (int jCoeff1=0; jCoeff1<_nCoeff; jCoeff1++)
                        _fCovar[jCondition][jc1][jc2]       += _contrastMatrix[jCondition][jc1][jCoeff]
                                * _XTWXm1_cc[jCoeff][jCoeff1] * _contrastMatrix[jCondition][jc2][jCoeff1];
                }
            }
        }
        // Make a copy for the inverse
        bool allZero=true;
        for (int jc1=0; jc1<_nContrasts[jCondition]; jc1++)
        {
            for (int jc2=0; jc2<_nContrasts[jCondition]; jc2++)
            {
                _fCovarInv[jCondition][jc1][jc2] = _fCovar[jCondition][jc1][jc2];
                allZero &= _fCovar[jCondition][jc1][jc2] == 0.;
            }
        }
        if ( allZero ) return; // don't bother trying to invert a zero matrix
        // Now invert the matrix.
        if ( ! utilMath::dInvertSquareMatrix(_fCovarInv[jCondition]) ) return;
    }
}

int GeneralGLM::getNumberCoefficientsInCondition(int iCondition)
{ // find the number of coefficients in the condition AND assign a vector matching them to the basis coefficient list
    int nCoefficientsInCondition = 0;
    // _nContrasts[iCondition] should be 1 contrast for a T test, but more for an F test
    for ( int jContrast=0; jContrast<_nContrasts[iCondition]; jContrast++ )
    {
        for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++ )
        {
            if ( _contrastMatrix[iCondition][jContrast][jCoeff] != 0. )
                nCoefficientsInCondition++;
        }
    }
    return nCoefficientsInCondition;
}

int GeneralGLM::getOnlyCoefficientInCondition(int iCondition)
{
    int iCoeff = -1;
    if ( getNumberCoefficientsInCondition(iCondition) == 1 )
    {
        // should be 1 contrast for a T test, but more for an F test
        for ( int jContrast=0; jContrast<_nContrasts[iCondition]; jContrast++ )
        {
            for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++ )
            {
                if ( _contrastMatrix[iCondition][jContrast][jCoeff] != 0. )
                    iCoeff = jCoeff;
            }
        }
    }
    return iCoeff;
}

int GeneralGLM::getFirstCoefficientInCondition(int iCondition)
{
    int iCoeff = -1;
    // should be 1 contrast for a T test, but more for an F test
    for ( int jContrast=0; jContrast<_nContrasts[iCondition] || iCoeff >=0; jContrast++ )
    {
        for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++ )
        {
            if ( _contrastMatrix[iCondition][jContrast][jCoeff] != 0. )
            {
                iCoeff = jCoeff;
                break;
            }
        }
    }
    return iCoeff;
}

void GeneralGLM::setNumberContrasts(int iCondition, int nContrasts )
{
    // This function needs to be called if either the condition changes or the # ncoefficients changes (call through defineConditions)
    if ( iCondition >= _conditionList.size() )
    {
        decodeConditions(false);  // this excludes bad conditions in the string and reset the string to valid conditions
        decodeConditions(true);
    }
    _nContrasts[iCondition] = nContrasts;
    _contrastMatrix[iCondition].resize(_nContrasts[iCondition]);
    for (int jContrast=0; jContrast<_nContrasts[iCondition]; jContrast++)
        _contrastMatrix[iCondition][jContrast].fill(0.,_nCoeff);
}

void GeneralGLM::evaluateCurrentCondition()
{
    int currentCondition = getCurrentCondition();
    if ( getNumberConditions() != 0 && currentCondition >= 0 )
        calculateFStat(currentCondition);
}

void GeneralGLM::calculateFStat (int lCondition)
{
    if ( _sigma2 == 0. )
    {
        _cnr = _effectSize = _effectStDev = 0.;
        _pRaw = 1.;
        return;
    }

    // Make sure the condition is within range.
    if ( lCondition < 0 || lCondition >= _conditionList.size() )
    {
        decodeConditions(false);  // this excludes bad conditions in the string and reset the string to valid conditions
        decodeConditions(true);
    }

    dVector CBeta;

    CBeta.fill(0.,_nContrasts[lCondition]);
    int lContrastMax = 0;  double CBetaMax  = 0.;
    // Calculate the product CBeta
    for (int jc1=0; jc1<_nContrasts[lCondition]; jc1++)
    {
        for (int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
            CBeta[jc1] += _contrastMatrix[lCondition][jc1][jCoeff] * _beta_c[jCoeff];
        if ( qAbs(CBeta[jc1]) > qAbs(CBetaMax) )
        {
            CBetaMax     = CBeta[jc1];
            lContrastMax = jc1;
        }
    }
    // Determine the effect size and standard deviation.
    _effectSize = CBetaMax;
    _effectStDev  = qSqrt(_fCovar[lCondition][lContrastMax][lContrastMax] * _sigma2);

    // Compute the numerator for F
    double numerator = 0.;
    for (int jc1=0; jc1<_nContrasts[lCondition]; jc1++)
    {
        for (int jc2=0; jc2<_nContrasts[lCondition]; jc2++)
            numerator += CBeta[jc1] * _fCovarInv[lCondition][jc1][jc2] * CBeta[jc2];
    }

    // Degrees of freedom
    double dof1 = _dof;
    double dof2 = _nContrasts[lCondition];
    // potentially override _dof by using an effective DOF due to regularization/smoothing
    if ( _dofEffective != 0. ) dof1 = _dofEffective;
    ///////////////////////////////////////////////
    // The covariance matrix, _fCovarInv, is formed from (_XTWXm1_cc)^-1
    // As such, _sigma2 can be included 2 ways:
    // 1) as long as the weights are normalized, _sigma2 can be determined post-hoc from data vs. fit
    // 2) in some cases (e.g., 2nd order GLM), weights are known explicitly (1/sigma2), so set _sigma2=1 in formula below
    // Note that in any case,  "weights/_sigma2" should = 1/sigma2_actual, so weight normalization is NOT arbitrary
    ///////////////////////////////////////////////
    _fStatistic = numerator / _sigma2 / dof2;
    _cnr        = qSqrt(_fStatistic);
    if ( CBetaMax < 0. ) _cnr *= -1.;

    double x = dof1/(dof1+dof2*_fStatistic);
    if (x < 0. || x > 1.0) x = 1.;
    double beta = utilMath::ibeta(0.5*dof1,0.5*dof2,x); // cephes
    if ( dof2 == 1. )
        _pRaw = beta;
    else
    {
        _pRaw = 2.*beta;
        if ( _pRaw == 2. ) _pRaw = 1.;
        if ( _pRaw > 1.0 ) _pRaw = 2. - _pRaw;
    }

    double tooSmall = exp(-100.);
    if ( _pRaw < tooSmall ) _pRaw = tooSmall;

    return;
}

void GeneralGLM::removeCoefficient(int iCoeff, dVector &data)
{
    removeCoefficient(iCoeff, data, _fit_t);
}

void GeneralGLM::removeCoefficient(int iCoeff, dVector &data, dVector &fit)
{
    if ( _nTime != data.size() )
    {
        qWarning() << "Error: the # time points (" << data.size() << ") does not match the GLM (" << _nTime << ").";
        exit(1);
    }
    for (int jt=0; jt<_nTime; jt++)
    {
        // remove coefficients
        fit[jt] -= _X_tc[jt][iCoeff] * _beta_c[iCoeff];
        data[jt]-= _X_tc[jt][iCoeff] * _beta_c[iCoeff];
    }
}

void GeneralGLM::orthogonalize(int iCoeff, dVector &data)
{ // orthogonalize data wrt basis[iCoeff]
//    return;
    double dotProduct = 0.;  double mag2 = 0.;
    for (int jt=0; jt<_nTime; jt++)
    {
        dotProduct += data[jt] * _X_tc[jt][iCoeff];
        mag2       += _X_tc[jt][iCoeff] * _X_tc[jt][iCoeff];
    }
    // Subtract the portion of the input vector that is parallel to the given coefficient vector
    for (int jt=0; jt<_nTime; jt++)
        data[jt] -= dotProduct / mag2 * _X_tc[jt][iCoeff];
}

void GeneralGLM::orthogonalizeExistingBasisFunction(int iCoeff, int iCoeffMod)
{ // orthogonalize basis[iCoeffMod] wrt coefficient iCoeff
//    return;
    if ( iCoeff == iCoeffMod ) return;  // do nothing if they are the same coefficient
    double dotProduct = 0.;  double mag2 = 0.;
    for (int jt=0; jt<_nTime; jt++)
    {
        dotProduct += _X_tc[jt][iCoeffMod] * _X_tc[jt][iCoeff];
        mag2       += _X_tc[jt][iCoeff]    * _X_tc[jt][iCoeff];
    }
    // Subtract the portion of the input vector that is parallel to the given coefficient vector
    for (int jt=0; jt<_nTime; jt++)
        _X_tc[jt][iCoeffMod] -= dotProduct / mag2 * _X_tc[jt][iCoeff];
}

void GeneralGLM::setOLS()
{
    if ( _nTime == 0 )
        qFatal("Error: GeneralGLM not initialized (nTime=0) in function setOLS");
    _weight_t.fill(1.,_nTime);
}


void GeneralGLM::init(int ntime , int nCoeff)
{ // this should be called first to initialize the GLM

    if ( ntime != _nTime )
    {
        _nTime  = _nIncluded = ntime;
        _fit_t.fill(0.,_nTime);
        _fitErr_t.fill(0.,_nTime);
        _X_tc.resize(_nTime);
    }

    _contrastMatrix.resize(MaxConditions);
    _fCovar.resize(MaxConditions);
    _fCovarInv.resize(MaxConditions);

    if ( nCoeff != _nCoeff )
    {
        _nCoeff = nCoeff;
        _XTWXm1_cc.resize(_nCoeff);
        _XTWXm1XTW_ct.resize(_nCoeff);
        for (int jc=0; jc<_nCoeff; jc++)
        {
            _XTWXm1_cc[jc].fill(0.,_nCoeff);
            _XTWXm1XTW_ct[jc].fill(0.,_nTime);
        }
    }

    _dof = _nIncluded - _nCoeff;   // this should be changed if time points are ignored
    for ( int jt=0; jt<_nTime; jt++ )
        _X_tc[jt].resize(_nCoeff);

    _initialized = true;
}

void GeneralGLM::setVariances( dVector variance_t)
{
    int nTime = variance_t.size();
    dVector weight_t; weight_t.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        if ( variance_t[jt] > 0. )
            weight_t[jt] = 1. / variance_t[jt];
        else
            weight_t[jt] = 1.;
    }
    /*
    // Normalize weights
    double average=0.;
    for (int jt=0; jt<nTime; jt++)
        average += weight_t[jt];
    average /= static_cast<double>(nTime);
    for (int jt=0; jt<nTime; jt++)
        weight_t[jt] /= average;
        */
    setWeights(weight_t);
}

void GeneralGLM::setWeights( dVector weight_t)
{
    if ( weight_t.size() != _nTime )
    {
        qWarning() << "Error: weight size " << weight_t.size() << " not = time size " << _nTime;
        exit(1);
    }
    _nIncluded = _nTime;
    _weight_t.resize(_nTime);
    for ( int jt=0; jt<_nTime; jt++ )
    {
        _weight_t[jt] = weight_t[jt];
        if ( _weight_t[jt] == 0. ) _nIncluded--;
    }
    _dof = _nIncluded - _nCoeff;
    _prepared = false;
}

bool GeneralGLM::basisFunctionsChanged()
{  // This is a tricky function: calling it once changes the state, even in a debug/print statement, so take care
    if ( _basisFunctionsChanged )
    {
        _basisFunctionsChanged=false;  // reset for next time
        return true;
    }
    return false;
}

void GeneralGLM::replaceBasisFunction( int iBasis, dVector X_t)
{
    if ( X_t.size() != _nTime )
    {
        qWarning() << "Error: time dimension " << X_t.size() << " not equal to GLM time dimension " << _nTime;
        exit(1);
    }
    else if ( iBasis < 0 || iBasis >= _nCoeff )
    {
        qWarning() << "Error: basis function index " << iBasis << " out of range in replaceBasisFunction";
        exit(1);
    }
    for ( int jt=0; jt<_nTime; jt++ )
    {
        _X_tc[jt][iBasis] = X_t[jt];
//        if ( _weight_t[jt] == 0. ) _X_tc[jt][iBasis] = 0.;
    }
    _basisFunctionsChanged=true;
    _prepared = false;
}

void GeneralGLM::addBasisFunction( dVector X_t)
{
    if ( X_t.size() != _nTime )
    {
        qWarning() << "Error: time dimension " << X_t.size() << " not equal to GLM time dimension " << _nTime;
        exit(1);
    }
    else
    {
        int nTime = X_t.size();
        init(nTime,_nCoeff+1);
        for (int jt=0; jt<X_t.size(); jt++)
        {
            _X_tc[jt][_nCoeff-1] = X_t[jt];
//            if ( _weight_t[jt] == 0. ) _X_tc[jt][_nCoeff-1] = 0.;
        }
    }
    _basisFunctionsChanged=true;
    _prepared = false;
}

void GeneralGLM::addOrInsertBasisFunction( int iBasis, dVector X_t)
{
    if ( X_t.size() != _nTime )
    {
        qWarning() << "Error: time dimension " << X_t.size() << " not equal to GLM time dimension " << _nTime;
        exit(1);
    }
    else if ( iBasis < 0 || iBasis > _nCoeff )
    {
        qWarning() << "Error: basis function index " << iBasis << " out of range in addBasisFunction";
        exit(1);
    }
    else if ( iBasis == _nCoeff )
    { // add basis function
        int nTime = X_t.size();
        init(nTime,_nCoeff+1);
    }
    for (int jt=0; jt<X_t.size(); jt++)
    {
        _X_tc[jt][iBasis] = X_t[jt];
//        if ( _weight_t[jt] == 0. ) _X_tc[jt][iBasis] = 0.;
    }
    _basisFunctionsChanged=true;
    _prepared = false;
}
void GeneralGLM::fitWLS(dVector data, bool computeSigma2)
{ // computeSigma2: if true, find sigma2 from residue; else assume sigma2 is known and used in weighting matrix
    // If any sigma values are zero, ignore the voxel.
//    FUNC_ENTER << _nCoeff;
    bool allZeroes=true;
    for (int jt=0; jt<_nTime; jt++)
        allZeroes &= (data[jt] == 0.);
    if ( allZeroes )
    {
        // Set everything to zero
        _average = _sigma2 = 0.;
        _beta_c.fill(0.,_nCoeff);   _sem_c.fill(0.,_nCoeff);
        _fit_t.fill(0.,_nTime);     _fitErr_t.fill(0.,_nTime);
        return;
    }

    if ( !_prepared )
        calculatePseudoInverse();
    // Calculate Beta_c = (X_transpose * W * X)^-1 * X_transpose * W * Y
    for (int jc=0; jc<_nCoeff; jc++)
    {
        _beta_c[jc] = 0.;
        // Calculate Beta, summing over time
        for (int jt=0; jt<_nTime; jt++)
            _beta_c[jc] += _XTWXm1XTW_ct[jc][jt] * data[jt];
    }
    // Compute the fit: Y_t = X_tc * beta_c
    for (int jt=0; jt<_nTime; jt++)
    {
        _fit_t[jt] = _fitErr_t[jt] = 0.;
        for (int jc=0; jc<_nCoeff; jc++)
            _fit_t[jt] += _X_tc[jt][jc] * _beta_c[jc];
    }

    if ( computeSigma2 )
    {
        // Compute sigma
        _sigma2 = _average = 0;
        for (int jt=0; jt<_nTime; jt++)
        {
            if ( _weight_t[jt] != 0. )
            {
                _sigma2 += (data[jt]-_fit_t[jt]) * (data[jt]-_fit_t[jt]);
                _average += data[jt];
            }
        }
        _sigma2  /= _dof;
        _average /= _dof;
        if ( _nTime == _nCoeff )
        {
            _sigma2 = 0.;
            _sem_c.fill(0.,_nCoeff);
        }
        else
        {
            // Compute the SEM for each coefficient
            for (int jc=0; jc<_nCoeff; jc++)
                _sem_c[jc] = sqrt(_sigma2*_XTWXm1_cc[jc][jc]);
        }
    }
    else
    {
        // Set _sigma2=1 in GeneralGLM. This variable is only used in post-hoc determination of variance.
        // In 2nd-order GLM, _sigma2 is known/set explicitly in the weighting matrix.
        _sigma2 = 1.;
        if ( _nTime == _nCoeff )
            _sem_c.fill(0.,_nCoeff);
        else
        {
            // Compute the SEM for each coefficient
            for (int jc=0; jc<_nCoeff; jc++)
                _sem_c[jc] = sqrt(_XTWXm1_cc[jc][jc]);
        }
    }
    // Compute the fit error
    _fitErr_t.fill(0.,_nTime);
    for (int jt=0; jt<_nTime; jt++)
    {
        for (int jc=0; jc<_nCoeff; jc++)
            _fitErr_t[jt] += _X_tc[jt][jc] * _sem_c[jc];
    }
//    FUNC_EXIT;
}

void GeneralGLM::defineConditions(QString conditionString )
{ // must be called AFTER GLM_define_basis_functions()
    //  DefineContrastMatrix:
    //    - Generally will be true,
    //    - except when one wants to interpret the conditionString to create a list of files without knowing about basis functions
    //  Set up the conditional tests based upon "condition_string".
    //  This string can contain more than 1 condition.
    //  Example: "1 2 3 4 1234 12-34 12-33 1-3,2-4
    //  yields 7 conditions, including 7 T tests:
    //    1-4) event 1 versus baseline, 2 versus baseline, ...
    //    5) event 1+2+3+4 versus baseline
    //    6) event 1+2-3-4 versus baseline
    //    7) event 1+2-3-3 versus baseline (e.g., [1 1 -2 ..] in contrast vector)
    //  ... and 1 F test
    //    1) events 1-3 or 2-4

    _conditionList = conditionString.split(QRegExp("[ ]"), QString::SkipEmptyParts);
    decodeConditions(false);  // this excludes bad conditions in the string and reset the string to valid conditions
    decodeConditions(true);
}

void GeneralGLM::decodeConditions( bool defineMatrices )
{ // generally call this function 2ce: 1) to exclude bad conditions with argument false, then define conditions with argument true
    if ( _conditionList.isEmpty() )
        return;
    else if ( defineMatrices )
    {
        _nContrasts.fill(1,_conditionList.size());
        _conditionShape.resize(_conditionList.size());
        _indexCoeffInCondition.resize(_conditionList.size());
    }

    QStringList validConditionList;
    int nConditions = _conditionList.size();
    for (int jCondition=0; jCondition<nConditions; jCondition++)
    {
        QString condition = _conditionList.at(jCondition);
        // allocate contrast matrix and set to 0
        if ( defineMatrices )
        {
            _indexCoeffInCondition[jCondition].resize(0);
            setNumberContrasts(jCondition,condition.count(',') + 1); // # of contrasts = # commas + 1 in the condition string
        }
        // Go through the condition string and decode it
        bool allGoodCoefficientsInCondition = true;
        int iContrast=0;
        int sign=1;
        for ( int jChar=0; jChar<condition.size(); jChar++)
        {
            QChar eventID = condition.at(jChar);
            if ( eventID == ',' )
            {
                iContrast++;
                sign = 1;
            }
            else if ( eventID == '-')
                sign = -1;
            else
            {
                iVector iCoeffVector = getEventCoefficients(eventID);
                for (int jCoeffInVector=0; jCoeffInVector<iCoeffVector.size(); jCoeffInVector++)
                {
                    allGoodCoefficientsInCondition &= (iCoeffVector[jCoeffInVector] >= 0);
                    if ( iCoeffVector[jCoeffInVector] >= 0 && defineMatrices )
                    {
                        _indexCoeffInCondition[jCondition].append(iCoeffVector[jCoeffInVector]);
                        _conditionShape[jCondition] = _basisShape[iCoeffVector[jCoeffInVector]];
                        if ( _conditionShape[jCondition] == Shape_FIR || _conditionShape[jCondition] == Shape_FR )
                        {
                            if ( _FIRContrast[jCoeffInVector] == 1 )
                                _contrastMatrix[jCondition][iContrast][iCoeffVector[jCoeffInVector]] += sign;
                        }
                        else if ( _conditionShape[jCondition] == Shape_RSFIR )
                            {
                                if ( _RSFIRContrast[jCoeffInVector] != 0. )
                                    _contrastMatrix[jCondition][iContrast][iCoeffVector[jCoeffInVector]] += sign * _RSFIRContrast[jCoeffInVector];
                            }
                        else
                            _contrastMatrix[jCondition][iContrast][iCoeffVector[jCoeffInVector]] += sign;
                    }
                } // jCoeffInVector
            } // apply sign
        } // jChar
        if ( allGoodCoefficientsInCondition && ! validConditionList.contains(condition) ) validConditionList.append(condition);
    }
    // reset the condition list
    _conditionList = validConditionList;

    // Go back through the contrast matrix and make the sum of all positives = 1 within a contrast (and same for the negatives)
    // Thus, the contrast ab-c = "(a+b)-c") becomes "(a+b)/2 - c"
    if ( defineMatrices )
    {
        for ( int jCondition=0; jCondition<nConditions; jCondition++)
        {
            if ( _conditionShape[jCondition] != Shape_RSFIR )
            {
                for ( int jContrast=0; jContrast<_nContrasts[jCondition]; jContrast++ )
                {
                    double sumPos=0;  double sumNeg=0;
                    for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
                    {
                        double cValue = _contrastMatrix[jCondition][jContrast][jCoeff];
                        if ( cValue > 0. ) sumPos += cValue;
                        if ( cValue < 0. ) sumNeg += qAbs(cValue);
                    }
                    if ( sumPos != 0. )
                    {
                        for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
                        {
                            double cValue = _contrastMatrix[jCondition][jContrast][jCoeff];
                            if ( cValue > 0. ) _contrastMatrix[jCondition][jContrast][jCoeff] /= sumPos;
                        }
                    }
                    if ( sumNeg != 0. )
                    {
                        for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
                        {
                            double cValue = _contrastMatrix[jCondition][jContrast][jCoeff];
                            if ( cValue < 0. ) _contrastMatrix[jCondition][jContrast][jCoeff] /= sumNeg;
                        }
                    }
                } // jContrast
            } // ! Shape_RSFIR
        } // jContrast
    } // jCondition

    // The conditions have been defined, so calculate the covariance matrix.
    if ( defineMatrices ) calculateCovarianceMatrix();
}

QChar GeneralGLM::getEventChar(int indexEvent)
{
    // indexEvent   char
    //  0-8         1-9
    //  9-34        a-z
    // 35-60        A-Z

    int iAscii;
    if ( indexEvent < 0)
        iAscii = 32;               // should only occur for erroneous entry
    else if ( indexEvent < 9 )
        iAscii = indexEvent + 49;  // event index 0 = '1' = ascii 49
    else if ( indexEvent < 35)
        iAscii = indexEvent + 88;  // event index 9 = 'a' = ascii 97
    else
        iAscii = indexEvent + 30;  // event index 35 = 'A' = ascii 65
    QChar qAscii = iAscii;
    return qAscii;
}

int GeneralGLM::getEventIndex(QChar eventID)
{
    int iAscii = eventID.toLatin1();
    int indexEvent;
    if ( iAscii >= 49 && iAscii <= 57 )
        indexEvent = iAscii - 49;  // '1' is 49 and goes into index 0, 9->8
    else if ( iAscii >= 97 && iAscii <= 122 )
        indexEvent = iAscii - 88;  // 'a' is 97 and goes into index 9
    else if ( iAscii >= 65 && iAscii <= 90 )
        indexEvent = iAscii - 30;  // 'A' is 65 and goes into index 35
    else
        indexEvent=-1; // invalid ID
    return indexEvent;
}

int GeneralGLM::getEventCoefficient(QChar eventID)
{
    for (int jCoeff=0; jCoeff<_basisID.size(); jCoeff++)
    {
        if ( _basisID[jCoeff] == eventID )
            return(jCoeff);
    }
    return -1;
}
int GeneralGLM::getEventCoefficient(QChar eventID, int iOccurrence)
{
    iVector eventCoefficients = getEventCoefficients(eventID);
    if ( iOccurrence >= eventCoefficients.size() )
        return eventCoefficients[0];
    else
        return eventCoefficients[iOccurrence];
}
iVector GeneralGLM::getEventCoefficients(QChar eventID )
{
    int nCoeffPerEvent=0;
    for (int jCoeff=0; jCoeff<_basisID.size(); jCoeff++)
        if ( _basisID[jCoeff] == eventID ) nCoeffPerEvent++;
    iVector eventCoefficients;
    if ( nCoeffPerEvent == 0 )
    {
        eventCoefficients.resize(1);
        eventCoefficients[0]=-1;
    }
    else
    {
        eventCoefficients.resize(nCoeffPerEvent);
        nCoeffPerEvent=0;
        for (int jCoeff=0; jCoeff<_basisID.size(); jCoeff++)
        {
            if ( _basisID[jCoeff] == eventID )
            {
                eventCoefficients[nCoeffPerEvent] = jCoeff;
                nCoeffPerEvent++;
            }
        }
    }
    return eventCoefficients;
}

void GeneralGLM::calculatePseudoInverse()
{
    // Zero all data
    _beta_c.fill(0.,_nCoeff);
    _sem_c.fill(0.,_nCoeff);
    _XTWXm1_cc.resize(_nCoeff);
    _XTWXm1XTW_ct.resize(_nCoeff);
    for (int jc=0; jc<_nCoeff; jc++)
    {
        _XTWXm1_cc[jc].fill(0.,_nCoeff);
        _XTWXm1XTW_ct[jc].fill(0.,_nTime);
    }

    // Calculate (X_transpose * W * X) by summing over the time index.
    for (int jc=0; jc<_nCoeff; jc++)
    {
        for (int jc1=0; jc1<_nCoeff; jc1++)
        {
            for (int jt=0; jt<_nTime; jt++)
                _XTWXm1_cc[jc][jc1] += _X_tc[jt][jc] * _weight_t[jt] * _X_tc[jt][jc1];
        }
    }

    // Invert to get (X_transpose * W * X)^-1
    if ( ! utilMath::dInvertSquareMatrix(_XTWXm1_cc) ) return;

    // Calculate (X_transpose * W * X)^-1 * X_transpose
    for (int jc=0; jc<_nCoeff; jc++)
    {
        for (int jt=0; jt<_nTime; jt++)
        {
            for (int jc1=0; jc1<_nCoeff; jc1++)
                _XTWXm1XTW_ct[jc][jt] += _XTWXm1_cc[jc][jc1] * _X_tc[jt][jc1]* _weight_t[jt];
        }
    }
    _prepared = true;
}

void PolynomialGLM::define(int nPoly, int nCoeff, int nTime)
{
//    FUNC_ENTER << nPoly << nCoeff << nTime;
    dVector xVector;  xVector.resize(nTime);
    double duration = qMax(1.,static_cast<double>(nTime-1));
    // convert to a vector on the range -1,1
    for ( int jt=0; jt<nTime; jt++ )
        xVector[jt] = 2.*static_cast<double>(jt)/duration - 1.;  // range=(-1,1) with center point at x=0 or x=1/(n-1)
    define(nPoly, nCoeff, xVector);
}

void PolynomialGLM::define(int nPoly, int nCoeff, dVector xVector)
{
    int iCenter = _xInputVector.size()/2;
    define(nPoly, nCoeff, xVector, iCenter);
}

void PolynomialGLM::define(int nPoly, int nCoeff, dVector xVector, int iCenter)
{ // This defines a simple polynomial GLM for fitting a single scan
    _nPoly = nPoly;  // usually this is the same as _nCoeff in generalGLM, but it could differ if additional basis functions are added
    _xInputVector = xVector;
    _centralPoint = iCenter;
    int nTime = _xInputVector.size();
    if ( nTime < nCoeff )
      {
        qInfo() << "Error: the # of points (" << nTime << ") must be >= # coefficients" << nCoeff;
        return;
      }
    else if ( nTime < nCoeff )
        nCoeff = nTime;
    nCoeff = qMin(nCoeff,6);
    // The following sets a flag; failure above will leave the flag unset
    init(nTime,nCoeff);

    // Create the basis functions.
    dVector X_t;
    X_t.resize(nTime);
    for (int jPoly=0; jPoly<nCoeff; jPoly++)
    {
        for ( int jt=0; jt<nTime; jt++ )
            X_t[jt] = baselineBasisFunction(jPoly,getNormalizedTime(jt));
//        FUNC_INFO << "jPoly" << jPoly << "X_t" << X_t;
        addOrInsertBasisFunction(jPoly,X_t);
    }
    // set uniform weights as default
    dVector weights;  weights.fill(1.,nTime);
    setWeights(weights);
    calculatePseudoInverse();
}
double PolynomialGLM::getNormalizedTime(double time)
{ // convert time (via the index into _xInputVector) to the range (-1,1)
    double minusOneToOne = 0.;
    if ( _centralPoint < 0 )
    {
        double duration = _xInputVector[_xInputVector.size()-1] - _xInputVector[0];  // duration = max - min
        if ( duration != 0. )
        {
            double scale = 2./duration;  double offset = 2.*_xInputVector[0]/duration - 1.;
            minusOneToOne = scale * time + offset;
        }
    }
    else
    { // LOESS; generally #pts should be odd, so there is a central point
        double x0 = _xInputVector[_centralPoint]; double xLow = _xInputVector[0]; double xHigh = _xInputVector[_xInputVector.size()-1];
        double maxHalfWidth = qMax(qAbs(x0-xLow),qAbs(xHigh-x0));
        if ( maxHalfWidth == 0.)
            minusOneToOne = 0.;
        else
            minusOneToOne = (time-x0)/maxHalfWidth;
//        FUNC_INFO << "(x0, xLow, xHigh), maxHalfWidth minusOneToOne" << "(" << xLow << x0 << xHigh << ")" << maxHalfWidth << time << minusOneToOne;
    }
    return minusOneToOne;
}
double PolynomialGLM::baselineBasisFunction(int iPoly, double x)
{ // return Legendre polynomial of x with order iPoly; use x in the range (-1,1)
    double value;
    if ( iPoly == 0 )
        value = 1.;
    else if ( iPoly == 1 )
        value = x;
    else if ( iPoly == 2 )
        value = 0.5 * (3*x*x -1.);
    else if ( iPoly == 3 )
        value = 0.5 * (5.*x*x*x - 3.*x);
    else if ( iPoly == 4 )
        value = 0.125 * (35.*x*x*x*x - 30.*x*x + 3.);
    else // if ( iPoly == 5 )
        value = 0.125 * (63.*x*x*x*x*x - 70.*x*x*x + 15.*x);
    return value;
}
double PolynomialGLM::getBasisDerivative(int iPoly, double x)
{ // return Legendre polynomial of x with order iPoly
    double value;
    if ( iPoly == 0 )
        value = 0.;
    else if ( iPoly == 1 )
        value = 1.;
    else if ( iPoly == 2 )
        value = 3*x;
    else if ( iPoly == 3 )
        value = 1.5 * (5.*x*x - 1.);
    else if ( iPoly == 4 )
        value = 2.5 * (7.*x*x*x - 3.*x);
    else // if ( iPoly == 5 )
        value = 1.875 * (21.*x*x*x*x - 14.*x*x + 1.);
    return value;
}
double PolynomialGLM::getFitInterpolation(double time)
{ // if _nCoeff > _nPoly, the extra terms will need to be computed outside of this class based upon the extra terms parametric definition
    double x = getNormalizedTime(time);
    FUNC_ENTER << _nPoly << _nCoeff << time << x;
    double fitValue=0.;
    for (int jPoly=0; jPoly<_nPoly; jPoly++)
        fitValue += baselineBasisFunction(jPoly, x) * getBeta(jPoly);
    FUNC_EXIT << fitValue;
    return fitValue;
}

double PolynomialGLM::getDerivative(int iTime)
{ // dBasis/dt = dBasis/dt' * dt'/dt
    double derivative = 0.;
    for (int jPoly=0; jPoly<_nCoeff; jPoly++)
    {
        double polyDeriv = getBasisDerivative(jPoly,getNormalizedTime(iTime));
        derivative += polyDeriv;
    }
    double duration = _xInputVector[_xInputVector.size()-1] - _xInputVector[0];  // duration = max - min
    return derivative * 2./duration;  // dBasis/dt = dBasis/dt' * dt'/dt
}

double PolynomialGLM::getDerivative(double time)
{ // dBasis/dt = dBasis/dt' * dt'/dt
    double derivative = 0.;
    for (int jPoly=0; jPoly<_nCoeff; jPoly++)
    {
        double polyDeriv = getBasisDerivative(jPoly,getNormalizedTime(time));
        derivative += polyDeriv;
    }
    double duration = _xInputVector[_xInputVector.size()-1] - _xInputVector[0];  // duration = max - min
    return derivative * 2./duration;  // dBasis/dt = dBasis/dt' * dt'/dt
}

void LOESS::define(dVector timeVector, double smoothingScaleMin)
{
    FUNC_ENTER << timeVector;
    FUNC_INFO << smoothingScaleMin;
    _xPoints = timeVector;
    _smoothingScaleMin = smoothingScaleMin;
    // Define local quadratic smoothing for each this run: (one function for each time point)
    int nTime = timeVector.size();
    _polyLOESS.resize(nTime);
    for ( int jt=0; jt<nTime; jt++)
    {
        FUNC_INFO << "jt" << jt << "\n";
        double time0 = timeVector[jt];
        int iLowSide=0;  double time = time0;  double width=0.;
        while ( jt - iLowSide >= 0 && width < _smoothingScaleMin )
        {
            time = timeVector[jt-iLowSide];
            width = qAbs(time-time0);
            iLowSide++;
        }
        iLowSide--;
        int iHighSide=0;  time = time0;  width = 0.;
        while ( jt + iHighSide < nTime && width < _smoothingScaleMin )
        {
            time = timeVector[jt+iHighSide];
            width = qAbs(time-time0);
            iHighSide++;
        }
        iHighSide--;
        // Set weights
        dVector weight;  weight.clear();
        int iHalfWidth = qMax(1,qMin(iLowSide,iHighSide));
        dVector localTime;
        int iLocal=0;
        int iCenter=0;
        for ( int jRel=-iLowSide; jRel<=iHighSide; jRel++)
        {
            int iTime = jt + jRel;
            if ( iTime >= 0 && iTime < nTime )
            {
                time = timeVector[iTime];
                double distance = qAbs(time-time0) / _smoothingScaleMin;
                if ( distance < 0.99 )
                {
                    localTime.append(time-time0);
                    weight.append(qPow(1-qPow(distance,3),3));  // use traditional tri-cube weight function
                    if ( jRel == 0 ) iCenter = iLocal;
                    iLocal++;
                }
                FUNC_INFO << "jRel, iTime, time, time0" << jRel << iTime << time << time0;
                FUNC_INFO << "_smoothingScaleMin, distance, weight" << _smoothingScaleMin << distance << weight.last();
            }
        }
        FUNC_INFO << "localTime" << localTime;
        FUNC_INFO << "define poly[" << jt << "] using central point" << iCenter << "and weights" << weight;
        //             define       (nCoeff,          nTime)
        _polyLOESS[jt].define(qMin(weight.size(),3),localTime,iCenter);  // qMax: don't define 3 terms with only 2 points
//        _polyLOESS[jt].define(qMin(weight.size(),3),weight.size());  // qMax: don't define 3 terms with only 2 points
        FUNC_EXIT << "poly[" << jt << "] = " << _polyLOESS[jt].getNumberTimePoints();
        _polyLOESS[jt].setWeights(weight);
    }
}
void LOESS::defineAndFit(dVector timeVector, dVector yVector, double smoothingScale, bool asymmetric, bool linkWidthToTime)
{
//    FUNC_ENTER << timeVector;
//    FUNC_INFO << smoothingScale;
    _xPoints = timeVector;
    _smoothingScaleMin    = smoothingScale;
    double smoothingWidth = smoothingScale;
    // Define local quadratic smoothing for each this run: (one function for each time point)
    int nTime = timeVector.size();
    _polyLOESS.resize(nTime);
    for ( int jt=0; jt<nTime; jt++)
    {
//        FUNC_INFO << "jt" << jt << "\n";
        double time0 = timeVector[jt];
        if ( linkWidthToTime ) smoothingWidth = time0 * smoothingScale;    // link smoothing width to time, e.g., smoothingScale=0.5
        int iLowSide=0;  double time = time0;  double width=0.;
        while ( jt - iLowSide >= 0 && width < smoothingWidth )
        {
            time = timeVector[jt-iLowSide];
            width = qAbs(time-time0);
            iLowSide++;
        }
        iLowSide--;
        int iHighSide=0;  time = time0;  width = 0.;
        while ( jt + iHighSide < nTime && width < smoothingWidth )
        {
            time = timeVector[jt+iHighSide];
            width = qAbs(time-time0);
            iHighSide++;
        }
        iHighSide--;
        // Set weights
        dVector weight;  weight.clear();
        int iHalfWidth = qMax(1,qMin(iLowSide,iHighSide));
        dVector localTime, localYData;
        int iLocal=0;
        int iCenter=0;
        if ( !asymmetric ) iLowSide = iHighSide = iHalfWidth;
        for ( int jRel=-iLowSide; jRel<=iHighSide; jRel++)
        {
            int iTime = jt + jRel;
            if ( iTime >= 0 && iTime < nTime )
            {
                time = timeVector[iTime];
                double distance = qAbs(time-time0) / smoothingWidth;
                if ( distance < 0.99 )
                {
                    localTime.append(time);
                    localYData.append(yVector[iTime]);
                    weight.append(qPow(1-qPow(distance,3),3));  // use traditional tri-cube weight function
                    if ( jRel == 0 ) iCenter = iLocal;
                    iLocal++;
                }
//                FUNC_INFO << "jRel, iTime, time, time0" << jRel << iTime << time << time0;
//                FUNC_INFO << "smoothingWidth, distance, weight" << smoothingWidth << distance << weight.last();
            }
        }
//        FUNC_INFO << "localTime[" << jt << "] =" << localTime;
//        FUNC_INFO << "localYData" << localYData;
//        FUNC_INFO << "define poly[" << jt << "] using central point" << iCenter << "and weights" << weight;
        //             define       (nCoeff,          nTime)
        _polyLOESS[jt].define(qMin(weight.size(),3),localTime,iCenter);  // qMax: don't define 3 terms with only 2 points
        _polyLOESS[jt].fitWLS(localYData,false);
//        FUNC_EXIT << "poly[" << jt << "] = " << _polyLOESS[jt].getNumberTimePoints();
        _polyLOESS[jt].setWeights(weight);
    }
}

double LOESS::getFitInterpolation(double xBetween)
{
    double minDistance=1.e10;  int iClosest = -1;
    for ( int jX=0; jX<_polyLOESS.size(); jX++)
    {
        double x = _xPoints[jX];
        double distance = qAbs(x-xBetween);
        if ( distance < minDistance)
        {
            minDistance = distance;
            iClosest = jX;
        }
    }
    FUNC_INFO << "for xBetween =" << xBetween << "closest =" << iClosest << "with x =" << _xPoints[iClosest];
    return _polyLOESS[iClosest].getFitInterpolation(xBetween);
}

dVector LOESS::fit(dVector data)
{
    dVector fit = data;
    for (int jt=0; jt<data.size(); jt++)
    {
        int nLocal  = _polyLOESS[jt].getNumberTimePoints();
        int iHalfWidth = nLocal/2;
        dVector localData;  localData.resize(nLocal);
        int jLocal=0;
        for ( int jRel=-iHalfWidth; jRel<=iHalfWidth; jRel++, jLocal++ )
        {
            int iTime = jt + jRel;
            if ( iTime >= 0 && iTime < data.size() )
                localData[jLocal] = data[iTime];
            else
                localData[jLocal] = 0.;  // should also have weight=0.
        }
        _polyLOESS[jt].fitWLS(localData,true);
        fit[jt] = getFitAtCentralPoint(jt);
    }
    return fit;
}

double LOESS::getFirstDerivative(int iTime)
{ // presume define() and fit() preceeded this
//    FUNC_ENTER << iTime;
    double deriv = 0.;
    if ( iTime > 0 && iTime < _polyLOESS.size()-1 )
        deriv = (getFitAtCentralPoint(iTime+1) - getFitAtCentralPoint(iTime-1)) / (_xPoints[iTime+1] - _xPoints[iTime-1]);
    else if ( iTime == 0 )
        deriv = (getFitAtCentralPoint(1) - getFitAtCentralPoint(0)) / (_xPoints[1] - _xPoints[0]);
    else // iTime == _polyLOESS.size()-1
    {
        int it = _polyLOESS.size()-1;
        deriv= (getFitAtCentralPoint(it) - getFitAtCentralPoint(it-1)) / (_xPoints[it] - _xPoints[it-1]);
    }
    return deriv;
}

double LOESS::getSecondDerivative(int iTime)
{ // presume define() and fit() preceeded this
//    FUNC_ENTER << iTime;
    double deriv = 0.;
    if ( iTime > 0 && iTime < _polyLOESS.size()-1 )
    {
        double x1 = _xPoints[iTime-1];  double y1 = getFitAtCentralPoint(iTime-1);
        double x2 = _xPoints[iTime];    double y2 = getFitAtCentralPoint(iTime);
        double x3 = _xPoints[iTime+1];  double y3 = getFitAtCentralPoint(iTime+1);
        double x21 = x2 - x1;           double x31 = x3 - x1;         double x32 = x3 - x2;
        deriv = 2.*y1/x21/x31 - 2.*y2/x32/x21 + 2.*y3/x32/x31;
    }
    else if ( iTime == 0 )
        deriv = getSecondDerivative(1);
    else
        deriv = getSecondDerivative(_polyLOESS.size()-2);
    return deriv;
}
