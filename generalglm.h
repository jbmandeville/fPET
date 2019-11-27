#ifndef GENERALGLM_H
#define GENERALGLM_H

#include "io.h"

enum fMRIEventShapes // used both in fMRI and timePage classes
{
    Shape_none,
    Shape_Gamma,
    Shape_Sigmoid,
    Shape_Square,
    Shape_RampUp,
    Shape_RampDown,
    Shape_Table,
    Shape_Nuisance,
    Shape_ASL,
    Shape_Sine,
    Shape_Cosine,
    Shape_Baseline0,
    Shape_Baseline1,
    Shape_Baseline2,
    Shape_Baseline3,
    Shape_Baseline4,
    Shape_Baseline5,
    Shape_Interaction,
    Shape_FIR,
    Shape_FR,
    Shape_RSFIR
};
enum PETWeightingModels
{
    Weights_Uniform,
    Weights_11C_Noiseless,
    Weights_11C,
    Weights_noUptake,
    Weights_18F_Noiseless,
    Weights_18F
};
enum PETEventTypes
{
    Type_R1,
    Type_k2,
    Type_k2a,
    Type_dCrdt,
    Type_challenge
};
class GeneralGLM
{
#define MaxConditions 100

private:
    int _nTime=0;
    int _nCoeff=0;
    int _currentCondition=0;
    iVector _nContrasts;     // [nConditions]; = 1 for T test or >= 1 for F test
    QVector<dMatrix> _fCovar;          // [nConditions][nContrasts][nContrasts]
    QVector<dMatrix> _fCovarInv;       // [nConditions][nContrasts][nContrasts]

    bool _initialized=false; // has the GLM been initialized successfully in constructor?
    bool _prepared=false;    // has the pseudo-inverse has been created?
    bool _sigmaREKnown=false;

    bool _basisFunctionsChanged=false;  // use this in a set/get fashion to determine when basis functions have changed
    dMatrix _XTWXm1_cc;      // used for variance for OLS
    int _nIncluded;          // sum of non-zero weights
    dVector _weight_t;
    dVector _fit_t;
    dVector _fitErr_t;
    dVector _beta_c;         // coefficient attached to a basis function
    dVector _sem_c;          // SEM of coefficient
    double _dof=0.;          // DOF for fit
    double _dofEffective=0.; // potentially override _dof by using non-zero _dofEffective (e.g, smoothing to boost dof alla Worsley 2002)
    double _average=0.;      // average across time course

protected:
    dMatrix _X_tc;           // Design matrix
    dMatrix _XTWXm1XTW_ct;   // used for Beta = (Xt*W*X)^-1 * Xt * W for either OLS or WLS (weights=1)
    double _sigma2=0.;       // SOS for fit; for 2nd-order GLM, sigma2 is incorporated into weights (=1/sigma2), so set _sigma2=1;
                             // otherwise, determine _sigma2 post-hoc from data-fit, so that 1/weights = _sigma2 * normalized weights

public:
    // conditions
    QStringList _conditionList; // nConditions = _conditionList.size()
    iVector _conditionShape;    // keep track of what type of condition (R1, k2, k2a, challenge)

    // vector of length _nCoeff (= # events)
    cVector _basisID;           // [_nCoeff]; points to eventID for the respective basis function
    iVector _basisShape;        // [_nCoeff]; R1, k2, k2a, challenge, or shape for fMRI
    QVector<dMatrix> _contrastMatrix;  // [nConditions][_nContrasts][_nCoeff]
    iVector _FIRContrast;
    dVector _RSFIRContrast;
    double _sinusoidPeriod=10.;

    // conditions
    double _fStatistic;      // F statitics for condition
    double _cnr;             // CNR for condition
    double _pRaw;            // raw p-value for condition
    double _effectSize;      // magnitude for condition
    double _effectStDev;       // SEM for condition
    iMatrix _indexCoeffInCondition; // [nConditions][# coeff in condition]; keep track of which coefficients are a condition

    // functions
    void init(int ntime, int nCoeff);  // call this first
    void setWeights(dVector weight_t);
    void setVariances(dVector variance_t);
    void setOLS();
    void replaceBasisFunction( int iBasis, dVector X_t);   // define/replace existing basis function
    void addBasisFunction( dVector X_t);
    void addOrInsertBasisFunction( int iBasis, dVector X_t);
    void fitWLS(dVector data) { fitWLS(data, true);}  // a convenience function for most situations (error not known a-priori)
    void fitWLS(dVector data, bool computeSigma2);
    void calculatePseudoInverse();
    void removeCoefficient(int iCoeff, dVector &data);
    void removeCoefficient(int iCoeff, dVector &data, dVector &fit);
    void orthogonalize(int iCoeff, dVector &data);
    void orthogonalizeExistingBasisFunction(int iCoeff, int iCoeffMod);
    void calculateCovarianceMatrix();
    void calculateFStat (int lCondition);
    void defineConditions(QString conditionString);
    void decodeConditions(bool defineMatrices );
    void evaluateCurrentCondition();

    // setters
    void setNumberContrasts(int iCondition, int nContrasts );
    inline void setPrepared(bool state) {_prepared = state;}
    inline void setBasisFunctionsChanged() {_basisFunctionsChanged=true;}  // allows recoloring of basis functions
    inline void setCurrentCondition(int lCondition) {_currentCondition = lCondition;}
    inline void setDOFEffective(double dof) {_dofEffective=dof;}

    // getters
    QString getConditionString();
    inline QChar getBasisID(int iBasis) {return _basisID[iBasis];}
    inline int getNumberTimePoints()    {return _nTime;}
    inline double getSigma2()           {return _sigma2;}
    inline double getAverage()          {return _average;}
    inline double getBeta(int iCoeff)   {return _beta_c[iCoeff];}
    inline double getSEM(int iCoeff)    {return _sem_c[iCoeff];}
    inline double getVar(int iCoeff)    {return _sem_c[iCoeff] * _sem_c[iCoeff];}
    inline double getCovar(int iC1, int iC2) {return _sigma2 * _XTWXm1_cc[iC1][iC2];}
    inline double getDOF()              {return _dof;}
    inline double getFit(int iTime)     {return _fit_t[iTime];}
    inline double getFitErr(int iTime)  {return _fitErr_t[iTime];}
    inline double getRegressor(int iCoeff, int iTime) {return _beta_c[iCoeff] * _X_tc[iTime][iCoeff];}
    inline bool isReady()               {return _prepared;}
    inline int getNumberCoefficients()  {return _nCoeff;}
    inline int getCurrentCondition()    {return _currentCondition;}
    inline double getEffectSize()       {return _effectSize;}  // call "setCurrentCondition" and "evaluateCurrentCondition" prior to this
    inline double getEffectStDev()      {return _effectStDev;} // ...
    inline double getpRaw()             {return _pRaw;}        // ...
    inline double getCNR()              {return _cnr;}         // ...
    inline double getBasisPoint(int iBasis, int iTime) {return _X_tc[iTime][iBasis];}
    inline int getNumberConditions()                {return _conditionList.size();}
    inline int getNumberContrasts(int iCondition)   {return _nContrasts[iCondition];}
    inline QString getConditionName(int lCondition) {return _conditionList.at(lCondition);}
    inline int getCurrentConditionShape() {return _conditionShape[_currentCondition];}
    inline iVector getFIRContrast() {return _FIRContrast;}
    inline dVector getRSFIRContrast() {return _RSFIRContrast;}
    inline bool isValidID(QChar ID) {return getEventIndex(ID) >= 0;}

    QChar getEventChar(int indexEvent);
    int getEventIndex(QChar eventID);
    int getNumberCoefficientsInCondition(int iCondition);
    int getOnlyCoefficientInCondition(int iCondition);
    int getFirstCoefficientInCondition(int iCondition);
    bool basisFunctionsChanged();
    int getEventCoefficient(QChar eventID);
    int getEventCoefficient(QChar eventID, int iOccurrence);
    iVector getEventCoefficients(QChar eventID);
    double getWeight(int iTime);
    inline dVector getWeights() {return _weight_t;}
};

////////////////////////////////// Polynomial GLM //////////////////////////////////
class PolynomialGLM : public GeneralGLM
{
public:
    void define(int nCoeff, int nTime);         // # coefficiens in polynomial
    void define(int nCoeff, dVector xVector);   // for non-evenly spaced points
private:
};

#endif // GENERALGLM_H
