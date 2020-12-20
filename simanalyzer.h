#ifndef SIMANALYZER_H
#define SIMANALYZER_H

#include <QDebug>

#include "io.h"
#include "generalglm.h"

enum PETChallengeShapes // used both in PETRTM and timePage classes
{
    Challenge_none,
    Challenge_Constant,
    Challenge_Sigmoid,
    Challenge_Gamma,
    Challenge_Square,
    Challenge_RampUp,
    Challenge_RampDown,
    Challenge_Table
};

class simAnalyzer
{
private:
    // Input quantities ///////////////////////
    // Setup bins
    int _nBins = 90;
    int _nCoeff = 3;   // frtm3: BPnd, R1, k2Ref
    iVector _durationBinSec;     // [_nBins]; seconds
    iVector _numberSamplesPerBin;  // [_nBins]
    dVector _dtFine;          // [sum of _numberSamplesPerBin = # of fine bins]
    dVector _timeFine;        // [sum of _numberSamplesPerBin = # of fine bins]
    dVector _dtCoarse;        // [_nBins]
    dVector _timeCoarse;      // [_nBins]

    // Reference region
    double _k2Ref=1./2.5;
    double _k2RefFixed=1./2.5;
    // Target region: state
    double _tau4 = 0.; // SRTM
    double _R1=1.;
    double _BP0=4.;     // BP = k3/koff = kon*Bavail0/koff
    double _plasmaMag=0.;
    bool _includePlasma=false;
    // errors
    double _BP0Err;
    double _R1Err;
    double _k2RefErr;
    double _dBPndErr;

    // Challenge
    int _challengeShape;
    double _challengeTau;
    dVector _challengeOn;      // [_maxStimuli]
    dVector _challengeOff;     // [_maxStimuli]
    dVector _challengeCurve;   // [_dtFine.size()]
    dVector _challengeCurveBinned;
    double _uptakeTime=10.;
    double _dBPndChallenge=0.;  // deltaBPnd in absolute units
    double _dBPndUptake=0.;     // deltaBPnd during uptake period

    bool _RTM3=true;
    bool _fitk4 = false;
    dVector _parValues     = {_BP0,_R1,_k2Ref,_tau4,_dBPndChallenge, _plasmaMag, _dBPndUptake};
    dVector _parIncrements = {_BP0,_R1,_k2Ref,_tau4,_dBPndChallenge, _plasmaMag, _dBPndUptake};
    bVector _parAdjust     = {true, true, _RTM3, _fitk4, false, false};
    bool _optHigh = false;

    // concentrations: fine scale
    dVector _Cr;        // fine samples Cr (from fit)
    dVector _CrDot;     // derivative of _Cr
    dVector _Ct;        // computed value of Ct (fine scale)
    // concentrations: binned/framed
    dVector _CrBinned;  // the fit to Cr prior to interpolation
    dVector _CtBinned;  // simulator value of Ct
    dVector _CtData;    // data value of Ct, passed for comparison/cost function
    dVector _weights;

    // derived quantities //////////////////////////
    double _k2;    // k2 = R1 * k2_ref
    double _k3;    // k3 = BP * koff

    double _sigma2;

    void initializeBins();
    void updateFineSamples();
    void updateChallengeShape();
    dVector downSample(dVector original);
    void fitTACByGridSearch(int level);
    void fitTACByLineScan(int level); // toleranceCost is a fraction (e.g., 0.01 is 1%)
    void lineScan1D( int iPar, double &costRelative, double &incrementOpt );
    void calculateErrorsByBootstrapping(int nBootStrap);

    void generateTargetTACFRTM();
    void generateTargetTACSRTM();

public:
    simAnalyzer();
    void fit(double BP0, double R1, double k2Ref, double dBPnd, double tau4, bool fitk4, bool RTM3, bool challenge, int nBootStrap);
    double runAndCalculateCost();

    // setters
    inline void setFitk4(bool state)          {_fitk4 = state;}
    inline void setRTM3(bool RTM3)            {_RTM3 = RTM3;}
    inline void setNumberBins(int nBins)      {_nBins = nBins; initializeBins();}
    inline void setk2RefRTM3(double value)    {_k2Ref = value; FUNC_INFO << _k2Ref; }
    inline void setTau2RefRTM3(double value)  {_k2Ref = 1./value; FUNC_INFO << value;}
    inline void setk2RefRTM2(double value)    {_k2Ref = _k2RefFixed = value; FUNC_INFO << _k2RefFixed; }
    inline void setTau2RefRTM2(double value)  {_k2Ref = _k2RefFixed = 1./value; FUNC_INFO << value;}
    inline void setTau4(double value)         {FUNC_ENTER << value; _tau4=value;}
    inline void setR1(double value)           {_R1 = value;}
    inline void setBP0(double value)          {_BP0 = value;}
    inline void setk2(double value)           {_k2 = value;}
    inline void setk3(double value)           {_k3 = value;}
    inline void setChallengeShape(int iShape)   {_challengeShape = iShape;}
    inline void setChallengeTau(double tau)     {_challengeTau = tau;}
    inline void setChallengeAlpha(double alpha) {_challengeTau = alpha;}
    inline void setIncludePlasma(bool state)    {_includePlasma = state;}
    void defineChallenge(int shape, double tau, dVector onTime, dVector offTime);

    inline void setUptakeTime(double value)   {_uptakeTime     = value;}
    inline void setSamplesPerBin(int lBin, int nSamples)  {_numberSamplesPerBin[lBin] = nSamples; updateFineSamples();}
    inline void setDurationBin(int lBin, int duration) {_durationBinSec[lBin] = duration;         updateFineSamples();}
    inline void setDurationBins(iVector durationVector)   {_durationBinSec = durationVector;      updateFineSamples();}
    inline void setCtData(dVector data) {_CtData = data;}
    inline void setWeights(dVector weights) {FUNC_ENTER << weights; _weights = weights;}

    void setCrFit(dVector CrFitCoarse); // this sets (_CrBinned, _Cr, _CrDot)

    // getters
    inline int getNumberTimeBinsFine()   {return _dtFine.size();}
    inline int getNumberTimeBinsCoarse() {return _timeCoarse.size();}
    inline int getNumberBins()      {return _nBins;}
    inline double getk2Ref()        {if ( _RTM3 ) return _k2Ref; else return _k2RefFixed;}
    inline double getk2RefRTM3()    {FUNC_INFO << _k2Ref;    return _k2Ref;}
    inline double getTau2RefRTM3()  {FUNC_INFO << 1./_k2Ref; return 1./_k2Ref;}
    inline double getk2RefRTM2()    {FUNC_INFO << _k2RefFixed;    return _k2RefFixed;}
    inline double getTau2RefRTM2()  {FUNC_INFO << 1./_k2RefFixed; return 1./_k2RefFixed;}
    inline double getk4()           {return 1./_tau4;}
    inline double getTau4()         {return _tau4;}
    inline double getR1()           {return _R1;}
    inline double getBP0()          {return _BP0;}
    inline double getBPndChallenge(){return _BP0 - _dBPndChallenge;}
    inline double getBPndUptake()   {return _BP0 - _dBPndUptake;}
    inline double getPlasmaMag()    {return _plasmaMag;}
    inline double getk2()           {return _k2;}
    inline double getk3()           {return _k3;}
    inline double getk2a()          {return _k2/(1.+_BP0);}
    inline double getk2k3()         {return _k2 * _k3;}  // or k2 * k4 * BPnd
    inline int getSamplesPerBin(int lBin) {return _numberSamplesPerBin[lBin];}
    inline double getDurationPerBin(int lBin) {return static_cast<double>(_durationBinSec[lBin])*60.;}
    inline double getDurationPerBinSec(int lBin) {return _durationBinSec[lBin];}
    inline iVector getTimeBinVectorSec()    {return _durationBinSec;}
    inline double getdtFine(int iTime)      {return _dtFine[iTime];}
    inline double getTimeFine(int iTime)    {return _timeFine[iTime];}
    inline double getTimeCoarse(int iTime)  {return _timeCoarse[iTime];}
    inline double getdtCoarse(int iTime)    {return _dtCoarse[iTime];}
    inline double getBPndCoarse(int iTime)  {return _BP0 - _challengeCurveBinned[iTime] * _dBPndChallenge;}
    inline double getChallengeCoarse(int iTime) {return _challengeCurveBinned[iTime] * _dBPndChallenge;}
    inline double getSigma2() {return _sigma2;}
    inline double getAIC()    {return 2*_nCoeff + _CrBinned.size() * qLn(_sigma2);}
    inline bool getFitk4State()  {return _fitk4;}
    double getDurationScan();
    double getdk2a();
    double getdk2k3();

    inline double getUptakeTime()   {return _uptakeTime;}
    inline double getChallengeMag() {return _dBPndChallenge;}
    inline double getUptakeMag()    {return _dBPndUptake;}
    inline double getBP0Err()       {return _BP0Err;}
    inline double getChallengeErr() {return _dBPndErr;}
    inline double getk2RefErr()     {return _k2RefErr;}
    inline double getTau2RefErr()   {return _k2RefErr/SQR(_k2Ref);}     // Tau=1/k ; dTau=dk/Tau^2
    // concentations
    inline double getCr(int iTime)       {return _Cr[iTime];}
    inline double getCt(int iTime)       {return _Ct[iTime];}
    inline double getCrDot(int iTime)    {return _CrDot[iTime];}
    inline double getCrCoarse(int iTime) {return _CrBinned[iTime];}
    inline double getCtCoarse(int iTime) {return _CtBinned[iTime];}
    inline dVector getTimeCourse()       {return _timeCoarse;}
    inline dVector getCtFit()            {return _CtBinned;}

};

#endif // SIMANALYZER_H
