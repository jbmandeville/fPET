#ifndef SIMANALYZER_H
#define SIMANALYZER_H

#include <QDebug>

#include "io.h"
#include "generalglm.h"

class simAnalyzer
{
private:
    // Input quantities ///////////////////////
    // Setup bins
    int _nBins = 90;
    iVector _durationBinSec;     // [_nBins]; seconds
    iVector _numberSamplesPerBin;  // [_nBins]
    dVector _dtFine;          // [sum of _numberSamplesPerBin = # of fine bins]
    dVector _timeFine;        // [sum of _numberSamplesPerBin = # of fine bins]
    dVector _dtCoarse;        // [_nBins]
    dVector _timeCoarse;      // [_nBins]

    // Reference region
    double _k2Ref=1./2.5;
    // Target region: state
    double _tau4 = 0.; // SRTM
    double _R1=1.;
    double _BP0=4.;     // BP = k3/koff = kon*Bavail0/koff
    // Challenge
    double _challengeTime=40.;
    double _deltaBPPercent=0.;  // deltaBPnd in %
    dVector _parValuesFRTM3     = {_BP0,_R1,_k2Ref};
    dVector _parIncrementsFRTM3 = {_BP0,_R1,_k2Ref};
    bVector _parAdjust          = {true, true, true};
    bool _optHigh = false;

    // concentrations: fine scale
    dVector _Cr;        // fine samples Cr (from fit)
    dVector _CrDot;     // derivative of _Cr
    dVector _Ct;        // computed value of Ct (fine scale)
    // concentrations: binned/framed
    dVector _CrBinned;  // the fit to Cr prior to interpolation
    dVector _CtBinned;  // simulator value of Ct
    dVector _CtData;    // data value of Ct, passed for comparison/cost function

    // derived quantities //////////////////////////
    double _k2;    // k2 = R1 * k2_ref
    double _k3;    // k3 = BP * koff

    double _sigma2;

    void initializeBins();
    void updateFineSamples();
    dVector downSample(dVector original);
    void fitTACByGridSearch(int level);
    void fitTACByLineScan(int level); // toleranceCost is a fraction (e.g., 0.01 is 1%)
    void lineScan1D( int iPar, double &costRelative, double &incrementOpt );

    void generateTargetTACFRTM();
    void generateTargetTACSRTM();

public:
    simAnalyzer();
    void fit(double BP0, double R1, double k2Ref, bool fitAll);
    double runAndCalculateCost();

    // setters
    inline void setNumberBins(int nBins)      {_nBins = nBins; initializeBins();}
    inline void setk2Ref(double value)        {_k2Ref = value;}
    inline void setTau2Ref(double value)      {_k2Ref = 1./value;}
    inline void setTau4(double value)         {FUNC_ENTER << value; _tau4=value;}
    inline void setR1(double value)           {_R1 = value;}
    inline void setBP0(double value)          {_BP0 = value;}
    inline void setk2(double value)           {_k2 = value;}
    inline void setk3(double value)           {_k3 = value;}
    inline void setChallengeTime(double value){_challengeTime  = value;}
    inline void setChallengeMag(double value) {_deltaBPPercent = value;}
    inline void setSamplesPerBin(int lBin, int nSamples)  {_numberSamplesPerBin[lBin] = nSamples; updateFineSamples();}
    inline void setDurationBin(int lBin, int duration) {_durationBinSec[lBin] = duration;         updateFineSamples();}
    inline void setDurationBins(iVector durationVector)   {_durationBinSec = durationVector;      updateFineSamples();}
    inline void setCtData(dVector data) {_CtData = data;}

    void setCrFit(dVector CrFitCoarse); // this sets (_CrBinned, _Cr, _CrDot)

    // getters
    inline int getNumberTimeBinsFine()   {return _dtFine.size();}
    inline int getNumberTimeBinsCoarse() {return _timeCoarse.size();}
    inline int getNumberBins()      {return _nBins;}
    inline double getk2Ref()        {return _k2Ref;}
    inline double getTau2Ref()      {return 1./_k2Ref;}
    inline double getk4()           {return 1./_tau4;}
    inline double getTau4()         {return _tau4;}
    inline double getR1()           {return _R1;}
    inline double getBP0()          {return _BP0;}
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
    inline double getSigma2() {return _sigma2;}
    double getDurationScan();
    double getdk2a();
    double getdk2k3();

    inline double getChallengeTime(){return _challengeTime;}
    inline double getChallengeMag() {return _deltaBPPercent;}
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
