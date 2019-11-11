#ifndef PETRTM_H
#define PETRTM_H

#include <QTextStream>
#include "io.h"
#include "generalglm.h"

enum RTMModelTypes
{
    RTM_SRTM3,    // 3-parameter SRTM:  R1, k2, k2a (GLM)
    RTM_SRTM2,    // 2-parameter SRTM:      k2, k2a (GLM)
    RTM_SRTM2Fit, // 2-parameter SRTM:      k2, k2a (iterative) using fit value of tau2Ref
    RTM_rFRTM3,   // 3-parameter rFRTM: R1, k2, k2a (iterative) as published (using dCt/dt X E)
    RTM_rFRTM2,   // 2-parameter rFRTM:     k2, k2a (iterative) as published (using dCt/dt X E)
    RTM_rFRTM3New,// 3-parameter rFRTM: R1, k2, k2a (iterative) using fixed k4 to determine BPnd by formulation "k2a"=k2*k4*BPnd
    RTM_rFRTM2New// 2-parameter rFRTM:     k2, k2a (iterative) using fixed k4 to determine BPnd by formulation "k2a"=k2*k4*BPnd
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

class PETRTM : public GeneralGLM    // for multi-threading, use no pointers
{

private:
    bool _visibleInGUI=false;

    QStringList *_warningErrors;
    QStringList *_fatalErrors;

    int _referenceRegionTableColumn=0;  // >0 means take RR from table; otherwise, take RR from overlay

    // overall
    RTMModelTypes _modelRTM=RTM_SRTM3;
    int _PETWeightingModel=Weights_Uniform;
    double _tau4Default=10.;
    double _smoothingScale=0.;
    bool _dCrdtIncluded=false;
    bool _referenceRegionIsDefined=false;
    QString _petConditionsFromTimeModelFile; // a copy to be applied when the time model is read completely (so it won't be trimmed by incomplete info)

    // # of runs
    int _nRuns=0;
    iVector _nTimePerRun;       // [_nRuns]
    // per run
    dMatrix _weights;               // [_nRuns][_nTimePerRun]
    // IDs must be independent across types: e.g., a k2 event in one run cannot be a k2a event in a different run
    cVector _R1EventID;             // [_nRuns]; not required for RTM2
    cVector _k2EventID;             // [_nRuns]; always required
    cVector _k2aEventID;            // [_nRuns]; always required
    cVector _dCrdtEventID;          // [_nRuns]; always required
    dVector _tau4;                  // [_nRuns]
    dVector _tau2RefSRTMFixed;      // [_nRuns]
    dVector _tau2RefFRTMFixed;      // [_nRuns]
    dMatrix _tau2RefSRTMCal;        // [_nRuns][3], with 0,1,2 = tau2 fit using polynomial in BPnd to powers 0,1,2
    double _tau2RefSRTMCalOffset=0.9; // Fraction of SRTM2 value of 1/k2' relative to regularized SRTM3 value
    iVector _R1EventCoefficient;    // [_nRuns]
    iVector _k2EventCoefficient;    // [_nRuns]
    iVector _k2aEventCoefficient;   // [_nRuns]
    iVector _dCrdtEventCoefficient; // [_nRuns]
    sVector _ignoreString;          // [_nRuns]

    // challenges
    cVector _challengeEventID; // [_maxChallenges];
    iVector _challengeShape;   // [_maxChallenges]
    dVector _challengeTau;              // [_maxChallenges]
    dVector _challengeAlpha;            // [_maxChallenges]
    iMatrix _challengeRun;     // [_maxChallenges][_maxStimuli]
    dMatrix _challengeOn;      // [_maxChallenges][_maxStimuli]
    dMatrix _challengeOff;     // [_maxChallenges][_maxStimuli]

    // tissue vectors
    QVector<QStringList> _columnNames; // [_nRuns]
    dMatrix3 _table;                  // [_nRuns][_nTimeInRun][nColumns]; used for creating integrals
    dMatrix _dtBins;                  // [_nRuns][_nTimeInRun]; this should be _table column 0 converted to minutes
    dMatrix _refRegionRaw;            // [_nRuns][_nTimeInRun]; actual ref region data
    dMatrix _refRegion;               // [_nRuns][_nTimeInRun]; value used in analysis (either raw or fit)
    dMatrix _refRegionIntegral;       // [_nRuns][_nTimeInRun]; integral of raw or fit
    dMatrix _refRegionDeriv;          // [_nRuns][_nTimeInRun]
    dMatrix _tissRegionRaw;           // [_nRuns][_nTimeInRun]
    dMatrix _tissRegion;              // [_nRuns][_nTimeInRun]
    dMatrix _tissRegionDeriv;         // [_nRuns][_nTimeInRun]
    dMatrix _frtmConv_dCtdtERaw;      // [_nRuns][_nTimeInRun]; for use with modified basis functions
    dMatrix _frtmConv_dCtdtE;         // [_nRuns][_nTimeInRun];
    dMatrix _frtmConv_CtE;            // [_nRuns][_nTimeInRun];
    dMatrix _frtmConvDeWeightUptake;  // [_nRuns][_nTimeInRun]; for de-weighting the uptake period (1-mag(_frtmConv_dCtdtE)/max)

    // iterative methods: save BPnd immediately after fitting so one can switch between models without affecting BPnd extraction
    dVector _BPndForIterations;       // [_nRuns]; use this for 1) SRTM2Fit (iterative BPnd),
                                      //                        2) 1st pass rFRTMNew (BPnd from SRTM)
                                      //                        3) _fitk4UsingFixedBPnd
    int _nIterations=0;
    bool _fitk4UsingFixedBPnd=false;

    // vector of length nCoeff (= # events)
    iVector _challengeIndex;    // [nCoeff]; index into challenge vectors

    // These matrices keep track of matching k2 & k2a events in order to make BP into a meta-GLM parameter
    iMatrix _iCoeffInCondition;       // [nconditions][nCoeffInCondition]
    iMatrix3 _matchingk2InCondition;  // [nconditions][nCoeffInCondition]
    iMatrix3 _matchingk2aInCondition; // [nconditions][nCoeffInCondition]

    int _challengeForStimAv=-1; // # <0 => ignore; else provides event index
    int _nPreForChallAv=10;  // # time points prior to stimulus onset
    int _nPostForChallAv=10; // # time points following stimulus onset (offset is ignored)
    dVector _xForChallAv;    // [nPre + nPost + 1]
    dVector _yForChallAv;    // [nPre + nPost + 1]
    dVector _ySEMForChallAv; // [nPre + nPost + 1]
    dVector _yFitForChallAv; // [nPre + nPost + 1]

    QVector<QVector<PolynomialGLM>> _quadLOESS;  // [_nRuns][_nTimePerRun]

    void createEventBasisFunction(QChar eventID, dVector &basis);
    static double referenceRegionFitFunction(double time, double *p);
    void smoothRunData(dMatrix &runData);
    void fitLoessCurve(dMatrix &runData);
    double Gauss(double x, double fwhm);
    void updateReferenceRegion();
    void defineLOESSFitting(int iRun);
    void readGLMIgnoreBlock(QTextStream *in_stream, int iRun, QString inputString);
    void integrateByRun(dMatrix &basisFunction );
    void differentiateByRun(dMatrix &basisFunction );
    dVector makeVectorFromROIData(QVector<ROI_data> timeSeriesVector);

    void fitDataByGLM(dMatrix timeSeriesVector);
    void fitDataByGLM(QVector<ROI_data> timeSeriesVector);
    void fitDataByGLM(QVector<ROI_data> timeSeriesVector, dMatrix &yFit);
    void fitDataByGLMIterationForConsistency(QVector<ROI_data> timeSeriesVector, dMatrix &yFit);
    void fitDataByGLMPlusLineScan(QVector<ROI_data> timeSeriesVector, dMatrix &yFit);
    double lineScan1D( int iRun, dVector valueInRun, double &valueIncrement, dVector data );
    double lineScanTau4(int iRun, double &valueIncrement, dVector data );
    double updateLineScanFit(dVector data);
    void fitWLSForIteration(dVector &data);
    void createAllBasisFunctions();
    void createRunBasisFunction(QChar eventID, dVector &basis);
    void createChallengeBasisFunction(int iCoeff, dVector &basis);
    void createChallengeShape(int iRun, int indexChallenge, dVector &shape);


    int readGLMFileOldFormat(int iRun, QString fileName);
    int readGLMFileNewFormat(int iRun, QString fileName);

public:
    QString _refRegionName;
    QString _brainRegionName;
    int _RRLocationInList=-1;   // location of RR in overlay list
    int _brainLocationInList=-1;// location of brain in overlay list
    int _maxChallenges=61;   // fixed size; a-z,A-z,1-9
    int _maxStimuli=20;      // initial allocation, could be resized
    QVector<QString> _frameFiles; // [_nRuns]

    void prepare();
    void calculateFRTMConvolution();
    dVector convolveEquilibration(int iRun, dVector tissue, dVector equilibration);
    void calculateBPndOrK4ForIteration();

    void fitData(dMatrix timeSeriesVector);
    void fitData(dMatrix timeSeriesVector, dMatrix &yFit);
    void fitData(QVector<ROI_data> timeSeriesVector, dMatrix &yFit);
    int readTimeBinsFile(int iRun, QString fileName);
    int readTimeModel(QString fileName);
    int readGLMFile(int iRun, QString fileName);
    QString createConditions();
    void definePETConditions(QString conditionString);
    void updateConditions();
    inline bool isSRTM() {return _modelRTM == RTM_SRTM2  || _modelRTM == RTM_SRTM3  || _modelRTM == RTM_SRTM2Fit;}
    inline bool isFRTM() {return _modelRTM == RTM_rFRTM2     || _modelRTM == RTM_rFRTM3
                ||               _modelRTM == RTM_rFRTM3New || _modelRTM == RTM_rFRTM2New;}
    inline bool isFRTMOld() {return _modelRTM == RTM_rFRTM3 || _modelRTM == RTM_rFRTM2;}
    inline bool isFRTMNew() {return _modelRTM == RTM_rFRTM3New || _modelRTM == RTM_rFRTM2New;}
    inline bool isRTM3() {return _modelRTM == RTM_SRTM3  || _modelRTM == RTM_rFRTM3 || _modelRTM == RTM_rFRTM3New;}
    inline bool isRTM2() {return _modelRTM == RTM_SRTM2  || _modelRTM == RTM_rFRTM2 || _modelRTM == RTM_rFRTM2New || _modelRTM == RTM_SRTM2Fit;}
    inline bool isSRTMReg() {return _modelRTM == RTM_SRTM2Fit;}
    inline bool isFRTMFitk4() {return isFRTMNew() && _fitk4UsingFixedBPnd;}
    void saveTimeModelFiles(QString dirName, QStringList dataFileName);
    void writeGLMFile(int iRun, QString fileName, bool allSameTau4);
    void writeGLMFileOldFormat(int iRun, QString fileName);
    void averageStimuli(dMatrix yData, dMatrix yFit);
    int writeReferenceRegion(int iRun, QString fileName, QString RRName);

    // setters
    inline void setFitk4State(bool state) {_fitk4UsingFixedBPnd = state;}
    inline void setErrorLists(QStringList *warningErrors, QStringList *fatalErrors) {_warningErrors = warningErrors; _fatalErrors = fatalErrors;}
    inline void setEventIndexForAveraging(int indexAv) {_challengeForStimAv = indexAv;}
    inline void setNPreForAveraging(int nPre) {_nPreForChallAv = nPre;}
    inline void setNPostForAveraging(int nPost) {_nPostForChallAv = nPost;}
    inline void setGUIVisibility(bool state) {_visibleInGUI=state;}
    inline void setWeightingModel(int whichWeightingModel) {_PETWeightingModel = whichWeightingModel; setPrepared(false);}
    inline void defineTimeModelFileConditions() {definePETConditions(_petConditionsFromTimeModelFile);}
    void setIgnoredPoints(int iRun, bool resetWeights, QString ignoreString);
    void setTimePointsInRun(int iRun, int nTime);
    void setNumberRuns(int nFiles);
    void setChallengeShape(int indexChallenge, int iShape);
    void setChallengeRun(int indexChallenge, int indexStimulus, int iRun);
    void setChallengeOnset(int indexChallenge, int indexStimulus, double time);
    void setChallengeOffset(int indexChallenge, int indexStimulus, double time);
    void setChallengeTau(int indexChallenge, double tau);
    void setChallengeAlpha(int indexChallenge, double alpha);
    void setTimeBins(int iRun, dVector timeBins);
    void setReferenceRegion(dMatrix referenceRegionRaw);
    void setReferenceRegion(dMatrix timeBins, dMatrix referenceRegionRaw);
    void setReferenceRegionFromTableColumn(int iColumn);
    void setTissueVector(dMatrix tissueRegion);
    void setRTMModelType(RTMModelTypes model);
    void setRTMModelType(QString model);
    void setWeightsInRun(int iRun);
    void setSmoothingScale(double smoothingScale);
    inline void setTau4(int iRun, double tau4) {_tau4[iRun] = tau4; setPrepared(false);}
    inline void setTau2RefSRTM(int iRun, double tau2Ref) {_tau2RefSRTMFixed[iRun] = tau2Ref; setPrepared(false);}
    inline void setTau2RefFRTM(int iRun, double tau2Ref) {_tau2RefFRTMFixed[iRun] = tau2Ref; setPrepared(false);}
    inline void setTau2RefSRTMCal(int iRun, double poly0, double poly1, double poly2)
    {_tau2RefSRTMCal[iRun][0] = poly0; _tau2RefSRTMCal[iRun][1] = poly1; _tau2RefSRTMCal[iRun][2] = poly2; setPrepared(false);}
    inline void setTau2RefSRTMCalOffset(double value) {_tau2RefSRTMCalOffset=value; setPrepared(false);}
    void setTau2Ref(int iRun, double tau2Ref);

    void setR1EventID(int iRun, QChar ID);
    void setdCrdtEventID(int iRun, QChar ID);
    inline void setk2EventID(int iRun, QChar ID) {_k2EventID[iRun] = ID;   setPrepared(false); setBasisFunctionsChanged();}
    inline void setk2aEventID(int iRun, QChar ID) {_k2aEventID[iRun] = ID; setPrepared(false); setBasisFunctionsChanged();}
    inline void setInclusionOfdCrdt(bool state) {_dCrdtIncluded = state;   setPrepared(false); setBasisFunctionsChanged();}
    inline void setReferenceRegionName(QString name) {_refRegionName = name;}
    inline void setBrainRegionName(QString name) {_brainRegionName = name;}

    // getters
    inline bool getFitk4State() {return _fitk4UsingFixedBPnd;}
    inline int getNumberColumns(int iRun) {if (_columnNames.size() > 0) return _columnNames[iRun].count(); else return 0;}
    inline QString getColumnName(int iRun, int iColumn) {return _columnNames[iRun].at(iColumn);}
    inline bool useReferenceRegionFromOverlay() {return _referenceRegionTableColumn <= 0;}
    inline bool useReferenceRegionFromTable()   {return _referenceRegionTableColumn >  0;}
    inline int getReferenceRegionColumn() {return _referenceRegionTableColumn;}
    inline int getChallengeIndexForAveraging() {return _challengeForStimAv;}
    inline int getNPreForAveraging() {return _nPreForChallAv;}
    inline int getNPostForAveraging() {return _nPostForChallAv;}
    void getStimAveragingVectors(dVector &xForChallAv, dVector &yForChallAv, dVector &ySEMForChallAv, dVector &yFitForChallAv);
    inline bool isVisible() {return _visibleInGUI;}
    double getTimeInRun(int iRun, int iTime);
    double getTimeInRun(int iRun, double rTime);
    double interpolateVector(dVector inputVector, double rBin);
    inline int getNumberRuns() {return _nRuns;}
    inline int getNumberEvents() {return _basisID.size();}
    inline int getNumberTimePointsInRun(int iRun) {return _nTimePerRun[iRun];}
    inline int getChallengeRun(int indexChallenge, int indexStimulus) {return _challengeRun[indexChallenge][indexStimulus];}
    inline double getChallengeOnset(int indexChallenge, int indexStimulus) {return _challengeOn[indexChallenge][indexStimulus];}
    inline double getChallengeOffset(int indexChallenge, int indexStimulus) {return _challengeOff[indexChallenge][indexStimulus];}
    inline int getChallengeShape(int indexChallenge) {return _challengeShape[indexChallenge];}
    inline double getChallengeTau(int indexChallenge) {return _challengeTau[indexChallenge];}
    inline double getChallengeAlpha(int indexChallenge) {return _challengeAlpha[indexChallenge];}
    inline double getSmoothingScale() {return _smoothingScale;}
    inline double getTau4(int iRun) {return _tau4[iRun];}
    inline int getNumberIterations() {return _nIterations;}
    inline bool getInclusionOfdCrdt() {return _dCrdtIncluded;}
    inline QString getTimeFrameFileName(int iRun) {return _frameFiles[iRun];}
    inline double getTau2RefSRTMInRun(int iRun) {return _tau2RefSRTMFixed[iRun];}
    inline double getTau2RefFRTMInRun(int iRun) {return _tau2RefFRTMFixed[iRun];}
    inline double getTau2RefSRTMCalInRun(int iRun, double BPnd)  // return a quadratic polynomial in BPnd
    {return _tau2RefSRTMCal[iRun][0] + _tau2RefSRTMCal[iRun][1] * BPnd + _tau2RefSRTMCal[iRun][2] * BPnd * BPnd;}
    inline double getTau2RefSRTMCalOffset() {return _tau2RefSRTMCalOffset;}
    inline bool isReferenceRegionDefined() {return _referenceRegionIsDefined;}

    bool getFrameStatus();
    inline QChar getR1EventID(int iRun) {return _R1EventID[iRun];}
    inline QChar getk2EventID(int iRun) {return _k2EventID[iRun];}
    inline QChar getk2aEventID(int iRun) {return _k2aEventID[iRun];}
    inline QChar getdCrdtEventID(int iRun) {return _dCrdtEventID[iRun];}
    int getNumberChallenges();
    int getNumberChallengesInRun(int iRun);
    int getFirstGoodStimulus(QChar ID);
    int getFirstGoodStimulus(int indexChallenge);
    int getFirstGoodStimulusInRun(int indexChallenge, int iRun);
    int getFirstBadStimulus(int indexChallenge);
    QChar getFirstGoodChallenge();
    int getFirstGoodChallengeIndex();
    int getFirstGoodChallengeIndexInRun(int iRun);
    int getBasisType(int iBasis) {return _basisShape[iBasis];}
    QString getIgnoredString(int iRun);
    int getTotalTimeIndex(int iRun, int iTimeInRun);

    int countEvents();
    bool isGoodChallenge(QChar ID);
    bool isGoodChallenge(int indexChallenge);
    bool isGoodChallengeInRun( QChar ID, int iRun );
    bool isGoodChallengeInRun( int indexChallenge, int iRun );
    bool isGoodStimulus(int indexChallenge, int indexStimulus);
    bool isGoodStimulusInRun( int indexChallenge, int indexStimulus, int iRun);
    int getChallengeInfoRequired(int iShape);
    inline int getPETWeightingModel() {return _PETWeightingModel;}
    dVector getEquilibrationVector(int iFile);
    dVector getBPndVector(int iFile);
    dPoint2D getBPndVersusTime(int iFile, int iTime);
    dPoint2D getBPndInCurrentCondition();
    double getBP0InRun(int iRun);
    double getTau4InRun(int iRun);
    dPoint2D averageParameter(iVector iCoeffVector);
    dPoint2D getk2RefInCurrentCondition();
    dPoint2D getTau2RefInCurrentCondition();
    dPoint2D getTau2InCurrentCondition();
    double getTau2RefInRun(int iRun);
    dPoint2D getk2RefInRun(int iRun);
    dPoint2D getR1InRun(int iRun);
    dPoint2D getk2InRun(int iRun);
    dPoint2D getk2aInRun(int iRun);
    dPoint2D getdk2aInRun(QChar challengeID);
    dPoint2D getReferenceRegionTimesR1(int iRun, int iTime);
    dPoint2D getReferenceRegion(bool useFit, int iRun, int iTime);
    dPoint2D getTissueRegion(bool useFit, int iRun, int iTime);
    QString getCurrentConditionTypeString();
    inline bool currentconditionIsBPType() {return (getCurrentConditionShape() == Type_k2a)
                ||                                 (getCurrentConditionShape() == Type_challenge);}
    inline bool currentconditionIsk2Type() {return getCurrentConditionShape() == Type_k2;}
    inline bool currentconditionIsR1Type() {return getCurrentConditionShape() == Type_R1;}
    inline bool currentconditionIsdCrdtType() {return getCurrentConditionShape() == Type_dCrdt;}
    dPoint2D getValueAndErrorForCurrentCondition();
    inline double getFRTMConvolution(int iFile, int iTime) {if ( isFRTMNew() ) return _frtmConv_CtE[iFile][iTime]; else return _frtmConv_dCtdtERaw[iFile][iTime];}
    inline double getFRTMConvolutionFit(int iFile, int iTime) {return _frtmConv_dCtdtE[iFile][iTime];}
    inline double getCtMinusCr(int iFile, int iTime) {return _tissRegion[iFile][iTime] - _refRegion[iFile][iTime];}

    bool isValidID(int iRun, int iType, QChar eventID);
};

#endif // PETRTM_H
