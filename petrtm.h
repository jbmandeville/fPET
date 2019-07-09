#ifndef PETRTM_H
#define PETRTM_H

#include <QTextStream>
#include "io.h"
#include "generalglm.h"

enum RTMModelTypes
{
    RTM_SRTM3,    // 3-parameter SRTM:  R1, k2, k2a (GLM)
    RTM_SRTM2,    // 2-parameter SRTM:      k2, k2a (GLM)
    RTM_rFRTM3,   // 3-parameter rFRTM: R1, k2, k2a (iterative)
    RTM_rFRTM2,   // 2-parameter rFRTM:     k2, k2a (iterative)
    RTM_SRTM2_R1, // 2-parameter SRTM:  R1, k2      (GLM)
    RTM_SRTM2x2,  // 2-parameter iterative SRTM that alternates between RTM_SRTM2 and RTM_SRTM2_R1
    RTM_rFRTM2_R1 // 2-parameter rFRTM: R1, k2      (GLM)
};

class PETRTM : public GeneralGLM    // for multi-threading, use no pointers
{

private:
    bool _visibleInGUI=false;
    bool _oneIterationOnly = false;

    QStringList *_warningErrors;
    QStringList *_fatalErrors;

    int _referenceRegionTableColumn=0;  // >0 means take RR from table; otherwise, take RR from overlay

    // overall
    int _modelRTM=RTM_SRTM3;
    bool _calibratingTau4=false;
    int _PETWeightingModel=Weights_Uniform;
    double _tau4Default=10.;
    dVector _tau4;
    int _nIterations=0;
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
    dVector _tau2RefSRTM;           // [_nRuns]
    dVector _tau2RefFRTM;           // [_nRuns];
    dMatrix _BPnd;                  // [_nRuns][_nTimePerRun]
    iVector _R1EventCoefficient;    // [_nRuns]
    iVector _k2EventCoefficient;    // [_nRuns]
    iVector _k2aEventCoefficient;   // [_nRuns]
    iVector _dCrdtEventCoefficient; // [_nRuns]
    sVector _ignoreString;          // [_nRuns]

    // challenges
    cVector _challengeEventID; // [_maxChallenges];
    iVector _challengeShape;   // [_maxChallenges]
    dVector _tau;              // [_maxChallenges]
    dVector _alpha;            // [_maxChallenges]
    iMatrix _challengeRun;     // [_maxChallenges][_maxStimuli]
    dMatrix _challengeOn;      // [_maxChallenges][_maxStimuli]
    dMatrix _challengeOff;     // [_maxChallenges][_maxStimuli]

    // tissue vectors
    QVector<QStringList> _columnNames; // [_nRuns]
    dMatrix3 _table;                  // [_nRuns][_nTimeInRun][nColumns]; used for creating integrals
    dMatrix _refRegionRaw;            // [_nRuns][_nTimeInRun]; actual ref region data
    dMatrix _refRegion;               // [_nRuns][_nTimeInRun]; value used in analysis (either raw or fit)
    dMatrix _refRegionIntegral;       // [_nRuns][_nTimeInRun]; integral of raw or fit
    dMatrix _refRegionDeriv;          // [_nRuns][_nTimeInRun]
    dMatrix _tissRegionRaw;           // [_nRuns][_nTimeInRun]
    dMatrix _tissRegion;              // [_nRuns][_nTimeInRun]
    dMatrix _tissRegionDeriv;         // [_nRuns][_nTimeInRun]
    dMatrix _frtmConvolution;         // [_nRuns][_nTimeInRun]; potentially a fitted version of raw
    dMatrix _frtmConvolutionRaw;      // [_nRuns][_nTimeInRun]; for use with modified basis functions
    dMatrix _frtmConvolutionRelative; // [_nRuns][_nTimeInRun]; rescaled to max 1 for weighting

    // vector of length _nCoeff (= # events)
    iVector _challengeIndex;    // [_nCoeff]; index into challenge vectors

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

    void createAllBasisFunctions();
    void createRunBasisFunction(bool newTAC, QChar eventID, dVector &basis);
    void createChallengeBasisFunction(int iCoeff, dVector &basis);
    void createChallengeShape(int iRun, int indexChallenge, dVector &shape);

    void recreateAllFRTMBasisFunctions();

public:
    QString _refRegionName;
    QString _brainRegionName;
    int _RRLocationInList=-1;   // location of RR in overlay list
    int _brainLocationInList=-1;// location of brain in overlay list
    int _maxChallenges=61;   // fixed size; a-z,A-z,1-9
    int _maxStimuli=20;      // initial allocation, could be resized
    QVector<QString> _frameFiles; // [_nRuns]

    void prepare();
    void resetAndCalculateFRTMConvolution(bool calculateConvolution);

    void fitData(dMatrix timeSeriesVector, dMatrix &yFit);
    void fitData(QVector<ROI_data> timeSeriesVector, dMatrix &yFit);
    void fitDataByGLM(QVector<ROI_data> timeSeriesVector, dMatrix &yFit);
    void fitDataByFRTMBasisFunctions(QVector<ROI_data> timeSeriesVector, dMatrix &yFit);
    void fitDataByIterativeAlternating2ParametersModels(QVector<ROI_data> timeSeriesVector, dMatrix &yFit);
    void switchToAlternateModel_rFRTM2();
    int readTimeBinsFile(int iRun, QString fileName);
    QString createConditions();
    void definePETConditions(QString conditionString);
    void updateConditions();
    inline bool isSRTM() {return _modelRTM == RTM_SRTM2  || _modelRTM == RTM_SRTM3  || _modelRTM == RTM_SRTM2x2 || _modelRTM == RTM_SRTM2_R1;}
    inline bool isFRTM() {return _modelRTM == RTM_rFRTM2 || _modelRTM == RTM_rFRTM3 || _modelRTM == RTM_rFRTM2_R1;}
    inline bool isRTM2_k2a() {return _modelRTM == RTM_SRTM2     || _modelRTM == RTM_rFRTM2;}
    inline bool isRTM2_R1()  {return _modelRTM == RTM_SRTM2_R1  || _modelRTM == RTM_rFRTM2_R1;}
    inline bool isRTM2()     {return isRTM2_k2a() || isRTM2_R1();}
    void saveTimeModelFiles(QString dirName, QStringList dataFileName);
    void writeGLMFile(int iRun, QString fileName, bool allSameTau4);
    void writeGLMFileOldFormat(int iRun, QString fileName);
    void averageStimuli(dMatrix yData, dMatrix yFit);
    int writeReferenceRegion(int iRun, QString fileName, QString RRName);
    void calculateTau2PrimeForSRTM2ForBiasFreeBPND(int iRun);

    // setters
    inline void setErrorLists(QStringList *warningErrors, QStringList *fatalErrors) {_warningErrors = warningErrors; _fatalErrors = fatalErrors;}
    inline void setEventIndexForAveraging(int indexAv) {_challengeForStimAv = indexAv;}
    inline void setNPreForAveraging(int nPre) {_nPreForChallAv = nPre;}
    inline void setNPostForAveraging(int nPost) {_nPostForChallAv = nPost;}
    inline void setOneIteration(bool state) {_oneIterationOnly=state;}
    inline void setGUIVisibility(bool state) {_visibleInGUI=state;}
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
    void setReferenceRegion(dMatrix timeBins, dMatrix referenceRegionRaw);
    void setReferenceRegion(dMatrix referenceRegionRaw);
    void setReferenceRegionFromTableColumn(int iColumn);
    void setTissueVector(bool newTAC, dMatrix tissueRegion);
    void setRTMModelType(QString model);
    void setRTMModelType(RTMModelTypes model);
    void setWeightingModel(int whichWeightingModel);
    void setWeightsInRun(int iRun);
    void setSmoothingScale(double smoothingScale);
    inline void setTau4(int iRun, double tau4) {_tau4[iRun] = tau4; setPrepared(false);}
    inline void setBPnd(int iRun, int iTime, double BPnd) {_BPnd[iRun][iTime] = BPnd; setPrepared(false);}
    inline void setTau2RefSRTM(int iRun, double tau2Ref) {_tau2RefSRTM[iRun] = tau2Ref; setPrepared(false);}
    inline void setTau2RefFRTM(int iRun, double tau2Ref) {_tau2RefFRTM[iRun] = tau2Ref; setPrepared(false);}
    void setTau2Ref(int iRun, double tau2Ref);

    void setR1EventID(int iRun, QChar ID);
    void setdCrdtEventID(int iRun, QChar ID);
    inline void setk2EventID(int iRun, QChar ID) {_k2EventID[iRun] = ID;   setPrepared(false); setBasisFunctionsChanged();}
    inline void setk2aEventID(int iRun, QChar ID) {_k2aEventID[iRun] = ID; setPrepared(false); setBasisFunctionsChanged();}
    inline void setInclusionOfdCrdt(bool state) {_dCrdtIncluded = state;   setPrepared(false); setBasisFunctionsChanged();}
    inline void setTau4Calibration(bool state) {_calibratingTau4 = state;}
    inline void setReferenceRegionName(QString name) {_refRegionName = name;}
    inline void setBrainRegionName(QString name) {_brainRegionName = name;}

    // getters
    inline int getNumberColumns(int iRun) {if (_columnNames.size() > 0) return _columnNames[iRun].count(); else return 0;}
    inline QString getColumnName(int iRun, int iColumn) {return _columnNames[iRun].at(iColumn);}
    inline bool useReferenceRegionFromOverlay() {return _referenceRegionTableColumn <= 0;}
    inline bool useReferenceRegionFromTable()   {return _referenceRegionTableColumn >  0;}
    inline int getReferenceRegionColumn() {return _referenceRegionTableColumn;}
    inline int getChallengeIndexForAveraging() {return _challengeForStimAv;}
    inline int getNPreForAveraging() {return _nPreForChallAv;}
    inline int getNPostForAveraging() {return _nPostForChallAv;}
    void getStimAveragingVectors(dVector &xForChallAv, dVector &yForChallAv, dVector &ySEMForChallAv, dVector &yFitForChallAv);
    inline bool isOneIteration() {return _oneIterationOnly;}
    inline bool isVisible() {return _visibleInGUI;}
    double getTimeInRun(int iRun, int iTime);
    inline int getNumberRuns() {return _nRuns;}
    inline int getNumberEvents() {return _basisID.size();}
    inline int getNumberTimePointsInRun(int iRun) {return _nTimePerRun[iRun];}
    inline int getChallengeRun(int indexChallenge, int indexStimulus) {return _challengeRun[indexChallenge][indexStimulus];}
    inline double getChallengeOnset(int indexChallenge, int indexStimulus) {return _challengeOn[indexChallenge][indexStimulus];}
    inline double getChallengeOffset(int indexChallenge, int indexStimulus) {return _challengeOff[indexChallenge][indexStimulus];}
    inline int getChallengeShape(int indexChallenge) {return _challengeShape[indexChallenge];}
    inline double getChallengeTau(int indexChallenge) {return _tau[indexChallenge];}
    inline double getChallengeAlpha(int indexChallenge) {return _alpha[indexChallenge];}
    inline double getSmoothingScale() {return _smoothingScale;}
    inline double getTau4(int iRun) {return _tau4[iRun];}
    inline int getNumberFRTMIterations() {return _nIterations;}
    inline bool getInclusionOfdCrdt() {return _dCrdtIncluded;}
    inline QString getTimeFrameFileName(int iRun) {return _frameFiles[iRun];}
    inline double getTau2RefSRTMInRun(int iRun) {return _tau2RefSRTM[iRun];}
    inline double getTau2RefFRTMInRun(int iRun) {return _tau2RefSRTM[iRun];}
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
    dPoint2D averageParameter(iVector iCoeffVector);
    dPoint2D getk2RefInCurrentCondition();
    dPoint2D getTau2RefInCurrentCondition();
    dPoint2D getTau2InCurrentCondition();
    double getTau2RefInRun(int iRun);
    dPoint2D getk2RefInRun(int iRun);
    dPoint2D getR1InRun(int iRun);
    dPoint2D getk2InRun(int iRun);
    dPoint2D getk2aInRun(int iRun);
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
    inline double getFRTMConvolution(int iFile, int iTime) {return _frtmConvolutionRaw[iFile][iTime];}
    inline double getFRTMConvolutionFit(int iFile, int iTime) {return _frtmConvolution[iFile][iTime];}
    bool isValidID(int iRun, int iType, QChar eventID);
};

#endif // PETRTM_H
