#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <QMainWindow>
#include <QWidget>
#include "plot.h"
#include "simengine.h"
#include "petrtm.h"

class Fitter
{
private:
    dVector _RRTimeVector;
    dVector _RRVector;
    int _nPoly = 3;
    int _nPar  = 3;
    PolynomialGLM _polyPlusGammaForRR; // fit the RR as a polynomial plus a gamma function
    dVector _gammaParVal;
    dVector _gammaParInc;
    bVector _gammaParAdj;
    bool _optHigh=false;

    inline double getCostFunction() {return _polyPlusGammaForRR.getSigma2();}
    void lineScan1D( int iPar, double &costRelative, double &incrementOpt );
    double fitRRComputeCost();

public:
    void init(dVector RRTimeVector, dVector RRVector);
    void setPolynomial(int nPoly);
    void setFitStage(int stage);
    void fitTAC(double toleranceCost);
    void fitTau();
    void fitTauAndOnset();
    void fitTauAndAlpha();

    // getters
    inline double getParValue(int iPar) {return _gammaParVal.at(iPar);}
    inline int getNumberTimeBins() {return _RRTimeVector.size();}
    inline double getFit(int iTime) {return _polyPlusGammaForRR.getFit(iTime);}
};

class SimWindow : public QMainWindow
{
    Q_OBJECT
public:
    SimWindow();
private:
    simEngine _simulator;
    PETRTM _PETRTM;
    // Import real data as ROIs from a table file
    QStringList _dataColumnNames; // column names
    dMatrix _dataTable;           // data table [time][columns]
    iVector _dtBinsSec;           // dt steps in time bins in seconds
    dVector _timeBins;            // time for bins (center of bin)

    double _BPndLowValue  = 1.0;
    double _BPndHighValue = 5.0;
    double _BPndStepValue = 0.5;
    double _timeLowValue  = 20.;
    double _timeHighValue = 80.;
    int _nThreads = 1;
    int _numberSimulationsPerThread = 1;

    dVector _BP0Vector; // [# BP0 values]
    dMatrix _errBPndMatrix, _errChallMatrix, _tau2RefMatrix, _errTau4Matrix; // [# BP0 values][nSimulations]

    Fitter _fitRR;

    QStringList _validBinSizeName = {"dt","delta-time","delta_time","deltaTime",
                                     "bin-size","bin_size","binSize",
                                     "bin-time","bin_time","binTime",
                                     "frame-duration","frame_duration","frameDuration",
                                     "frame-size","frame_size","frameSize"};
    QStringList _validBinTimeName = {"t","time","time-point","time_point", "timePoint",
                                     "frame-time","frame_time","frameTime"};

    QTabWidget *_tabTimeSpace;
    QWidget *_setupPage;
    QWidget *_targetPage;
    QWidget *_sweepBPndPage;
    QWidget *_sweepTimePage;

    QComboBox *_threadsComboBox;
    QStatusBar *_statusBar;
    QProgressBar *_progressBar;
    QStatusBar *_plasmaStatusBar;
    QStatusBar *_RRStatusBar;
    QStatusBar *_TRStatusBar;

    // plots
    plotData *_plasmaPlot;   // setup page
    plotData *_RRPlot;       // setup page

    plotData *_basisPlot;    // target page
    plotData *_targetPlot;   // target page

    plotData *_errBPndPlot;  // vsBPnd page
    plotData *_errChallPlot; // vsBPnd page
    plotData *_tau2RefPlot;  // vsBPnd page
    plotData *_errk4Plot;

    plotData *_errBPndOrChallVsTimePlot;

    // setup page
    QComboBox *_whichPlasmaPlot;
    QVBoxLayout *_setupPlotLayout;
    QVBoxLayout *_targetPlotLayout;
    // timing
    QLineEdit *_numberTimeBins;
    QSpinBox *_binIndex;
    QLineEdit *_binDuration;
    QLineEdit *_timeStep;
    QLineEdit *_subSample;
    QCheckBox *_applyToAllBinDuration;
    // bolus in/out
    QLineEdit *_bolusMag;
    QLineEdit *_tauDecay;
    QLineEdit *_KBol;
    QLineEdit *_fastTau;
    QLineEdit *_slowTau;
    QLineEdit *_fastFraction;

    QCheckBox *_bolusMagCheckBox;
    QCheckBox *_tauDecayCheckBox;
    QCheckBox *_KBolCheckBox;
    QCheckBox *_fastTauCheckBox;
    QCheckBox *_slowTauCheckBox;
    QCheckBox *_fastFractionCheckBox;
    // reference region
    QLineEdit *_tau2Ref;
    QLineEdit *_tau1Ref;
    QLineEdit *_plasmaFracRef;
    QLineEdit *_noiseRef;
    QComboBox *_dataRefRegion;
    QRadioButton *_startWithPlasma;
    QRadioButton *_startWithRR;
    QPushButton *_calcRRMatch;
    QPushButton *_readROIFile;
    QLabel *_ROIFileName;

    // target region input
    QGroupBox *_targetSimulationGroupBox;
    QGroupBox *_targetDataGroupBox;
    QRadioButton *_analyzeRealData;
    QRadioButton *_analyzeSimulation;
    QComboBox *_dataTargetRegion;
    QLineEdit *_BPnd;
    QLineEdit *_R1;
    QLineEdit *_tau4;
    QLineEdit *_plasmaFracTar;
    QLineEdit *_noiseTar;
    QLineEdit *_challengeTime;
    QLineEdit *_challengeMag;
    // target region analysis
    QComboBox *_modelType;
    QComboBox *_weightType;
    QLineEdit *_ignoreString;
    QLineEdit *_tau2RefAnalysis;
    QLineEdit *_tau4Analysis;
    QLabel *_tau2RefAnalysisLabel;
    QLabel *_tau4AnalysisLabel;
    QCheckBox *_fitk4;
    QCheckBox *_fitChallenge;
    // errors
    QLabel *_errorBPnd;
    QLabel *_errork2;
    QLabel *_errork2a;
    QLabel *_errordk2a;
    QLabel *_errorR1;
    QLabel *_errorR1Label;
    QLabel *_errorTau2Ref;
    QLabel *_errorTau2RefLabel;
    QLabel *_errorChallenge;
    QLabel *_errorChallengeLabel;
    QLabel *_errork4;
    QLabel *_errork4Label;
    QLabel *_sigma;
    QLabel *_sigmaLabel;
    QLabel *_errordk2aLabel;

    // vsBPnd page
    QLineEdit *_BPndLow;
    QLineEdit *_BPndHigh;
    QLineEdit *_BPndStep;

    QLabel *_nSamplesBPndPerThreadLabel;
    QLineEdit *_nSamplesBPndPerThread;
    QLabel *_nSamplesBPndLabel;
    QLabel *_nSamplesBPnd;
    QPushButton *_calculateBPndCurves;
    QPushButton *_clearBPndCurves;
    QCheckBox *_checkBoxBPndErrGraph;
    QCheckBox *_checkBoxChallErrGraph;
    QCheckBox *_checkBoxk4ErrGraph;
    QCheckBox *_checkBoxTau2RefGraph;

    // vsTime page
    QLineEdit *_timeLow;
    QLineEdit *_timeHigh;
    QPushButton *_calculateTimeCurves;
    QPushButton *_clearTimeCurves;

    void createSetupPage();
    void createTargetPage();
    void createSweepBPndPage();
    void createVersusTimePage();

    void updatePlasmaGraph();
    void updateBasisGraph();
    void updateTargetGraph();
    void enablePlasmaMatching(bool state);
    void addSimulationCurveRR(plotData *whichPlot);
    void addDataCurveRR(plotData *whichPlot);
    void addSimulationCurveTarget();
    void addDataCurveTarget();
    void defineRTMModel();

    void scaleSimulationToDataAverage();
    void analyzeTAC();
    QString analyzeString(double truth, double guess);
    inline double percentageError(double guess, double truth) {return 100.*(guess/truth-1.);}
    double getChallengeMagFromAnalysis();
    double bestTau2RefForRTM2();
    void finishedLieDetectorAllThreads();
    void setThreadVisibility(bool state);
    double calculateMean(dVector vec);
    double calculateStDev(double mean, dVector vec);
    void enableComboBoxItem(QComboBox *comboBox, int itemNumber, bool enable);
    inline bool realDataAvailable() {return _dataTable.size() != 0 && _dataRefRegion->count() != 0;}
    inline bool analyzeRealData() {return _analyzeRealData->isChecked();}

private slots:
    void exitApp();
    void aboutApp();
    void aboutROI();
    void updateAllGraphs();
    void changedNumberThreads(int indexInBox);
    void showPlasmaRR();
    void showPlasma();
    void showRR();
    void showBasisTarget();
    void showBasis();
    void showTarget();
    inline void changedDataReferenceRegion() {updateAllGraphs();}
    inline void changedDataTargetRegion()    {updateAllGraphs();}
    void clickedAnalyzeStimulation(bool state);
    void clickedAnalyzeRealData(bool state);
    void updateReferenceGraph();

    void changedNumberBins();
    void changedBinIndex(int indexPlusOne);
    void changedBinDuration();
    void changedApplyToAllBinDuration(bool state);
    void changedSubSample();
    void changedBolusMag();
    void changedTauBolus();
    void changedInfusion();

    void changedFastElimination();
    void changedSlowElimination();
    void changedFastEliminationFraction();

    void changedTau2Ref();
    void changedTau1Ref();
    void changedPlasmaFracRef();
    void changedNoiseRef();

    void changedWeightType(int indexInBox);
    void changedModelType(int indexInBox);

    void changedBPND();
    void changedR1();
    void changedTau4();
    void changedChallengeTime();
    void changedChallengeMag();
    void changedPlasmaFracTar();
    void changedNoiseTar();

    void calculateRRMatch();

    void changedIgnoreString();
    void changedTau2RefAnalysis();
    void changedTau4Analysis();
    void changedCheckBoxFitk4(bool state);
    void changedCheckBoxChallenge(bool state);

    void changedBPndLow();
    void changedBPndHigh();
    void changedBPndStep();
    void changedNumberSimulationsBPnd();
    void calculateBPndCurves();
    void calculateBPndCurvesInThreads();
    void clearBPndCurves();
    void changedVersusBPndGraphs();

    void changedTimeLow();
    void changedTimeHigh();
    void calculateTimeCurves();
    void clearTimeCurves();

    void getTableDataFile();
    bool defineTimeBinsFromBinSize(QStringList validBinName, int &iColumn);
    bool defineTimeBinsFromTimePoints(QStringList validBinName, int &iColumn);

    void startWithRR(bool state);

public slots:
    void updateLieDetectorProgress(int iProgress);
    void finishedLieDetectorOneThread(dMatrix errBPnd, dMatrix errChall, dMatrix tau2Ref, dMatrix errTau4);
};

#endif // SIMWINDOW_H
