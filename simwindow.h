#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <QMainWindow>
#include <QWidget>
#include "plot.h"
#include "simengine.h"
#include "petrtm.h"

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
    dMatrix _dataTable;           // data table [columns][time]
    iVector _dtBinsSec;           // dt steps in time bins in seconds
    dVector _timeBins;            // time for bins (center of bin)

    double _BPndLowValue  = 1.0;
    double _BPndHighValue = 5.0;
    double _BPndStepValue = 0.5;
    double _timeLowValue  = 20.;
    double _timeHighValue = 80.;
    double _tau4HighValue = 30.;
    double _tau4StepValue = 1.;
    int _nThreads = 1;
    int _numberSimulationsPerThread = 1;

    dVector _BP0Vector; // [# BP0 values]
    dMatrix _errBPndMatrix, _errChallMatrix, _tau2RefMatrix, _errTau4Matrix, _errDVMatrix; // [# BP0 values][nSimulations]

    dVector _tau4Vector;
    dMatrix _AICMatrix;

    LOESS _quadLOESS;

    QStringList _validBinSizeName = {"dt","delta-time","delta_time","deltaTime",
                                     "bin-size","bin_size","binSize",
                                     "bin-time","bin_time","binTime",
                                     "frame-duration","frame_duration","frameDuration",
                                     "frame-size","frame_size","frameSize"};
    QStringList _validBinTimeName = {"t","time","time-point","time_point", "timePoint",
                                     "frame-time","frame_time","frameTime"};

    struct global_fit
    {
        dVector tau4Vector;
        dVector tau2RefVector;
        dVector DVVector;
        dVector AICVector1D;
        dMatrix3 AICMatrix;
        double  AICSRTM;
    } _globalFit;


    QTabWidget *_tabTimeSpace;
    QWidget *_setupPage;
    QWidget *_targetPage;
    QWidget *_sweepBPndPage;
    QWidget *_sweepTimePage;
    QWidget *_sweepTau4Page;

    QComboBox *_threadsComboBox;
    QStatusBar *_statusBar;
    QProgressBar *_progressBar;
    QStatusBar *_plasmaStatusBar;
    QStatusBar *_RRStatusBar;
    QStatusBar *_TRStatusBar;

    // plots
    plotData *_plotPlasma;   // setup page
    plotData *_plotRR;       // setup page

    plotData *_plotBasis;
    plotData *_plotTarget;
    plotData *_plotAIC;

    plotData *_plotErrBPndVsBPnd;  // vsBPnd page
    plotData *_plotErrChallVsBPnd; // vsBPnd page
    plotData *_plotTau2RefVsBPnd;  // vsBPnd page
    plotData *_plotErrk4VsBPnd;
    plotData *_plotErrDVVsBPnd;

    plotData *_plotErrBPndOrChallVsTime;  // vs time page

    plotData *_plotErrBPndVsTau4;  // vs tau4 page
    plotData *_plotErrChallVsTau4; // vs tau4 page
    plotData *_plotAICVsTau4;  // vs tau4 page

    // setup page
    QComboBox *_whichPlasmaPlot;
    QCheckBox *_clearPlasmaPlot;
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
    QLineEdit *_KBolDelay;
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
    QComboBox *_simulationStartingPoint;
    QPushButton *_calcRRMatch;
    QPushButton *_readROIFile;
    QLabel *_ROIFileName;

    // target region input
    QGroupBox *_targetSimulationGroupBox;
    QGroupBox *_targetDataGroupBox;
    QRadioButton *_radioShowAIC;
    QRadioButton *_analyzeRealData;
    QRadioButton *_analyzeSimulation;
    QComboBox *_dataTargetRegion;
    QLineEdit *_BPnd;
    QLineEdit *_R1;
    QLineEdit *_tau4;
    QLineEdit *_DV;
    QLineEdit *_plasmaFracTar;
    QLineEdit *_noiseTar;
    QLineEdit *_challengeTime;
    QLabel *_challengeMagLabel;
    QLineEdit *_challengeMagPercent;
    // target region analysis
    QComboBox *_modelType;
    QComboBox *_weightType;
    QLineEdit *_ignoreString;
    QLineEdit *_tau2RefAnalysis;
    QLineEdit *_tau4Analysis;
    QLineEdit *_DVAnalysis;
    QGroupBox *_analysisGroupBox;
    QLabel *_tau2RefAnalysisLabel;
    QLabel *_tau4AnalysisLabel;
    QLabel *_DVAnalysisLabel;
    QCheckBox *_fitChallenge;

    QRadioButton *_fit2ParRadioButton;
    QRadioButton *_fitk4RadioButton;
    QRadioButton *_fitk2RefRadioButton;
    QRadioButton *_fitk4k2RefRadioButton;
    QRadioButton *_fitDVRadioButton;
    QPushButton *_runSimulationAndAnalysis;

    QLineEdit *_parSearchStart;
    QLineEdit *_parSearchEnd;
    QLineEdit *_parSearchStep;
    // errors
    QLabel *_errorBPnd;
    QLabel *_errork2;
    QLabel *_errorR1;
    QLabel *_errorR1Label;
    QLabel *_errorDV;
    QLabel *_errorDVLabel;
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

    // vsk4 page
    QLineEdit *_tau4High;
    QLineEdit *_tau4Step;

    QLabel *_nSamplesBPndPerThreadLabel;
    QLineEdit *_nSamplesBPndPerThread;
    QLabel *_nSamplesBPndLabel;
    QLabel *_nSamplesBPnd;
    QLabel *_sigma2Label;
    QPushButton *_calculateBPndCurves;
    QPushButton *_clearBPndCurves;
    QCheckBox *_checkBoxBPndErrVsBPnd;
    QCheckBox *_checkBoxChallErrVsBPnd;
    QCheckBox *_checkBoxk4ErrGraph;
    QCheckBox *_checkBoxDVErrGraph;
    QCheckBox *_checkBoxTau2RefGraph;

    // vsTime page
    QLineEdit *_timeLow;
    QLineEdit *_timeHigh;
    QPushButton *_calculateTimeCurves;
    QPushButton *_clearTimeCurves;

    // vs k4 page
    QCheckBox *_checkBoxBPndErrVsTau4;
    QCheckBox *_checkBoxChallErrVsTau4;
    QCheckBox *_checkBoxAICVsTau4;
    QPushButton *_calculateTau4Curves;
    QPushButton *_clearTau4Curves;

    void createSetupPage();
    void createTargetPage();
    void createSweepBPndPage();
    void createSweepTimePage();
    void createSweepTau4Page();
    void clearErrorMatrices(bool sweepBP);

    void updatePlasmaGraph();
    void updateBasisGraph();
    void updateTargetGraph();
    void enablePlasmaMatching(bool state);
    void addSimulationRR(plotData *whichPlot);
    void addFitRR(plotData *whichPlot);
    void addDataCurveRR(plotData *whichPlot);
    void addSimulationTarget();
    void addDataCurveTarget();
    void defineRTMModel();
    void setReferenceRegionForAnalysis();

    void scaleSimulationToDataAverage();
    void analyzeTAC();
    QString analyzeString(double truth, double guess);
    inline double percentageError(double guess, double truth) {return 100.*(guess/truth-1.);}
    double bestTau2RefForRTM2();
    void finishedLieDetectorBPndAllThreads();
    void finishedLieDetectorTau4AllThreads();
    void setThreadVisibility(bool state);
    double calculateMean(dVector vec);
    double calculateStDev(double mean, dVector vec);
    void enableComboBoxItem(QComboBox *comboBox, int itemNumber, bool enable);
    inline bool realDataAvailable() {return _dataTable.size() != 0 && _dataRefRegion->count() != 0;}
    inline bool analyzeRealData() {return _analyzeRealData->isChecked();}
    void calculateForwardTau4andTau2Ref();

    // getters
    double getChallengeMagPercentFromAnalysis(double &deltaBPAbs);
    inline bool simStartsFromPlasma()    {return _simulationStartingPoint->currentIndex() == simStart_fromPlasma;}
    inline bool simStartsFromPlasmaFit() {return _simulationStartingPoint->currentIndex() == simStart_fromPlasmaFit;}
    inline bool simStartsFromDataFit()   {return _simulationStartingPoint->currentIndex() == simStart_fromDataFit;}

private slots:
    void exitApp();
    void aboutApp();
    void aboutROI();
    void updateAllGraphs();
    void changedNumberThreads(int indexInBox);
    void showPlasma();
    void showRR();
    void showBasis();
    void showTarget();
    void showAICPlot();
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
    void changedInfusionDelay();

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
    void changedDV();
    void changedChallengeTime();
    void changedChallengeMag();
    void changedPlasmaFracTar();
    void changedNoiseTar();

    void runSimulationAndAnalysis();
    void calculateRRMatch();

    void changedIgnoreString();
    void changedTau2RefAnalysis();
    void changedTau4Analysis();
    void changedDVAnalysis();
    void changedCheckBoxFitk4(bool state);
    void changedCheckBoxFitk2Ref(bool state);
    void changedCheckBoxFitDV(bool state);
    void changedCheckBoxChallenge(bool state);
    void changedParSearchStart();
    void changedParSearchEnd();
    void changedParSearchStep();

    void changedBPndLow();
    void changedBPndHigh();
    void changedBPndStep();
    void changedNumberSimulationsBPnd();
    void calculateBPndCurves();
    void calculateBPndCurvesInThreads();
    void clearBPndCurves();
    void changedVersusBPndGraphs();
    void changedVersusTau4Graphs();
    void updateAICvsTau4Graph();

    void changedTau4High();
    void changedTau4Step();

    void changedTimeLow();
    void changedTimeHigh();
    void calculateTimeCurves();
    void clearTimeCurves();
    inline void clearPlasmaPlot(bool isChecked) {if (isChecked) updateAllGraphs();}

    void getTableDataFile();
    bool defineTimeBinsFromBinSize(QStringList validBinName, int &iColumn);
    bool defineTimeBinsFromTimePoints(QStringList validBinName, int &iColumn);

    void changedSimulationStartingPoint();

    void calculateTau4Curves();
    void calculateTau4CurvesInThreads();
    void clearTau4Curves();

public slots:
    void updateLieDetectorProgress(int iProgress);
    void finishedLieDetectorBPndOneThread(dMatrix errBPnd, dMatrix errChall, dMatrix tau2Ref,
                                          dMatrix errTau4, dMatrix errDV, double sigma2);
    void finishedLieDetectorTau4OneThread(dMatrix errBPnd, dMatrix errChall, dMatrix AIC, double sigma2);
    void globalk4FitterUpdate(dPoint4D tau4Tau2RefDVAIC);
    void globalk4FitterFinished(int mode);
};

class globalk4Fitter : public QObject, public QRunnable
{
    Q_OBJECT   // this can create a problem of "missing vtable"; fix by "run qmake" from Qt Creator build menu
public:
    globalk4Fitter(double tau4, dVector tau2RefVector, dVector DVVector, dVector ROI, const PETRTM petRTM);
    void run(); // if public, it can be run directly in the same thread to test thread safety
private:
    double _tau4;
    dVector _tau2RefVector;
    dVector _DVVector;
    dVector _ROI;     // [nTime]
    PETRTM _petRTM;
signals:
    void update(dPoint4D tau4Tau2RefDVAIC);
    void finished(int mode);
public slots:
};


#endif // SIMWINDOW_H
