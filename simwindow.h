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
    QStringList _dataColumnNames;
    dMatrix _dataTable;
    bool _analyzeRealData = false;

    double _BPndLowValue  = 1.0;
    double _BPndHighValue = 5.0;
    double _BPndStepValue = 0.5;
    double _timeLowValue  = 20.;
    double _timeHighValue = 80.;
    int _nThreads = 1;
    int _numberSimulationsPerThread = 1;

    dVector _BP0Vector; // [# BP0 values]
    dMatrix _errBPndMatrix, _errChallMatrix, _tau2RefMatrix, _errTau4Matrix; // [# BP0 values][nSimulations]

    QTabWidget *_tabTimeSpace;
    QWidget *_setupPage;
    QWidget *_targetPage;
    QWidget *_sweepBPndPage;
    QWidget *_sweepTimePage;

    QComboBox *_threadsComboBox;
    QStatusBar *_statusBar;
    QProgressBar *_progressBar;
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
    QVBoxLayout *_setupPlotLayout;
    QVBoxLayout *_targetPlotLayout;
    // timing
    QLineEdit *_numberTimeBins;
    QSpinBox *_binIndex;
    QLineEdit *_binDuration;
    QLineEdit *_timeStep;
    QLineEdit *_subSample;
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
    QPushButton *_calcRRMatch;
    QPushButton *_readROIFile;
    QLabel *_ROIFileName;

    // target region input
    QRadioButton *_analyzeTarget;
    QRadioButton *_analyzeSimulation;
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
    void updateReferenceGraph();
    void updateBasisGraph();
    void updateTargetGraph();
    void updateAllGraphs();
    QString readTableFile(QString fileName, QStringList &columnNames, dMatrix &table);
    void enablePlasmaMatching(bool state);
    void addSimulationCurveRR();
    void addDataCurveRR();
    void defineRTMModel();

    void analyzeSimulatedTAC();
    QString analyzeString(double truth, double guess);
    inline double percentageError(double guess, double truth) {return 100.*(guess/truth-1.);}
    double getChallengeMagFromAnalysis();
    double bestTau2RefForRTM2();
    void finishedLieDetectorAllThreads();
    void setThreadVisibility(bool state);
    double calculateMean(dVector vec);
    double calculateStDev(double mean, dVector vec);
    void enableComboBoxItem(QComboBox *comboBox, int itemNumber, bool enable);
    inline bool realDataAvailable() {return _dataTable.size() != 0;}

private slots:
    void changedNumberThreads(int indexInBox);
    void showPlasmaRR();
    void showPlasma();
    void showRR();
    void showBasisTarget();
    void showBasis();
    void showTarget();

    void changedNumberBins();
    void changedBinDuration();
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
    void changedDataRefRegion(int indexInBox);

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

public slots:
    void updateLieDetectorProgress(int iProgress);
    void finishedLieDetectorOneThread(dMatrix errBPnd, dMatrix errChall, dMatrix tau2Ref, dMatrix errTau4);
};

#endif // SIMWINDOW_H
