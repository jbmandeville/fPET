#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <QMainWindow>
#include <QWidget>
#include "plot.h"
#include "simengine.h"
#include "fmriglm.h"

class SimWindow : public QMainWindow
{
    Q_OBJECT
public:
    SimWindow();
private:
    simEngine _simulator;
    fMRIGLM _fMRIGLM;
    // Import real data as ROIs from a table file
    QStringList _dataColumnNames; // column names
    dMatrix _dataTable;           // data table [columns][time]
    iVector _dtBinsSec;           // dt steps in time bins in seconds
    dVector _timeBins;            // time for bins (center of bin)

    int _nThreads = 1;
    int _numberSimulationsPerThread = 1;

    LOESS _quadLOESS;

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

    QComboBox *_threadsComboBox;
    QStatusBar *_statusBar;
    QProgressBar *_progressBar;
    QStatusBar *_plasmaStatusBar;
    QStatusBar *_TRStatusBar;

    // plots
    plotData *_plotPlasma;   // setup page
    // target page
    plotData *_plotTarget;
    plotData *_plotTarget1;  // plot on plasma page
    plotData *_plotBasis;

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
    QPushButton *_readROIFile;
    QLabel *_ROIFileName;

    // target region input
    QGroupBox *_targetSimulationGroupBox;
    QGroupBox *_targetDataGroupBox;
    QRadioButton *_analyzeRealData;
    QRadioButton *_analyzeSimulation;
    QComboBox *_dataTargetRegion;
    QLineEdit *_K1;
    QLineEdit *_k2;
    QLineEdit *_k3;
    QLineEdit *_plasmaFracTar;
    QLineEdit *_noiseTar;
    QLineEdit *_challengeOnsetTime1;
    QLineEdit *_challengeOffsetTime1;
    QLineEdit *_challengeOnsetTime2;
    QLineEdit *_challengeOffsetTime2;
    QLineEdit *_challengeMagPercent1;
    QLineEdit *_challengeMagPercent2;
    // target region analysis
    QComboBox *_modelType;
    QComboBox *_weightType;
    QComboBox *_baselineTerms;
    QLineEdit *_tauAnalysis;
    QLineEdit *_ignoreString;
    QGroupBox *_analysisGroupBox;
    QCheckBox *_fitChallenge;
    QCheckBox *_removeBaseline;

    QRadioButton *_signalRaw;
    QRadioButton *_signalDerivative;

    QPushButton *_runSimulationAndAnalysis;

    // errors
    QLabel *_errork3;
    QLabel *_errork2;
    QLabel *_errorK1;
    QLabel *_errorK1Label;
    QLabel *_errorChallenge;
    QLabel *_errorChallengeLabel;
    QLabel *_sigma;
    QLabel *_sigmaLabel;

    void createPlasmaPage();
    void createTargetPage();

    void updatePlasmaGraph();
    void updateBasisGraph();
    void updateTargetGraph(dMatrix yData, dMatrix yFit);

    void analyzeTAC();
    QString analyzeString(double truth, double guess);
    inline double percentageError(double guess, double truth) {return 100.*(guess/truth-1.);}
    void setThreadVisibility(bool state);
    double calculateMean(dVector vec);
    double calculateStDev(double mean, dVector vec);
    void enableComboBoxItem(QComboBox *comboBox, int itemNumber, bool enable);
    inline bool realDataAvailable() {return _dataTable.size() != 0;}
    inline bool analyzeRealData() {return _analyzeRealData->isChecked();}
    void differentiateSignal(dMatrix &signal);

    // getters

private slots:
    void exitApp();
    void aboutApp();
    void aboutROI();
    void updateAllGraphs();
    void changedNumberThreads(int indexInBox);
    void showBasis(bool state);
    void showTarget(bool state);
    inline void changedDataTargetRegion()    {updateAllGraphs();}
    void clickedAnalyzeStimulation(bool state);
    void clickedAnalyzeRealData(bool state);

    void changedNumberBins();
    void changedBinIndex(int indexPlusOne);
    void changedBinDuration();
    void changedSubSample();
    void changedBolusMag();
    void changedTauBolus();
    void changedInfusion();
    void changedInfusionDelay();

    void changedFastElimination();
    void changedSlowElimination();
    void changedFastEliminationFraction();

    void changedWeightType(int indexInBox);
    void changedModelType(int indexInBox);
    void changedBaselineTerms(int indexInBox);

    void changedK1();
    void changedk2();
    void changedk3();
    void changedChallengeOnsetTime1();
    void changedChallengeOffsetTime1();
    void changedChallengeOnsetTime2();
    void changedChallengeOffsetTime2();
    void changedChallengeMag1();
    void changedChallengeMag2();
    void changedPlasmaFracTar();
    void changedNoiseTar();

    void changedTauAnalysis();
    void changedIgnoreString();
    void changedCheckBoxChallenge(bool state);
    void changedCheckBoxBaseline(bool state);
    void changedDerivativeRadioButton(bool state);

    void getTableDataFile();
    bool defineTimeBinsFromBinSize(QStringList validBinName, int &iColumn);
    bool defineTimeBinsFromTimePoints(QStringList validBinName, int &iColumn);

public slots:
};
#endif // SIMWINDOW_H
