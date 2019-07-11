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
    double _BPndLowValue  = 1.0;
    double _BPndHighValue = 5.0;
    double _BPndStepValue = 0.5;
    double _timeLowValue  = 20.;
    double _timeHighValue = 80.;

    QTabWidget *_tabTimeSpace;
    QWidget *_setupPage;
    QWidget *_targetPage;
    QWidget *_vsBPndPage;
    QWidget *_vsTimePage;

    // plots
    plotData *_plasmaPlot;   // setup page
    plotData *_RRPlot;       // setup page

    plotData *_basisPlot;    // target page
    plotData *_targetPlot;   // target page

    plotData *_errBPndPlot;  // vsBPnd page
    plotData *_errChallPlot; // vsBPnd page
    plotData *_tau2RefPlot;  // vsBPnd page

    plotData *_errBPndOrChallVsTimePlot;

    // setup page

    QVBoxLayout *_setupPlotLayout;
    // timing
    QLineEdit *_timeDuration;
    QLineEdit *_timeStep;
    QLineEdit *_downSample;
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
    QLineEdit *_noiseRef;
    QComboBox *_dataRefRegion;
    QPushButton *_calcRRMatch;

    // target region input
    QLineEdit *_BPnd;
    QLineEdit *_R1;
    QLineEdit *_tau4;
    QLineEdit *_noiseTar;
    QLineEdit *_challengeTime;
    QLineEdit *_challengeMag;
    // target region analysis
    QComboBox *_modelType;
    QComboBox *_weightType;
    QLineEdit *_tau2RefAnalysis;
    QLineEdit *_tau4Analysis;
    QLabel *_tau2RefAnalysisLabel;
    QLabel *_tau4AnalysisLabel;
    QCheckBox *_includeChallenge;
    // errors
    QLabel *_errorBPnd;
    QLabel *_errork2;
    QLabel *_errork2a;
    QLabel *_errorR1;
    QLabel *_errorR1Label;
    QLabel *_errorTau2Ref;
    QLabel *_errorTau2RefLabel;
    QLabel *_errorChallenge;
    QLabel *_errorChallengeLabel;

    // vsBPnd page
    QLineEdit *_BPndLow;
    QLineEdit *_BPndHigh;
    QLineEdit *_BPndStep;
    QPushButton *_calculateBPndCurves;
    QPushButton *_clearBPndCurves;
    QCheckBox *_checkBoxBPndErrGraph;
    QCheckBox *_checkBoxChallErrGraph;
    QCheckBox *_checkBoxTau2RefGraph;

    // vsTime page
    QLineEdit *_timeLow;
    QLineEdit *_timeHigh;
    QPushButton *_calculateTimeCurves;
    QPushButton *_clearTimeCurves;

    void createSetupPage();
    void createTargetPage();
    void createVersusBPndPage();
    void createVersusTimePage();

    void updatePlasmaGraph();
    void updateReferenceGraph();
    void updateBasisGraph();
    void updateTargetGraph();
    void updateAllGraphs();
    void getTableDataFile();
    QString readTableFile(QString fileName, QStringList &columnNames, dMatrix &table);
    void enablePlasmaMatching(bool state);

    void analyzeSimulatedTAC();
    QString analyzeString(double truth, double guess);
    inline double percentageError(double truth, double guess) {return 100.*(guess/truth-1.);}
    double getChallengeMagFromAnalysis();
    double bestTau2RefForRTM2();

private slots:
    void changedGraphSizes(int iSelection);

    void changedTimeDuration();
    void changedTimeStep();
    void changedDownSample();
    void changedBolusMag();
    void changedTauBolus();
    void changedInfusion();

    void changedFastElimination();
    void changedSlowElimination();
    void changedFastEliminationFraction();

    void changedTau2Ref();
    void changedTau1Ref();
    void changedNoiseRef();
    void changedDataRefRegion(int indexInBox);

    void changedWeightType(int indexInBox);
    void changedModelType(int indexInBox);

    void changedBPND();
    void changedR1();
    void changedTau4();
    void changedChallengeTime();
    void changedChallengeMag();
    void changedNoiseTar();

    void calculateRRMatch();

    void changedTau2RefAnalysis();
    void changedTau4Analysis();
    void changedCheckBoxChallenge(bool state);

    void changedBPndLow();
    void changedBPndHigh();
    void changedBPndStep();
    void calculateBPndCurves();
    void clearBPndCurves();
    void changedVersusBPndGraphs();

    void changedTimeLow();
    void changedTimeHigh();
    void calculateTimeCurves();
    void clearTimeCurves();
};

#endif // SIMWINDOW_H
