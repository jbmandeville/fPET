#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <QMainWindow>
#include <QWidget>
#include "plot.h"
#include "simengine.h"

class SimWindow : public QMainWindow
{
    Q_OBJECT
public:
    SimWindow();
private:
    simEngine _simulator;
    QStringList _dataColumnNames;
    dMatrix _dataTable;

    QTabWidget *_tabTimeSpace;
    QWidget *_setupPage;
    QWidget *_targetPage;

    // setup page
    plotData *_plasmaPlot;
    plotData *_RRPlot;
    QVBoxLayout *_setupPlotLayout;
    // timing
    QLineEdit *_timeDuration;
    QLineEdit *_timeStep;
    QLineEdit *_downSample;
    // bolus in/out
    QLineEdit *_bolusMag;
    QLineEdit *_tauDecay;
    QLineEdit *_alphaDecay;
    QLineEdit *_infusion;
    QLineEdit *_fastTau;
    QLineEdit *_slowTau;
    QLineEdit *_fastFraction;
    // reference region
    QLineEdit *_tau2Ref;
    QLineEdit *_tau1Ref;
    QLineEdit *_noiseRef;
    QComboBox *_dataRefRegion;
    // target region
    QLineEdit *_BPnd;
    QLineEdit *_tau4;
    QLineEdit *_noiseTar;
    QLineEdit *_challengeTime;
    QLineEdit *_challengeMag;

    // targetPage
    plotData *_targetPlot;

    void createSetupPage();
    void createTargetPage();

    void updatePlasmaGraph();
    void updateReferenceGraph();
    void updateTargetGraph();
    void updateAllGraphs();
    void getTableDataFile();
    QString readTableFile(QString fileName, QStringList &columnNames, dMatrix &table);

private slots:
    void changedGraphSizes(int iSelection);

    void changedTimeDuration();
    void changedTimeStep();
    void changedDownSample();
    void changedBolusMag();
    void changedTauBolus();
    void changedAlphaBolus();
    void changedInfusion();

    void changedFastElimination();
    void changedSlowElimination();
    void changedFastEliminationFraction();

    void changedTau2Ref();
    void changedTau1Ref();
    void changedNoiseRef();
    void changedDataRefRegion(int indexInBox);

    void changedBPND();
    void changedTau4();
    void changedChallengeTime();
    void changedChallengeMag();
    void changedNoiseTar();

};

#endif // SIMWINDOW_H
