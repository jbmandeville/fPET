#include "simwindow.h"
#include <QtWidgets>
#include <QDebug>
#include <QVector>
#include <QObject>
#include <QFileDialog>
#include <QString>

SimWindow::SimWindow()
{
    qDebug() << "SimWindow::SimWindow enter";
    _tabTimeSpace = new QTabWidget();
    createSetupPage();
    createTargetPage();
    createDependenciesPage();
    _tabTimeSpace->addTab(_setupPage, tr("setup"));
    _tabTimeSpace->addTab(_targetPage, tr("target"));
    _tabTimeSpace->addTab(_dependenciesPage, tr("dependencies"));

    QWidget *centralWidget = new QWidget(this);
    this->setCentralWidget( centralWidget );
    auto *mainLayout = new QVBoxLayout( centralWidget );
    mainLayout->addWidget(_tabTimeSpace);

    // add a menu
    QMenuBar *menuBar = new QMenuBar;
    mainLayout->setMenuBar(menuBar);
    QMenu *openMenu  = new QMenu(tr("&Open"), this);
    menuBar->addMenu(openMenu);
    QAction *openDataFileAction = openMenu->addAction(tr("Open table file (data)"));
    openDataFileAction->setShortcut(QKeySequence::Open);
    connect(openDataFileAction,  &QAction::triggered, this, &SimWindow::getTableDataFile);

    QSize defaultWindowSize;
    QRect rec = QApplication::desktop()->screenGeometry();
//    QRect rec = QApplication::desktop()->screenGeometry();
    defaultWindowSize.setWidth(rec.width()*2/3);
    defaultWindowSize.setHeight(rec.height()*2/3);
    resize(defaultWindowSize);

    qDebug() << "SimWindow::SimWindow updateAllGraphs";
    updateAllGraphs();

    qDebug() << "SimWindow::SimWindow exit";
}

void SimWindow::getTableDataFile()
{
    QString fileName;
    QFileDialog fileDialog;

    fileName = fileDialog.getOpenFileName(this,
                                          "Select an overlay file to open",
                                          QDir::currentPath(),
                                          "Overlay list files (overlay-list.* *.dat *.list)");
    if ( fileName.isEmpty() )
        return;
    else
    {
        QString errorString = readTableFile(fileName, _dataColumnNames, _dataTable);
        if ( !errorString.isEmpty() )
        {
            QMessageBox msgBox;
            msgBox.setText(errorString);
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.exec();
        }
        else
        {
            _dataRefRegion->clear();
            for (int jColumn=1; jColumn<_dataColumnNames.count(); jColumn++)
                _dataRefRegion->addItem(_dataColumnNames.at(jColumn));
            _dataRefRegion->setCurrentIndex(_dataRefRegion->count()-1);
            enablePlasmaMatching(true);
        }
    }
}

void SimWindow::enablePlasmaMatching(bool state)
{
    _calcRRMatch->setEnabled(state);
    _bolusMagCheckBox->setVisible(state);
    _tauDecayCheckBox->setVisible(state);
    _KBolCheckBox->setVisible(state);
    _fastTauCheckBox->setVisible(state);
    _slowTauCheckBox->setVisible(state);
    _fastFractionCheckBox->setVisible(state);

}

void SimWindow::createSetupPage()
{
    _setupPage = new QWidget();

    _plasmaPlot = new plotData(0);
    _RRPlot     = new plotData(1);
    _setupPlotLayout = new QVBoxLayout();
    _setupPlotLayout->addWidget(_plasmaPlot->getPlotSurface());
    _setupPlotLayout->addWidget(_RRPlot->getPlotSurface());
    changedGraphSizes(0);
    QString numberString;
    int editTextSize=80;

    // setup (duration, step, down-sampling
    auto *setupGroupBox = new QGroupBox("Setup the simulation");
    QLabel *timeDurationLabel = new QLabel("Duration  (min)");
    QLabel *timeStepLabel     = new QLabel("Time step (min)");
    QLabel *downSampleLabel   = new QLabel("Downsampling (int)");
    _timeDuration = new QLineEdit();
    _timeStep     = new QLineEdit();
    _downSample   = new QLineEdit();
    _timeDuration->setText(numberString.setNum(_simulator.getDuration()));
    _timeStep->setText(numberString.setNum(_simulator.getStepSize()));
    _downSample->setText(numberString.setNum(_simulator.getDownSampling()));
    _timeDuration->setFixedWidth(editTextSize);
    _timeStep->setFixedWidth(editTextSize);
    _downSample->setFixedWidth(editTextSize);
    auto *setupLayout = new QGridLayout();
    setupLayout->addWidget(timeDurationLabel,0,0);
    setupLayout->addWidget(_timeDuration,0,1);
    setupLayout->addWidget(timeStepLabel,1,0);
    setupLayout->addWidget(_timeStep,1,1);
    setupLayout->addWidget(downSampleLabel,2,0);
    setupLayout->addWidget(_downSample,2,1);
    setupGroupBox->setLayout(setupLayout);
    setupLayout->setSpacing(0);
    connect(_timeDuration, SIGNAL(editingFinished()), this, SLOT(changedTimeDuration()));
    connect(_timeStep,     SIGNAL(editingFinished()), this, SLOT(changedTimeStep()));
    connect(_downSample,   SIGNAL(editingFinished()), this, SLOT(changedDownSample()));

    // Plasma Input
    auto *plasmaInGroupBox     = new QGroupBox("plasma: adjust TAC to match ref. region");
    QLabel *bolusMagLabel      = new QLabel("Bolus Magnitude");
    QLabel *tauDecayLabel      = new QLabel("Bolus shape (tau)");
    QLabel *infusionLabel      = new QLabel("'Kbol' for BI (min)");
    QLabel *fastTauLabel       = new QLabel("Fast elimination (tau)");
    QLabel *slowTauLabel       = new QLabel("Slow elimination (tau)");
    QLabel *fastFractionLabel  = new QLabel("Fast fraction (<=1.)");
    QLabel *calcLabel       = new QLabel("Match to data");
    _bolusMag     = new QLineEdit();
    _tauDecay     = new QLineEdit();
    _KBol         = new QLineEdit();
    _fastTau      = new QLineEdit();
    _slowTau      = new QLineEdit();
    _fastFraction = new QLineEdit();
    _bolusMagCheckBox     = new QCheckBox();    _bolusMagCheckBox->setVisible(false);       _bolusMagCheckBox->setChecked(true);
    _tauDecayCheckBox     = new QCheckBox();    _tauDecayCheckBox->setVisible(false);       _tauDecayCheckBox->setChecked(true);
    _KBolCheckBox         = new QCheckBox();    _KBolCheckBox->setVisible(false);           _KBolCheckBox->setChecked(false);
    _fastTauCheckBox      = new QCheckBox();    _fastTauCheckBox->setVisible(false);        _fastTauCheckBox->setChecked(true);
    _slowTauCheckBox      = new QCheckBox();    _slowTauCheckBox->setVisible(false);        _slowTauCheckBox->setChecked(true);
    _fastFractionCheckBox = new QCheckBox();    _fastFractionCheckBox->setVisible(false);   _fastFractionCheckBox->setChecked(true);
    _calcRRMatch   = new QPushButton();
    QPixmap pixmapCalculate(":/My-Icons/calculator.png");
    QIcon calculatorIcon(pixmapCalculate);
    _calcRRMatch->setIcon(calculatorIcon);
    _calcRRMatch->setEnabled(false);
    _bolusMag->setText(numberString.setNum(_simulator.getMagBolus()));
    _tauDecay->setText(numberString.setNum(_simulator.getTauBolus()));
    _KBol->setText(numberString.setNum(_simulator.getKBol()));
    _fastTau->setText(numberString.setNum(_simulator.getTauFastElim()));
    _slowTau->setText(numberString.setNum(_simulator.getTauSlowElim()));
    _fastFraction->setText(numberString.setNum(_simulator.getFastFraction()));
    _bolusMag->setFixedWidth(editTextSize);
    _tauDecay->setFixedWidth(editTextSize);
    _KBol->setFixedWidth(editTextSize);
    _fastTau->setFixedWidth(editTextSize);
    _slowTau->setFixedWidth(editTextSize);
    _fastFraction->setFixedWidth(editTextSize);
    _KBol->setToolTip("The time at which the integral of infusion matches the bolus;\n'Kbol' is a misnomer: this is really 'Tbol';\nSet to 0 for no BI");
    auto *plasmaInLayout = new QGridLayout();
    plasmaInLayout->addWidget(bolusMagLabel,0,0);
    plasmaInLayout->addWidget(_bolusMag,0,1);
    plasmaInLayout->addWidget(_bolusMagCheckBox,0,2);
    plasmaInLayout->addWidget(tauDecayLabel,1,0);
    plasmaInLayout->addWidget(_tauDecay,1,1);
    plasmaInLayout->addWidget(_tauDecayCheckBox,1,2);
    plasmaInLayout->addWidget(infusionLabel,2,0);
    plasmaInLayout->addWidget(_KBol,2,1);
    plasmaInLayout->addWidget(_KBolCheckBox,2,2);
    plasmaInLayout->addWidget(fastTauLabel,3,0);
    plasmaInLayout->addWidget(_fastTau,3,1);
    plasmaInLayout->addWidget(_fastTauCheckBox,3,2);
    plasmaInLayout->addWidget(slowTauLabel,4,0);
    plasmaInLayout->addWidget(_slowTau,4,1);
    plasmaInLayout->addWidget(_slowTauCheckBox,4,2);
    plasmaInLayout->addWidget(fastFractionLabel,5,0);
    plasmaInLayout->addWidget(_fastFraction,5,1);
    plasmaInLayout->addWidget(_fastFractionCheckBox,5,2);
    plasmaInLayout->addWidget(calcLabel,6,0);
    plasmaInLayout->addWidget(_calcRRMatch,6,1);
    plasmaInGroupBox->setLayout(plasmaInLayout);
    plasmaInLayout->setSpacing(0);
    connect(_bolusMag, SIGNAL(editingFinished()), this, SLOT(changedBolusMag()));
    connect(_tauDecay, SIGNAL(editingFinished()), this, SLOT(changedTauBolus()));
    connect(_KBol, SIGNAL(editingFinished()), this, SLOT(changedInfusion()));
    connect(_fastTau, SIGNAL(editingFinished()), this, SLOT(changedFastElimination()));
    connect(_slowTau, SIGNAL(editingFinished()), this, SLOT(changedSlowElimination()));
    connect(_fastFraction, SIGNAL(editingFinished()), this, SLOT(changedFastEliminationFraction()));
    connect(_calcRRMatch, SIGNAL(pressed()), this, SLOT(calculateRRMatch()));

    // reference region
    auto *RRGroupBox        = new QGroupBox("Reference Region");
    QLabel *tau2RefLabel    = new QLabel("1/k2 (min)");
    QLabel *tau1RefLabel    = new QLabel("1/K1 (ml/cm^3/min)");
    QLabel *noiseLabel      = new QLabel("Noise");
    QLabel *dataLabel       = new QLabel("Data reference");
    _tau2Ref = new QLineEdit();
    _tau1Ref = new QLineEdit();
    _noiseRef   = new QLineEdit();
    _tau2Ref->setText(numberString.setNum(_simulator.getTau2Ref()));
    _tau1Ref->setText(numberString.setNum(_simulator.getTau1Ref()));
    _noiseRef->setText(numberString.setNum(_simulator.getNoiseTar()));
    _tau2Ref->setFixedWidth(editTextSize);
    _tau1Ref->setFixedWidth(editTextSize);
    _noiseRef->setFixedWidth(editTextSize);
    _dataRefRegion = new QComboBox();
    auto *RRLayout = new QGridLayout();
    RRLayout->addWidget(tau2RefLabel,0,0);
    RRLayout->addWidget(_tau2Ref,0,1);
    RRLayout->addWidget(tau1RefLabel,1,0);
    RRLayout->addWidget(_tau1Ref,1,1);
    RRLayout->addWidget(noiseLabel,2,0);
    RRLayout->addWidget(_noiseRef,2,1);
    RRLayout->addWidget(dataLabel,3,0);
    RRLayout->addWidget(_dataRefRegion,3,1);
    RRGroupBox->setLayout(RRLayout);
    RRLayout->setSpacing(0);
    connect(_tau2Ref,  SIGNAL(editingFinished()), this, SLOT(changedTau2Ref()));
    connect(_tau1Ref,  SIGNAL(editingFinished()), this, SLOT(changedTau1Ref()));
    connect(_noiseRef, SIGNAL(editingFinished()), this, SLOT(changedNoiseRef()));
    connect(_dataRefRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(changedDataRefRegion(int)));

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(setupGroupBox);
    rightLayout->addWidget(plasmaInGroupBox);
    rightLayout->addWidget(RRGroupBox);

    QHBoxLayout *fullLayout = new QHBoxLayout();
    fullLayout->addLayout(_setupPlotLayout);
    fullLayout->addLayout(rightLayout);
    fullLayout->setStretch(0,100);
    fullLayout->setStretch(1,1);

    /////////////// Plotting tool bars //////////////////////
    const QIcon *dragX = new QIcon(":/My-Icons/dragX.png");
    auto *dragXAction = new QAction(*dragX, tr("drag/zoom X axis"), this);
    const QIcon *dragY = new QIcon(":/My-Icons/dragY.png");
    auto *dragYAction = new QAction(*dragY, tr("drag/zoom Y axis"), this);
    const QIcon *rescaleXY = new QIcon(":/My-Icons/rescaleGraph.png");
    auto *rescaleXYAction = new QAction(*rescaleXY, tr("Auto-scale X and Y ranges"), this);
    dragXAction->setCheckable(true);
    dragYAction->setCheckable(true);
    rescaleXYAction->setCheckable(true);
    rescaleXYAction->setChecked(true);
    QActionGroup *graphButtons = new QActionGroup(this);
    graphButtons->addAction(dragXAction);
    graphButtons->addAction(dragYAction);
    graphButtons->addAction(rescaleXYAction);

    auto *graphSizes = new QComboBox();
    graphSizes->addItem("plasma/RR");
    graphSizes->addItem("plasma");
    graphSizes->addItem("RR");
    graphSizes->setMaximumWidth(150);  // important when adding atlas regions, which have very long names

    QToolBar *graphToolBar = new QToolBar("time tool bar");

    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addAction(rescaleXYAction);
    graphToolBar->addWidget(graphSizes);
    fullLayout->setMenuBar(graphToolBar);

    // Make toolbar connections
    connect(dragXAction,     SIGNAL(triggered(bool)), _plasmaPlot, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plasmaPlot, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plasmaPlot, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plasmaPlot, SLOT(setSelectPoints()));
    connect(dragXAction,     SIGNAL(triggered(bool)), _RRPlot, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _RRPlot, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _RRPlot, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _RRPlot, SLOT(setSelectPoints()));
    connect(graphSizes,      SIGNAL(currentIndexChanged(int)), this, SLOT(changedGraphSizes(int)));

    _setupPage->setLayout(fullLayout);
}

void SimWindow::calculateRRMatch()
{
    bVector includeParameter;
    includeParameter.append(_bolusMagCheckBox->checkState());
    includeParameter.append(_tauDecayCheckBox->checkState());
    includeParameter.append(_KBolCheckBox->checkState());
    includeParameter.append(_fastTauCheckBox->checkState());
    includeParameter.append(_slowTauCheckBox->checkState());
    includeParameter.append(_fastFractionCheckBox->checkState());
//    _simulator.calculateRRMatch(includeParamter, )
}

void SimWindow::changedDataRefRegion(int indexInBox)
{
    updateAllGraphs();
}

void SimWindow::changedWeightType(int indexInBox)
{
    _PETRTM.setWeightingModel(indexInBox);
    _PETRTM.setPrepared(false);
    updateAllGraphs();
}

void SimWindow::changedModelType(int indexInBox)
{
    if ( indexInBox == RTM_SRTM2 || indexInBox == RTM_rFRTM2 )
    {
        _tau2RefAnalysisLabel->setVisible(true);
        _tau2RefAnalysis->setVisible(true);
        _errorR1Label->setVisible(false);
        _errorR1->setVisible(false);
        _errorTau2RefLabel->setVisible(false);
        _errorTau2Ref->setVisible(false);
    }
    else
    {
        _tau2RefAnalysisLabel->setVisible(false);
        _tau2RefAnalysis->setVisible(false);
        _errorR1Label->setVisible(true);
        _errorR1->setVisible(true);
        _errorTau2RefLabel->setVisible(true);
        _errorTau2Ref->setVisible(true);
    }
    if ( indexInBox == RTM_rFRTM3 || indexInBox == RTM_rFRTM2 )
    {
        _tau4AnalysisLabel->setVisible(true);
        _tau4Analysis->setVisible(true);
    }
    else
    {
        _tau4AnalysisLabel->setVisible(false);
        _tau4Analysis->setVisible(false);
    }
    if ( indexInBox == RTM_SRTM3 )
        _PETRTM.setRTMModelType("SRTM3");
    else if ( indexInBox == RTM_SRTM2 )
        _PETRTM.setRTMModelType("SRTM2");
    else if ( indexInBox == RTM_rFRTM3 )
        _PETRTM.setRTMModelType("rFRTM3");
    else if ( indexInBox == RTM_rFRTM2 )
        _PETRTM.setRTMModelType("rFRTM2");
    _PETRTM.setPrepared(false);
    updateAllGraphs();
}

void SimWindow::createTargetPage()
{
    qDebug() << "SimWindow::createTargetPage enter";
    _targetPage = new QWidget();

    _basisPlot  = new plotData(2);
    _targetPlot = new plotData(3);
    auto *plotLayout = new QVBoxLayout();
    plotLayout->addWidget(_basisPlot->getPlotSurface());
    plotLayout->addWidget(_targetPlot->getPlotSurface());
    QString numberString;
    int editTextSize=80;

    // target region: input
    auto *targetGroupBox = new QGroupBox("Target Region: Input");
    QLabel *BPndLabel    = new QLabel("BPnd");
    QLabel *R1Label      = new QLabel("R1");
    QLabel *tau4Label    = new QLabel("1/k4 (min)");
    QLabel *challengeTimeLabel = new QLabel("Challenge time");
    QLabel *challengeMagLabel  = new QLabel("Delta_BPnd (%)");
    QLabel *noiseLabel   = new QLabel("Noise");
    _BPnd          = new QLineEdit();
    _R1            = new QLineEdit();
    _tau4          = new QLineEdit();
    _challengeTime = new QLineEdit();
    _challengeMag  = new QLineEdit();
    _noiseTar = new QLineEdit();
    _BPnd->setText(numberString.setNum(_simulator.getBP0()));
    _R1->setText(numberString.setNum(_simulator.getR1()));
    _tau4->setText(numberString.setNum(_simulator.getTau4()));
    _challengeTime->setText(numberString.setNum(_simulator.getChallengeTime()));
    _challengeMag->setText(numberString.setNum(_simulator.getChallengeMag()));
    _noiseTar->setText(numberString.setNum(_simulator.getNoiseTar()));
    _BPnd->setFixedWidth(editTextSize);
    _R1->setFixedWidth(editTextSize);
    _tau4->setFixedWidth(editTextSize);
    _challengeTime->setFixedWidth(editTextSize);
    _challengeMag->setFixedWidth(editTextSize);
    _noiseTar->setFixedWidth(editTextSize);
    auto *targetLayout = new QGridLayout();
    targetLayout->addWidget(BPndLabel,0,0);
    targetLayout->addWidget(_BPnd,0,1);
    targetLayout->addWidget(R1Label,1,0);
    targetLayout->addWidget(_R1,1,1);
    targetLayout->addWidget(tau4Label,2,0);
    targetLayout->addWidget(_tau4,2,1);
    targetLayout->addWidget(challengeTimeLabel,3,0);
    targetLayout->addWidget(_challengeTime,3,1);
    targetLayout->addWidget(challengeMagLabel,4,0);
    targetLayout->addWidget(_challengeMag,4,1);
    targetLayout->addWidget(noiseLabel,5,0);
    targetLayout->addWidget(_noiseTar,5,1);
    targetGroupBox->setLayout(targetLayout);
    targetLayout->setSpacing(0);
    connect(_BPnd,          SIGNAL(editingFinished()), this, SLOT(changedBPND()));
    connect(_R1,            SIGNAL(editingFinished()), this, SLOT(changedR1()));
    connect(_tau4,          SIGNAL(editingFinished()), this, SLOT(changedTau4()));
    connect(_challengeTime, SIGNAL(editingFinished()), this, SLOT(changedChallengeTime()));
    connect(_challengeMag,  SIGNAL(editingFinished()), this, SLOT(changedChallengeMag()));
    connect(_noiseTar,      SIGNAL(editingFinished()), this, SLOT(changedNoiseTar()));

    // target region: analysis
    auto *analysisGroupBox = new QGroupBox("Target Region: Analysis");
    QLabel *modelLabel      = new QLabel("RTM type");
    QLabel *weightLabel     = new QLabel("Weighting scheme");
    QLabel *challengeLabel  = new QLabel("Challenge?   ");
    _tau2RefAnalysisLabel   = new QLabel("Analysis Tau2'");
    _tau4AnalysisLabel      = new QLabel("Analysis Tau4");
    _modelType = new QComboBox();
    _modelType->addItem("SRTM3");
    _modelType->addItem("SRTM2");
    _modelType->addItem("rFRTM3");
    _modelType->addItem("rFRTM2");
    _modelType->setCurrentIndex(0);
    _weightType  = new QComboBox();
    _weightType->addItem("Uniform weights");
    _weightType->addItem("11C-Noiseless");
    _weightType->addItem("11C");
    _weightType->addItem("Custom");
    _weightType->setToolTip("Set a weighting scheme for WLS");
    _tau2RefAnalysis = new QLineEdit();
    _tau4Analysis    = new QLineEdit();
    _tau2RefAnalysis->setText(numberString.setNum(_simulator.getTau2Ref()));
    _tau4Analysis->setText(numberString.setNum(_simulator.getTau4()));
    _tau2RefAnalysis->setFixedWidth(editTextSize);
    _tau4Analysis->setFixedWidth(editTextSize);
    _tau2RefAnalysisLabel->setVisible(false);
    _tau4AnalysisLabel->setVisible(false);
    _tau2RefAnalysis->setVisible(false);
    _tau4Analysis->setVisible(false);
    connect(_modelType,  SIGNAL(currentIndexChanged(int)), this, SLOT(changedModelType(int)));
    connect(_weightType, SIGNAL(currentIndexChanged(int)), this, SLOT(changedWeightType(int)));
    _includeChallenge = new QCheckBox();
    _includeChallenge->setChecked(false);

    auto *analysisLayout = new QGridLayout();
    analysisLayout->addWidget(modelLabel,0,0);
    analysisLayout->addWidget(_modelType,0,1);
    analysisLayout->addWidget(weightLabel,1,0);
    analysisLayout->addWidget(_weightType,1,1);
    analysisLayout->addWidget(challengeLabel,2,0);
    analysisLayout->addWidget(_includeChallenge,2,1);
    analysisLayout->addWidget(_tau2RefAnalysisLabel,3,0);
    analysisLayout->addWidget(_tau2RefAnalysis,3,1);
    analysisLayout->addWidget(_tau4AnalysisLabel,4,0);
    analysisLayout->addWidget(_tau4Analysis,4,1);
    analysisGroupBox->setLayout(analysisLayout);
    analysisLayout->setSpacing(0);
    connect(_tau2RefAnalysis,  SIGNAL(editingFinished()), this, SLOT(changedTau2RefAnalysis()));
    connect(_tau4Analysis,     SIGNAL(editingFinished()), this, SLOT(changedTau4Analysis()));
    connect(_includeChallenge, SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxChallenge(bool)));
    // target region: errors
    auto *errorGroupBox = new QGroupBox("Target Region: Percentage Errors");
    QLabel *errorBPndLabel  = new QLabel("BPnd ");
    QLabel *errork2Label    = new QLabel("k2   ");
    QLabel *errork2aLabel   = new QLabel("k2a  ");
    _errorR1Label           = new QLabel("R1   ");
    _errorTau2RefLabel      = new QLabel("Tau2Ref");
    _errorChallengeLabel    = new QLabel("Challenge");
    _errorBPnd  = new QLabel();
    _errork2    = new QLabel();
    _errork2a   = new QLabel();
    _errorR1    = new QLabel();
    _errorTau2Ref   = new QLabel();
    _errorChallenge = new QLabel();
    _errorChallengeLabel->setVisible(false);
    _errorChallenge->setVisible(false);
    auto *errorLayout = new QGridLayout();
    errorLayout->addWidget(errorBPndLabel,0,0);
    errorLayout->addWidget(_errorBPnd,0,1);
    errorLayout->addWidget(errork2Label,1,0);
    errorLayout->addWidget(_errork2,1,1);
    errorLayout->addWidget(errork2aLabel,2,0);
    errorLayout->addWidget(_errork2a,2,1);
    errorLayout->addWidget(_errorR1Label,3,0);
    errorLayout->addWidget(_errorR1,3,1);
    errorLayout->addWidget(_errorTau2RefLabel,4,0);
    errorLayout->addWidget(_errorTau2Ref,4,1);
    errorLayout->addWidget(_errorChallengeLabel,5,0);
    errorLayout->addWidget(_errorChallenge,5,1);
    errorGroupBox->setLayout(errorLayout);
    errorLayout->setSpacing(0);

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(targetGroupBox);
    rightLayout->addWidget(analysisGroupBox);
    rightLayout->addWidget(errorGroupBox);
    rightLayout->setSpacing(0);

    QHBoxLayout *fullLayout = new QHBoxLayout();
    fullLayout->addLayout(plotLayout);
    fullLayout->addLayout(rightLayout);
    fullLayout->setStretch(0,5);
    fullLayout->setStretch(1,1);

    /////////////// Plotting tool bars //////////////////////
    const QIcon *dragX = new QIcon(":/My-Icons/dragX.png");
    auto *dragXAction = new QAction(*dragX, tr("drag/zoom X axis"), this);
    const QIcon *dragY = new QIcon(":/My-Icons/dragY.png");
    auto *dragYAction = new QAction(*dragY, tr("drag/zoom Y axis"), this);
    const QIcon *rescaleXY = new QIcon(":/My-Icons/rescaleGraph.png");
    auto *rescaleXYAction = new QAction(*rescaleXY, tr("Auto-scale X and Y ranges"), this);
    dragXAction->setCheckable(true);
    dragYAction->setCheckable(true);
    rescaleXYAction->setCheckable(true);
    rescaleXYAction->setChecked(true);
    QActionGroup *graphButtons = new QActionGroup(this);
    graphButtons->addAction(dragXAction);
    graphButtons->addAction(dragYAction);
    graphButtons->addAction(rescaleXYAction);

    QToolBar *graphToolBar = new QToolBar("time tool bar");
    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addAction(rescaleXYAction);
    fullLayout->setMenuBar(graphToolBar);

    // Make toolbar connections to the main plot, not the basis plot
    connect(dragXAction,     SIGNAL(triggered(bool)), _basisPlot,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _basisPlot,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _basisPlot,  SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _targetPlot, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _targetPlot, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _targetPlot, SLOT(autoScale(bool)));

    _targetPage->setLayout(fullLayout);
    qDebug() << "SimWindow::createTargetPage exit";
}

void SimWindow::createDependenciesPage()
{
    qDebug() << "SimWindow::createDependenciesPage enter";
    // For RTM3:
    // 1) BPnd_err vs. BPnd
    // 2) Challenge error vs. BPnd
    // 3) Tau2Ref  vs. BPnd

    // For RTM2:
    // 1) BPnd_err vs. BPnd
    // 2) Challenge error vs. BPnd
    // 3) Tau2Ref (no bias) vs. BPnd (requires root finding)

    _dependenciesPage = new QWidget();

    _errBPndPlot  = new plotData(4);
    _errChallPlot = new plotData(5);
    _tau2RefPlot   = new plotData(6);
    auto *plotLayout = new QVBoxLayout();
    plotLayout->addWidget(_errBPndPlot->getPlotSurface());
    plotLayout->addWidget(_errChallPlot->getPlotSurface());
    plotLayout->addWidget(_tau2RefPlot->getPlotSurface());
    QString numberString;
    int editTextSize=80;

    auto *BPndGroupBox    = new QGroupBox("BPnd graph range and step");
    QLabel *BPndLowLabel  = new QLabel("BPnd low");
    QLabel *BPndHighLabel = new QLabel("BPnd high");
    QLabel *BPndStepLabel = new QLabel("BPnd step");
    _BPndLow          = new QLineEdit();
    _BPndHigh         = new QLineEdit();
    _BPndStep         = new QLineEdit();
    _BPndLow->setFixedWidth(editTextSize);
    _BPndHigh->setFixedWidth(editTextSize);
    _BPndStep->setFixedWidth(editTextSize);
    _BPndLow->setText(numberString.setNum(_BPndLowValue));
    _BPndHigh->setText(numberString.setNum(_BPndHighValue));
    _BPndStep->setText(numberString.setNum(_BPndStepValue));
    connect(_BPndLow,   SIGNAL(editingFinished()), this, SLOT(changedBPndLow()));
    connect(_BPndHigh,  SIGNAL(editingFinished()), this, SLOT(changedBPndHigh()));
    connect(_BPndStep,  SIGNAL(editingFinished()), this, SLOT(changedBPndStep()));
    auto *BPndLayout = new QGridLayout();
    BPndLayout->addWidget(BPndLowLabel,0,0);
    BPndLayout->addWidget(_BPndLow,0,1);
    BPndLayout->addWidget(BPndHighLabel,1,0);
    BPndLayout->addWidget(_BPndHigh,1,1);
    BPndLayout->addWidget(BPndStepLabel,2,0);
    BPndLayout->addWidget(_BPndStep,2,1);
    BPndGroupBox->setLayout(BPndLayout);
    BPndLayout->setSpacing(0);

    auto *calcGroupBox = new QGroupBox("Calculate and clear curves");
    QLabel *calcLabel  = new QLabel("Calculate curves");
    QLabel *clearLabel = new QLabel("Clear curves");
    _calculateCurves   = new QPushButton();
    QPixmap pixmapCalculate(":/My-Icons/calculator.png");
    QIcon calculatorIcon(pixmapCalculate);
    _calculateCurves->setIcon(calculatorIcon);
    _clearCurves  = new QPushButton();
    QPixmap eraser(":/My-Icons/eraser.png");
    QIcon eraserIcon(eraser);
    _clearCurves->setIcon(eraserIcon);
    auto *calcLayout = new QGridLayout();
    calcLayout->addWidget(calcLabel,0,0);
    calcLayout->addWidget(_calculateCurves,0,1);
    calcLayout->addWidget(clearLabel,1,0);
    calcLayout->addWidget(_clearCurves,1,1);
    calcGroupBox->setLayout(calcLayout);
    calcLayout->setSpacing(0);
    connect(_calculateCurves, SIGNAL(pressed()), this, SLOT(calculateBPndCurves()));
    connect(_clearCurves, SIGNAL(pressed()),     this, SLOT(clearBPndCurves()));

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(BPndGroupBox);
    rightLayout->addWidget(calcGroupBox);

    QHBoxLayout *fullLayout = new QHBoxLayout();
    fullLayout->addLayout(plotLayout);
    fullLayout->addLayout(rightLayout);
    fullLayout->setStretch(0,100);
    fullLayout->setStretch(1,1);

    /////////////// Plotting tool bars //////////////////////
    const QIcon *dragX = new QIcon(":/My-Icons/dragX.png");
    auto *dragXAction = new QAction(*dragX, tr("drag/zoom X axis"), this);
    const QIcon *dragY = new QIcon(":/My-Icons/dragY.png");
    auto *dragYAction = new QAction(*dragY, tr("drag/zoom Y axis"), this);
    const QIcon *rescaleXY = new QIcon(":/My-Icons/rescaleGraph.png");
    auto *rescaleXYAction = new QAction(*rescaleXY, tr("Auto-scale X and Y ranges"), this);
    dragXAction->setCheckable(true);
    dragYAction->setCheckable(true);
    rescaleXYAction->setCheckable(true);
    rescaleXYAction->setChecked(true);
    QActionGroup *graphButtons = new QActionGroup(this);
    graphButtons->addAction(dragXAction);
    graphButtons->addAction(dragYAction);
    graphButtons->addAction(rescaleXYAction);

    QToolBar *graphToolBar = new QToolBar("time tool bar");
    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addAction(rescaleXYAction);
    fullLayout->setMenuBar(graphToolBar);

    // Make toolbar connections to the main plot, not the basis plot
    connect(dragXAction,     SIGNAL(triggered(bool)), _errBPndPlot,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _errBPndPlot,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _errBPndPlot,  SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _errChallPlot, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _errChallPlot, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _errChallPlot, SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _tau2RefPlot,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _tau2RefPlot,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _tau2RefPlot,  SLOT(autoScale(bool)));

    clearBPndCurves();

    _dependenciesPage->setLayout(fullLayout);
    qDebug() << "SimWindow::createTargetPage exit";
}

void SimWindow::changedGraphSizes(int iSelection)
{
    if ( iSelection == 0 )
    { // show both
        _plasmaPlot->getPlotSurface()->setVisible(true);
        _RRPlot->getPlotSurface()->setVisible(true);
        _setupPlotLayout->setStretch(0,1);     // basis plot
        _setupPlotLayout->setStretch(1,1);     // time plot
    }
    else if ( iSelection == 1 )
    { // data only
        _plasmaPlot->getPlotSurface()->setVisible(true);
        _RRPlot->getPlotSurface()->setVisible(false);
        _setupPlotLayout->setStretch(0,1);     // basis plot
        _setupPlotLayout->setStretch(1,1);     // time plot
    }
    else
    { // basis only
        _plasmaPlot->getPlotSurface()->setVisible(false);
        _RRPlot->getPlotSurface()->setVisible(true);
        _setupPlotLayout->setStretch(0,1);     // basis plot
        _setupPlotLayout->setStretch(1,1);     // time plot
    }
}
void SimWindow::updatePlasmaGraph()
{
    qDebug() << "SimWindow::updatePlasmaGraph enter";
    // run the simulation
    _simulator.generatePlasmaTAC();
    _simulator.generateReferenceTAC();

    // update the plot
    _plasmaPlot->init();
    _plasmaPlot->setLegendOn(true);
    _plasmaPlot->addCurve(0,"plasma");
    _plasmaPlot->setColor(Qt::red);
    double duration = _simulator.getDuration();
    double stepSize = _simulator.getStepSize();
    int nTime = static_cast<int>(duration/stepSize);
    dVector xTime; xTime.resize(nTime);
    dVector yTAC;  yTAC.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt]   = jt * stepSize;
        yTAC[jt] = _simulator.getCp(jt);
    }
    _plasmaPlot->setData(xTime,yTAC);

    // update the plot: RR
    _plasmaPlot->addCurve(0,"RR");
    for (int jt=0; jt<nTime; jt++)
        yTAC[jt] = _simulator.getCr(jt);
    _plasmaPlot->setData(xTime,yTAC);

    _plasmaPlot->conclude(0,true);
    _plasmaPlot->plotDataAndFit(true);
    qDebug() << "SimWindow::updatePlasmaGraph exit";
}
void SimWindow::updateReferenceGraph()
{
    qDebug() << "SimWindow::updateReferenceGraph enter";
    // run the simulation
    _simulator.generateReferenceTAC();

    // update the plot: RR
    _RRPlot->init();
    _RRPlot->setLegendOn(true);
    _RRPlot->addCurve(0,"RR");
    _RRPlot->setColor(Qt::red);
    _RRPlot->setPointSize(5);
    double duration = _simulator.getDuration();
    double stepSize = _simulator.getStepSize();
    int lDownSample = _simulator.getDownSampling();
    int nTime = static_cast<int>(duration/stepSize) / lDownSample;
    dVector xTime;   xTime.resize(nTime);
    dVector yTAC;    yTAC.resize(nTime);
    dMatrix refRegion;  refRegion.resize(1);  refRegion[0].resize(nTime);
    dMatrix timeBins;   timeBins.resize(1);   timeBins[0].resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = jt * lDownSample * stepSize;
        yTAC[jt]  = refRegion[0][jt] = _simulator.getCrDown(jt);
        timeBins[0][jt] = stepSize * lDownSample;
    }
    _RRPlot->setData(xTime,yTAC);
    static int firstTime=true;
    if ( firstTime || nTime != _PETRTM.getNumberTimePoints() )
    {
        // Start with SRTM and no challenge
        _PETRTM.setReferenceRegion(timeBins,refRegion);
        // Assign these IDs:
        _PETRTM.setR1EventID(0,'R');
        _PETRTM.setk2EventID(0,'k');
        _PETRTM.setk2aEventID(0,'a');
        int indexChallenge = _PETRTM.getEventIndex('c');
        _PETRTM.setChallengeShape(indexChallenge,Challenge_Sigmoid);
        _PETRTM.setChallengeTau(indexChallenge,1.);
        _PETRTM.setChallengeOnset(indexChallenge,0,_simulator.getChallengeTime());
        _PETRTM.setTau2RefSRTM(0,_simulator.getTau2Ref());
        _PETRTM.setTau2RefFRTM(0,_simulator.getTau2Ref());
        _PETRTM.setTau4(0,_simulator.getTau4());
        firstTime = false;
    }
    _PETRTM.setReferenceRegion(refRegion);

    if ( _dataRefRegion->count() != 0 )
    {
        int indexInBox = _dataRefRegion->currentIndex();
        int indexData  = indexInBox + 1;  // data index 0 is time
        if ( indexData >= _dataTable.size() )
            qFatal("Fatal Error: the combo-box index exceeds the table index.");
        _RRPlot->addCurve(0,"real data");
        _RRPlot->setData(_dataTable[0],_dataTable[indexData]);
    }

    _RRPlot->conclude(0,true);
    _RRPlot->plotDataAndFit(true);
    qDebug() << "SimWindow::updateReferenceGraph exit";
}
void SimWindow::updateBasisGraph()
{
    _basisPlot->init();
    _basisPlot->setLegendOn(true);

    _basisPlot->addCurve(0,"weights");
    _basisPlot->setPointSize(3);
    double averageY=0;
    dVector xData;  xData.clear();
    dVector yData;  yData.clear();
    for (int jt=0; jt<_PETRTM.getNumberTimePoints(); jt++)
    {
        xData.append(jt);
        yData.append(_PETRTM.getWeight(jt));
        averageY += _PETRTM.getWeight(jt);
    }
    averageY /= static_cast<double>(_PETRTM.getNumberTimePoints());
    _basisPlot->setData(xData, yData);

    int nBasis = _PETRTM.getNumberCoefficients();
    qDebug() << "SimWindow::updateBasisGraph nBasis" << nBasis;
    for ( int jBasis=0; jBasis<nBasis; jBasis++)
    {
        qDebug() << "SimWindow::updateBasisGraph jBasis" << jBasis << _PETRTM.getBasisType(jBasis);
        if ( _PETRTM.getBasisType(jBasis) == Type_R1 )
        {
            _basisPlot->addCurve(0,"R1");
            _basisPlot->setColor(Qt::red);
        }
        else if ( _PETRTM.getBasisType(jBasis) == Type_k2 )
        {
            _basisPlot->addCurve(0,"k2");
            _basisPlot->setColor(Qt::blue);
        }
        else if ( _PETRTM.getBasisType(jBasis) == Type_k2a )
        {
            _basisPlot->addCurve(0,"k2a");
            _basisPlot->setColor(Qt::green);
        }
        else if ( _PETRTM.getBasisType(jBasis) == Type_dCrdt )
        {
            _basisPlot->addCurve(0,"dCr/dt");
            _basisPlot->setColor(Qt::yellow);
        }
        else
        {
            _basisPlot->addCurve(0,"challenge");
            _basisPlot->setColor(Qt::magenta);
        }
        // add data
        yData.clear();
        for (int jt=0; jt<_PETRTM.getNumberTimePoints(); jt++)
            yData.append(_PETRTM.getBasisPoint(jBasis,jt));
        _basisPlot->setData(xData, yData);
        // rescale weights
        if ( _PETRTM.getBasisType(jBasis) == Type_k2 )
        {
            double yFraction = 0.67;
            plotCurve *targetCurve = _basisPlot->getThisCurvePointer();  // this
            plotCurve sourceCurve  = _basisPlot->_listOfCurves[0];       // weights
            double sourceMax = _basisPlot->getMaxY(&sourceCurve) / yFraction;
            double targetMax = _basisPlot->getMaxY(targetCurve);
            if ( targetMax != 0. && sourceMax != 0. )
            {
                _basisPlot->_yAxis2Ratio = sourceMax / targetMax;
                _basisPlot->_listOfCurves[0].scaleFactor = 1. / _basisPlot->_yAxis2Ratio;
                for (int jt=0; jt<_PETRTM.getNumberTimePoints(); jt++)
                    _basisPlot->_listOfCurves[0].yData[jt] /= _basisPlot->_yAxis2Ratio;
            }
        }
    }

    _basisPlot->conclude(0,true);
    _basisPlot->plotDataAndFit(true);

}
void SimWindow::updateTargetGraph()
{
    qDebug() << "SimWindow::updateTargetGraph enter";
    // run the simulation
    _simulator.generateTargetTAC();
    qDebug() << "SimWindow::updateTargetGraph 1";

    // update the plot
    _targetPlot->init();
    _targetPlot->setLegendOn(true);
    _targetPlot->addCurve(0,"target");
    double duration = _simulator.getDuration();
    double stepSize = _simulator.getStepSize();
    int lDownSample = _simulator.getDownSampling();
    int nTime = static_cast<int>(duration/stepSize) / lDownSample;
    dVector xTime;   xTime.resize(nTime);
    dVector yTAC;    yTAC.resize(nTime);
    dMatrix tissueVector;  tissueVector.resize(1);  tissueVector[0].resize(nTime);
    dMatrix fitVector;     fitVector.resize(1);     fitVector[0].resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = jt * lDownSample * stepSize;
        yTAC[jt]  = tissueVector[0][jt] = _simulator.getCtDown(jt);
    }
    _targetPlot->setData(xTime,yTAC);

    // update the RTM model
    _PETRTM.setTissueVector(true,tissueVector);
    _PETRTM.definePETConditions("a c"); // don't define R1, which is not valid for RTM2
    _PETRTM.prepare();
    _PETRTM.fitData(tissueVector,fitVector);
    qDebug() << "SimWindow::updateTargetGraph 3";
    analyzeSimulatedTAC();
    qDebug() << "SimWindow::updateTargetGraph 4";

    // add fit
    _targetPlot->addCurve(0,"fit");
    _targetPlot->setColor(Qt::red);
    for (int jt=0; jt<nTime; jt++)
        yTAC[jt] = _PETRTM.getFit(jt);
    _targetPlot->setData(xTime,yTAC);

    // add RR
    _targetPlot->addCurve(0,"RR");
    _targetPlot->setColor(Qt::gray);
    for (int jt=0; jt<nTime; jt++)
        yTAC[jt]  = _simulator.getCrDown(jt);
    _targetPlot->setData(xTime,yTAC);

    _targetPlot->conclude(0,true);
    _targetPlot->plotDataAndFit(true);
    qDebug() << "SimWindow::updateTargetGraph exit";
}
void SimWindow::updateAllGraphs()
{
    updatePlasmaGraph();
    updateReferenceGraph();
    updateTargetGraph();
    updateBasisGraph();  // update basis graph AFTER target graph, because basis functions use target curve
}
void SimWindow::analyzeSimulatedTAC()
{
    qDebug() << "SimWindow::analyzeSimulatedTAC enter";
    // BPnd
    double truth = _simulator.getBP0();
    double guess = _PETRTM.getBP0InRun(0);
    _errorBPnd->setText(analyzeString(truth,guess));
    // k2
    truth = _simulator.getk2();
    guess = _PETRTM.getk2InRun(0).x;
    _errork2->setText(analyzeString(truth,guess));
    // k2a
    truth = _simulator.getk2a();
    guess = _PETRTM.getk2aInRun(0).x;
    _errork2a->setText(analyzeString(truth,guess));
    // R1
    truth = _simulator.getR1();
    guess = _PETRTM.getR1InRun(0).x;
    _errorR1->setText(analyzeString(truth,guess));
    // tau2Ref
    truth = _simulator.getTau2Ref();
    guess = _PETRTM.getTau2RefInRun(0);
    _errorTau2Ref->setText(analyzeString(truth,guess));
    // challenge
    truth = _simulator.getChallengeMag();  // delta_BPnd abs
    guess = getChallengeMagFromAnalysis();
    QString valueString, diffString;
    valueString.setNum(guess,'g',2);
    if ( guess == 0. )
        diffString = "----";
    else
        diffString.setNum(guess-truth,'g',2);
    _errorChallenge->setText(valueString + " %, diff = " + diffString + " %");
}
double SimWindow::getChallengeMagFromAnalysis()
{
    _PETRTM.updateConditions();
    int nConditions = _PETRTM.getNumberConditions();
    if ( nConditions == 2 )
    {
        _PETRTM.setCurrentCondition(1);
        _PETRTM.evaluateCurrentCondition();
        double percentChange = _PETRTM.getBPndInCurrentCondition().x;
        percentChange *= 100./_PETRTM.getBP0InRun(0);
        return percentChange;
    }
    else
        return 0.;
}
QString SimWindow::analyzeString(double truth, double guess)
{
    QString percentString, valueString;
    percentString.setNum(percentageError(truth,guess),'g',2);
    valueString.setNum(guess,'g',3);
    return percentString + " % " + QString("(abs = %1)").arg(valueString);
}

///////////////////////////////////////
// Slots
///////////////////////////////////////
void SimWindow::changedTimeDuration()
{
    QString stringEntered = _timeDuration->text();
    bool ok;
    double duration = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setDuration(duration);
        double stepSize = _simulator.getStepSize();
        int lDownSample = _simulator.getDownSampling();
        int nTime = static_cast<int>(duration/stepSize) / lDownSample;
        _PETRTM.setTimePointsInRun(0,nTime);
        updateAllGraphs();
    }
    else
        _timeDuration->setText(stringEntered.setNum(_simulator.getDuration()));
}
void SimWindow::changedTimeStep()
{
    QString stringEntered = _timeStep->text();
    bool ok;
    double stepSize = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setStepSize(stepSize);
        double duration = _simulator.getDuration();
        int lDownSample = _simulator.getDownSampling();
        int nTime = static_cast<int>(duration/stepSize) / lDownSample;
        _PETRTM.setTimePointsInRun(0,nTime);
        updateAllGraphs();
    }
    else
        _timeStep->setText(stringEntered.setNum(_simulator.getStepSize()));
}
void SimWindow::changedDownSample()
{
    QString stringEntered = _downSample->text();
    bool ok;
    int lDownSample = stringEntered.toInt(&ok);
    if ( ok )
    {
        _simulator.setDownSampling(lDownSample);
        double duration = _simulator.getDuration();
        double stepSize = _simulator.getStepSize();
        int lDownSample = _simulator.getDownSampling();
        int nTime = static_cast<int>(duration/stepSize) / lDownSample;
        _PETRTM.setTimePointsInRun(0,nTime);
        updateAllGraphs();
    }
    else
        _downSample->setText(stringEntered.setNum(_simulator.getDownSampling()));
}
void SimWindow::changedBolusMag()
{
    QString stringEntered = _bolusMag->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setMagBolus(value);
        updateAllGraphs();
    }
    else
        _bolusMag->setText(stringEntered.setNum(_simulator.getMagBolus()));
}
void SimWindow::changedTauBolus()
{
    QString stringEntered = _tauDecay->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setTauBolus(value);
        updateAllGraphs();
    }
    else
        _tauDecay->setText(stringEntered.setNum(_simulator.getTauBolus()));
}
void SimWindow::changedInfusion()
{
    QString stringEntered = _KBol->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setKBol(value);
        updateAllGraphs();
    }
    else
        _KBol->setText(stringEntered.setNum(_simulator.getKBol()));
}
void SimWindow::changedTau2Ref()
{
    QString stringEntered = _tau2Ref->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setTau2Ref(value);
        updateAllGraphs();
    }
    else
        _tau2Ref->setText(stringEntered.setNum(_simulator.getTau2Ref()));
}
void SimWindow::changedTau1Ref()
{
    QString stringEntered = _tau1Ref->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setTau1Ref(value);
        updateAllGraphs();
    }
    else
        _tau1Ref->setText(stringEntered.setNum(_simulator.getTau1Ref()));
}
void SimWindow::changedNoiseRef()
{
    QString stringEntered = _noiseRef->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setNoiseRef(value);
        updateAllGraphs();
    }
    else
        _noiseRef->setText(stringEntered.setNum(_simulator.getNoiseRef()));
}
void SimWindow::changedFastElimination()
{
    QString stringEntered = _fastTau->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setTauFastElim(value);
        updateAllGraphs();
    }
    else
        _fastTau->setText(stringEntered.setNum(_simulator.getTauFastElim()));
}
void SimWindow::changedSlowElimination()
{
    QString stringEntered = _slowTau->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setTauSlowElim(value);
        updateAllGraphs();
    }
    else
        _slowTau->setText(stringEntered.setNum(_simulator.getTauSlowElim()));
}
void SimWindow::changedFastEliminationFraction()
{
    QString stringEntered = _fastFraction->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setFastFraction(value);
        updateAllGraphs();
    }
    else
        _fastFraction->setText(stringEntered.setNum(_simulator.getFastFraction()));
}
void SimWindow::changedBPND()
{
    QString stringEntered = _BPnd->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setBP0(value);
        updateAllGraphs();
    }
    else
        _BPnd->setText(stringEntered.setNum(_simulator.getBP0()));
}
void SimWindow::changedR1()
{
    QString stringEntered = _R1->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setR1(value);
        updateAllGraphs();
    }
    else
        _R1->setText(stringEntered.setNum(_simulator.getR1()));
}
void SimWindow::changedTau4()
{
    QString stringEntered = _tau4->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        if ( value == 0. )
            _simulator.setSRTM();
        else
        {
            _simulator.setFRTM();
            _simulator.setTau4(value);
        }
        updateAllGraphs();
    }
    else
        _tau4->setText(stringEntered.setNum(_simulator.getTau4()));
}
void SimWindow::changedChallengeTime()
{
    QString stringEntered = _challengeTime->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeTime(value);
        updateAllGraphs();
    }
    else
        _challengeTime->setText(stringEntered.setNum(_simulator.getChallengeTime()));
}
void SimWindow::changedChallengeMag()
{
    QString stringEntered = _challengeMag->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeMag(value);
        updateAllGraphs();
    }
    else
        _challengeMag->setText(stringEntered.setNum(_simulator.getChallengeMag()));
}
void SimWindow::changedNoiseTar()
{
    QString stringEntered = _noiseTar->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setNoiseTar(value);
        updateAllGraphs();
    }
    else
        _noiseTar->setText(stringEntered.setNum(_simulator.getNoiseTar()));
}
void SimWindow::changedTau2RefAnalysis()
{
    QString stringEntered = _tau2RefAnalysis->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        qDebug() << "changed tau2Ref (analysis) to" << value;
        _PETRTM.setTau2Ref(0,value);
        updateAllGraphs();
    }
    else
        _tau2RefAnalysis->setText(stringEntered.setNum(_PETRTM.getTau2RefInRun(0)));
}
void SimWindow::changedTau4Analysis()
{
    QString stringEntered = _tau4Analysis->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _PETRTM.setTau4(0,value);
        updateAllGraphs();
    }
    else
        _tau4Analysis->setText(stringEntered.setNum(_PETRTM.getTau4(0)));
}
void SimWindow::changedCheckBoxChallenge(bool state)
{
    qDebug() << "SimWindow::changedCheckBoxChallenge" << state;
    int indexChallenge = _PETRTM.getEventIndex('c');
    if ( state )
    {
        _PETRTM.setChallengeRun(indexChallenge,0,0);
        _PETRTM.prepare();  // prepare basis functions so that challenge condition does not get rejected
        _PETRTM.definePETConditions("a c"); // don't define R1, which is not valid for RTM2
    }
    else
    {
        _PETRTM.setChallengeRun(indexChallenge,0,-1);
        _PETRTM.prepare();
        _PETRTM.definePETConditions("a"); // don't define R1, which is not valid for RTM2
    }
    _errorChallengeLabel->setVisible(state);
    _errorChallenge->setVisible(state);
    _PETRTM.setPrepared(false);
    updateAllGraphs();
}

QString SimWindow::readTableFile(QString fileName, QStringList &columnNames, dMatrix &table)
{
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QString errorString = "Error attempting to open file " + fileName;
        return errorString;
    }
    QTextStream in_stream(&file);
    QString line = in_stream.readLine();
    QString unCommented = line.left(line.indexOf("#"));
    QRegExp rx("[,\\s]");// match a comma or a space
    columnNames = unCommented.split(rx, QString::SkipEmptyParts);
    int nColumns = columnNames.size();
    table.resize(nColumns);

    int iTime = 0;
    while ( !in_stream.atEnd() )
    {
        QString line = in_stream.readLine();
        QString unCommented = line.left(line.indexOf("#"));
        if ( !unCommented.isEmpty() )
        {
            QStringList valueList = unCommented.split(QRegExp("[,\\s]"), QString::SkipEmptyParts);
            if ( valueList.size() != 0 )
            {
                if ( valueList.size() < nColumns )
                {
                    QString errorString = QString("Error: # columns = %1 but expected # = %2 on line %3.").
                            arg(valueList.size()).arg(nColumns).arg(iTime);
                    return errorString;
                }
                else
                {
                    for ( int jColumn=0; jColumn<nColumns; jColumn++)
                    {
                        QString valueString = valueList.at(jColumn);
                        bool ok;
                        double value = valueString.toDouble(&ok);
                        if ( ok )
                            table[jColumn].append(value);
                    } // jColumn
                } // < nColumns
            } // valueSize != 0
        } // !empty
        iTime++;
    } // new line
    file.close();
    return "";
}

void SimWindow::changedBPndLow()
{
    QString stringEntered = _BPndLow->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
        _BPndLowValue = value;
    else
        _BPndLow->setText(stringEntered.setNum(_BPndLowValue));
}
void SimWindow::changedBPndHigh()
{
    QString stringEntered = _BPndHigh->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
        _BPndHighValue = value;
    else
        _BPndHigh->setText(stringEntered.setNum(_BPndHighValue));
}
void SimWindow::changedBPndStep()
{
    QString stringEntered = _BPndStep->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
        _BPndStepValue = value;
    else
        _BPndStep->setText(stringEntered.setNum(_BPndStepValue));
}

void SimWindow::clearBPndCurves()
{
    QCPRange xRange;
    xRange.lower = _BPndLowValue  - _BPndStepValue;
    xRange.upper = _BPndHighValue + _BPndStepValue;

    _errBPndPlot->init();
//    _errBPndPlot->setLegendOn(true);
    _errBPndPlot->setLabelXAxis("BPnd");
    _errBPndPlot->setLabelYAxis("BPnd % err");
    _errBPndPlot->plotDataAndFit(true);
    _errBPndPlot->setXRange(xRange);

    _errChallPlot->init();
//    _errChallPlot->setLegendOn(true);
    _errChallPlot->setLabelXAxis("BPnd");
    _errChallPlot->setLabelYAxis("dBPnd abs err");
    _errChallPlot->plotDataAndFit(true);
    _errChallPlot->setXRange(xRange);

    _tau2RefPlot->init();
//    _tau2RefPlot->setLegendOn(true);
    _tau2RefPlot->setLabelXAxis("BPnd");
    _tau2RefPlot->setLabelYAxis("Tau2'");
    _tau2RefPlot->plotDataAndFit(true);
    _tau2RefPlot->setXRange(xRange);
}

void SimWindow::calculateBPndCurves()
{
    QColor colors[10] = {Qt::black, Qt::red, Qt::blue, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};
    int nCurves = _errBPndPlot->getNumberCurves();
    int iColor  = nCurves%10;

    _errBPndPlot->addCurve(0,"error_BPnd");
    _errBPndPlot->setPointSize(5);
    _errBPndPlot->setColor(colors[iColor]);

    _errChallPlot->addCurve(0,"error_Challenge");
    _errChallPlot->setPointSize(5);
    _errChallPlot->setColor(colors[iColor]);

    _tau2RefPlot->addCurve(0,"error_Tau2Ref");
    _tau2RefPlot->setPointSize(5);
    _tau2RefPlot->setColor(colors[iColor]);

    dVector xVector;
    for (double BP0=_BPndLowValue; BP0<_BPndHighValue; BP0 += _BPndStepValue)
        xVector.append(BP0);
    dVector errBPnd, errChall, tau2Ref, errBPndRTM2, errChallRTM2;
    for (double BP0=_BPndLowValue; BP0<_BPndHighValue; BP0 += _BPndStepValue)
    {
        QString numberString;  numberString.setNum(BP0);
        _BPnd->setText(numberString);
        changedBPND();

        double truth = _simulator.getBP0();
        double guess = _PETRTM.getBP0InRun(0);
        errBPnd.append(percentageError(BP0,guess));
        qDebug() << "BPnd err = " << percentageError(BP0,guess) << guess << BP0;

        truth = _simulator.getChallengeMag();
        guess = getChallengeMagFromAnalysis();
        errChall.append(guess - truth);
        qDebug() << "challenge err = " << guess - truth << guess << truth;

        if ( !_PETRTM.isRTM2() )
        {
            guess = _PETRTM.getTau2RefInRun(0);
            tau2Ref.append(guess);
            qDebug() << "tau2Ref = " << guess;
        }
        else
        {
            double bestTau2Ref = bestTau2RefForRTM2();
            tau2Ref.append(bestTau2Ref);
            numberString.setNum(bestTau2Ref);  _tau2RefAnalysis->setText(numberString);  changedTau2RefAnalysis();
            // Update error vectors for optimized RTM2
            truth = _simulator.getBP0();
            guess = _PETRTM.getBP0InRun(0);
            errBPndRTM2.append(percentageError(BP0,guess));
            truth = _simulator.getChallengeMag();
            guess = getChallengeMagFromAnalysis();
            errChallRTM2.append(guess - truth);
        }
    }
    _errBPndPlot->setData(xVector,errBPnd);
    _errChallPlot->setData(xVector,errChall);
    _tau2RefPlot->setData(xVector,tau2Ref);

    if ( _PETRTM.isRTM2() )
    { // add a new curves to the two error plots
        _errBPndPlot->addCurve(0,"error_BPnd");
        _errBPndPlot->setPointSize(1);
        _errBPndPlot->setColor(colors[iColor]);
        _errBPndPlot->setData(xVector,errBPndRTM2);
        _errChallPlot->addCurve(0,"error_Challenge");
        _errChallPlot->setPointSize(1);
        _errChallPlot->setColor(colors[iColor]);
        _errChallPlot->setData(xVector,errChallRTM2);
    }

    _errBPndPlot->conclude(0,true);
    _errChallPlot->conclude(0,true);
    _tau2RefPlot->conclude(0,true);

    _errBPndPlot->plotDataAndFit(true);
    _errChallPlot->plotDataAndFit(true);
    _tau2RefPlot->plotDataAndFit(true);
}

double SimWindow::bestTau2RefForRTM2()
{ // search for root
    double tau2Ref = _PETRTM.getTau2RefInRun(0);
    double trueBP0 = _simulator.getBP0();
    double estBP0  = _PETRTM.getBP0InRun(0);
    double error   = estBP0 - trueBP0;
    QString numberString;
    if ( error > 0. )
    {
        while ( error > 0. )
        {
            tau2Ref -= 0.1;
            numberString.setNum(tau2Ref);  _tau2RefAnalysis->setText(numberString);  changedTau2RefAnalysis();
            estBP0 = _PETRTM.getBP0InRun(0);
            error = estBP0 - trueBP0;
        }
    }
    else
    {
        while ( error < 0. )
        {
            tau2Ref += 0.1;
            numberString.setNum(tau2Ref);  _tau2RefAnalysis->setText(numberString);  changedTau2RefAnalysis();
            estBP0 = _PETRTM.getBP0InRun(0);
            error = estBP0 - trueBP0;
        }
    }
    return tau2Ref;
}
