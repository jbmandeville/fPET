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
    _tabTimeSpace->addTab(_setupPage, tr("setup"));
    _tabTimeSpace->addTab(_targetPage, tr("target"));

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
            _dataRefRegion->setCurrentIndex(0);
        }
    }
}

void SimWindow::createSetupPage()
{
    _setupPage = new QWidget();
    QHBoxLayout *fullLayout = new QHBoxLayout();

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
    QLabel *alphaDecayLabel    = new QLabel("Bolus shape (alpha)");
    QLabel *infusionLabel      = new QLabel("Constant Infusion mag");
    QLabel *fastTauLabel       = new QLabel("Fast elimination (tau)");
    QLabel *slowTauLabel       = new QLabel("Slow elimination (tau)");
    QLabel *fastFractionLabel  = new QLabel("Fast fraction (<=1.)");
    _bolusMag     = new QLineEdit();
    _tauDecay     = new QLineEdit();
    _alphaDecay   = new QLineEdit();
    _infusion     = new QLineEdit();
    _fastTau      = new QLineEdit();
    _slowTau      = new QLineEdit();
    _fastFraction = new QLineEdit();
    _bolusMag->setText(numberString.setNum(_simulator.getMagBolus()));
    _tauDecay->setText(numberString.setNum(_simulator.getTauBolus()));
    _alphaDecay->setText(numberString.setNum(_simulator.getAlphaBolus()));
    _infusion->setText(numberString.setNum(_simulator.getMagInfusion()));
    _fastTau->setText(numberString.setNum(_simulator.getKFastElim()));
    _slowTau->setText(numberString.setNum(_simulator.getKSlowElim()));
    _fastFraction->setText(numberString.setNum(_simulator.getFastFraction()));
    _bolusMag->setFixedWidth(editTextSize);
    _tauDecay->setFixedWidth(editTextSize);
    _alphaDecay->setFixedWidth(editTextSize);
    _infusion->setFixedWidth(editTextSize);
    _fastTau->setFixedWidth(editTextSize);
    _slowTau->setFixedWidth(editTextSize);
    _fastFraction->setFixedWidth(editTextSize);
    auto *plasmaInLayout = new QGridLayout();
    plasmaInLayout->addWidget(bolusMagLabel,0,0);
    plasmaInLayout->addWidget(_bolusMag,0,1);
    plasmaInLayout->addWidget(tauDecayLabel,1,0);
    plasmaInLayout->addWidget(_tauDecay,1,1);
    plasmaInLayout->addWidget(alphaDecayLabel,2,0);
    plasmaInLayout->addWidget(_alphaDecay,2,1);
    plasmaInLayout->addWidget(infusionLabel,3,0);
    plasmaInLayout->addWidget(_infusion,3,1);
    plasmaInLayout->addWidget(fastTauLabel,4,0);
    plasmaInLayout->addWidget(_fastTau,4,1);
    plasmaInLayout->addWidget(slowTauLabel,5,0);
    plasmaInLayout->addWidget(_slowTau,5,1);
    plasmaInLayout->addWidget(fastFractionLabel,6,0);
    plasmaInLayout->addWidget(_fastFraction,6,1);
    plasmaInGroupBox->setLayout(plasmaInLayout);
    plasmaInLayout->setSpacing(0);
    connect(_bolusMag, SIGNAL(editingFinished()), this, SLOT(changedBolusMag()));
    connect(_tauDecay, SIGNAL(editingFinished()), this, SLOT(changedTauBolus()));
    connect(_alphaDecay, SIGNAL(editingFinished()), this, SLOT(changedAlphaBolus()));
    connect(_infusion, SIGNAL(editingFinished()), this, SLOT(changedInfusion()));
    connect(_fastTau, SIGNAL(editingFinished()), this, SLOT(changedFastElimination()));
    connect(_slowTau, SIGNAL(editingFinished()), this, SLOT(changedSlowElimination()));
    connect(_fastFraction, SIGNAL(editingFinished()), this, SLOT(changedFastEliminationFraction()));

    // reference region
    auto *RRGroupBox        = new QGroupBox("Reference Region");
    QLabel *tau2RefLabel    = new QLabel("1/k2 (min)");
    QLabel *tau1RefLabel    = new QLabel("1/K1 (min)");
    QLabel *noiseLabel      = new QLabel("Noise");
    QLabel *dataLabel       = new QLabel("Data matching");
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

void SimWindow::changedDataRefRegion(int indexInBox)
{
    updateReferenceGraph();
}

void SimWindow::createTargetPage()
{
    _targetPage = new QWidget();
    QHBoxLayout *fullLayout = new QHBoxLayout();

    _targetPlot = new plotData(0);
    auto *plotLayout = new QVBoxLayout();
    plotLayout->addWidget(_targetPlot->getPlotSurface());
    QString numberString;
    int editTextSize=80;

    // target region
    auto *targetGroupBox = new QGroupBox("Target Region");
    QLabel *BPndLabel    = new QLabel("BPnd");
    QLabel *tau4Label    = new QLabel("1/k4 (min)");
    QLabel *challengeTimeLabel = new QLabel("Challenge time");
    QLabel *challengeMagLabel  = new QLabel("Delta_BPnd (%)");
    QLabel *noiseLabel   = new QLabel("Noise");
    _BPnd          = new QLineEdit();
    _tau4          = new QLineEdit();
    _challengeTime = new QLineEdit();
    _challengeMag  = new QLineEdit();
    _noiseTar = new QLineEdit();
    _BPnd->setText(numberString.setNum(_simulator.getBP0()));
    _tau4->setText(numberString.setNum(_simulator.getTau4()));
    _challengeTime->setText(numberString.setNum(_simulator.getChallengeTime()));
    _challengeMag->setText(numberString.setNum(_simulator.getChallengeMag()));
    _noiseTar->setText(numberString.setNum(_simulator.getNoiseTar()));
    _BPnd->setFixedWidth(editTextSize);
    _tau4->setFixedWidth(editTextSize);
    _challengeTime->setFixedWidth(editTextSize);
    _challengeMag->setFixedWidth(editTextSize);
    _noiseTar->setFixedWidth(editTextSize);
    auto *targetLayout = new QGridLayout();
    targetLayout->addWidget(BPndLabel,0,0);
    targetLayout->addWidget(_BPnd,0,1);
    targetLayout->addWidget(challengeTimeLabel,1,0);
    targetLayout->addWidget(_challengeTime,1,1);
    targetLayout->addWidget(challengeMagLabel,2,0);
    targetLayout->addWidget(_challengeMag,2,1);
    targetLayout->addWidget(tau4Label,3,0);
    targetLayout->addWidget(_tau4,3,1);
    targetLayout->addWidget(noiseLabel,4,0);
    targetLayout->addWidget(_noiseTar,4,1);
    targetGroupBox->setLayout(targetLayout);
    targetLayout->setSpacing(0);

    connect(_BPnd,          SIGNAL(editingFinished()), this, SLOT(changedBPND()));
    connect(_tau4,          SIGNAL(editingFinished()), this, SLOT(changedTau4()));
    connect(_challengeTime, SIGNAL(editingFinished()), this, SLOT(changedChallengeTime()));
    connect(_challengeMag,  SIGNAL(editingFinished()), this, SLOT(changedChallengeMag()));
    connect(_noiseTar,      SIGNAL(editingFinished()), this, SLOT(changedNoiseTar()));

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(targetGroupBox);
    rightLayout->setSpacing(0);

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
    connect(dragXAction,     SIGNAL(triggered(bool)), _targetPlot, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _targetPlot, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _targetPlot, SLOT(autoScale(bool)));

    _targetPage->setLayout(fullLayout);
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
}
void SimWindow::updateReferenceGraph()
{
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
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = jt * lDownSample * stepSize;
        yTAC[jt]  = _simulator.getCrDown(jt);
    }
    _RRPlot->setData(xTime,yTAC);

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
}
void SimWindow::updateTargetGraph()
{
    // run the simulation
    _simulator.generateTargetTAC();

    // update the plot
    _targetPlot->init();
    _targetPlot->setLegendOn(true);
    _targetPlot->addCurve(0,"target");
    _targetPlot->setColor(Qt::red);
    double duration = _simulator.getDuration();
    double stepSize = _simulator.getStepSize();
    int lDownSample = _simulator.getDownSampling();
    int nTime = static_cast<int>(duration/stepSize) / lDownSample;
    dVector xTime;   xTime.resize(nTime);
    dVector yTAC;    yTAC.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = jt * lDownSample * stepSize;
        yTAC[jt]  = _simulator.getCtDown(jt);
    }
    _targetPlot->setData(xTime,yTAC);

    // add RR
    _targetPlot->addCurve(0,"RR");
    for (int jt=0; jt<nTime; jt++)
        yTAC[jt]  = _simulator.getCrDown(jt);
    _targetPlot->setData(xTime,yTAC);

    _targetPlot->conclude(0,true);
    _targetPlot->plotDataAndFit(true);
}
void SimWindow::updateAllGraphs()
{
    updatePlasmaGraph();
    updateReferenceGraph();
    updateTargetGraph();
}

///////////////////////////////////////
// Slots
///////////////////////////////////////
void SimWindow::changedTimeDuration()
{
    QString stringEntered = _timeDuration->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setDuration(value);
        updateAllGraphs();
    }
    else
        _timeDuration->setText(stringEntered.setNum(_simulator.getDuration()));
}
void SimWindow::changedTimeStep()
{
    QString stringEntered = _timeStep->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setStepSize(value);
        updateAllGraphs();
    }
    else
        _timeStep->setText(stringEntered.setNum(_simulator.getStepSize()));
}
void SimWindow::changedDownSample()
{
    QString stringEntered = _downSample->text();
    bool ok;
    int value = stringEntered.toInt(&ok);
    if ( ok )
    {
        _simulator.setDownSampling(value);
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
void SimWindow::changedAlphaBolus()
{
    QString stringEntered = _alphaDecay->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setAlphaBolus(value);
        updateAllGraphs();
    }
    else
        _alphaDecay->setText(stringEntered.setNum(_simulator.getAlphaBolus()));
}
void SimWindow::changedInfusion()
{
    QString stringEntered = _infusion->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setMagInfusion(value);
        updateAllGraphs();
    }
    else
        _infusion->setText(stringEntered.setNum(_simulator.getMagInfusion()));
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
void SimWindow::changedTau4()
{
    QString stringEntered = _tau4->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setTau4(value);
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
