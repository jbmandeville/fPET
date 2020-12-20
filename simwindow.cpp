#include <QtWidgets>
#include <QDebug>
#include <QVector>
#include <QObject>
#include <QFileDialog>
#include <QString>
#include "simwindow.h"
#include "liedetector.h"

SimWindow::SimWindow()
{
    _threadsComboBox = new QComboBox();
    int maxThreads = QThread::idealThreadCount();
    for (int jThread=0; jThread<maxThreads; jThread++)
        _threadsComboBox->addItem(QString("%1 threads").arg(jThread+1));
    _nThreads = maxThreads - 2;
    _threadsComboBox->setCurrentIndex(_nThreads-1);
    _threadsComboBox->setVisible(false);
    connect(_threadsComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(changedNumberThreads(int)));

    _tabTimeSpace = new QTabWidget();
    createSetupPage();
    createTargetPage();
    createSweepBPndPage();
    createVersusTimePage();
    _tabTimeSpace->addTab(_setupPage, tr("Ref Region"));
    _tabTimeSpace->addTab(_targetPage, tr("Target Region"));
    _tabTimeSpace->addTab(_sweepBPndPage, tr("Sweep BPnd"));
    _tabTimeSpace->addTab(_sweepTimePage, tr("Sweep time"));

    QWidget *centralWidget = new QWidget(this);
    this->setCentralWidget( centralWidget );
    auto *mainLayout = new QVBoxLayout( centralWidget );
    mainLayout->addWidget(_tabTimeSpace);

    _statusBar = this->statusBar();
    _statusBar->setStyleSheet("color:Darkred");

    _progressBar = new QProgressBar(_statusBar);
    _progressBar->setAlignment(Qt::AlignRight);
    _progressBar->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    _progressBar->setRange(0,1);

    QHBoxLayout *statusBarLayout = new QHBoxLayout();
    // A page-independent help region or status bar can be added here if desired
    statusBarLayout->addWidget(_threadsComboBox);
    statusBarLayout->addWidget(_statusBar);
    statusBarLayout->addWidget(_progressBar);
    mainLayout->addLayout(statusBarLayout);

    // add a menu
    QMenuBar *menuBar = new QMenuBar;
    mainLayout->setMenuBar(menuBar);
    QMenu *openMenu  = new QMenu(tr("&Menu"), this);
    QMenu *helpMenu  = new QMenu(tr("Help"), this);
    menuBar->addMenu(openMenu);
    menuBar->addMenu(helpMenu);

    QAction *aboutAppAct = helpMenu->addAction(tr("About Simulator"));
    QAction *aboutROIAct = helpMenu->addAction(tr("ROI file format"));
    connect(aboutAppAct,     &QAction::triggered, this, &SimWindow::aboutApp);
    connect(aboutROIAct,     &QAction::triggered, this, &SimWindow::aboutROI);

    QAction *openDataFileAction = openMenu->addAction(tr("Open table file (data)"));
    openDataFileAction->setShortcut(QKeySequence::Open);
    connect(openDataFileAction,  &QAction::triggered, this, &SimWindow::getTableDataFile);

    QAction *quitAction = openMenu->addAction(tr("&Quit"));
    // short-cuts and tooltips
    quitAction->setShortcut(Qt::ControlModifier + Qt::Key_Q);
    connect(quitAction, &QAction::triggered, this, &SimWindow::exitApp);

    QSize defaultWindowSize;
    QRect rec = QApplication::desktop()->screenGeometry();
//    QRect rec = QApplication::desktop()->screenGeometry();
    defaultWindowSize.setWidth(rec.width()*2/3);
    defaultWindowSize.setHeight(rec.height()*2/3);
    resize(defaultWindowSize);

    setThreadVisibility(false);
    updateAllGraphs();
}

void SimWindow::aboutApp()
{
    QMessageBox msgBox;
    QString version = qVersion();
    QString text = "Qt Version " + version;
    msgBox.setText(text);
    msgBox.setIcon(QMessageBox::Information);
    msgBox.exec();

    QMessageBox::information(nullptr, QGuiApplication::applicationDisplayName(),
                             QGuiApplication::applicationDisplayName() + ' '
                             + QCoreApplication::applicationVersion() + " , by Joe Mandeville;\n" +
                             "Request bug fixes by email to\njbm@nmr.mgh.harvard.edu\nwith subject line 'simulator bug'.");
}
void SimWindow::aboutROI()
{
    QMessageBox msgBox;
    QString ROIFile = "ROI file format (use extension .dat, .roi, .list, or .table)\n\n";
    ROIFile += "label1 label2 label3 ...\n";
    ROIFile += "value1(1) value2(1) value3(1)...\n";
    ROIFile += "value1(2) value2(2) value3(2)...\n";
    ROIFile += "value1(3) value2(3) value3(3)...\n";
    ROIFile += "...\n\n";
    QString binList  = _validBinSizeName.join(", ");
    QString timeList = _validBinTimeName.join(", ");
    QString text = ROIFile + "Valid labels for frame durations (units=sec):\n" + binList;
    text += "\n\nValid labels for frame time points (units=min):\n" + timeList;
    msgBox.setText(text);
//    msgBox.setIcon(QMessageBox::Information);
    msgBox.exec();
}
void SimWindow::exitApp()
{
    exit(0);
}

void SimWindow::setThreadVisibility(bool state)
{
    _threadsComboBox->setVisible(state);
    _nSamplesBPndPerThreadLabel->setVisible(state);
    _nSamplesBPndPerThread->setVisible(state);
    _nSamplesBPndLabel->setVisible(state);
    _nSamplesBPnd->setVisible(state);
}

void SimWindow::getTableDataFile()
{
    QString fileName;
    QFileDialog fileDialog;

    fileName = fileDialog.getOpenFileName(this,
                                          "Select an overlay file to open",
                                          QDir::currentPath(),
                                          "Overlay list files (*.dat *.roi *.list *.table)");
    if ( fileName.isEmpty() )
        return;
    else
    {
        QString errorString = utilIO::readTimeTableFile(fileName, _dataColumnNames, _dataTable);
        if ( !errorString.isEmpty() )
        {
            QMessageBox msgBox;
            msgBox.setText(errorString);
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.exec();
        }
        else
        {
            int iColumnBinSize=-1;  int iColumnBinTime=-1;
            if ( !defineTimeBinsFromBinSize(_validBinSizeName, iColumnBinSize) )
            {
                if ( !defineTimeBinsFromTimePoints(_validBinTimeName, iColumnBinTime) )
                {
                    QString binList = _validBinSizeName.join(", ");
                    QString binTime = _validBinTimeName.join(", ");
                    QMessageBox msgBox;
                    errorString = QString("The table must include either frame durations or time points.\n\n");
                    errorString += QString("Frame durations should use seconds and one of these column labels:\n%1\n\n").arg(binList);
                    errorString += QString("Time points should use minutes and one of these column labels:\n%1\n").arg(binTime);
                    msgBox.setText(errorString);
                    msgBox.setIcon(QMessageBox::Critical);
                    msgBox.exec();
                    return;
                }
            }

            // update the GUI to refect the new information
            disconnect(_dataRefRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(changedDataReferenceRegion()));
            disconnect(_dataTargetRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(changedDataTargetRegion()));
            _dataRefRegion->clear();
            _dataTargetRegion->clear();
            for (int jColumn=0; jColumn<_dataColumnNames.count(); jColumn++)
            {
                _dataRefRegion->addItem(_dataColumnNames.at(jColumn));
                _dataTargetRegion->addItem(_dataColumnNames.at(jColumn));
            }
            _dataRefRegion->setCurrentIndex(_dataRefRegion->count()-1);
            // set the default target region to the first region after the bin info
            int iColumnFirstRegion = qMin(_dataColumnNames.count()-1,qMax(iColumnBinSize,iColumnBinTime)+1);
            _dataTargetRegion->setCurrentIndex(iColumnFirstRegion);
            connect(_dataRefRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(changedDataReferenceRegion()));
            connect(_dataTargetRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(changedDataTargetRegion()));

            enablePlasmaMatching(true);
            _ROIFileName->setToolTip(fileName);
            _ROIFileName->setText(utilString::getFileNameWithoutDirectory(fileName));
            // Resample the simulation to match the binning of the real data
            _simulator.setDurationBins(_dtBinsSec);
            QString text;
            _numberTimeBins->setText(text.setNum(_dtBinsSec.size())); changedNumberBins();  // this will update _PETRTM to match _simulator
            // update the bin duration
            int iBin = _binIndex->value() - 1;
            _binDuration->setText(text.setNum(_simulator.getDurationPerBinSec(iBin)));  // this will update graphs

            // enable fitting of real data
            enableComboBoxItem(_simulationStartingPoint,simStart_fromDataFit,true);
            _simulationStartingPoint->setCurrentIndex(simStart_fromDataFit);
        }
    }
}

bool SimWindow::defineTimeBinsFromBinSize(QStringList validBinName, int &iColumn)
{ // read time bin sizes
    // Find bin sizes
    iColumn = -1;
    for (int jName=0; jName<validBinName.count(); jName++)
    {
        iColumn = _dataColumnNames.indexOf(QRegExp(validBinName[jName], Qt::CaseInsensitive));
        if ( iColumn >= 0 ) break;
    }
    if ( iColumn < 0 )
        return false;
    // Convert selected column to time bins
    _dtBinsSec.clear();  _timeBins.clear();
    double dt = _dataTable[iColumn][0]; // bin width in seconds
    double time = dt/60./2.;  // center of 1st bin in min
    for ( int jt=0; jt<_dataTable[iColumn].size(); jt++)
    {
        dt = _dataTable[iColumn][jt]; // input data is seconds
        _dtBinsSec.append(static_cast<int>(dt));
        _timeBins.append(time);
        time += dt * 60.;
    }
    return true;
}

bool SimWindow::defineTimeBinsFromTimePoints(QStringList validBinName, int &iColumn)
{ // convert time (MINUTES) to dt (SECONDS)
    // Find bin sizes
    iColumn = -1;
    for (int jName=0; jName<validBinName.count(); jName++)
    {
        iColumn = _dataColumnNames.indexOf(QRegExp(validBinName[jName], Qt::CaseInsensitive));
        if ( iColumn >= 0 ) break;
    }
    if ( iColumn < 0 )
        return false;
    // Convert selected column to time bins
    _dtBinsSec.clear();  _timeBins.clear();
    for ( int jt=0; jt<_dataTable[iColumn].size(); jt++)
    {
        double time = _dataTable[iColumn][jt];  // input data is time
        double dt;
        if ( jt == 0 )
            dt = 2.*time*60.;
        else
        {
            double timeSec = time * 60.;
            double lowerEdgeSec = _timeBins.last()*60. + static_cast<double>(_dtBinsSec.last())/2.;
            dt = 2.*(timeSec-lowerEdgeSec);
        }
        _dtBinsSec.append(qRound(dt));
        _timeBins.append(time);
    }
    return true;
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

    //////// The setup page status bar
    _plasmaStatusBar = new QStatusBar;  // must be global so it doesn't go out of scope
    _plasmaStatusBar->setStyleSheet("color:blue");
    _plasmaPlot->setQCStatusBar(_plasmaStatusBar);
    _RRStatusBar = new QStatusBar;  // must be global so it doesn't go out of scope
    _RRStatusBar->setStyleSheet("color:blue");
    _RRPlot->setQCStatusBar(_RRStatusBar);
//    _RRPlot->setMainPage(this);

    _setupPlotLayout = new QVBoxLayout();
    _setupPlotLayout->addWidget(_plasmaPlot->getPlotSurface());
    _setupPlotLayout->addWidget(_plasmaStatusBar);
    _setupPlotLayout->addWidget(_RRPlot->getPlotSurface());
    _setupPlotLayout->addWidget(_RRStatusBar);
    _setupPlotLayout->setStretch(0,10);
    _setupPlotLayout->setStretch(1,1);
    _setupPlotLayout->setStretch(2,10);
    _setupPlotLayout->setStretch(3,1);
    QString numberString;
    int editTextSize=80;

    // setup (duration, step, down-sampling
    _binIndex = new QSpinBox();
    auto *setupGroupBox = new QGroupBox("Setup the simulation");
    QLabel *numberTimeBinsLabel = new QLabel("# time bins");
    QLabel *binIndexLabel = new QLabel("bin index");
    QLabel *binDurationLabel = new QLabel("Bin duration (sec)");
    QLabel *subSampleLabel   = new QLabel("Downsampling (int)");
    _numberTimeBins = new QLineEdit();
    _binDuration  = new QLineEdit();
    _subSample    = new QLineEdit();
    _binDuration->setFixedWidth(editTextSize);
    _subSample->setFixedWidth(editTextSize);
    QString text;
    _numberTimeBins->setText(text.setNum(_simulator.getNumberBins()));
    _binDuration->setText(text.setNum(_simulator.getDurationPerBinSec(0)));
    _subSample->setText(numberString.setNum(_simulator.getSamplesPerBin(0)));
    _binIndex->setRange(1,_simulator.getNumberBins());
    _applyToAllBinDuration = new QCheckBox("all");
    _applyToAllBinDuration->setToolTip("apply this bin duration to all bins");

    auto *setupLayout = new QGridLayout();
    setupLayout->addWidget(numberTimeBinsLabel,0,0);
    setupLayout->addWidget(_numberTimeBins,0,1);
    setupLayout->addWidget(binIndexLabel,1,0);
    setupLayout->addWidget(_binIndex,1,1);
    setupLayout->addWidget(binDurationLabel,2,0);
    setupLayout->addWidget(_binDuration,2,1);
    setupLayout->addWidget(_applyToAllBinDuration,2,2);
    setupLayout->addWidget(subSampleLabel,3,0);
    setupLayout->addWidget(_subSample,3,1);
    setupGroupBox->setLayout(setupLayout);
    setupLayout->setSpacing(0);
    connect(_numberTimeBins, SIGNAL(editingFinished()), this, SLOT(changedNumberBins()));
    connect(_binDuration, SIGNAL(editingFinished()), this, SLOT(changedBinDuration()));
    connect(_subSample,   SIGNAL(editingFinished()), this, SLOT(changedSubSample()));
    connect(_binIndex, SIGNAL(valueChanged(int)), this, SLOT(changedBinIndex(int)));
    connect(_applyToAllBinDuration, SIGNAL(toggled(bool)), this, SLOT(changedApplyToAllBinDuration(bool)));

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
    auto *RRGroupBox        = new QGroupBox("Reference Region: truth");
    QLabel *tau2RefLabel    = new QLabel("1/k2' (min)");
    QLabel *tau1RefLabel    = new QLabel("1/K1' (ml/cm^3/min)^-1");
    QLabel *plasmaFracLabel = new QLabel("Plasma %");
    QLabel *noiseLabel      = new QLabel("Noise");
    _tau2Ref = new QLineEdit();
    _tau1Ref = new QLineEdit();
    _plasmaFracRef = new QLineEdit();
    _noiseRef   = new QLineEdit();
    _tau2Ref->setText(numberString.setNum(_simulator.getTau2Ref()));
    _tau1Ref->setText(numberString.setNum(_simulator.getTau1Ref()));
    _plasmaFracRef->setText(numberString.setNum(_simulator.getPlasmaPercentRef()));
    _noiseRef->setText(numberString.setNum(_simulator.getNoiseTar()));
    _tau2Ref->setFixedWidth(editTextSize);
    _tau1Ref->setFixedWidth(editTextSize);
    _plasmaFracRef->setFixedWidth(editTextSize);
    _noiseRef->setFixedWidth(editTextSize);
    auto *RRLayout = new QGridLayout();
    RRLayout->addWidget(tau2RefLabel,0,0);
    RRLayout->addWidget(_tau2Ref,0,1);
    RRLayout->addWidget(tau1RefLabel,1,0);
    RRLayout->addWidget(_tau1Ref,1,1);
    RRLayout->addWidget(plasmaFracLabel,2,0);
    RRLayout->addWidget(_plasmaFracRef,2,1);
    RRLayout->addWidget(noiseLabel,3,0);
    RRLayout->addWidget(_noiseRef,3,1);
    RRGroupBox->setLayout(RRLayout);
    RRLayout->setSpacing(0);

    QPixmap pixmapOpenAlignment(":/My-Icons/openFile.png");
    QIcon openIcon(pixmapOpenAlignment);
    _readROIFile = new QPushButton(openIcon,"open",_setupPage);
    _ROIFileName = new QLabel("No real data",_setupPage);
    auto *realDataGroupBox = new QGroupBox("Real data (imported ROIs from table)");
    QLabel *dataLabel = new QLabel("Ref Region ROI");
    _dataRefRegion = new QComboBox();
    QLabel *startingPointLabel = new QLabel("Simulation starting point");
    _simulationStartingPoint = new QComboBox();
    _simulationStartingPoint->addItem("Simulation from plasma");
    _simulationStartingPoint->addItem("Fit to simulated RR");
    _simulationStartingPoint->addItem("Fit to real data");
    enableComboBoxItem(_simulationStartingPoint,simStart_fromDataFit,false);
    connect(_simulationStartingPoint, SIGNAL(currentIndexChanged(int)), this, SLOT(changedSimulationStartingPoint()));

    auto *realDataLayout = new QGridLayout();
    realDataLayout->addWidget(_readROIFile,0,0);
    realDataLayout->addWidget(_ROIFileName,0,1);
    realDataLayout->addWidget(dataLabel,1,0);
    realDataLayout->addWidget(_dataRefRegion,1,1);

    auto *startingPointLayout = new QHBoxLayout();
    startingPointLayout->addWidget(startingPointLabel);
    startingPointLayout->addWidget(_simulationStartingPoint);
    startingPointLayout->setEnabled(false);

    auto *realDataAndStartingLayout = new QVBoxLayout();
    realDataAndStartingLayout->addLayout(realDataLayout);
    realDataAndStartingLayout->addLayout(startingPointLayout);
    realDataGroupBox->setLayout(realDataAndStartingLayout);

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(setupGroupBox);
    rightLayout->addWidget(plasmaInGroupBox);
    rightLayout->addWidget(RRGroupBox);
    rightLayout->addWidget(realDataGroupBox);

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
    const QIcon *crossIcon = new QIcon(":/My-Icons/cursor-cross.png");
    auto *crossCursorAct = new QAction(*crossIcon, tr("Select time points"), this);
    dragXAction->setCheckable(true);
    dragYAction->setCheckable(true);
    rescaleXYAction->setCheckable(true);
    rescaleXYAction->setChecked(true);
    crossCursorAct->setCheckable(true);
    crossCursorAct->setChecked(true);
    QActionGroup *graphButtons = new QActionGroup(this);
    graphButtons->addAction(dragXAction);
    graphButtons->addAction(dragYAction);
    graphButtons->addAction(rescaleXYAction);

    QFrame* separator1 = new QFrame();
    separator1->setFrameShape(QFrame::VLine);
    separator1->setLineWidth(3);
    separator1->setFrameShadow(QFrame::Raised);
    QFrame* separator2 = new QFrame();
    separator2->setFrameShape(QFrame::VLine);
    separator2->setLineWidth(3);
    separator2->setFrameShadow(QFrame::Raised);
    QLabel *showLabel = new QLabel("show:",_setupPage);
    auto *radioShowPlasma    = new QRadioButton("All",_setupPage);
    auto *radioShowRR        = new QRadioButton("Ref Region",_setupPage);
    radioShowRR->setChecked(true);

    _whichPlasmaPlot  = new QComboBox();
    _whichPlasmaPlot->addItem("Cr (coarse) "); // 0
    _whichPlasmaPlot->addItem("Cp (coarse) "); // 1
    _whichPlasmaPlot->addItem("Cr (fine)   "); // 2
    _whichPlasmaPlot->addItem("Cp (fine)   "); // 3
    _whichPlasmaPlot->addItem("dt (coarse) "); // 4
    _whichPlasmaPlot->addItem("dt (fine)   "); // 5
    _whichPlasmaPlot->addItem("Cr_fit (coarse)");   // 6
    _whichPlasmaPlot->addItem("Cr_fit (fine)");     // 7
    _whichPlasmaPlot->addItem("dCr_fit/dt (fine)"); // 8
    _whichPlasmaPlot->addItem("dCr/dt (fine)");     // 9
    _whichPlasmaPlot->setCurrentIndex(0);
    _whichPlasmaPlot->setVisible(false);
    enableComboBoxItem(_whichPlasmaPlot,6,false);
    enableComboBoxItem(_whichPlasmaPlot,7,false);
    enableComboBoxItem(_whichPlasmaPlot,8,false);

    _clearPlasmaPlot = new QCheckBox("refresh plot");
    _clearPlasmaPlot->setChecked(true);
    _clearPlasmaPlot->setVisible(false);
    connect(_clearPlasmaPlot, SIGNAL(toggled(bool)), this, SLOT(clearPlasmaPlot(bool)));

    auto *showWidget = new QWidget();
    auto *showHBoxLayout = new QHBoxLayout();
    showHBoxLayout->addWidget(showLabel);
    showHBoxLayout->addWidget(_whichPlasmaPlot);
    showHBoxLayout->addWidget(_clearPlasmaPlot);
    showHBoxLayout->addWidget(radioShowPlasma);
    showHBoxLayout->addWidget(radioShowRR);
    showWidget->setLayout(showHBoxLayout);
    showWidget->setStyleSheet("color:Darkred");

    QToolBar *graphToolBar = new QToolBar("graph tool bar");
    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addAction(rescaleXYAction);
//    graphToolBar->addAction(crossCursorAct);
    graphToolBar->addWidget(separator1);
    graphToolBar->addWidget(showLabel);
    graphToolBar->addWidget(showWidget);
    graphToolBar->addWidget(separator2);

    fullLayout->setMenuBar(graphToolBar);

    // Make toolbar connections
    connect(dragXAction,     SIGNAL(triggered(bool)), _plasmaPlot, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plasmaPlot, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plasmaPlot, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plasmaPlot, SLOT(setSelectPoints()));
//    connect(crossCursorAct,  SIGNAL(triggered(bool)), _plasmaPlot, SLOT(setSelectPoints()));

    connect(dragXAction,     SIGNAL(triggered(bool)), _RRPlot, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _RRPlot, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _RRPlot, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _RRPlot, SLOT(setSelectPoints()));
//    connect(crossCursorAct,  SIGNAL(triggered(bool)), _RRPlot, SLOT(setSelectPoints()));

    connect(radioShowPlasma,   SIGNAL(clicked(bool)), this, SLOT(showPlasma()));
    connect(radioShowRR,       SIGNAL(clicked(bool)), this, SLOT(showRR()));

    connect(_tau2Ref,  SIGNAL(editingFinished()), this, SLOT(changedTau2Ref()));
    connect(_tau1Ref,  SIGNAL(editingFinished()), this, SLOT(changedTau1Ref()));
    connect(_plasmaFracRef, SIGNAL(editingFinished()), this, SLOT(changedPlasmaFracRef()));
    connect(_noiseRef, SIGNAL(editingFinished()), this, SLOT(changedNoiseRef()));

    connect(_readROIFile, SIGNAL(pressed()), this, SLOT(getTableDataFile()));
    connect(_dataRefRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(changedDataReferenceRegion()));

    connect(_whichPlasmaPlot, SIGNAL(currentIndexChanged(int)), this, SLOT(updateAllGraphs()));

    showRR();
    // showPlasma();

    _setupPage->setLayout(fullLayout);
}

void SimWindow::changedSimulationStartingPoint()
{ // 0) from plasma  1) fit plasma  2) fit data
    FUNC_ENTER;
    if ( !simStartsFromPlasma() )
    {
        dVector RRTimeVector = _simulator.getTimeCourse();
        dVector RRVector;
        if ( simStartsFromDataFit() )
            RRVector = _dataTable[_dataRefRegion->currentIndex()];
        else if ( simStartsFromPlasmaFit() )
            RRVector = _simulator.getCrCoarse();
        // LOESS: coarse (binned) RR
        double smoothingScale = 0.5;
        _quadLOESS.defineAndFit(RRTimeVector,RRVector, smoothingScale, true, true);
        dVector CrFitCoarse;
        for (int jTime=0; jTime<_simulator.getNumberTimeBinsCoarse(); jTime++)
            CrFitCoarse.append(_quadLOESS(jTime));
        _simulator.setCrFit(CrFitCoarse);

        enableComboBoxItem(_whichPlasmaPlot,6,true);
        enableComboBoxItem(_whichPlasmaPlot,7,true);
        enableComboBoxItem(_whichPlasmaPlot,8,true);
    }
    simulationStartingPoint startingPoint = static_cast<simulationStartingPoint>(_simulationStartingPoint->currentIndex());
    _simulator.setSimulationStartingPoint(startingPoint);
    updateAllGraphs();
}

void SimWindow::enableComboBoxItem(QComboBox *comboBox, int itemNumber, bool enable)
{
    QStandardItemModel* model = qobject_cast<QStandardItemModel*>(comboBox->model());
    QStandardItem* item= model->item(itemNumber);
    int nItems = comboBox->count();
    QVariant backgroundColor;  // setting default background color only works if at least 1 item is enabled!!
    for (int j=0; j<nItems; j++)
    {
        QStandardItem *thisItem = model->item(j);
        if ( thisItem->isEnabled() ) backgroundColor = comboBox->itemData(j,Qt::BackgroundRole);
    }
    if ( enable )
    {
        item->setFlags(item->flags() | Qt::ItemIsEnabled);   // enable
        comboBox->setItemData(itemNumber,backgroundColor,Qt::BackgroundRole);
    }
    else
    {
        item->setFlags(item->flags() & ~Qt::ItemIsEnabled);  // disable
        comboBox->setItemData(itemNumber,QColor(Qt::lightGray),Qt::BackgroundRole);
    }
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

void SimWindow::changedWeightType(int indexInBox)
{
    _PETRTM.setWeightingModel(indexInBox);
    _PETRTM.setPrepared(false);
    updateAllGraphs();
}

void SimWindow::changedModelType(int indexInBox)
{
    if ( indexInBox == RTM_SRTM3 )
        _PETRTM.setRTMModelType("SRTM3");
    else if ( indexInBox == RTM_SRTM2 )
        _PETRTM.setRTMModelType("SRTM2");
    else if ( indexInBox == RTM_rFRTM3 )
        _PETRTM.setRTMModelType("rFRTM3");
    else if ( indexInBox == RTM_rFRTM2 )
        _PETRTM.setRTMModelType("rFRTM2");
    else if ( indexInBox == RTM_SRTM2Fit )
        _PETRTM.setRTMModelType("SRTM2Fit");
    else if ( indexInBox == RTM_rFRTM3New )
        _PETRTM.setRTMModelType("rFRTM3New");
    else if ( indexInBox == RTM_rFRTM2New )
        _PETRTM.setRTMModelType("rFRTM2New");
    else if ( indexInBox == RTM_fmFRTM3 )
        _PETRTM.setRTMModelType("fmFRTM3");
    else if ( indexInBox == RTM_fmFRTM2 )
        _PETRTM.setRTMModelType("fmFRTM2");
    _PETRTM.setPrepared(false);
    updateAllGraphs();

    _radioShowAIC->setEnabled(_PETRTM.isForwardModel());
    _tau2RefAnalysisLabel->setVisible(_PETRTM.isRTM2());
    _tau2RefAnalysis->setVisible(_PETRTM.isRTM2());
    _errorR1Label->setVisible(_PETRTM.isRTM3());
    _errorR1->setVisible(_PETRTM.isRTM3());
    _errorTau2RefLabel->setVisible(_PETRTM.isRTM3());
    _errorTau2Ref->setVisible(_PETRTM.isRTM3());
    if ( _PETRTM.isRTM2() )
    {
        bool noisy = _simulator.getNoiseRef() != 0. || _simulator.getNoiseTar() != 0.;
        if ( indexInBox == RTM_SRTM2 && !noisy )
            _checkBoxTau2RefGraph->setText("1/k2' (for unbiased BPnd);\nUse numerical methods");
        else
            _checkBoxTau2RefGraph->setVisible(false);
        _checkBoxTau2RefGraph->setChecked(false);
    }
    else
    {
        _checkBoxTau2RefGraph->setVisible(true);
        _checkBoxTau2RefGraph->setText("1/k2' (=R1/k2)");
    }
    _tau4AnalysisLabel->setVisible( _PETRTM.isFRTM() || _PETRTM.isForwardModel() );
    _tau4Analysis->setVisible(      _PETRTM.isFRTM() || _PETRTM.isForwardModel() );
}

void SimWindow::createTargetPage()
{
    _targetPage = new QWidget();

    _basisPlot  = new plotData(2);
    _targetPlot = new plotData(3);

    //////// The setup page status bar
    _TRStatusBar = new QStatusBar;  // must be global so it doesn't go out of scope
    _TRStatusBar->setStyleSheet("color:blue");
    _targetPlot->setQCStatusBar(_TRStatusBar);
//    _targetPlot->setMainPage(this);

    connect(_targetPlot, SIGNAL(changedPointFromGraph(int,int,int)), _basisPlot,   SLOT(changePoint(int,int,int)));
    connect(_targetPlot, SIGNAL(changedPointFromGraph(int,int,int)), _RRPlot,      SLOT(changePoint(int,int,int)));
    connect(_basisPlot,  SIGNAL(changedPointFromGraph(int,int,int)), _targetPlot,  SLOT(changePoint(int,int,int)));
    connect(_basisPlot,  SIGNAL(changedPointFromGraph(int,int,int)), _RRPlot,      SLOT(changePoint(int,int,int)));
    connect(_RRPlot,     SIGNAL(changedPointFromGraph(int,int,int)), _basisPlot,   SLOT(changePoint(int,int,int)));
    connect(_RRPlot,     SIGNAL(changedPointFromGraph(int,int,int)), _targetPlot,  SLOT(changePoint(int,int,int)));

    _targetPlotLayout = new QVBoxLayout();
    _targetPlotLayout->addWidget(_basisPlot->getPlotSurface());
    _targetPlotLayout->addWidget(_targetPlot->getPlotSurface());
    _targetPlotLayout->addWidget(_TRStatusBar);
    _targetPlotLayout->setStretch(0,1);     // basis plot
    _targetPlotLayout->setStretch(1,1);     // target plot
    QString numberString;
    int editTextSize=80;

    // target region: input
    _targetSimulationGroupBox = new QGroupBox("Target Simulation: truth");
    QLabel *BPndLabel    = new QLabel("BPnd");
    QLabel *R1Label      = new QLabel("R1");
    QLabel *tau4Label    = new QLabel("1/k4 (min)");
    QLabel *challengeTimeLabel = new QLabel("Challenge time");
    QChar delta = QChar(0x0394);
    QLabel *challengeMagLabel  = new QLabel(QString("%1BPnd (%)").arg(delta));
    QLabel *plasmaFraclabel    = new QLabel("Plasma %");
    QLabel *noiseLabel   = new QLabel("Noise");
    _BPnd          = new QLineEdit();
    _R1            = new QLineEdit();
    _tau4          = new QLineEdit();
    _challengeTime = new QLineEdit();
    _challengeMag  = new QLineEdit();
    _plasmaFracTar = new QLineEdit();
    _noiseTar      = new QLineEdit();
    _BPnd->setText(numberString.setNum(_simulator.getBP0()));
    _R1->setText(numberString.setNum(_simulator.getR1()));
    _tau4->setText(numberString.setNum(_simulator.getTau4()));
    _challengeTime->setText(numberString.setNum(_simulator.getChallengeTime()));
    _challengeMag->setText(numberString.setNum(_simulator.getChallengeMag()));
    _noiseTar->setText(numberString.setNum(_simulator.getNoiseTar()));
    _plasmaFracTar->setText(numberString.setNum(_simulator.getPlasmaPercentTar()));
    _BPnd->setFixedWidth(editTextSize);
    _R1->setFixedWidth(editTextSize);
    _tau4->setFixedWidth(editTextSize);
    _challengeTime->setFixedWidth(editTextSize);
    _challengeMag->setFixedWidth(editTextSize);
    _noiseTar->setFixedWidth(editTextSize);
    _plasmaFracTar->setFixedWidth(editTextSize);
    auto *targetSimulationLayout = new QGridLayout();
    targetSimulationLayout->addWidget(BPndLabel,0,0);
    targetSimulationLayout->addWidget(_BPnd,0,1);
    targetSimulationLayout->addWidget(R1Label,1,0);
    targetSimulationLayout->addWidget(_R1,1,1);
    targetSimulationLayout->addWidget(tau4Label,2,0);
    targetSimulationLayout->addWidget(_tau4,2,1);
    targetSimulationLayout->addWidget(challengeTimeLabel,3,0);
    targetSimulationLayout->addWidget(_challengeTime,3,1);
    targetSimulationLayout->addWidget(challengeMagLabel,4,0);
    targetSimulationLayout->addWidget(_challengeMag,4,1);
    targetSimulationLayout->addWidget(plasmaFraclabel,5,0);
    targetSimulationLayout->addWidget(_plasmaFracTar,5,1);
    targetSimulationLayout->addWidget(noiseLabel,6,0);
    targetSimulationLayout->addWidget(_noiseTar,6,1);
    _targetSimulationGroupBox->setLayout(targetSimulationLayout);
    targetSimulationLayout->setSpacing(0);
    connect(_BPnd,          SIGNAL(editingFinished()), this, SLOT(changedBPND()));
    connect(_R1,            SIGNAL(editingFinished()), this, SLOT(changedR1()));
    connect(_tau4,          SIGNAL(editingFinished()), this, SLOT(changedTau4()));
    connect(_challengeTime, SIGNAL(editingFinished()), this, SLOT(changedChallengeTime()));
    connect(_challengeMag,  SIGNAL(editingFinished()), this, SLOT(changedChallengeMag()));
    connect(_plasmaFracTar, SIGNAL(editingFinished()), this, SLOT(changedPlasmaFracTar()));
    connect(_noiseTar,      SIGNAL(editingFinished()), this, SLOT(changedNoiseTar()));

    _targetDataGroupBox = new QGroupBox("Target ROI from imported table");
    auto *targetDataLayout = new QGridLayout();
    auto *targetROILabel = new QLabel("ROI");
    _dataTargetRegion = new QComboBox();
    targetDataLayout->addWidget(targetROILabel,0,0);
    targetDataLayout->addWidget(_dataTargetRegion,0,1);
    targetDataLayout->setSpacing(0);
    _targetDataGroupBox->setLayout(targetDataLayout);
    _targetDataGroupBox->setVisible(false);
    connect(_dataTargetRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(changedDataTargetRegion()));

    // target region: analysis
    _analysisGroupBox = new QGroupBox("Target Region: analysis");
    QLabel *modelLabel      = new QLabel("RTM type");
    QLabel *weightLabel     = new QLabel("Weighting scheme");
    QLabel *ignoreLabel     = new QLabel("Ignore points");
    _tau2RefAnalysisLabel   = new QLabel("Analysis 1/k2'");
    _tau4AnalysisLabel      = new QLabel("Analysis 1/k4");
    _modelType = new QComboBox();
    _modelType->addItem("SRTM3");
    _modelType->addItem("SRTM2");
    _modelType->addItem("SRTM2Fit");
    _modelType->addItem("rFRTM3");
    _modelType->addItem("rFRTM2");
    _modelType->addItem("rFRTM3New");
    _modelType->addItem("rFRTM2New");
    _modelType->addItem("fmFRTM3");
    _modelType->addItem("fmFRTM2");
    enableComboBoxItem(_modelType,2,false);
//    enableComboBoxItem(_modelType,5,false);  // xxx tmp
//    enableComboBoxItem(_modelType,6,false);  // xxx tmp
    _modelType->setCurrentIndex(0);
    _weightType  = new QComboBox();
    _weightType->addItem("Uniform weights");
    _weightType->addItem("11C-Noiseless");
    _weightType->addItem("11C");
    _weightType->addItem("Custom");
    _weightType->setToolTip("Set a weighting scheme for WLS");
    _ignoreString     = new QLineEdit();
    _tau2RefAnalysis = new QLineEdit();
    _tau4Analysis    = new QLineEdit();
    _ignoreString->setText("");
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
    _fitk4CheckBox = new QCheckBox("fit k4?");
    _fitk4CheckBox->setChecked(false);
//    _fitk4CheckBox->setEnabled(false);  // xxx tmp
    _fitChallenge = new QCheckBox("fit challenge?");
    _fitChallenge->setChecked(false);

    auto *analysisLayout = new QGridLayout();
    analysisLayout->addWidget(modelLabel,0,0);
    analysisLayout->addWidget(_modelType,0,1);
    analysisLayout->addWidget(_fitChallenge,1,0);
    analysisLayout->addWidget(_fitk4CheckBox,1,1);
    analysisLayout->addWidget(weightLabel,2,0);
    analysisLayout->addWidget(_weightType,2,1);
    analysisLayout->addWidget(ignoreLabel,3,0);
    analysisLayout->addWidget(_ignoreString,3,1);
    analysisLayout->addWidget(_tau2RefAnalysisLabel,4,0);
    analysisLayout->addWidget(_tau2RefAnalysis,4,1);
    analysisLayout->addWidget(_tau4AnalysisLabel,5,0);
    analysisLayout->addWidget(_tau4Analysis,5,1);
    _analysisGroupBox->setLayout(analysisLayout);
    analysisLayout->setSpacing(0);
    connect(_ignoreString,     SIGNAL(editingFinished()), this, SLOT(changedIgnoreString()));
    connect(_tau2RefAnalysis,  SIGNAL(editingFinished()), this, SLOT(changedTau2RefAnalysis()));
    connect(_tau4Analysis,     SIGNAL(editingFinished()), this, SLOT(changedTau4Analysis()));
    connect(_fitChallenge,     SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxChallenge(bool)));
    connect(_fitk4CheckBox,    SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxFitk4(bool)));
    // target region: errors
    auto *errorGroupBox = new QGroupBox("Target Region: Percentage Errors");
    QLabel *errorBPndLabel  = new QLabel("BPnd ");
    QLabel *errork2Label    = new QLabel("k2   ");
    QLabel *errork2aLabel   = new QLabel("k2a  ");
    _errordk2aLabel         = new QLabel("dk2a ");
    _errorR1Label           = new QLabel("R1   ");
    _errorTau2RefLabel      = new QLabel("1/k2'");
    _errorChallengeLabel    = new QLabel("Challenge");
    _errork4Label           = new QLabel("1/k4");
    _sigmaLabel             = new QLabel("sigma");
    _errorBPnd      = new QLabel();
    _errork2        = new QLabel();
    _errork2a       = new QLabel();
    _errordk2a      = new QLabel();
    _errorR1        = new QLabel();
    _errorTau2Ref   = new QLabel();
    _errorChallenge = new QLabel();
    _errork4        = new QLabel();
    _sigma          = new QLabel();
    _errorChallengeLabel->setVisible(false);
    _errorChallenge->setVisible(false);
    _errork4Label->setVisible(false);
    _errork4->setVisible(false);
    _errordk2a->setVisible(false);
    _errordk2aLabel->setVisible(false);
    auto *errorLayout = new QGridLayout();
    errorLayout->addWidget(errorBPndLabel,0,0);
    errorLayout->addWidget(_errorBPnd,0,1);
    errorLayout->addWidget(errork2Label,1,0);
    errorLayout->addWidget(_errork2,1,1);
    errorLayout->addWidget(errork2aLabel,2,0);
    errorLayout->addWidget(_errork2a,2,1);
    errorLayout->addWidget(_errordk2aLabel,3,0);
    errorLayout->addWidget(_errordk2a,3,1);
    errorLayout->addWidget(_errorR1Label,4,0);
    errorLayout->addWidget(_errorR1,4,1);
    errorLayout->addWidget(_errorTau2RefLabel,5,0);
    errorLayout->addWidget(_errorTau2Ref,5,1);
    errorLayout->addWidget(_errorChallengeLabel,6,0);
    errorLayout->addWidget(_errorChallenge,6,1);
    errorLayout->addWidget(_errork4Label,7,0);
    errorLayout->addWidget(_errork4,7,1);
    errorLayout->addWidget(_sigmaLabel,8,0);
    errorLayout->addWidget(_sigma,8,1);
    errorGroupBox->setLayout(errorLayout);
    errorLayout->setSpacing(0);

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(_targetSimulationGroupBox);
    rightLayout->addWidget(_targetDataGroupBox);
    rightLayout->addWidget(_analysisGroupBox);
    rightLayout->addWidget(errorGroupBox);
    rightLayout->setSpacing(0);

    QHBoxLayout *fullLayout = new QHBoxLayout();
    fullLayout->addLayout(_targetPlotLayout);
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

    QFrame* separator1 = new QFrame(_targetPage);
    separator1->setFrameShape(QFrame::VLine);
    separator1->setLineWidth(3);
    separator1->setFrameShadow(QFrame::Raised);
    QFrame* separator2 = new QFrame(_targetPage);
    separator2->setFrameShape(QFrame::VLine);
    separator2->setLineWidth(3);
    separator2->setFrameShadow(QFrame::Raised);
    QFrame* separator3 = new QFrame(_targetPage);
    separator3->setFrameShape(QFrame::VLine);
    separator3->setLineWidth(3);
    separator3->setFrameShadow(QFrame::Raised);

    QLabel *showLabel = new QLabel("show:",_targetPage);
    auto *radioShowBasisTarget = new QRadioButton("Basis/Target",_targetPage);
    auto *radioShowBasis       = new QRadioButton("Basis",_targetPage);
    auto *radioShowTarget      = new QRadioButton("Target",_targetPage);
    _radioShowAIC              = new QRadioButton("AIC",_targetPage);
    _radioShowAIC->setEnabled(false);
    QButtonGroup *showGroup = new QButtonGroup(_targetPage);
    showGroup->addButton(radioShowBasisTarget);
    showGroup->addButton(radioShowBasis);
    showGroup->addButton(radioShowTarget);
    showGroup->addButton(_radioShowAIC);
    radioShowTarget->setChecked(true);
    _basisPlot->getPlotSurface()->setVisible(false);
    _targetPlot->getPlotSurface()->setVisible(true);

    QLabel *analyzeLabel = new QLabel("analyze:",_targetPage);
    _analyzeSimulation = new QRadioButton("Simulation(s)",_targetPage);
    _analyzeRealData = new QRadioButton("ROI Data",_targetPage);
    QButtonGroup *analyzeGroup = new QButtonGroup(_targetPage);
    analyzeGroup->addButton(_analyzeSimulation);
    analyzeGroup->addButton(_analyzeRealData);
    _analyzeSimulation->setChecked(true);

    auto *showWidget = new QWidget();
    auto *showHBoxLayout = new QHBoxLayout();
    showHBoxLayout->addWidget(showLabel);
    showHBoxLayout->addWidget(radioShowBasisTarget);
    showHBoxLayout->addWidget(radioShowBasis);
    showHBoxLayout->addWidget(radioShowTarget);
    showHBoxLayout->addWidget(_radioShowAIC);
    showWidget->setLayout(showHBoxLayout);
    showWidget->setStyleSheet("color:Darkred");

    auto *analyzeWidget = new QWidget();
    auto *analyzeHBoxLayout = new QHBoxLayout();
    analyzeHBoxLayout->addWidget(analyzeLabel);
    analyzeHBoxLayout->addWidget(_analyzeSimulation);
    analyzeHBoxLayout->addWidget(_analyzeRealData);
    analyzeWidget->setLayout(analyzeHBoxLayout);
    analyzeWidget->setStyleSheet("color:Darkgreen");

    QToolBar *graphToolBar = new QToolBar("graph tool bar");
    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addAction(rescaleXYAction);
    graphToolBar->addWidget(separator1);
    graphToolBar->addWidget(showWidget);
    graphToolBar->addWidget(separator2);
    graphToolBar->addWidget(analyzeWidget);
    graphToolBar->addWidget(separator3);
    fullLayout->setMenuBar(graphToolBar);

    // Make toolbar connections to the main plot, not the basis plot
    connect(dragXAction,     SIGNAL(triggered(bool)), _basisPlot,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _basisPlot,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _basisPlot,  SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _targetPlot, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _targetPlot, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _targetPlot, SLOT(autoScale(bool)));
    connect(radioShowBasisTarget, SIGNAL(clicked(bool)), this, SLOT(showBasisTarget()));
    connect(radioShowBasis,  SIGNAL(clicked(bool)), this, SLOT(showBasis()));
    connect(radioShowTarget, SIGNAL(clicked(bool)), this, SLOT(showTarget()));
    connect(_radioShowAIC,   SIGNAL(clicked(bool)), this, SLOT(showAICvsTau4()));
    connect(_analyzeSimulation, SIGNAL(clicked(bool)), this, SLOT(clickedAnalyzeStimulation(bool)));
    connect(_analyzeRealData,   SIGNAL(clicked(bool)), this, SLOT(clickedAnalyzeRealData(bool)));

    _targetPage->setLayout(fullLayout);
}

void SimWindow::clickedAnalyzeStimulation(bool state)
{
    _calculateBPndCurves->setEnabled(state);
    _calculateTimeCurves->setEnabled(state);
    updateAllGraphs();
}
void SimWindow::clickedAnalyzeRealData(bool state)
{
    _calculateBPndCurves->setEnabled(!state);
    _calculateTimeCurves->setEnabled(!state);
    updateAllGraphs();
}

void SimWindow::createSweepBPndPage()
{
    // For RTM3:
    // 1) BPnd_err vs. BPnd
    // 2) Challenge error vs. BPnd
    // 3) Tau2Ref  vs. BPnd

    // For RTM2:
    // 1) BPnd_err vs. BPnd
    // 2) Challenge error vs. BPnd
    // 3) Tau2Ref (no bias) vs. BPnd (requires root finding)

    _sweepBPndPage = new QWidget();

    _errBPndPlot  = new plotData(4);
    _errChallPlot = new plotData(5);
    _tau2RefPlot  = new plotData(6);
    _errk4Plot    = new plotData(7);
    auto *plotLayout = new QVBoxLayout();
    plotLayout->addWidget(_errBPndPlot->getPlotSurface());
    plotLayout->addWidget(_errChallPlot->getPlotSurface());
    plotLayout->addWidget(_tau2RefPlot->getPlotSurface());
    plotLayout->addWidget(_errk4Plot->getPlotSurface());
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

    auto *calcGroupBox   = new QGroupBox("Calculate or clear curves");
    _nSamplesBPndPerThreadLabel = new QLabel("# samples/thread");
    _nSamplesBPndLabel = new QLabel("# samples");
    QLabel *calcLabel    = new QLabel("Calculate curves");
    QLabel *clearLabel   = new QLabel("Clear curves");
    _nSamplesBPndPerThread = new QLineEdit();
    _nSamplesBPndPerThread->setText(numberString.setNum(_numberSimulationsPerThread));
    _nSamplesBPnd = new QLabel(numberString.setNum(_nThreads*_numberSimulationsPerThread));
    _nSamplesBPnd->setFixedWidth(editTextSize);
    _nSamplesBPnd->setAlignment(Qt::AlignCenter);
    _calculateBPndCurves = new QPushButton();
    QPixmap pixmapCalculate(":/My-Icons/calculator.png");
    QIcon calculatorIcon(pixmapCalculate);
    _calculateBPndCurves->setIcon(calculatorIcon);
    _clearBPndCurves  = new QPushButton();
    QPixmap eraser(":/My-Icons/eraser.png");
    QIcon eraserIcon(eraser);
    _clearBPndCurves->setIcon(eraserIcon);
    auto *calcLayout = new QGridLayout();
    calcLayout->addWidget(_nSamplesBPndPerThreadLabel,0,0);
    calcLayout->addWidget(_nSamplesBPndPerThread,0,1);
    calcLayout->addWidget(_nSamplesBPndLabel,1,0);
    calcLayout->addWidget(_nSamplesBPnd,1,1);
    calcLayout->addWidget(calcLabel,2,0);
    calcLayout->addWidget(_calculateBPndCurves,2,1);
    calcLayout->addWidget(clearLabel,3,0);
    calcLayout->addWidget(_clearBPndCurves,3,1);
    calcGroupBox->setLayout(calcLayout);
    calcLayout->setSpacing(0);
    connect(_nSamplesBPndPerThread, SIGNAL(editingFinished()), this, SLOT(changedNumberSimulationsBPnd()));
    connect(_calculateBPndCurves, SIGNAL(pressed()), this, SLOT(calculateBPndCurves()));
    connect(_clearBPndCurves,     SIGNAL(pressed()), this, SLOT(clearBPndCurves()));

    _checkBoxBPndErrGraph  = new QCheckBox("BPnd error");
    _checkBoxChallErrGraph = new QCheckBox("Challenge error");
    _checkBoxTau2RefGraph  = new QCheckBox("1/k2' (=R1/k2)");
    _checkBoxk4ErrGraph    = new QCheckBox("1/k4");
    _sigma2Label = new QLabel("0");
    _checkBoxBPndErrGraph->setChecked(true);
    _checkBoxChallErrGraph->setChecked(false);
    _checkBoxChallErrGraph->setVisible(false);
    _checkBoxk4ErrGraph->setVisible(false);
    _checkBoxTau2RefGraph->setChecked(true);
    auto *checkGroupBox = new QGroupBox("Show or hide graphs");
    auto *checkLayout = new QVBoxLayout();
    checkLayout->addWidget(_checkBoxBPndErrGraph);
    checkLayout->addWidget(_checkBoxChallErrGraph);
    checkLayout->addWidget(_checkBoxTau2RefGraph);
    checkLayout->addWidget(_checkBoxk4ErrGraph);
    checkLayout->addWidget(_sigma2Label);
    checkGroupBox->setLayout(checkLayout);
    connect(_checkBoxBPndErrGraph,  SIGNAL(toggled(bool)), this, SLOT(changedVersusBPndGraphs()));
    connect(_checkBoxChallErrGraph, SIGNAL(toggled(bool)), this, SLOT(changedVersusBPndGraphs()));
    connect(_checkBoxTau2RefGraph,  SIGNAL(toggled(bool)), this, SLOT(changedVersusBPndGraphs()));
    connect(_checkBoxk4ErrGraph,    SIGNAL(toggled(bool)), this, SLOT(changedVersusBPndGraphs()));

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(BPndGroupBox);
    rightLayout->addWidget(calcGroupBox);
    rightLayout->addWidget(checkGroupBox);

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
    connect(dragXAction,     SIGNAL(triggered(bool)), _errk4Plot,    SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _errk4Plot,    SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _errk4Plot,    SLOT(autoScale(bool)));

    clearBPndCurves();
    changedVersusBPndGraphs();

    _sweepBPndPage->setLayout(fullLayout);
}

void SimWindow::setReferenceRegionForAnalysis()
{ // xxx this needs modification
    if ( analyzeRealData() )
    {
        FUNC_INFO << "defineRTMModel real data";
        // matrices: [nRuns][nTime]
        // set reference region to data RR
        iMatrix timeBinsSec;  timeBinsSec.resize(1);
        dMatrix refRegion; refRegion.resize(1);
        timeBinsSec[0] = _dtBinsSec;          // _dataTable[nColums]
        refRegion[0] = _dataTable[_dataRefRegion->currentIndex()];  // getYData(iCurve)
        _PETRTM.setReferenceRegion(timeBinsSec,refRegion);
    }
    else
    {
        FUNC_INFO << "defineRTMModel simulation";
        int nTime = _simulator.getNumberTimeBinsCoarse();
        iMatrix timeBinsSec;  timeBinsSec.resize(1);
        timeBinsSec[0] = _simulator.getTimeBinVectorSec();
        dMatrix refRegion; refRegion.resize(1);
        refRegion[0].resize(nTime);
        if ( simStartsFromPlasma() )
        {
            for (int jt=0; jt<nTime; jt++)
                refRegion[0][jt] = _simulator.getCrCoarse(jt);
        }
        else
        {
            for (int jt=0; jt<nTime; jt++)
                refRegion[0][jt] = _simulator.getCrFitBinned(jt);
        }
        _PETRTM.setReferenceRegion(timeBinsSec,refRegion);
    }
}

void SimWindow::createVersusTimePage()
{
    // 1) BPnd_err vs. time
    // 2) Challenge error vs. time

    _sweepTimePage = new QWidget();

    _errBPndOrChallVsTimePlot = new plotData(6);
    auto *plotLayout = new QVBoxLayout();
    plotLayout->addWidget(_errBPndOrChallVsTimePlot->getPlotSurface());
    QString numberString;
    int editTextSize=80;

    auto *timeGroupBox    = new QGroupBox("challenge time lower & upper limits");
    QLabel *timeLowLabel  = new QLabel("lowest time");
    QLabel *timeHighLabel = new QLabel("highest time");
    _timeLow   = new QLineEdit();
    _timeHigh  = new QLineEdit();
    _timeLow->setFixedWidth(editTextSize);
    _timeHigh->setFixedWidth(editTextSize);
    _timeLow->setText(numberString.setNum(_timeLowValue));
    _timeHigh->setText(numberString.setNum(_timeHighValue));
    connect(_timeLow,   SIGNAL(editingFinished()), this, SLOT(changedTimeLow()));
    connect(_timeHigh,  SIGNAL(editingFinished()), this, SLOT(changedTimeHigh()));
    auto *timeLayout = new QGridLayout();
    timeLayout->addWidget(timeLowLabel,0,0);
    timeLayout->addWidget(_timeLow,0,1);
    timeLayout->addWidget(timeHighLabel,1,0);
    timeLayout->addWidget(_timeHigh,1,1);
    timeGroupBox->setLayout(timeLayout);
    timeLayout->setSpacing(0);

    auto *calcGroupBox = new QGroupBox("Calculate or clear curves");
    QLabel *calcLabel  = new QLabel("Calculate curves");
    QLabel *clearLabel = new QLabel("Clear curves");
    _calculateTimeCurves   = new QPushButton();
    QPixmap pixmapCalculate(":/My-Icons/calculator.png");
    QIcon calculatorIcon(pixmapCalculate);
    _calculateTimeCurves->setIcon(calculatorIcon);
    _clearTimeCurves  = new QPushButton();
    QPixmap eraser(":/My-Icons/eraser.png");
    QIcon eraserIcon(eraser);
    _clearTimeCurves->setIcon(eraserIcon);
    auto *calcLayout = new QGridLayout();
    calcLayout->addWidget(calcLabel,0,0);
    calcLayout->addWidget(_calculateTimeCurves,0,1);
    calcLayout->addWidget(clearLabel,1,0);
    calcLayout->addWidget(_clearTimeCurves,1,1);
    calcGroupBox->setLayout(calcLayout);
    calcLayout->setSpacing(0);
    connect(_calculateTimeCurves, SIGNAL(pressed()), this, SLOT(calculateTimeCurves()));
    connect(_clearTimeCurves, SIGNAL(pressed()),     this, SLOT(clearTimeCurves()));

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(timeGroupBox);
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
    connect(dragXAction,     SIGNAL(triggered(bool)), _errBPndOrChallVsTimePlot,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _errBPndOrChallVsTimePlot,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _errBPndOrChallVsTimePlot,  SLOT(autoScale(bool)));

    clearTimeCurves();

    _sweepTimePage->setLayout(fullLayout);

}

void SimWindow::showPlasma()
{
    _plasmaPlot->getPlotSurface()->setVisible(true);
    _plasmaStatusBar->setVisible(true);
    _RRPlot->getPlotSurface()->setVisible(false);
    _RRStatusBar->setVisible(false);
    _whichPlasmaPlot->setVisible(true);
    _clearPlasmaPlot->setVisible(true);
}
void SimWindow::showRR()
{
    _plasmaPlot->getPlotSurface()->setVisible(false);
    _plasmaStatusBar->setVisible(false);
    _RRPlot->getPlotSurface()->setVisible(true);
    _RRStatusBar->setVisible(true);
    _whichPlasmaPlot->setVisible(false);
    _clearPlasmaPlot->setVisible(false);
}
void SimWindow::showBasisTarget()
{
    _basisPlot->getPlotSurface()->setVisible(true);
    _targetPlot->getPlotSurface()->setVisible(true);
}
void SimWindow::showBasis()
{
    _basisPlot->getPlotSurface()->setVisible(true);
    _targetPlot->getPlotSurface()->setVisible(false);
}
void SimWindow::showTarget()
{
    _basisPlot->getPlotSurface()->setVisible(false);
    _targetPlot->getPlotSurface()->setVisible(true);
    updateAllGraphs();
}
void SimWindow::showAICvsTau4()
{
    _basisPlot->getPlotSurface()->setVisible(false);
    _targetPlot->getPlotSurface()->setVisible(true);
    updateAICvsTau4Graph();
}

void SimWindow::updatePlasmaGraph()
{
    FUNC_ENTER;
    static int whichColor=0;
    QColor colors[10] = {Qt::black, Qt::red, Qt::blue, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};

    int index = _whichPlasmaPlot->currentIndex();
    bool fineScale = (index == 2) || (index == 3) || (index == 5) || (index == 7) || (index == 8) || (index==9);

    // update the plot: RR
    if ( _clearPlasmaPlot->isChecked() )
        _plasmaPlot->init();
    _plasmaPlot->setLegendOn(true);
    if ( index != 0 )
        _plasmaPlot->addCurve(0,_whichPlasmaPlot->currentText());

    if ( index == 0 )
    { // Cr coarse
        if ( analyzeRealData() )
            addDataCurveRR(_plasmaPlot);
        else
            addSimulationRR(_plasmaPlot);
    }
    else if ( fineScale )
    {
        int nTime = _simulator.getNumberTimeBinsFine();
        dVector xTime; xTime.resize(nTime);
        dVector yTAC;  yTAC.resize(nTime);
        for (int jt=0; jt<nTime; jt++)
        {
            xTime[jt] = _simulator.getTimeFine(jt);
            if ( index == 2 )
                yTAC[jt]  = _simulator.getCr(jt);
            else if ( index == 3 )
                yTAC[jt]  = _simulator.getCp(jt);
            else if ( index == 5 )
                yTAC[jt]  = _simulator.getdtFine(jt);
            else if ( index == 7 )
                yTAC[jt]  = _simulator.getCrFit(jt);
            else if ( index == 8 )
                yTAC[jt]  = _simulator.getCrFitDot(jt);
            else if ( index == 9 )
                yTAC[jt]  = _simulator.getCrDot(jt);
        }
        _plasmaPlot->setData(xTime,yTAC);
    }
    else
    { // Cp or dt in coarse units
        index = _whichPlasmaPlot->currentIndex();
        int nTime = _simulator.getNumberTimeBinsCoarse();
        dVector xTime; xTime.resize(nTime);
        dVector yTAC;  yTAC.resize(nTime);
        for (int jt=0; jt<nTime; jt++)
        {
            xTime[jt] = _simulator.getTimeCoarse(jt);
            if ( index == 1 )
                yTAC[jt] = _simulator.getCpCoarse(jt);
            else if ( index == 4 )
                yTAC[jt]  = _simulator.getdtCoarse(jt);
            else if ( index == 6 )
                yTAC[jt]  = _simulator.getCrFitBinned(jt);
        }
        _plasmaPlot->setData(xTime,yTAC);
    }

    if ( _clearPlasmaPlot->isChecked() || whichColor > 9 )
        whichColor = 0;
    FUNC_INFO << "which color" << whichColor;
    _plasmaPlot->setColor(colors[whichColor]);
    whichColor++;
    if ( fineScale )
        _plasmaPlot->setPointSize(2);
    else
        _plasmaPlot->setPointSize(5);
    _plasmaPlot->conclude(0,true);
    _plasmaPlot->plotDataAndFit(true);
    FUNC_EXIT;
}

void SimWindow::updateReferenceGraph()
{
    FUNC_ENTER;
    // update the plot: RR
    _RRPlot->init();
    _RRPlot->setLegendOn(true);

    if ( simStartsFromPlasmaFit() || simStartsFromDataFit() )
        addFitRR(_RRPlot);
    if ( simStartsFromPlasma() || simStartsFromPlasmaFit() )
        addSimulationRR(_RRPlot);
    else
        addDataCurveRR(_RRPlot);

    _RRPlot->conclude(0,true);
    _RRPlot->plotDataAndFit(true);
    FUNC_EXIT;
}

void SimWindow::addSimulationRR(plotData *whichPlot)
{
    FUNC_ENTER;
    whichPlot->addCurve(0,"RR: simulation");
    if ( whichPlot->getNumberCurves() == 1 )
        whichPlot->setColor(Qt::red); // 1st curve
    else
        whichPlot->setColor(Qt::gray);
    whichPlot->setPointSize(5);
    int nTimeSim = _simulator.getNumberTimeBinsCoarse();
    dVector xTime;   xTime.resize(nTimeSim);
    dVector yTAC;    yTAC.resize(nTimeSim);
    for (int jt=0; jt<nTimeSim; jt++)
    {
        xTime[jt] = _simulator.getTimeCoarse(jt);
        yTAC[jt] = _simulator.getCrCoarse(jt);
    }
    whichPlot->setData(xTime,yTAC);
    FUNC_EXIT;
}
void SimWindow::addFitRR(plotData *whichPlot)
{
    whichPlot->addCurve(0,"RR: fit");
    if ( whichPlot->getNumberCurves() == 1 )
        whichPlot->setColor(Qt::red); // 1st curve is red
    else
        whichPlot->setColor(Qt::gray); // 2nd curve is gray
    whichPlot->setPointSize(2);
    int nTimeSim = _simulator.getNumberTimeBinsCoarse();
    dVector xTime;   xTime.resize(nTimeSim);
    dVector yTAC;    yTAC.resize(nTimeSim);
    for (int jt=0; jt<nTimeSim; jt++)
    {
        xTime[jt] = _simulator.getTimeCoarse(jt);
        yTAC[jt]  = _quadLOESS(jt);
    }
    whichPlot->setData(xTime,yTAC);
}
void SimWindow::addDataCurveRR(plotData *whichPlot)
{
    if ( !realDataAvailable() ) return;
    int indexInBox = _dataRefRegion->currentIndex();
    if ( indexInBox >= _dataTable.size() )
        qFatal("Fatal Error: the combo-box index exceeds the table index.");
    whichPlot->addCurve(0,"RR: ROI data");
    if ( whichPlot->getNumberCurves() == 1 )
        whichPlot->setColor(Qt::red); // 1st curve
    else
        whichPlot->setColor(Qt::gray);
    whichPlot->setPointSize(5);
    whichPlot->setData(_timeBins,_dataTable[indexInBox]);
    // enable buttons
    _analyzeSimulation->setEnabled(true);
    _analyzeRealData->setEnabled(true);
    _targetDataGroupBox->setVisible(true);
}
void SimWindow::addSimulationTarget()
{
    FUNC_ENTER;
    if ( realDataAvailable() )
        _targetPlot->addCurve(0,"target: simulation");
    else
        _targetPlot->addCurve(0,"target");
    if ( _targetPlot->getNumberCurves() == 1 )
        _targetPlot->setColor(Qt::red); // 1st curve
    else
        _targetPlot->setColor(Qt::gray);
    _targetPlot->setPointSize(3);
    int nTime = _simulator.getNumberTimeBinsCoarse();
    dVector xTime;   xTime.resize(nTime);
    dVector yTAC;    yTAC.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = _simulator.getTimeCoarse(jt);
        yTAC[jt]  = _simulator.getCtCoarse(jt);
    }
    _targetPlot->setData(xTime,yTAC);
    FUNC_EXIT;
}
void SimWindow::addDataCurveTarget()
{
    if ( realDataAvailable() )
    {
        _targetPlot->addCurve(0,"target: ROI data");
        if ( _targetPlot->getNumberCurves() == 1 )
            _targetPlot->setColor(Qt::red); // 1st curve
        else
            _targetPlot->setColor(Qt::gray);
        _targetPlot->setPointSize(5);
        int indexTarget = _dataTargetRegion->currentIndex();
        if ( indexTarget >= _dataTable.size() )
            qFatal("Fatal Error: the combo-box index exceeds the table index.");
        else if ( indexTarget < 0 ) return;
        _targetPlot->setData(_timeBins,_dataTable[indexTarget]);
    }
}

void SimWindow::updateTargetGraph()
{
    FUNC_ENTER;
    // update the plot
    _targetPlot->init();
    _targetPlot->setLegendOn(true);

    if ( analyzeRealData() )
    {
        addDataCurveTarget();
        addSimulationTarget();
    }
    else
    {
        addSimulationTarget();
        addDataCurveTarget();
    }

    int nTime = _PETRTM.getNumberTimePointsInRun(0);
    dVector xTime;   xTime.resize(nTime);
    dVector yTAC;    yTAC.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
        xTime[jt] = _PETRTM.getReferenceRegion(false,0,jt).x;

    // add fit
    FUNC_INFO << "add fit";
    _targetPlot->addCurve(0,"fit");
    _targetPlot->setLineThickness(2);
    _targetPlot->setColor(Qt::blue);
    if ( _PETRTM.isForwardModel() )
    {
        for (int jt=0; jt<nTime; jt++)
            yTAC[jt] = _PETRTM.getSimulationFit(0,jt);
    }
    else
    {
        for (int jt=0; jt<nTime; jt++)
            yTAC[jt] = _PETRTM.getFit(jt);
    }
    _targetPlot->setData(xTime,yTAC);

    // add RR
    _targetPlot->addCurve(0,"RR");
    _targetPlot->setColor(Qt::gray);
    for (int jt=0; jt<nTime; jt++)
        yTAC[jt] = _PETRTM.getReferenceRegion(false,0,jt).y;
    _targetPlot->setData(xTime,yTAC);

    // add FRTM convolution
    if ( _PETRTM.isFRTM() )
    {
        if ( _PETRTM.isFRTMNew() )
//            _targetPlot->addCurve(0,"BP0*convolution");
            _targetPlot->addCurve(0,"convolution");
        else
            _targetPlot->addCurve(0,"convolution");
        _targetPlot->setColor(Qt::magenta);
        for (int jt=0; jt<nTime; jt++)
            yTAC[jt]  = _simulator.getBP0() * _PETRTM.getFRTMConvolution(0,jt);
//            yTAC[jt]  = _PETRTM.getFRTMConvolution(0,jt);
        _targetPlot->setData(xTime,yTAC);
    }
    if ( _PETRTM.isFRTMNew() )
    {
        _targetPlot->addCurve(0,"Ct-Cr");
        _targetPlot->setColor(Qt::cyan);
        for (int jt=0; jt<nTime; jt++)
            yTAC[jt]  = _PETRTM.getCtMinusCr(0,jt);
        _targetPlot->setData(xTime,yTAC);
    }

    _targetPlot->conclude(0,true);
    _targetPlot->plotDataAndFit(true);
    FUNC_EXIT;
}

void SimWindow::defineRTMModel()
{
    static int firstTime=true;
    setReferenceRegionForAnalysis();
    if ( firstTime )
    {
        // Assign these IDs:
        _PETRTM.setR1EventID(0,'R');
        _PETRTM.setk2EventID(0,'k');
        _PETRTM.setk2aEventID(0,'a');
        int indexChallenge = _PETRTM.getEventIndex('c');
        _PETRTM.setChallengeShape(indexChallenge,Challenge_Sigmoid);
        _PETRTM.setChallengeTau(indexChallenge,0.1);
        _PETRTM.setChallengeOnset(indexChallenge,0,_simulator.getChallengeTime());
        _PETRTM.setTau2RefSRTM(0,_simulator.getTau2Ref());
        _PETRTM.setTau2RefFRTM(0,_simulator.getTau2Ref());
        _PETRTM.setTau2Ref(0,_simulator.getTau2Ref());
        _PETRTM.setTau4(0,_simulator.getTau4());
    }
    firstTime = false;
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
    for ( int jBasis=0; jBasis<nBasis; jBasis++)
    {
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
            double sourceMax = _basisPlot->getMaxYAbs(&sourceCurve) / yFraction;
            double targetMax = _basisPlot->getMaxYAbs(targetCurve);
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

void SimWindow::updateAICvsTau4Graph()
{
    dVector tau4Vector, AICVector, tau2RefVector;
    tau4Vector.clear();  AICVector.clear();  tau2RefVector.clear();
    bool ok;  double bestTau4 = _tau4Analysis->text().toDouble(&ok);
    if ( bestTau4 < 10. ) bestTau4 = 10.;
    if ( _PETRTM.isRTM3() && !_PETRTM.getFitk4State() )
    {
        double bestTau2Ref;
        _PETRTM.calculateTau4andTau2Ref(0,AICVector, tau4Vector, tau2RefVector, bestTau4, bestTau2Ref);
        _tau4Analysis->setText(QString::number(bestTau4,'g',3));
        _tau2RefAnalysis->setText(QString::number(bestTau2Ref));
    }
    else if ( !_PETRTM.getFitk4State() )
    {
        _PETRTM.calculateTau4atFixedTau2Ref(0,AICVector, tau4Vector, bestTau4);
        _tau4Analysis->setText(QString::number(bestTau4,'g',3));
    }
    if ( isnan(bestTau4) ) bestTau4 = 10.;
    _PETRTM.setTau4Nominal(bestTau4);

    _targetPlot->init();
    _targetPlot->setLegendOn(false);
    _targetPlot->addCurve(0,"AIC vs 1/k4");
    _targetPlot->setData(tau4Vector, AICVector);
    _targetPlot->conclude(0,true);
    _targetPlot->plotDataAndFit(true);
    _targetPlot->setPositionTracer(0,bestTau4);
}

void SimWindow::updateAllGraphs()
{
    // run the simulation
    _simulator.run();

    updatePlasmaGraph();
    setReferenceRegionForAnalysis();
    updateReferenceGraph();
    analyzeTAC();
    if ( _radioShowAIC->isChecked() )
        updateAICvsTau4Graph();
    else
        updateTargetGraph();
    updateBasisGraph();  // update basis graph AFTER target graph, because basis functions use target curve
    double AIC;
    if ( _PETRTM.isForwardModel() )
        AIC = _PETRTM.getSimulationAIC(0);
    else
        AIC = _PETRTM.getAIC();
    _analysisGroupBox->setTitle(QString("Target Region: analysis, AIC = %1").arg(AIC));
}
void SimWindow::analyzeTAC()
{
    // define analysis
    defineRTMModel();
    // update the RTM model

    int nTime;
    dMatrix tissueVector; tissueVector.resize(1);
    if ( analyzeRealData() )
    {
        nTime = _dataTable[0].size();
        int indexTarget = _dataTargetRegion->currentIndex();
        tissueVector[0] = _dataTable[indexTarget];
    }
    else
    {
        nTime = _simulator.getNumberTimeBinsCoarse();
        tissueVector[0].resize(nTime);
        for (int jt=0; jt<nTime; jt++)
            tissueVector[0][jt] = _simulator.getCtCoarse(jt);
    }
    _PETRTM.setTissueVector(tissueVector);
    _PETRTM.definePETConditions("a c"); // don't define R1, which is not valid for RTM2
    _PETRTM.prepare();
    dMatrix fitVector;     fitVector.resize(1);     fitVector[0].resize(nTime);
    _PETRTM.fitData(tissueVector,fitVector);

    // BPnd
    double truth = _simulator.getBP0();
    double guess = _PETRTM.getBP0InRun(0);
    _errorBPnd->setText(analyzeString(truth,guess));
    // k2
    truth = _simulator.getk2();
    guess = _PETRTM.getk2InRun(0).x;
    _errork2->setText(analyzeString(truth,guess));
    // k2a
    if ( _PETRTM.isFRTMNew() )
    { // get the "k2a" parameter, which really is k4 * k2 * BPND
        double k2 = _simulator.getk2();
        double BP0 = _simulator.getBP0();
        truth = k2 * BP0;
    }
    else
        truth = _simulator.getk2a();
    guess = _PETRTM.getk2aInRun(0).x;
    _errork2a->setText(analyzeString(truth,guess));
    QString valueString;
    if ( _fitChallenge->isChecked() )
    {
        if ( _PETRTM.isFRTMNew() )
        { // dk2k3
            truth = _simulator.getdk2k3();
            guess = _PETRTM.getdk2aInRun('c').x;
            _errordk2a->setText(analyzeString(truth,guess));
        }
        else
        { // dka2
            truth = _simulator.getdk2a();
            guess = _PETRTM.getdk2aInRun('c').x;
            _errordk2a->setText(analyzeString(truth,guess));
        }
    }

    // R1
    truth = _simulator.getR1();
    guess = _PETRTM.getR1InRun(0).x;
    _errorR1->setText(analyzeString(truth,guess));
    // tau2Ref
    truth = _simulator.getTau2Ref();
    guess = _PETRTM.getTau2RefInRun(0);
    _errorTau2Ref->setText(analyzeString(truth,guess));
    dPoint2D k2Ref = _PETRTM.getk2RefInRun(0);

    // challenge
    if ( _fitChallenge->isChecked() )
    {
        truth = _simulator.getChallengeMag();  // delta_BPnd abs
        guess = getChallengeMagFromAnalysis();
        QString diffString;
        valueString.setNum(guess,'g',2);
        if ( guess == 0. )
            diffString = "----";
        else
            diffString.setNum(guess-truth,'g',2);
        if ( analyzeRealData() )
            _errorChallenge->setText(valueString + " %");
        else
            _errorChallenge->setText(valueString + " %, diff = " + diffString + " %");
    }

    // k4
    if ( _fitk4CheckBox->isChecked() )
    {
        truth = _simulator.getTau4();  // delta_BPnd abs
        guess = _PETRTM.getTau4InRun(0);
        valueString.setNum(guess,'g',2);
        if ( guess == 0. )
            _errork4->setText("----");
        else
            _errork4->setText(analyzeString(truth,guess));
    }

    // sigma2
    double sigma;
    if ( _PETRTM.isForwardModel() )
        sigma = qSqrt(_PETRTM.getSimulationSigma2(0));
    else
        sigma = qSqrt(_PETRTM.getSigma2());
    valueString.setNum(sigma,'g',3);
    _sigma->setText(valueString);

    if ( analyzeRealData() )
    {
//        if ( changedAnalysisRegion )
//        {
            double BPnd = _PETRTM.getBP0InRun(0);
            double R1   = _PETRTM.getR1InRun(0).x;
            _simulator.setBP0(BPnd);  _simulator.setR1(R1);
            QString text; _BPnd->setText(text.setNum(BPnd));  _R1->setText(text.setNum(R1));
            _simulator.generateTargetTAC();
//        }
        addSimulationTarget();
        updateTargetGraph();
        _targetPlot->conclude(0,true);
        _targetPlot->plotDataAndFit(true);
    }

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
    if ( analyzeRealData() )
    {
        QString valueString;
        valueString.setNum(guess,'g',3);
        return QString("%1").arg(valueString);
    }
    else
    {
        QString percentString, valueString;
        percentString.setNum(percentageError(guess,truth),'g',2);
        valueString.setNum(guess,'g',3);
        return percentString + " % " + QString("(abs = %1)").arg(valueString);
    }
}

void SimWindow::scaleSimulationToDataAverage()
{
    if ( !realDataAvailable() ) return;

    // update the simulation
    _simulator.run();

    dVector data = _PETRTM.getReferenceRegionVector(0);
    double dataAverage = 0.;
    for (int jBin=0; jBin<data.size(); jBin++)
        dataAverage += data[jBin];
    dataAverage /= static_cast<double>(data.size());

    dVector sim; sim.resize(data.size());
    double simAverage = 0.;
    for (int jBin=0; jBin<data.size(); jBin++)
        simAverage += _simulator.getCrCoarse(jBin);
    simAverage /= static_cast<double>(data.size());

    double ratio = dataAverage / simAverage;

    bool ok;  QString bolusText;
    double bolusMag = _bolusMag->text().toDouble(&ok);
    bolusMag *= ratio;
    _bolusMag->setText(bolusText.setNum(bolusMag,'g',4));  changedBolusMag();

}

///////////////////////////////////////
// Slots
///////////////////////////////////////
void SimWindow::changedNumberThreads(int indexInBox)
{
    _nThreads = indexInBox + 1;
    QString numberString;
    _nSamplesBPnd->setText((numberString.setNum(_nThreads*_numberSimulationsPerThread)));
}

void SimWindow::changedNumberBins()
{
    bool ok;
    int numberBins = _numberTimeBins->text().toInt(&ok);
    if ( ok && numberBins != _simulator.getNumberBins() )
    {
        _simulator.setNumberBins(numberBins);
        if ( _PETRTM.getNumberTimePointsInRun(0) != numberBins ) _PETRTM.setTimePointsInRun(0,numberBins);
        _PETRTM.setTimeBinsSec(0,_simulator.getTimeBinVectorSec());
        _binIndex->setRange(1,numberBins);
        updateAllGraphs();
    }
    else
    {
        QString text;
        _numberTimeBins->setText(text.setNum(_simulator.getNumberTimeBinsCoarse()));
    }
}
void SimWindow::changedBinIndex(int indexPlusOne)
{
    int iBin = indexPlusOne-1;
    double duration = _simulator.getDurationPerBinSec(iBin);
    QString text; _binDuration->setText(text.setNum(duration));
}
void SimWindow::changedBinDuration()
{
    bool ok;
    int duration = _binDuration->text().toInt(&ok);
    int iBin = _binIndex->value() - 1;
    if ( ok )
    {
        if ( _applyToAllBinDuration->isChecked() )
        {
            int numberBins = _numberTimeBins->text().toInt(&ok);
            for (int jBin=0; jBin<numberBins; jBin++)
                _simulator.setDurationBin(jBin,duration);  // convert sec to min
            _PETRTM.setTimeBinsSec(0,_simulator.getTimeBinVectorSec());
        }
        else
        {
            _simulator.setDurationBin(iBin,duration);  // convert sec to min
            _PETRTM.setTimeBinsSec(0,_simulator.getTimeBinVectorSec());
        }
        updateAllGraphs();
    }
    else
    {
        QString text;
        _binDuration->setText(text.setNum(_simulator.getDurationPerBinSec(iBin)));
    }
}
void SimWindow::changedApplyToAllBinDuration(bool state)
{
    if ( state )
    {
        bool ok;
        int numberBins = _numberTimeBins->text().toInt(&ok);
        int duration = _binDuration->text().toInt(&ok);
        for (int jBin=0; jBin<numberBins; jBin++)
            _simulator.setDurationBin(jBin,duration);  // convert sec to min
        _PETRTM.setTimeBinsSec(0,_simulator.getTimeBinVectorSec());
        updateAllGraphs();
    }
}
void SimWindow::changedSubSample()
{
    QString stringEntered = _subSample->text();
    bool ok;
    int lSubSample = stringEntered.toInt(&ok);
    int iBin = _binIndex->value() - 1;
    if ( ok )
    {
        _simulator.setSamplesPerBin(iBin,lSubSample);
        _PETRTM.setTimeBinsSec(0,_simulator.getTimeBinVectorSec());
        updateAllGraphs();
    }
    else
        _subSample->setText(stringEntered.setNum(_simulator.getSamplesPerBin(iBin)));
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
        _simulator.setTauBolus(value); scaleSimulationToDataAverage();
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
        _simulator.setKBol(value);  scaleSimulationToDataAverage();
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
        _simulator.setTau2Ref(value);  scaleSimulationToDataAverage();
        updateAllGraphs();
    }
    else
        _tau2Ref->setText(stringEntered.setNum(_simulator.getTau2Ref()));
}
void SimWindow::changedTau1Ref()
{
    FUNC_ENTER;
    QString stringEntered = _tau1Ref->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setTau1Ref(value);  scaleSimulationToDataAverage();
        updateAllGraphs();
    }
    else
        _tau1Ref->setText(stringEntered.setNum(_simulator.getTau1Ref()));
    FUNC_EXIT << value;
}

void SimWindow::changedPlasmaFracRef()
{
    QString utilityString = _plasmaFracRef->text();
    bool ok;
    double value = utilityString.toDouble(&ok);
    if ( ok )
    {
        _simulator.setPlasmaPercentRef(value);  scaleSimulationToDataAverage();
        updateAllGraphs();
    }
    else
        _plasmaFracRef->setText(utilityString.setNum(_simulator.getPlasmaPercentRef()));
}
void SimWindow::changedPlasmaFracTar()
{
    QString utilityString = _plasmaFracTar->text();
    bool ok;
    double value = utilityString.toDouble(&ok);
    if ( ok )
    {
        _simulator.setPlasmaPercentTar(value);
        updateAllGraphs();
    }
    else
        _plasmaFracTar->setText(utilityString.setNum(_simulator.getPlasmaPercentTar()));
}
void SimWindow::changedNoiseRef()
{
    QString utilityString = _noiseRef->text();
    bool ok;
    double value = utilityString.toDouble(&ok);
    if ( ok )
    {
        _simulator.setNoiseRef(value);
        updateAllGraphs();
    }
    else
        _noiseRef->setText(utilityString.setNum(_simulator.getNoiseRef()));
    bool noisy = _simulator.getNoiseRef() != 0. || _simulator.getNoiseTar() != 0.;
    setThreadVisibility(noisy);
}
void SimWindow::changedFastElimination()
{
    QString stringEntered = _fastTau->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setTauFastElim(value);  scaleSimulationToDataAverage();
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
        _simulator.setTauSlowElim(value);  scaleSimulationToDataAverage();
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
        _simulator.setFastFraction(value);  scaleSimulationToDataAverage();
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
        int indexChallenge = _PETRTM.getEventIndex('c');
        _PETRTM.setChallengeOnset(indexChallenge,0,_simulator.getChallengeTime());
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
    QString utilityString = _noiseTar->text();
    bool ok;
    double value = utilityString.toDouble(&ok);
    if ( ok )
    {
        _simulator.setNoiseTar(value);
        updateAllGraphs();
    }
    else
        _noiseTar->setText(utilityString.setNum(_simulator.getNoiseTar()));
    bool noisy = _simulator.getNoiseRef() != 0. || _simulator.getNoiseTar() != 0.;
    setThreadVisibility(noisy);
}
void SimWindow::changedIgnoreString()
{
    QString stringEntered = _ignoreString->text();
    _PETRTM.setIgnoredPoints(0,true,stringEntered);
    updateAllGraphs();
}
void SimWindow::changedTau2RefAnalysis()
{
    QString stringEntered = _tau2RefAnalysis->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
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
        _PETRTM.setTau4Nominal(value);
        updateAllGraphs();
    }
    else
        _tau4Analysis->setText(stringEntered.setNum(_PETRTM.getTau4(0)));
}
void SimWindow::changedCheckBoxFitk4(bool state)
{
    _errork4Label->setVisible(state);
    _errork4->setVisible(state);
    _checkBoxk4ErrGraph->setChecked(state);
    _checkBoxk4ErrGraph->setVisible(state);
//    clearTimeCurves();
    _PETRTM.setFitk4State(state);
    _PETRTM.setPrepared(false);
    updateAllGraphs();
}
void SimWindow::changedCheckBoxChallenge(bool state)
{
    _errorChallengeLabel->setVisible(state);
    _errorChallenge->setVisible(state);
    _checkBoxChallErrGraph->setChecked(state);
    _checkBoxChallErrGraph->setVisible(state);
    _errordk2aLabel->setVisible(state);
    _errordk2a->setVisible(state);
    clearTimeCurves();

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
    _PETRTM.setPrepared(false);
    updateAllGraphs();
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
void SimWindow::changedNumberSimulationsBPnd()
{
    QString utilityString = _nSamplesBPndPerThread->text();
    bool ok;
    int value = utilityString.toInt(&ok);
    if ( ok )
        _numberSimulationsPerThread = value;
    else
        _nSamplesBPndPerThread->setText(utilityString.setNum(_numberSimulationsPerThread));
    _nSamplesBPnd->setText((utilityString.setNum(_nThreads*_numberSimulationsPerThread)));
}
void SimWindow::changedTimeLow()
{
    QString stringEntered = _timeLow->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
        _timeLowValue = value;
    else
        _timeLow->setText(stringEntered.setNum(_timeLowValue));
}
void SimWindow::changedTimeHigh()
{
    QString stringEntered = _timeHigh->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
        _timeHighValue = value;
    else
        _timeHigh->setText(stringEntered.setNum(_timeHighValue));
}
void SimWindow::changedVersusBPndGraphs()
{
    _errBPndPlot->getPlotSurface()->setVisible(_checkBoxBPndErrGraph->isChecked());
    _errChallPlot->getPlotSurface()->setVisible(_checkBoxChallErrGraph->isChecked());
    _tau2RefPlot->getPlotSurface()->setVisible(_checkBoxTau2RefGraph->isChecked());
    _errk4Plot->getPlotSurface()->setVisible(_checkBoxk4ErrGraph->isChecked());
}
void SimWindow::clearBPndCurves()
{
    QCPRange xRange;
    xRange.lower = _BPndLowValue  - _BPndStepValue;
    xRange.upper = _BPndHighValue + _BPndStepValue;

    _errBPndPlot->init();
    _errBPndPlot->setLabelXAxis("BPnd");
//    _errBPndPlot->setLabelXAxis("");
    _errBPndPlot->setLabelYAxis("BPnd % err");
    _errBPndPlot->plotDataAndFit(true);
    _errBPndPlot->setXRange(xRange);

    _errChallPlot->init();
    _errChallPlot->setLabelXAxis("BPnd");
    _errChallPlot->setLabelYAxis(QString("%1BPnd abs err").arg(QChar(0x0394)));
    _errChallPlot->plotDataAndFit(true);
    _errChallPlot->setXRange(xRange);

    _tau2RefPlot->init();
    _tau2RefPlot->setLabelXAxis("BPnd");
    _tau2RefPlot->setLabelYAxis("1/k2'");
    _tau2RefPlot->plotDataAndFit(true);
    _tau2RefPlot->setXRange(xRange);

    _errk4Plot->init();
    _errk4Plot->setLabelXAxis("BPnd");
    _errk4Plot->setLabelYAxis("1/k4");
    _errk4Plot->plotDataAndFit(true);
    _errk4Plot->setXRange(xRange);
}

void SimWindow::calculateBPndCurves()
{
    bool noisy = _simulator.getNoiseRef() != 0. || _simulator.getNoiseTar() != 0.;
    if ( noisy || _PETRTM.isForwardFitk4() )
    {
        calculateBPndCurvesInThreads();
        return;
    }

    bool calculateTau2Ref = !_PETRTM.isRTM2() || _checkBoxTau2RefGraph->isChecked();
    double saveBPnd = _simulator.getBP0(); // BPnd will be changed, so save the value for later restoration
    QString numberString;
    int nTime = _simulator.getNumberTimeBinsCoarse();
    dMatrix refRegion;      refRegion.resize(1);    refRegion[0].resize(nTime);
    dMatrix tissueVector;   tissueVector.resize(1); tissueVector[0].resize(nTime);
    dMatrix fitVector;      fitVector.resize(1);    fitVector[0].resize(nTime);

    dVector xVector, errBPnd, errChall, tau2Ref, errTau4, errBPndRTM2, errChallRTM2;
    double sigma2Sum=0.;
    double tau4Guess = _PETRTM.getTau4InRun(0);

    for (double BP0=_BPndLowValue; BP0<=_BPndHighValue; BP0 += _BPndStepValue)
    {
        xVector.append(BP0);
        _simulator.setBP0(BP0);  // set BP0 and run a series of samples

        _simulator.run();  // run the simulations with randomized noise
        // Perform the analysis
        for (int jt=0; jt<nTime; jt++)
        {
            refRegion[0][jt]    = _simulator.getCrCoarse(jt);
            tissueVector[0][jt] = _simulator.getCtCoarse(jt);
        }
        _PETRTM.setReferenceRegion(refRegion);
        _PETRTM.setTissueVector(tissueVector);
        if ( _PETRTM.isForwardFitk4() )
            _PETRTM.setTau4(0,tau4Guess);

        _PETRTM.prepare();
        _PETRTM.fitData(tissueVector,fitVector);

        if ( _PETRTM.isForwardModel() )
            sigma2Sum += _PETRTM.getSimulationSigma2(0);
        else
            sigma2Sum += _PETRTM.getSigma2();

        // update the BP error
        double guess = _PETRTM.getBP0InRun(0);
//        errBPnd.append(percentageError(guess,BP0));
        errBPnd.append(guess-BP0);
        // update the challenge error
        double truth = _simulator.getChallengeMag();
        guess = getChallengeMagFromAnalysis();
        errChall.append(guess - truth);
        // update the 1/k4 error
        if ( _PETRTM.isForwardFitk4() )
        {
            truth = _simulator.getTau4();
            guess = _PETRTM.getTau4InRun(0);
            errTau4.append(percentageError(guess,truth));
        }

        // update the tau2Ref value
        if ( !_PETRTM.isRTM2() )
        { // 3-parameter RTM, so find k2'
            guess = _PETRTM.getTau2RefInRun(0);
            tau2Ref.append(guess);
        }
        else if ( _checkBoxTau2RefGraph->isChecked() )
        {
            double saveTau2Ref = _PETRTM.getTau2RefInRun(0); // tau2Ref will be changed, so save the value for later restoration
            double bestTau2Ref = bestTau2RefForRTM2();
            tau2Ref.append(bestTau2Ref);
            _tau2RefAnalysis->setText(numberString.setNum(bestTau2Ref));  changedTau2RefAnalysis();
            // Update error vectors for optimized RTM2
            truth = _simulator.getBP0();
            guess = _PETRTM.getBP0InRun(0);
            errBPndRTM2.append(percentageError(guess,BP0));
            truth = _simulator.getChallengeMag();
            guess = getChallengeMagFromAnalysis();
            errChallRTM2.append(guess - truth);
            // restore the value of tau2Ref
            _tau2RefAnalysis->setText(numberString.setNum(saveTau2Ref));  changedTau2RefAnalysis();
        }
    } // loop over BP0 values
    sigma2Sum /= static_cast<double>(xVector.size());
    _sigma2Label->setText(QString("sigma2 = %1").arg(sigma2Sum));

    // restore the value of BPnd
    _BPnd->setText(numberString.setNum(saveBPnd));  changedBPND();

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

    _tau2RefPlot->addCurve(0,"1/k2'");
    _tau2RefPlot->setPointSize(5);
    _tau2RefPlot->setColor(colors[iColor]);
    if ( _PETRTM.isRTM2() )
        _tau2RefPlot->setPointStyle(QCPScatterStyle::ssCross);

    _errk4Plot->addCurve(0,"1/k4");
    _errk4Plot->setPointSize(5);
    _errk4Plot->setColor(colors[iColor]);

    _errBPndPlot->setData(xVector,errBPnd);
    _errChallPlot->setData(xVector,errChall);
    if ( calculateTau2Ref )
        _tau2RefPlot->setData(xVector,tau2Ref);
    if ( _PETRTM.getFitk4State() )
    {
//        if ( _PETRTM.isForwardModel() && _PETRTM.isRTM3() )
//            _errk4Plot->setData(_tau4Vector,_AICVector);
//        else
            _errk4Plot->setData(xVector,errTau4);
    }

    if ( _PETRTM.isRTM2() && _checkBoxTau2RefGraph->isChecked() )
    { // add a new curves to the two error plots
        _errBPndPlot->addCurve(0,"error_BPnd");
        _errBPndPlot->setPointSize(5);
        _errBPndPlot->setPointStyle(QCPScatterStyle::ssCross);
        _errBPndPlot->setColor(colors[iColor]);
        _errBPndPlot->setData(xVector,errBPndRTM2);
        _errChallPlot->addCurve(0,"error_Challenge");
        _errChallPlot->setPointSize(5);
        _errChallPlot->setPointStyle(QCPScatterStyle::ssCross);
        _errChallPlot->setColor(colors[iColor]);
        _errChallPlot->setData(xVector,errChallRTM2);
    }

    _errBPndPlot->conclude(0,true);
    _errChallPlot->conclude(0,true);
    _tau2RefPlot->conclude(0,true);
    _errk4Plot->conclude(0,true);

    _errBPndPlot->plotDataAndFit(true);
    _errChallPlot->plotDataAndFit(true);
    _tau2RefPlot->plotDataAndFit(true);
    _errk4Plot->plotDataAndFit(true);
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
            tau2Ref -= 0.05;
            _tau2RefAnalysis->setText(numberString.setNum(tau2Ref));  changedTau2RefAnalysis();
            estBP0 = _PETRTM.getBP0InRun(0);
            error = estBP0 - trueBP0;
        }
    }
    else
    {
        while ( error < 0. )
        {
            tau2Ref += 0.05;
            _tau2RefAnalysis->setText(numberString.setNum(tau2Ref));  changedTau2RefAnalysis();
            estBP0 = _PETRTM.getBP0InRun(0);
            error = estBP0 - trueBP0;
        }
    }
    return tau2Ref;
}

void SimWindow::clearTimeCurves()
{
    QCPRange xRange;
    xRange.lower = _timeLowValue;
    xRange.upper = _timeHighValue;

    _errBPndOrChallVsTimePlot->init();
    _errBPndOrChallVsTimePlot->setLabelXAxis("time");
    if ( _fitChallenge->isChecked() )
        _errBPndOrChallVsTimePlot->setLabelYAxis("Challenge err (abs)");
    else
        _errBPndOrChallVsTimePlot->setLabelYAxis("BPnd % err");
    _errBPndOrChallVsTimePlot->plotDataAndFit(true);
    _errBPndOrChallVsTimePlot->setXRange(xRange);
}
/*
void SimWindow::calculateTimeCurves()
{
    QMessageBox msgBox;
    msgBox.setText("Function calculateTimeCurves() currently disabled");
    msgBox.setIcon(QMessageBox::Critical);
    msgBox.exec();
}
*/
void SimWindow::calculateTimeCurves()
{
    QColor colors[10] = {Qt::black, Qt::red, Qt::blue, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};
    int nCurves = _errBPndOrChallVsTimePlot->getNumberCurves();
    int iColor  = nCurves%10;

    iVector saveTimeBinVector = _simulator.getTimeBinVectorSec();
    double saveChallengeTime = _simulator.getChallengeTime();
    dVector xVector, errVector;
    for (int jBin=0; jBin<_simulator.getNumberTimeBinsCoarse(); jBin++)
    {
        double time = _simulator.getTimeCoarse(jBin);
        if ( time >= _timeLowValue && time <= _timeHighValue)
        {
            xVector.append(time);
            if ( _fitChallenge->isChecked() )
            { // calculate the error in the challenge magnitude
                QString numberString;  numberString.setNum(time);
                _challengeTime->setText(numberString);  changedChallengeTime();
                double truth = _simulator.getChallengeMag();
                double guess = getChallengeMagFromAnalysis();
                errVector.append(guess - truth);
            }
            else
            { // calculate the percent error in BPnd
                QString numberString;  numberString.setNum(jBin);
                _numberTimeBins->setText(numberString);  changedNumberBins();

                double truth = _simulator.getBP0();
                double guess = _PETRTM.getBP0InRun(0);
                errVector.append(percentageError(guess,truth));

                _simulator.setNumberBins(saveTimeBinVector.length());
                _simulator.setDurationBins(saveTimeBinVector);
                _PETRTM.setTimeBinsSec(0,_simulator.getTimeBinVectorSec());
            }
        }
    }

    // restore the time duration and challenge time
    int nBins = saveTimeBinVector.length();
    _simulator.setNumberBins(nBins);
    _simulator.setDurationBins(saveTimeBinVector);
    QString numberString;  _numberTimeBins->setText(numberString.setNum(nBins));  changedNumberBins();

    _challengeTime->setText(numberString.setNum(saveChallengeTime));
    changedChallengeTime();  // this will update the simulation and also the challenge time in analysis

    _errBPndOrChallVsTimePlot->addCurve(0,"error_BPnd");
    _errBPndOrChallVsTimePlot->setPointSize(5);
    _errBPndOrChallVsTimePlot->setColor(colors[iColor]);
    _errBPndOrChallVsTimePlot->setData(xVector,errVector);
    _errBPndOrChallVsTimePlot->conclude(0,true);
    _errBPndOrChallVsTimePlot->plotDataAndFit(true);
}

void SimWindow::calculateBPndCurvesInThreads()
{
    updateLieDetectorProgress(10*_nThreads);

    _BP0Vector.clear();
    for (double BP0=_BPndLowValue; BP0<=_BPndHighValue; BP0 += _BPndStepValue)
        _BP0Vector.append(BP0);
    _errBPndMatrix.clear(); _errChallMatrix.clear(); _tau2RefMatrix.clear(); _errTau4Matrix.clear();

    /////////////////////////////////////////////////////////
    // Create the average volume.
    qRegisterMetaType<dVector>("dMatrix");
    QVector<lieDetector *> simSegment;
    simSegment.resize(_nThreads);
    for (int jThread=0; jThread<_nThreads; jThread++)
    {
        simSegment[jThread] = new lieDetector(_numberSimulationsPerThread, _BP0Vector, _simulator, _PETRTM);
        connect(simSegment[jThread], SIGNAL(progressLieDetector(int)), this, SLOT(updateLieDetectorProgress(int)));
        connect(simSegment[jThread], SIGNAL(finishedLieDetector(dMatrix,dMatrix,dMatrix,dMatrix,double)),
                this,SLOT(finishedLieDetectorOneThread(dMatrix,dMatrix,dMatrix,dMatrix,double)));
        QThreadPool::globalInstance()->start(simSegment[jThread]);
    }
}

void SimWindow::updateLieDetectorProgress(int iProgress)
{
    static int progressCounter=0;
    if ( iProgress > 0 )  // this gives the range
    { // initialize with nVolumes = iProgress
        progressCounter = 0;
        _progressBar->setMaximum(iProgress);
    }
    else
    {
        progressCounter++;
        _progressBar->setValue(qMin(progressCounter,_progressBar->maximum()));
    }
}

void SimWindow::finishedLieDetectorOneThread(dMatrix errBPnd, dMatrix errChall,
                                             dMatrix tau2Ref, dMatrix errTau4,
                                             double sigma2)
{
    static int nThreadsFinished=0;
    static double sigma2Sum=0.;

    if ( errBPnd.size() != _BP0Vector.size() )
        qFatal("Error: mismatched vector sizes in SimWindow::finishedLieDetectorOneThread");

    int nBP0Values = _BP0Vector.size();
    int nSamples   = errBPnd[0].size();

    QMutex mutex;
    mutex.lock();
    if ( _errBPndMatrix.size() == 0 )
    {
        _errBPndMatrix.resize(nBP0Values); _errChallMatrix.resize(nBP0Values);
        _tau2RefMatrix.resize(nBP0Values); _errTau4Matrix.resize(nBP0Values);
    }
    for (int jBP=0; jBP<nBP0Values; jBP++)
    {
        for ( int jSample=0; jSample<nSamples; jSample++)
        {
            _errBPndMatrix[jBP].append(errBPnd[jBP][jSample]);
            _errChallMatrix[jBP].append(errChall[jBP][jSample]);
            _tau2RefMatrix[jBP].append(tau2Ref[jBP][jSample]);
            _errTau4Matrix[jBP].append(errTau4[jBP][jSample]);
        }
    }
    mutex.unlock();

    nThreadsFinished++;
    sigma2Sum += sigma2;
    if ( nThreadsFinished == _nThreads )
    {
        _progressBar->reset();
        _progressBar->setMaximum(1);  // indicates that the bar is available
        sigma2Sum /= static_cast<double>(_nThreads);
        _sigma2Label->setText(QString("sigma2 = %1").arg(sigma2Sum));
        nThreadsFinished = 0;  // reset for next time
        sigma2Sum = 0.;
        finishedLieDetectorAllThreads();
    }
}

double SimWindow::calculateMean(dVector vec)
{
    double mean = 0.;
    for (int j=0; j<vec.size(); j++)
        mean += vec[j];
    mean /= static_cast<double>(vec.size());
    return mean;
}
double SimWindow::calculateStDev(double mean, dVector vec)
{
    double stdev = 0.;
    for (int j=0; j<vec.size(); j++)
        stdev += SQR(vec[j] - mean);
    stdev /= static_cast<double>(vec.size()-1);
    stdev = qSqrt(stdev);
    return stdev;
}

void SimWindow::finishedLieDetectorAllThreads()
{
    int nBP0Values = _BP0Vector.size();
    int nSamples   = _nThreads * _numberSimulationsPerThread;

    dVector errBP, errChall, tau2Ref, errTau4, errBPSEM, errChallSEM, tau2RefSEM, errTau4SEM;
    // determine the means
    for (int jBP=0; jBP<nBP0Values; jBP++)
    {
        errBP.append(calculateMean(_errBPndMatrix[jBP]));
        errChall.append(calculateMean(_errChallMatrix[jBP]));
        tau2Ref.append(calculateMean(_tau2RefMatrix[jBP]));
        errTau4.append(calculateMean(_errTau4Matrix[jBP]));

        errBPSEM.append(calculateStDev(errBP[jBP],       _errBPndMatrix[jBP])  / qSqrt(nSamples));
        errChallSEM.append(calculateStDev(errChall[jBP], _errChallMatrix[jBP]) / qSqrt(nSamples));
        tau2RefSEM.append(calculateStDev(tau2Ref[jBP],   _tau2RefMatrix[jBP])  / qSqrt(nSamples));
        errTau4SEM.append(calculateStDev(errTau4[jBP],   _errTau4Matrix[jBP])  / qSqrt(nSamples));
    }

    // restore the value of BPnd in simulator, plus regenerate graphs
    changedBPND();  // this will use the value saved in the text field of _BPnd, which should have changed during the simulation

    QColor colors[10] = {Qt::black, Qt::red, Qt::blue, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};
    int nCurves = _errBPndPlot->getNumberCurves();
    int iColor  = nCurves%10;

    _errBPndPlot->addCurve(0,"error_BPnd");
    _errBPndPlot->setPointSize(5);
    _errBPndPlot->setColor(colors[iColor]);
    _errBPndPlot->setErrorBars(1);

    _errChallPlot->addCurve(0,"error_Challenge");
    _errChallPlot->setPointSize(5);
    _errChallPlot->setColor(colors[iColor]);
    _errChallPlot->setErrorBars(1);

    _tau2RefPlot->addCurve(0,"Tau2Ref");
    _tau2RefPlot->setPointSize(5);
    _tau2RefPlot->setColor(colors[iColor]);
    _tau2RefPlot->setErrorBars(1);
    if ( _PETRTM.isRTM2() )
        _tau2RefPlot->setPointStyle(QCPScatterStyle::ssCross);

    _errk4Plot->addCurve(0,"1/k4");
    _errk4Plot->setPointSize(5);
    _errk4Plot->setColor(colors[iColor]);
    _errk4Plot->setErrorBars(1);

    dVector xVector;
    for (double BP0=_BPndLowValue; BP0<=_BPndHighValue; BP0 += _BPndStepValue)
        xVector.append(BP0);

    _errBPndPlot->setData(_BP0Vector,errBP,errBPSEM);
    _errChallPlot->setData(_BP0Vector,errChall,errChallSEM);
    if ( !_PETRTM.isRTM2() )
        _tau2RefPlot->setData(_BP0Vector,tau2Ref,tau2RefSEM);
    if ( _PETRTM.getFitk4State() )
        _errk4Plot->setData(_BP0Vector,errTau4,errTau4SEM);

    _errBPndPlot->conclude(0,true);
    _errChallPlot->conclude(0,true);
    _tau2RefPlot->conclude(0,true);
    _errk4Plot->conclude(0,true);

    _errBPndPlot->plotDataAndFit(true);
    _errChallPlot->plotDataAndFit(true);
    _tau2RefPlot->plotDataAndFit(true);
    _errk4Plot->plotDataAndFit(true);

    _progressBar->reset();
}

void Fitter::init(dVector RRTimeVector, dVector RRVector)
{
    FUNC_ENTER << RRTimeVector;
    _RRTimeVector   = RRTimeVector;
    _RRVector       = RRVector;
    // define a polynomial based upon the RR bins, which may be unevenly spaced
    _polyPlusGammaForRR.define(_nPoly, _nPoly+1,_RRTimeVector); // extra basis function is gamma-variate
    int nTime = _RRTimeVector.size();
    // Now add a gamma basis function, with variable parameters 1) onset time, 2) alpha, 3) tau
    // find starting point for tau
    double rrMaxTime = 0.;   double rrMaxValue = 0.;
    for (int jt=0; jt<nTime; jt++)
    {
        if ( _RRVector[jt] > rrMaxValue )
        {
            rrMaxValue = _RRVector[jt];
            rrMaxTime  = _RRTimeVector[jt];
        }
    }
    double tau = rrMaxTime;  double onsetTime = 0;  double alpha = 1.;
    FUNC_INFO << "set tau guess to " << tau;
    _gammaParVal.resize(_nPar);     _gammaParInc.resize(_nPar);     _gammaParAdj.resize(_nPar);
    _gammaParVal[0] = onsetTime;    _gammaParVal[1] = alpha;        _gammaParVal[2] = tau;
    _gammaParInc[0] = 0.5;          _gammaParInc[1] = alpha/10.;    _gammaParInc[2] = tau/10.;
    _gammaParAdj[0] = true;         _gammaParAdj[1] = true;         _gammaParAdj[2] = true;
}

void Fitter::setPolynomial(int nPoly)
{ // only call this this function AFTER initializing the RR
    _nPoly = nPoly;
    // define a polynomial based upon the RR bins, which may be unevenly spaced
    _polyPlusGammaForRR.define(_nPoly,_nPoly+1,_RRTimeVector);   // extra basis function is gamma-variate
    _gammaParInc[0] = 0.5;
    _gammaParInc[1] = _gammaParVal[1]/10.;
    _gammaParInc[2] = _gammaParVal[2]/10.;
}

double Fitter::computeGammaFunction(double time)
{
    // gamma(x) = x^alpha * exp*(alpha(1-x))
    double onsetTime = _gammaParVal.at(0);
    double alpha     = _gammaParVal.at(1);
    double tau       = _gammaParVal.at(2);
    double t = (time - onsetTime) / tau;
    double gammaValue=0.;
    if ( t > 0. )
        gammaValue = qPow(t,alpha) * qExp(alpha*(1.-t));
    return gammaValue;
}
double Fitter::computeGammaFunctionDerivative(double time)
{
    double onsetTime = _gammaParVal.at(0);
    double alpha     = _gammaParVal.at(1);
    double tau       = _gammaParVal.at(2);
    double t = (time - onsetTime) / tau;
    double gammaDerivative=0.;
    if ( t > 0. )
        gammaDerivative = alpha/tau * qPow(t,alpha-1.) * qExp(alpha*(1.-t)) * (1.-t);
    return gammaDerivative;
}
double Fitter::fitRRComputeCost()
{ // fit RR by GLM and compute sigma2
    // get the RR time vector (x, not y)
    dVector gammaBasis;  gammaBasis.fill(0.,_RRTimeVector.size());
    for (int jt=0; jt<_RRTimeVector.size(); jt++)
        gammaBasis[jt] = computeGammaFunction(jt);
    _polyPlusGammaForRR.addOrInsertBasisFunction(_nPoly,gammaBasis);
    // The basis functions are complete
    _polyPlusGammaForRR.fitWLS(_RRVector,true);
    double cost = getCostFunction();
    return cost;
}
void Fitter::fitTAC(double toleranceCost) // toleranceCost is a fraction (e.g., 0.01 is 1%)
{
    int iCount = 0;
    // Find the best set of parameters at this resolution level.
    bool converged = false;
    int MAXIT = 50;
    while ( !converged && iCount < MAXIT )
    {
        double costInitial = getCostFunction();
        ////////////////////////
        dVector costPar; costPar.resize(_nPar);
        dVector incPar;  incPar.resize(_nPar);
        // Test each parameter by increasing and decreasing its value by the increment.
        for (int jPar=0; jPar<_nPar; jPar++)
        {
            if ( _gammaParAdj[jPar] && _gammaParInc[jPar] != 0. )
            { // an increment of zero eliminates the parameter from the search
                // find the best value for this parameter.
                lineScan1D(jPar, costPar[jPar], incPar[jPar]);
            }
            else
                costPar[jPar] = incPar[jPar]  = 0.;
        }

        // Normalize the costPar array.
        double mag=0.;
        for (int jPar=0; jPar<_nPar; jPar++)
            mag += SQR(costPar[jPar]);
        mag = qSqrt(mag);
        // If the magnitude of the costPar vector is 0, then no values are lower than the
        // original value.
        if ( mag == 0. )
            break;
        for (int jPar=0; jPar<_nPar; jPar++)
            costPar[jPar] /= mag;

        // Calculate the projection of the incPar along the normalized costPar vector.
        for (int jPar=0; jPar<_nPar; jPar++)
        {
            incPar[jPar] = incPar[jPar] * costPar[jPar];
            if ( _gammaParAdj[jPar] ) _gammaParVal[jPar] += incPar[jPar];
        }
        qDebug() << "iteration" << iCount << "pars:" << _gammaParVal[0] << _gammaParVal[1] << _gammaParVal[2];
        qDebug() << "iteration" << iCount << "incs:" << _gammaParInc[0] << _gammaParInc[1] << _gammaParInc[2];
        qDebug() << "iteration" << iCount << "rCost:" << getCostFunction()/costInitial;
        ////////////////////////
        double costFinal = getCostFunction();
        converged = qAbs(costFinal/costInitial - 1.) < toleranceCost;
        converged &= iCount > 10;
        qDebug() << "converged?" << costInitial << costFinal << qAbs(costFinal/costInitial - 1.) << toleranceCost;
        iCount++;
    }
    qDebug() << "Fitter::fitTAC exit";
}

void Fitter::lineScan1D( int iPar, double &costRelative, double &incrementOpt )
{ // return value is true if a change was made
    double valueInitial = _gammaParVal[iPar];
    double costInitial  = fitRRComputeCost();  // refit with all initial parameter set
    // Save the initial cost function.
    double xParabola[3], yParabola[3];
    xParabola[1] = valueInitial;
    yParabola[1] = costInitial;

    // Now test + and - directions
    // test the - direction
    double value0 = valueInitial - _gammaParInc[iPar];
    _gammaParVal[iPar] = value0;
    double cost0 = fitRRComputeCost();
    xParabola[0] = value0;
    yParabola[0] = cost0;

    // test the + direction
    double value2 = valueInitial + _gammaParInc[iPar];
    _gammaParVal[iPar] = value2;
    double cost2 = fitRRComputeCost();
    xParabola[2] = value2;
    yParabola[2] = cost2;

    // allocate values that might be required if we need to keep going in one direction to find the max
    bool capturedMinOrMax;
    if ( _optHigh )
        capturedMinOrMax = costInitial > cost0 && costInitial > cost2;  // costInitial is higher than either side, so good
    else
        capturedMinOrMax = costInitial < cost0 && costInitial < cost2;  // costInitial is lower than either side, so good
    if ( capturedMinOrMax )
    { // the increment range contains the maximum/minimum of the cost function, so interpolate to get the best estimate
        double xMax, yMax;
        if ( utilMath::ParabolicInterpolation(xParabola, yParabola, xMax, yMax) )
        { // set the final increment to the interpolated value; half the increment range for next time
            if ( _optHigh )
                costRelative = yMax - costInitial;  // cost function should be growing
            else
                costRelative = costInitial - yMax;  // cost function should be shrinking
            incrementOpt = xMax - valueInitial;
            _gammaParInc[iPar] /= 2.;
        }
        else
            // This should never happen.
            costRelative = incrementOpt = 0.;
    }
    else if ( qFuzzyCompare(costInitial,cost0) || qFuzzyCompare(costInitial,cost2) )
        // not enough information
        costRelative = incrementOpt = 0.;
    else
    { // set the increment at the best edge
        if ( cost2 > cost0 )  // cost2 > costInitial > cost0
        {
            if ( _optHigh )
            {
                incrementOpt = _gammaParInc[iPar];
                costRelative = cost2 - costInitial;
            }
            else
            {
                incrementOpt = - _gammaParInc[iPar];
                costRelative = costInitial - cost0;
            }
        }
        else // cost0 > costInitial > cost2
        {
            if ( _optHigh )
            {
                incrementOpt = -_gammaParInc[iPar];
                costRelative = cost0 - costInitial;
            }
            else
            {
                incrementOpt = _gammaParInc[iPar];
                costRelative = costInitial - cost2;
            }
        }
    }
    _gammaParVal[iPar] = valueInitial;  // reset parameter value
    return;
}

void Fitter::fitGammaFunctionByGridSearch(double widthRatio)
{
    FUNC_ENTER;
    // Fit a gamma function plus a polynomial by repeated GLM using a grid search for the best gamma parameters (3)
    double onset = _gammaParVal.at(0);
    double alpha = _gammaParVal.at(1);
    double tau   = _gammaParVal.at(2);

    double tauWidth = tau / widthRatio;
    double alphaWidth = alpha / widthRatio;
    double onsetWidth = qMax(onset/widthRatio,tauWidth);
    FUNC_INFO << "tau, alpha, onset" << tauWidth << alphaWidth << onsetWidth;

    double onSetStart = onset - onsetWidth; double onSetStop = onset + onsetWidth;  double onSetStep = onsetWidth/5.;
    double tauStart   = tau - tauWidth;     double tauStop   = tau + tauWidth;      double tauStep   = tauWidth/5.;
    double alphaStart = alpha - alphaWidth; double alphaStop = alpha + alphaWidth;  double alphaStep = alphaWidth/5.;
    double bestSigma2 = fitRRComputeCost(); dVector bestGammaParVal = _gammaParVal;

    for ( onset=onSetStart; onset<onSetStop; onset+=onSetStep)
    {
        _gammaParVal[0] = onset;
        for ( alpha=alphaStart; alpha<alphaStop; alpha+=alphaStep)
        {
            _gammaParVal[1] = qMax(alpha,0.01);
            for ( tau=qMax(tauStart,0.1); tau<tauStop; tau+=tauStep)
            {
                _gammaParVal[2] = tau;
                double sigma2 = fitRRComputeCost();
                if ( sigma2 < bestSigma2 )
                {
                    bestSigma2 = sigma2;
                    bestGammaParVal = _gammaParVal;
                }
            }
        }
    }
    _gammaParVal = bestGammaParVal;
    fitRRComputeCost();
    FUNC_INFO << "gammaParVal" << _gammaParVal;
}
