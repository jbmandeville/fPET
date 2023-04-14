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
    FUNC_ENTER;
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
    createSweepTimePage();
    createSweepTau4Page();
    _tabTimeSpace->addTab(_setupPage, tr("Ref Region"));
    _tabTimeSpace->addTab(_targetPage, tr("Target Region"));
    _tabTimeSpace->addTab(_sweepBPndPage, tr("Sweep BPnd"));
    _tabTimeSpace->addTab(_sweepTimePage, tr("Sweep time"));
    _tabTimeSpace->addTab(_sweepTau4Page, tr("Sweep k4"));
    _tabTimeSpace->setTabEnabled(4,false);  // usually start in SRTM (no k4)

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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
    QCoreApplication::exit(0);
}

void SimWindow::setThreadVisibility(bool state)
{
    FUNC_ENTER;
    _threadsComboBox->setVisible(state);
    _nSamplesBPndPerThreadLabel->setVisible(state);
    _nSamplesBPndPerThread->setVisible(state);
    _nSamplesBPndLabel->setVisible(state);
    _nSamplesBPnd->setVisible(state);
}

void SimWindow::getTableDataFile()
{
    FUNC_ENTER;
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
        FUNC_INFO << 1;
        FUNC_INFO << 1.25 << _dataTable.size();
        QString errorString = utilIO::readTimeTableFile(fileName, _dataColumnNames, _dataTable);
        FUNC_INFO << 1.5;
        if ( !errorString.isEmpty() )
        {
            QMessageBox msgBox;
            msgBox.setText(errorString);
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.exec();
        }
        else
        {
            FUNC_INFO << 2;
            int iColumnBinSize=-1;  int iColumnBinTime=-1;
            if ( !defineTimeBinsFromBinSize(_validBinSizeName, iColumnBinSize) )
            {
                FUNC_INFO << 3;
                if ( !defineTimeBinsFromTimePoints(_validBinTimeName, iColumnBinTime) )
                {
                    FUNC_INFO << 4;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
    _setupPage = new QWidget();

    _plotPlasma = new plotData(0);
    _plotRR     = new plotData(1);

    //////// The setup page status bar
    _plasmaStatusBar = new QStatusBar;  // must be global so it doesn't go out of scope
    _plasmaStatusBar->setStyleSheet("color:blue");
    _plotPlasma->setQCStatusBar(_plasmaStatusBar);
    _RRStatusBar = new QStatusBar;  // must be global so it doesn't go out of scope
    _RRStatusBar->setStyleSheet("color:blue");
    _plotRR->setQCStatusBar(_RRStatusBar);
//    _plotRR->setMainPage(this);

    auto *setupPlotLayout = new QVBoxLayout();
    setupPlotLayout->addWidget(_plotPlasma->getPlotSurface());
    setupPlotLayout->addWidget(_plasmaStatusBar);
    setupPlotLayout->addWidget(_plotRR->getPlotSurface());
    setupPlotLayout->addWidget(_RRStatusBar);
    setupPlotLayout->setStretch(0,10);
    setupPlotLayout->setStretch(1,1);
    setupPlotLayout->setStretch(2,10);
    setupPlotLayout->setStretch(3,1);
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
    QLabel *delayLabel         = new QLabel("'Kbol' delay (min)");
    QLabel *fastTauLabel       = new QLabel("Fast elimination (tau)");
    QLabel *slowTauLabel       = new QLabel("Slow elimination (tau)");
    QLabel *fastFractionLabel  = new QLabel("Fast fraction (<=1.)");
    QLabel *calcLabel       = new QLabel("Match to data");
    _bolusMag     = new QLineEdit();
    _tauDecay     = new QLineEdit();
    _KBol         = new QLineEdit();
    _KBolDelay    = new QLineEdit();
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
    _KBolDelay->setText(numberString.setNum(_simulator.getKBolDelay()));
    _fastTau->setText(numberString.setNum(_simulator.getTauFastElim()));
    _slowTau->setText(numberString.setNum(_simulator.getTauSlowElim()));
    _fastFraction->setText(numberString.setNum(_simulator.getFastFraction()));
    _bolusMag->setFixedWidth(editTextSize);
    _tauDecay->setFixedWidth(editTextSize);
    _KBol->setFixedWidth(editTextSize);
    _KBolDelay->setFixedWidth(editTextSize);
    _fastTau->setFixedWidth(editTextSize);
    _slowTau->setFixedWidth(editTextSize);
    _fastFraction->setFixedWidth(editTextSize);
    _KBol->setToolTip("The time at which the integral of infusion matches the bolus;\n'Kbol' is a misnomer: this is really 'Tbol';\nSet to 0 for no BI");
    _KBolDelay->setToolTip("The time at which the infusion starts");
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
    plasmaInLayout->addWidget(delayLabel,3,0);
    plasmaInLayout->addWidget(_KBolDelay,3,1);

    plasmaInLayout->addWidget(fastTauLabel,4,0);
    plasmaInLayout->addWidget(_fastTau,4,1);
    plasmaInLayout->addWidget(_fastTauCheckBox,4,2);
    plasmaInLayout->addWidget(slowTauLabel,5,0);
    plasmaInLayout->addWidget(_slowTau,5,1);
    plasmaInLayout->addWidget(_slowTauCheckBox,5,2);
    plasmaInLayout->addWidget(fastFractionLabel,6,0);
    plasmaInLayout->addWidget(_fastFraction,6,1);
    plasmaInLayout->addWidget(_fastFractionCheckBox,6,2);
    plasmaInLayout->addWidget(calcLabel,7,0);
    plasmaInLayout->addWidget(_calcRRMatch,7,1);
    plasmaInGroupBox->setLayout(plasmaInLayout);
    plasmaInLayout->setSpacing(0);
    connect(_bolusMag, SIGNAL(editingFinished()), this, SLOT(changedBolusMag()));
    connect(_tauDecay, SIGNAL(editingFinished()), this, SLOT(changedTauBolus()));
    connect(_KBol, SIGNAL(editingFinished()), this, SLOT(changedInfusion()));
    connect(_KBolDelay, SIGNAL(editingFinished()), this, SLOT(changedInfusionDelay()));
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
    fullLayout->addLayout(setupPlotLayout);
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
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotPlasma, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotPlasma, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotPlasma, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotPlasma, SLOT(setSelectPoints(bool)));
//    connect(crossCursorAct,  SIGNAL(triggered(bool)), _plotPlasma, SLOT(setSelectPoints()));

    connect(dragXAction,     SIGNAL(triggered(bool)), _plotRR, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotRR, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotRR, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotRR, SLOT(setSelectPoints(bool)));
//    connect(crossCursorAct,  SIGNAL(triggered(bool)), _plotRR, SLOT(setSelectPoints()));

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

    _tau1Ref->setEnabled( !simStartsFromDataFit() );
    _plasmaFracRef->setEnabled( !simStartsFromDataFit() );
    _noiseRef->setEnabled( !simStartsFromDataFit() );

    updateAllGraphs();
}

void SimWindow::enableComboBoxItem(QComboBox *comboBox, int itemNumber, bool enable)
{
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
    _PETRTM.setWeightingModel(indexInBox);
    _PETRTM.setPrepared(false);
    updateAllGraphs();
}

void SimWindow::changedModelType(int indexInBox)
{
    FUNC_ENTER;
    _PETRTM.setRTMModelType(indexInBox);

    bool forward = _PETRTM.isForwardModel();
    _tabTimeSpace->setTabEnabled(4,forward);  // only enable with k4 (not SRTM)
    _radioShowAIC->setVisible(forward);

    _PETRTM.setPrepared(false);
    if ( _PETRTM.isSRTM() ) updateAllGraphs();

    _fit2ParRadioButton->setVisible(forward);
    _fitk4RadioButton->setVisible(forward);
    _fitk2RefRadioButton->setVisible(forward);
    _fitk4k2RefRadioButton->setVisible(forward);
    _fitDVRadioButton->setVisible(forward);
    _runSimulationAndAnalysis->setVisible(forward);

    _tau2RefAnalysisLabel->setVisible(_PETRTM.isRTM2());
    _tau2RefAnalysis->setVisible(_PETRTM.isRTM2());
    _errorR1Label->setVisible(_PETRTM.isRTM3() || forward);
    _errorR1->setVisible(_PETRTM.isRTM3() || forward);
    _errorDVLabel->setVisible(forward && _PETRTM.getRTMModel() == RTM_fFRTM3DV);
    _errorDV->setVisible(forward && _PETRTM.getRTMModel() == RTM_fFRTM3DV);
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
    _tau4AnalysisLabel->setVisible( forward );
    _tau4Analysis->setVisible(      forward );
    _DVAnalysisLabel->setVisible(   forward );
    _DVAnalysis->setVisible(        forward );
}

void SimWindow::createTargetPage()
{
    FUNC_ENTER;
    _targetPage = new QWidget();

    _plotBasis  = new plotData(2);
    _plotTarget = new plotData(3);
    _plotAIC    = new plotData(4);

    //////// The setup page status bar
    _TRStatusBar = new QStatusBar;  // must be global so it doesn't go out of scope
    _TRStatusBar->setStyleSheet("color:blue");
    _plotTarget->setQCStatusBar(_TRStatusBar);
    _plotAIC->setQCStatusBar(_TRStatusBar);
//    _plotTarget->setMainPage(this);

    auto *plotTargetLayout = new QVBoxLayout();
    plotTargetLayout->addWidget(_plotBasis->getPlotSurface());
    plotTargetLayout->addWidget(_plotTarget->getPlotSurface());
    plotTargetLayout->addWidget(_plotAIC->getPlotSurface());
    plotTargetLayout->addWidget(_TRStatusBar);
    plotTargetLayout->setStretch(0,1);
    plotTargetLayout->setStretch(1,1);
    plotTargetLayout->setStretch(2,1);
//    plotTargetLayout->setStretch(2,1);

    QString numberString;
    int editTextSize=80;

    // target region: input
    _targetSimulationGroupBox = new QGroupBox("Target Simulation: truth");
    QLabel *BPndLabel    = new QLabel("BPnd");
    QLabel *R1Label      = new QLabel("R1");
    QLabel *tau4Label    = new QLabel("1/k4 (min)");
    QLabel *DVLabel      = new QLabel("DV");
    QLabel *challengeTimeLabel = new QLabel("Challenge time");
    QChar delta = QChar(0x0394);
    QLabel *plasmaFraclabel    = new QLabel("Plasma %");
    QLabel *noiseLabel = new QLabel("Noise");
    _challengeMagLabel = new QLabel(QString("%1BPnd in % (abs = %2)").arg(delta).arg(_simulator.getChallengeMagAbs()));
    _BPnd          = new QLineEdit();
    _R1            = new QLineEdit();
    _tau4          = new QLineEdit();
    _DV            = new QLineEdit();
    _challengeTime = new QLineEdit();
    _challengeMagPercent = new QLineEdit();
    _plasmaFracTar = new QLineEdit();
    _noiseTar      = new QLineEdit();
    _BPnd->setText(numberString.setNum(_simulator.getBP0()));
    _R1->setText(numberString.setNum(_simulator.getR1()));
    _tau4->setText(numberString.setNum(_simulator.getTau4()));
    _DV->setText(numberString.setNum(_simulator.getDV()));
    _challengeTime->setText(numberString.setNum(_simulator.getChallengeTime()));
    _challengeMagPercent->setText(numberString.setNum(_simulator.getChallengeMagPercent()));
    _noiseTar->setText(numberString.setNum(_simulator.getNoiseTar()));
    _plasmaFracTar->setText(numberString.setNum(_simulator.getPlasmaPercentTar()));
    _BPnd->setFixedWidth(editTextSize);
    _R1->setFixedWidth(editTextSize);
    _tau4->setFixedWidth(editTextSize);
    _DV->setFixedWidth(editTextSize);
    _challengeTime->setFixedWidth(editTextSize);
    _challengeMagPercent->setFixedWidth(editTextSize);
    _noiseTar->setFixedWidth(editTextSize);
    _plasmaFracTar->setFixedWidth(editTextSize);

    _runSimulationAndAnalysis = new QPushButton();
    QPixmap pixmapCalculate(":/My-Icons/calculator.png");
    QIcon calculatorIcon(pixmapCalculate);
    _runSimulationAndAnalysis->setIcon(calculatorIcon);
    connect(_runSimulationAndAnalysis, SIGNAL(pressed()), this, SLOT(runSimulationAndAnalysis()));

    auto *targetSimulationLayout = new QGridLayout();
    targetSimulationLayout->addWidget(BPndLabel,0,0);
    targetSimulationLayout->addWidget(_BPnd,0,1);
    targetSimulationLayout->addWidget(R1Label,1,0);
    targetSimulationLayout->addWidget(_R1,1,1);
    targetSimulationLayout->addWidget(tau4Label,2,0);
    targetSimulationLayout->addWidget(_tau4,2,1);
    targetSimulationLayout->addWidget(DVLabel,3,0);
    targetSimulationLayout->addWidget(_DV,3,1);
    targetSimulationLayout->addWidget(challengeTimeLabel,4,0);
    targetSimulationLayout->addWidget(_challengeTime,4,1);
    targetSimulationLayout->addWidget(_challengeMagLabel,5,0);
    targetSimulationLayout->addWidget(_challengeMagPercent,5,1);
    targetSimulationLayout->addWidget(plasmaFraclabel,6,0);
    targetSimulationLayout->addWidget(_plasmaFracTar,6,1);
    targetSimulationLayout->addWidget(noiseLabel,7,0);
    targetSimulationLayout->addWidget(_noiseTar,7,1);
    targetSimulationLayout->addWidget(_runSimulationAndAnalysis,8,1);

    _targetSimulationGroupBox->setLayout(targetSimulationLayout);
    targetSimulationLayout->setSpacing(0);
    connect(_BPnd,          SIGNAL(editingFinished()), this, SLOT(changedBPND()));
    connect(_R1,            SIGNAL(editingFinished()), this, SLOT(changedR1()));
    connect(_tau4,          SIGNAL(editingFinished()), this, SLOT(changedTau4()));
    connect(_DV,            SIGNAL(editingFinished()), this, SLOT(changedDV()));
    connect(_challengeTime, SIGNAL(editingFinished()), this, SLOT(changedChallengeTime()));
    connect(_challengeMagPercent,  SIGNAL(editingFinished()), this, SLOT(changedChallengeMag()));
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
    _DVAnalysisLabel        = new QLabel("Analysis DV");
    _modelType = new QComboBox();
    _modelType->addItem("SRTM3");
    _modelType->addItem("SRTM2");
    _modelType->addItem("SRTM2Fit");
    _modelType->addItem("fFRTM");
    _modelType->addItem("fFRTM+k4");
    _modelType->addItem("fFRTM+DV");
    _modelType->addItem("fFRTM+k2Ref");
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
    _tau2RefAnalysis  = new QLineEdit();
    _tau4Analysis     = new QLineEdit();
    _DVAnalysis       = new QLineEdit();
    _ignoreString->setText("");
    _tau2RefAnalysis->setText(numberString.setNum(_simulator.getTau2Ref()));
    _tau4Analysis->setText(numberString.setNum(_simulator.getTau4()));
    _DVAnalysis->setText(numberString.setNum(_simulator.getDV()));
    _tau2RefAnalysis->setFixedWidth(editTextSize);
    _tau4Analysis->setFixedWidth(editTextSize);
    _DVAnalysis->setFixedWidth(editTextSize);
    _tau2RefAnalysisLabel->setVisible(false);
    _tau4AnalysisLabel->setVisible(false);
    _DVAnalysisLabel->setVisible(false);
    _tau2RefAnalysis->setVisible(false);
    _tau4Analysis->setVisible(false);
    _DVAnalysis->setVisible(false);
    connect(_modelType,  SIGNAL(currentIndexChanged(int)), this, SLOT(changedModelType(int)));
    connect(_weightType, SIGNAL(currentIndexChanged(int)), this, SLOT(changedWeightType(int)));
    _fitChallenge = new QCheckBox("fit challenge?");
    _fitChallenge->setChecked(false);

    _fit2ParRadioButton = new QRadioButton("2 pars");
    _fitk4RadioButton = new QRadioButton("k4");
    _fitk2RefRadioButton = new QRadioButton("k2Ref");
    _fitk4k2RefRadioButton = new QRadioButton("k4/k2Ref");
    _fitDVRadioButton = new QRadioButton("DV");

    _fit2ParRadioButton->setVisible(false);
    _fitk4RadioButton->setVisible(false);
    _fitk2RefRadioButton->setVisible(false);
    _fitk4k2RefRadioButton->setVisible(false);
    _fitDVRadioButton->setVisible(false);

    auto *buttonGroup = new QButtonGroup();
    buttonGroup->addButton(_fit2ParRadioButton);
    buttonGroup->addButton(_fitk4RadioButton);
    buttonGroup->addButton(_fitk2RefRadioButton);
    buttonGroup->addButton(_fitk4k2RefRadioButton);
    buttonGroup->addButton(_fitDVRadioButton);
    _fit2ParRadioButton->setChecked(true);
    _fitk4RadioButton->setChecked(false);
    _fitk2RefRadioButton->setChecked(false);
    _fitk4k2RefRadioButton->setChecked(false);
    _fitDVRadioButton->setChecked(false);

    _parSearchStart = new QLineEdit("");
    _parSearchEnd   = new QLineEdit("");
    _parSearchStep  = new QLineEdit("");
    _parSearchStart->setVisible(false);
    _parSearchEnd->setVisible(false);
    _parSearchStep->setVisible(false);
    connect(_parSearchStart,  SIGNAL(editingFinished()), this, SLOT(changedParSearchStart()));
    connect(_parSearchEnd,    SIGNAL(editingFinished()), this, SLOT(changedParSearchEnd()));
    connect(_parSearchStep,   SIGNAL(editingFinished()), this, SLOT(changedParSearchStep()));

    auto *analysisLayout = new QGridLayout();
    analysisLayout->addWidget(modelLabel,0,0);
    analysisLayout->addWidget(_modelType,0,1);
    analysisLayout->addWidget(_fitChallenge,0,2);

    analysisLayout->addWidget(weightLabel,1,0);
    analysisLayout->addWidget(_weightType,1,1);
    analysisLayout->addWidget(_fit2ParRadioButton,1,2);
    analysisLayout->addWidget(ignoreLabel,2,0);
    analysisLayout->addWidget(_ignoreString,2,1);
    analysisLayout->addWidget(_fitk4k2RefRadioButton,2,2);
    analysisLayout->addWidget(_tau2RefAnalysisLabel,3,0);
    analysisLayout->addWidget(_tau2RefAnalysis,3,1);
    analysisLayout->addWidget(_fitk2RefRadioButton,3,2);
    analysisLayout->addWidget(_tau4AnalysisLabel,4,0);
    analysisLayout->addWidget(_tau4Analysis,4,1);
    analysisLayout->addWidget(_fitk4RadioButton,4,2);
    analysisLayout->addWidget(_DVAnalysisLabel,5,0);
    analysisLayout->addWidget(_DVAnalysis,5,1);
    analysisLayout->addWidget(_fitDVRadioButton,5,2);

    analysisLayout->addWidget(_parSearchStart,6,0);
    analysisLayout->addWidget(_parSearchEnd,6,1);
    analysisLayout->addWidget(_parSearchStep,6,2);

    _analysisGroupBox->setLayout(analysisLayout);
    analysisLayout->setSpacing(0);
    connect(_ignoreString,     SIGNAL(editingFinished()), this, SLOT(changedIgnoreString()));
    connect(_tau2RefAnalysis,  SIGNAL(editingFinished()), this, SLOT(changedTau2RefAnalysis()));
    connect(_tau4Analysis,     SIGNAL(editingFinished()), this, SLOT(changedTau4Analysis()));
    connect(_DVAnalysis,       SIGNAL(editingFinished()), this, SLOT(changedDVAnalysis()));
    connect(_fitChallenge,     SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxChallenge(bool)));
    connect(_fitk4RadioButton,    SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxFitk4(bool)));
    connect(_fitk2RefRadioButton, SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxFitk2Ref(bool)));
    connect(_fitDVRadioButton,    SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxFitDV(bool)));
    // target region: errors
    auto *errorGroupBox = new QGroupBox("Target Region: Percentage Errors");
    QLabel *errorBPndLabel  = new QLabel("BPnd ");
    QLabel *errork2Label    = new QLabel("k2   ");
    _errorR1Label           = new QLabel("R1   ");
    _errorDVLabel           = new QLabel("DV   ");
    _errorTau2RefLabel      = new QLabel("1/k2'");
    _errorChallengeLabel    = new QLabel("Challenge");
    _errork4Label           = new QLabel("1/k4");
    _sigmaLabel             = new QLabel("sigma");
    _errorBPnd      = new QLabel();
    _errork2        = new QLabel();
    _errorR1        = new QLabel();
    _errorDV        = new QLabel();
    _errorTau2Ref   = new QLabel();
    _errorChallenge = new QLabel();
    _errork4        = new QLabel();
    _sigma          = new QLabel();
    _errorChallengeLabel->setVisible(false);
    _errorChallenge->setVisible(false);
    _errork4Label->setVisible(false);
    _errork4->setVisible(false);
    auto *errorLayout = new QGridLayout();
    errorLayout->addWidget(errorBPndLabel,0,0);
    errorLayout->addWidget(_errorBPnd,0,1);
    errorLayout->addWidget(errork2Label,1,0);
    errorLayout->addWidget(_errork2,1,1);
    errorLayout->addWidget(_errorR1Label,2,0);
    errorLayout->addWidget(_errorR1,2,1);
    errorLayout->addWidget(_errorTau2RefLabel,3,0);
    errorLayout->addWidget(_errorTau2Ref,3,1);
    errorLayout->addWidget(_errorDVLabel,4,0);
    errorLayout->addWidget(_errorDV,4,1);
    errorLayout->addWidget(_errorChallengeLabel,5,0);
    errorLayout->addWidget(_errorChallenge,5,1);
    errorLayout->addWidget(_errork4Label,6,0);
    errorLayout->addWidget(_errork4,6,1);
    errorLayout->addWidget(_sigmaLabel,7,0);
    errorLayout->addWidget(_sigma,7,1);
    errorGroupBox->setLayout(errorLayout);
    errorLayout->setSpacing(0);

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(_targetSimulationGroupBox);
    rightLayout->addWidget(_targetDataGroupBox);
    rightLayout->addWidget(_analysisGroupBox);
    rightLayout->addWidget(errorGroupBox);
    rightLayout->setSpacing(0);

    QHBoxLayout *fullLayout = new QHBoxLayout();
    fullLayout->addLayout(plotTargetLayout);
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
    auto *radioShowBasis       = new QRadioButton("Basis",_targetPage);
    auto *radioShowTarget      = new QRadioButton("Target",_targetPage);
    _radioShowAIC              = new QRadioButton("AIC",_targetPage);
    _radioShowAIC->setVisible(false);
    QButtonGroup *showGroup = new QButtonGroup(_targetPage);
    showGroup->addButton(radioShowBasis);
    showGroup->addButton(radioShowTarget);
    showGroup->addButton(_radioShowAIC);
    radioShowTarget->setChecked(true);
    _plotBasis->getPlotSurface()->setVisible(false);
    _plotAIC->getPlotSurface()->setVisible(false);
    _plotTarget->getPlotSurface()->setVisible(true);

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
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotBasis,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotBasis,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotBasis,  SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotTarget, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotTarget, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotTarget, SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotAIC, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotAIC, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotAIC, SLOT(autoScale(bool)));
    connect(radioShowBasis,  SIGNAL(clicked(bool)), this, SLOT(showBasis()));
    connect(radioShowTarget, SIGNAL(clicked(bool)), this, SLOT(showTarget()));
    connect(_radioShowAIC,   SIGNAL(clicked(bool)), this, SLOT(showAICPlot()));
    connect(_analyzeSimulation, SIGNAL(clicked(bool)), this, SLOT(clickedAnalyzeStimulation(bool)));
    connect(_analyzeRealData,   SIGNAL(clicked(bool)), this, SLOT(clickedAnalyzeRealData(bool)));

    _targetPage->setLayout(fullLayout);
}

void SimWindow::clickedAnalyzeStimulation(bool state)
{
    FUNC_ENTER;
    _calculateBPndCurves->setEnabled(state);
    _calculateTimeCurves->setEnabled(state);
    _calculateTau4Curves->setEnabled(state);
    updateAllGraphs();
}
void SimWindow::clickedAnalyzeRealData(bool state)
{
    FUNC_ENTER;
    _calculateBPndCurves->setEnabled(!state);
    _calculateTimeCurves->setEnabled(!state);
    _calculateTau4Curves->setEnabled(!state);
    updateAllGraphs();
}

void SimWindow::createSweepBPndPage()
{
    FUNC_ENTER;
    // For RTM3:
    // 1) BPnd_err vs. BPnd
    // 2) Challenge error vs. BPnd
    // 3) Tau2Ref  vs. BPnd

    // For RTM2:
    // 1) BPnd_err vs. BPnd
    // 2) Challenge error vs. BPnd
    // 3) Tau2Ref (no bias) vs. BPnd (requires root finding)

    _sweepBPndPage = new QWidget();

    _plotErrBPndVsBPnd  = new plotData(4);
    _plotErrChallVsBPnd = new plotData(5);
    _plotTau2RefVsBPnd  = new plotData(6);
    _plotErrk4VsBPnd    = new plotData(7);
    _plotErrDVVsBPnd    = new plotData(8);
    auto *plotLayout = new QVBoxLayout();
    plotLayout->addWidget(_plotErrBPndVsBPnd->getPlotSurface());
    plotLayout->addWidget(_plotErrChallVsBPnd->getPlotSurface());
    plotLayout->addWidget(_plotTau2RefVsBPnd->getPlotSurface());
    plotLayout->addWidget(_plotErrk4VsBPnd->getPlotSurface());
    plotLayout->addWidget(_plotErrDVVsBPnd->getPlotSurface());
    QString numberString;
    int editTextSize=80;

    auto *BPndGroupBox    = new QGroupBox("BPnd graph range and step");
    auto *BPndLowLabel  = new QLabel("BPnd low");
    auto *BPndHighLabel = new QLabel("BPnd high");
    auto *BPndStepLabel = new QLabel("BPnd step");
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

    _checkBoxBPndErrVsBPnd  = new QCheckBox("BPnd error");
    _checkBoxChallErrVsBPnd = new QCheckBox("Challenge error");
    _checkBoxTau2RefGraph  = new QCheckBox("1/k2' (=R1/k2)");
    _checkBoxk4ErrGraph    = new QCheckBox("1/k4");
    _checkBoxDVErrGraph    = new QCheckBox("DV");
    _sigma2Label = new QLabel("sigma2 = 0");
    _checkBoxBPndErrVsBPnd->setChecked(true);
    _checkBoxChallErrVsBPnd->setChecked(false);
    _checkBoxChallErrVsBPnd->setVisible(false);
    _checkBoxk4ErrGraph->setVisible(false);
    _checkBoxDVErrGraph->setVisible(false);
    _checkBoxTau2RefGraph->setChecked(true);
    auto *checkGroupBox = new QGroupBox("Show or hide graphs");
    auto *checkLayout = new QVBoxLayout();
    checkLayout->addWidget(_checkBoxBPndErrVsBPnd);
    checkLayout->addWidget(_checkBoxChallErrVsBPnd);
    checkLayout->addWidget(_checkBoxTau2RefGraph);
    checkLayout->addWidget(_checkBoxk4ErrGraph);
    checkLayout->addWidget(_checkBoxDVErrGraph);
    checkLayout->addWidget(_sigma2Label);
    checkGroupBox->setLayout(checkLayout);
    connect(_checkBoxBPndErrVsBPnd,  SIGNAL(toggled(bool)), this, SLOT(changedVersusBPndGraphs()));
    connect(_checkBoxChallErrVsBPnd, SIGNAL(toggled(bool)), this, SLOT(changedVersusBPndGraphs()));
    connect(_checkBoxTau2RefGraph,  SIGNAL(toggled(bool)), this, SLOT(changedVersusBPndGraphs()));
    connect(_checkBoxk4ErrGraph,    SIGNAL(toggled(bool)), this, SLOT(changedVersusBPndGraphs()));
    connect(_checkBoxDVErrGraph,    SIGNAL(toggled(bool)), this, SLOT(changedVersusBPndGraphs()));

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
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotErrBPndVsBPnd,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotErrBPndVsBPnd,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotErrBPndVsBPnd,  SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotErrChallVsBPnd, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotErrChallVsBPnd, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotErrChallVsBPnd, SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotTau2RefVsBPnd,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotTau2RefVsBPnd,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotTau2RefVsBPnd,  SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotErrk4VsBPnd,    SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotErrk4VsBPnd,    SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotErrk4VsBPnd,    SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotErrDVVsBPnd,    SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotErrDVVsBPnd,    SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotErrDVVsBPnd,    SLOT(autoScale(bool)));

    clearBPndCurves();
    changedVersusBPndGraphs();

    _sweepBPndPage->setLayout(fullLayout);
}

void SimWindow::setReferenceRegionForAnalysis()
{ // xxx this needs modification
    FUNC_ENTER;
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

    if ( _PETRTM.isForwardModel() )
        _PETRTM.setSmoothingScale(0.5);
    else
        _PETRTM.setSmoothingScale(0.01);
}

void SimWindow::createSweepTau4Page()
{
    FUNC_ENTER;
    // 1) BPnd_err vs. time
    // 2) Challenge error vs. time

    _sweepTau4Page = new QWidget();
//    _sweepTau4Page->setVisible(false);  // start in SRTM

    _plotErrBPndVsTau4  = new plotData(8);
    _plotErrChallVsTau4 = new plotData(9);
    _plotAICVsTau4      = new plotData(10);
    auto *plotLayout = new QVBoxLayout();
    plotLayout->addWidget(_plotErrBPndVsTau4->getPlotSurface());
    plotLayout->addWidget(_plotErrChallVsTau4->getPlotSurface());
    plotLayout->addWidget(_plotAICVsTau4->getPlotSurface());
    QString numberString;

    int editTextSize=80;
    auto *tau4GroupBox    = new QGroupBox("1/k4 graph range and step");
    auto *tau4HighLabel = new QLabel("1/k4 high (min)");
    auto *tau4StepLabel = new QLabel("1/k4 step (min)");
    _tau4High         = new QLineEdit();
    _tau4Step         = new QLineEdit();
    _tau4High->setFixedWidth(editTextSize);
    _tau4Step->setFixedWidth(editTextSize);
    _tau4High->setText(numberString.setNum(_tau4HighValue));
    _tau4Step->setText(numberString.setNum(_tau4StepValue));
    connect(_tau4High,  SIGNAL(editingFinished()), this, SLOT(changedTau4High()));
    connect(_tau4Step,  SIGNAL(editingFinished()), this, SLOT(changedTau4Step()));

    auto *tau4Layout = new QGridLayout();
    tau4Layout->addWidget(tau4HighLabel,1,0);
    tau4Layout->addWidget(_tau4High,1,1);
    tau4Layout->addWidget(tau4StepLabel,2,0);
    tau4Layout->addWidget(_tau4Step,2,1);
    tau4GroupBox->setLayout(tau4Layout);
    tau4Layout->setSpacing(0);

    auto *calcGroupBox = new QGroupBox("Calculate or clear curves");
    QLabel *calcLabel  = new QLabel("Calculate curves");
    QLabel *clearLabel = new QLabel("Clear curves");
    _calculateTau4Curves   = new QPushButton();
    QPixmap pixmapCalculate(":/My-Icons/calculator.png");
    QIcon calculatorIcon(pixmapCalculate);
    _calculateTau4Curves->setIcon(calculatorIcon);
    _clearTau4Curves  = new QPushButton();
    QPixmap eraser(":/My-Icons/eraser.png");
    QIcon eraserIcon(eraser);
    _clearTau4Curves->setIcon(eraserIcon);
    auto *calcLayout = new QGridLayout();
    calcLayout->addWidget(calcLabel,0,0);
    calcLayout->addWidget(_calculateTau4Curves,0,1);
    calcLayout->addWidget(clearLabel,1,0);
    calcLayout->addWidget(_clearTau4Curves,1,1);
    calcGroupBox->setLayout(calcLayout);
    calcLayout->setSpacing(0);
    connect(_calculateTau4Curves, SIGNAL(pressed()), this, SLOT(calculateTau4Curves()));
    connect(_clearTau4Curves, SIGNAL(pressed()),     this, SLOT(clearTau4Curves()));

    _checkBoxBPndErrVsTau4  = new QCheckBox("BPnd error");
    _checkBoxChallErrVsTau4 = new QCheckBox("Challenge error");
    _checkBoxAICVsTau4      = new QCheckBox("AIC");
    _checkBoxBPndErrVsTau4->setChecked(true);
    _checkBoxChallErrVsTau4->setChecked(false);
    _checkBoxAICVsTau4->setChecked(true);
    auto *checkGroupBox = new QGroupBox("Show or hide graphs");
    auto *checkLayout = new QVBoxLayout();
    checkLayout->addWidget(_checkBoxBPndErrVsTau4);
    checkLayout->addWidget(_checkBoxChallErrVsTau4);
    checkLayout->addWidget(_checkBoxAICVsTau4);
    checkGroupBox->setLayout(checkLayout);
    connect(_checkBoxBPndErrVsTau4,  SIGNAL(toggled(bool)), this, SLOT(changedVersusTau4Graphs()));
    connect(_checkBoxChallErrVsTau4, SIGNAL(toggled(bool)), this, SLOT(changedVersusTau4Graphs()));
    connect(_checkBoxAICVsTau4,      SIGNAL(toggled(bool)), this, SLOT(changedVersusTau4Graphs()));

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(tau4GroupBox);
    rightLayout->addWidget(calcGroupBox);
    rightLayout->addWidget(checkGroupBox);

    QHBoxLayout *fullLayout = new QHBoxLayout();
    fullLayout->addLayout(plotLayout);
    fullLayout->addLayout(rightLayout);
    fullLayout->setStretch(0,100);
    fullLayout->setStretch(1,1);

    _sweepTau4Page->setLayout(fullLayout);

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
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotErrBPndVsTau4,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotErrBPndVsTau4,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotErrBPndVsTau4,  SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotErrChallVsTau4, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotErrChallVsTau4, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotErrChallVsTau4, SLOT(autoScale(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotAICVsTau4,      SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotAICVsTau4,      SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotAICVsTau4,      SLOT(autoScale(bool)));

    clearTau4Curves();
    changedVersusTau4Graphs();
}

void SimWindow::createSweepTimePage()
{
    FUNC_ENTER;
    // 1) BPnd_err vs. time
    // 2) Challenge error vs. time

    _sweepTimePage = new QWidget();

    _plotErrBPndOrChallVsTime = new plotData(6);
    auto *plotLayout = new QVBoxLayout();
    plotLayout->addWidget(_plotErrBPndOrChallVsTime->getPlotSurface());
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
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotErrBPndOrChallVsTime,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotErrBPndOrChallVsTime,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotErrBPndOrChallVsTime,  SLOT(autoScale(bool)));

    clearTimeCurves();

    _sweepTimePage->setLayout(fullLayout);

}

void SimWindow::showPlasma()
{
    FUNC_ENTER;
    _plotPlasma->getPlotSurface()->setVisible(true);
    _plasmaStatusBar->setVisible(true);
    _plotRR->getPlotSurface()->setVisible(false);
    _RRStatusBar->setVisible(false);
    _whichPlasmaPlot->setVisible(true);
    _clearPlasmaPlot->setVisible(true);
}
void SimWindow::showRR()
{
    FUNC_ENTER;
    _plotPlasma->getPlotSurface()->setVisible(false);
    _plasmaStatusBar->setVisible(false);
    _plotRR->getPlotSurface()->setVisible(true);
    _RRStatusBar->setVisible(true);
    _whichPlasmaPlot->setVisible(false);
    _clearPlasmaPlot->setVisible(false);
}
void SimWindow::showBasis()
{
    FUNC_ENTER;
    _plotBasis->getPlotSurface()->setVisible(true);
    _plotTarget->getPlotSurface()->setVisible(false);
    _plotAIC->getPlotSurface()->setVisible(false);
}
void SimWindow::showTarget()
{
    qInfo() << "showTarget";
    FUNC_ENTER;
    _plotBasis->getPlotSurface()->setVisible(false);
    _plotTarget->getPlotSurface()->setVisible(true);
    _plotAIC->getPlotSurface()->setVisible(false);
//    updateAllGraphs();
}
void SimWindow::showAICPlot()
{
    qInfo() << "showAICPlot";
    FUNC_ENTER;
    _plotBasis->getPlotSurface()->setVisible(false);
    _plotTarget->getPlotSurface()->setVisible(false);
    _plotAIC->getPlotSurface()->setVisible(true);
    _runSimulationAndAnalysis->setVisible(true);
}

void SimWindow::updatePlasmaGraph()
{
    FUNC_ENTER;
    static int whichColor=0;
    QColor colors[10] = {Qt::black, Qt::blue, Qt::red, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};

    int index = _whichPlasmaPlot->currentIndex();
    bool fineScale = (index == 2) || (index == 3) || (index == 5) || (index == 7) || (index == 8) || (index==9);

    // update the plot: RR
    if ( _clearPlasmaPlot->isChecked() )
        _plotPlasma->init();
    _plotPlasma->setLegendOn(true);
    if ( index != 0 )
        _plotPlasma->addCurve(0,_whichPlasmaPlot->currentText());

    if ( index == 0 )
    { // Cr coarse
        if ( analyzeRealData() )
            addDataCurveRR(_plotPlasma);
        else
            addSimulationRR(_plotPlasma);
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
        _plotPlasma->setData(xTime,yTAC);
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
        _plotPlasma->setData(xTime,yTAC);
    }

    if ( _clearPlasmaPlot->isChecked() || whichColor > 9 )
        whichColor = 0;
    FUNC_INFO << "which color" << whichColor;
    _plotPlasma->setColor(colors[whichColor]);
    whichColor++;
    if ( fineScale )
        _plotPlasma->setPointSize(2);
    else
        _plotPlasma->setPointSize(5);
    _plotPlasma->conclude(0,true);
    _plotPlasma->plotDataAndFit(true);
    FUNC_EXIT;
}

void SimWindow::updateReferenceGraph()
{
    FUNC_ENTER;
    // update the plot: RR
    _plotRR->init();
    _plotRR->setLegendOn(true);
    _plotRR->setLabelYAxis("TAC (e.g., kBq/ml)");

    if ( simStartsFromPlasmaFit() || simStartsFromDataFit() )
        addFitRR(_plotRR);
    if ( simStartsFromPlasma() || simStartsFromPlasmaFit() )
        addSimulationRR(_plotRR);
    else
        addDataCurveRR(_plotRR);

    _plotRR->conclude(0,true);
    _plotRR->plotDataAndFit(true);
    FUNC_EXIT;
}

void SimWindow::addSimulationRR(plotData *whichPlot)
{
    FUNC_ENTER;
    whichPlot->addCurve(0,"RR: simulation");
    if ( whichPlot->getNumberCurves() == 1 )
        whichPlot->setColor(Qt::blue); // 1st curve
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
    FUNC_ENTER;
    whichPlot->addCurve(0,"RR: fit");
    if ( whichPlot->getNumberCurves() == 1 )
        whichPlot->setColor(Qt::blue); // 1st curve
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
    FUNC_ENTER;
    if ( !realDataAvailable() ) return;
    int indexInBox = _dataRefRegion->currentIndex();
    if ( indexInBox >= _dataTable.size() )
        qFatal("Fatal Error: the combo-box index exceeds the table index.");
    whichPlot->addCurve(0,"RR: ROI data");
    if ( whichPlot->getNumberCurves() == 1 )
        whichPlot->setColor(Qt::blue); // 1st curve
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
        _plotTarget->addCurve(0,"target: simulation");
    else
        _plotTarget->addCurve(0,"target");
    if ( _plotTarget->getNumberCurves() == 1 )
        _plotTarget->setColor(Qt::blue); // 1st curve
    else
        _plotTarget->setColor(Qt::gray);
    _plotTarget->setLabelYAxis("TAC (e.g., kBq/ml)");
    _plotTarget->setPointSize(3);
    int nTime = _simulator.getNumberTimeBinsCoarse();
    dVector xTime;   xTime.resize(nTime);
    dVector yTAC;    yTAC.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = _simulator.getTimeCoarse(jt);
        yTAC[jt]  = _simulator.getCtCoarse(jt);
    }
    _plotTarget->setData(xTime,yTAC);

    _plotTarget->addCurve(0,"exclusions");
    _plotTarget->setColor(Qt::black);
    _plotTarget->setLineThickness(0);
    _plotTarget->setPointSize(20);
    _plotTarget->setPointStyle(QCPScatterStyle::ssStar);
    _plotTarget->setExport(false);
    // set time vectors
    dVector xDataIgnore, yDataIgnore;
    bVector ignorePoint;
    for (int jt=0; jt<nTime; jt++) // loop over ALL points
    {
        if ( !(_PETRTM.getWeight(jt) != 0.) )
        {
            xDataIgnore.append(xTime[jt]);
            yDataIgnore.append(yTAC[jt]);
            ignorePoint.append(true);
        }
    }
    _plotTarget->setData(xDataIgnore, yDataIgnore, ignorePoint);

    FUNC_EXIT;
}
void SimWindow::addDataCurveTarget()
{
    FUNC_ENTER;
    if ( realDataAvailable() )
    {
        _plotTarget->addCurve(0,"target: ROI data");
        if ( _plotTarget->getNumberCurves() == 1 )
            _plotTarget->setColor(Qt::blue); // 1st curve
        else
            _plotTarget->setColor(Qt::gray);
        _plotTarget->setPointSize(5);
        int indexTarget = _dataTargetRegion->currentIndex();
        if ( indexTarget >= _dataTable.size() )
            qFatal("Fatal Error: the combo-box index exceeds the table index.");
        else if ( indexTarget < 0 ) return;
        _plotTarget->setData(_timeBins,_dataTable[indexTarget]);
    }
}

void SimWindow::updateTargetGraph()
{
    FUNC_ENTER;
    // update the plot
    _plotTarget->init();
    _plotTarget->setLegendOn(true);

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
    _plotTarget->addCurve(0,"fit");
    _plotTarget->setLineThickness(2);
    _plotTarget->setColor(Qt::red);
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
    _plotTarget->setData(xTime,yTAC);

    // add RR
    _plotTarget->addCurve(0,"RR");
    _plotTarget->setColor(Qt::gray);
    for (int jt=0; jt<nTime; jt++)
        yTAC[jt] = _PETRTM.getReferenceRegion(false,0,jt).y;
    _plotTarget->setData(xTime,yTAC);

    // add RR fit for fFRTM
    if ( _PETRTM.isForwardModel() )
    {
        _plotTarget->addCurve(0,"RR fit");
        _plotTarget->setColor(Qt::green);
        for (int jt=0; jt<nTime; jt++)
            yTAC[jt] = _PETRTM.getReferenceRegion(true,0,jt).y;
        _plotTarget->setData(xTime,yTAC);
    }

    _plotTarget->conclude(0,true);
    _plotTarget->plotDataAndFit(true);
    FUNC_EXIT;
}

void SimWindow::defineRTMModel()
{
    FUNC_ENTER;
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
    FUNC_ENTER;
    _plotBasis->init();
    _plotBasis->setLegendOn(true);

    _plotBasis->addCurve(0,"weights");
    _plotBasis->setPointSize(3);
    dVector xData;  xData.clear();
    dVector yData;  yData.clear();
    for (int jt=0; jt<_PETRTM.getNumberTimePoints(); jt++)
    {
        xData.append(jt);
        yData.append(_PETRTM.getWeight(jt));
    }
    _plotBasis->setData(xData, yData);

    int nBasis = _PETRTM.getNumberCoefficients();
    for ( int jBasis=0; jBasis<nBasis; jBasis++)
    {
        if ( _PETRTM.getBasisType(jBasis) == Type_R1 )
        {
            _plotBasis->addCurve(0,"R1");
            _plotBasis->setColor(Qt::red);
        }
        else if ( _PETRTM.getBasisType(jBasis) == Type_k2 )
        {
            _plotBasis->addCurve(0,"k2");
            _plotBasis->setColor(Qt::blue);
        }
        else if ( _PETRTM.getBasisType(jBasis) == Type_k2a )
        {
            _plotBasis->addCurve(0,"k2a");
            _plotBasis->setColor(Qt::green);
        }
        else if ( _PETRTM.getBasisType(jBasis) == Type_dCrdt )
        {
            _plotBasis->addCurve(0,"dCr/dt");
            _plotBasis->setColor(Qt::yellow);
        }
        else
        {
            _plotBasis->addCurve(0,"challenge");
            _plotBasis->setColor(Qt::magenta);
        }
        // add data
        yData.clear();
        for (int jt=0; jt<_PETRTM.getNumberTimePoints(); jt++)
            yData.append(_PETRTM.getBasisPoint(jBasis,jt));
        _plotBasis->setData(xData, yData);
        // rescale weights
        if ( _PETRTM.getBasisType(jBasis) == Type_k2 )
        {
            double yFraction = 0.67;
            plotCurve *targetCurve = _plotBasis->getThisCurvePointer();  // this
            plotCurve sourceCurve  = _plotBasis->_listOfCurves[0];       // weights
            double sourceMax = _plotBasis->getMaxYAbs(&sourceCurve) / yFraction;
            double targetMax = _plotBasis->getMaxYAbs(targetCurve);
            if ( targetMax != 0. && sourceMax != 0. )
            {
                _plotBasis->setYAxisRatio(sourceMax / targetMax);
                _plotBasis->setScaleFactorY(0,1. / _plotBasis->getYAxisRatio());
                for (int jt=0; jt<_PETRTM.getNumberTimePoints(); jt++)
                    _plotBasis->_listOfCurves[0].yData[jt] /= _plotBasis->getYAxisRatio();
            }
        }
    }

    _plotBasis->conclude(0,true);
    _plotBasis->plotDataAndFit(true);
}

void SimWindow::updateAICvsTau4Graph()
{
    qInfo() << "updateAICvsTau4Graph enter";
    FUNC_ENTER;
    bool sweep = _fitk4RadioButton->isChecked() ||
            _fitk2RefRadioButton->isChecked()   ||
            _fitDVRadioButton->isChecked()      ||
            _fitk4k2RefRadioButton->isChecked();
    if ( !sweep ) return;

    dVector valueVector, AICVector;
    double bestValue;

    if ( _fitk4k2RefRadioButton->isChecked() )
    {
        valueVector = _globalFit.tau4Vector;
        AICVector   = _globalFit.AICVector1D;
        bestValue   = _PETRTM.getTau4InRun(0);
    }
    else
    {
        valueVector = _PETRTM.getAICXVector();
        AICVector   = _PETRTM.getAICYVector();
        bestValue   = _PETRTM.getAICXBestValue();
    }

    _plotAIC->init();
    _plotAIC->setLegendOn(false);
    _plotAIC->addCurve(0,"AIC vs 1/k4");
    if ( valueVector.size() > 0 )
        _plotAIC->setData(valueVector, AICVector);
    _plotAIC->conclude(0,true);
    _plotAIC->plotDataAndFit(true);
    _plotAIC->setPositionTracer(0,bestValue);
    qInfo() << "updateAICvsTau4Graph exit";
}

void SimWindow::runSimulationAndAnalysis()
{
    if ( _fitk4k2RefRadioButton->isChecked() )
        calculateForwardTau4andTau2Ref();
    else
        updateAllGraphs();
}

void SimWindow::updateAllGraphs()
{
    FUNC_ENTER;
    // run the simulation
    qInfo() << "run simulation";
    _simulator.run();
    qInfo() << "run simulation done";

    updatePlasmaGraph();
    setReferenceRegionForAnalysis();
    updateReferenceGraph();
    qInfo() << "analyze TAC";
    analyzeTAC();
    qInfo() << "analyze TAC done";
    updateAICvsTau4Graph();
    updateTargetGraph();
    updateBasisGraph();  // update basis graph AFTER target graph, because basis functions use target curve
    double AIC;
    if ( _PETRTM.isForwardModel() )
        AIC = _PETRTM.getSimulationAIC(0);
    else
        AIC = _PETRTM.getAIC();
    _analysisGroupBox->setTitle(QString("Target Region: analysis, AIC = %1").arg(AIC));
    if ( _radioShowAIC->isChecked() )
        showAICPlot();
    else
        showTarget();
}
void SimWindow::analyzeTAC()
{
    FUNC_ENTER;
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
    _PETRTM.setTissueRegion(tissueVector);
    _PETRTM.definePETConditions("a c"); // don't define R1, which is not valid for RTM2
    _PETRTM.prepare();
    dMatrix fitVector;     fitVector.resize(1);     fitVector[0].resize(nTime);
    _PETRTM.fitData(tissueVector,fitVector);

    // BPnd
    double truth = _simulator.getBP0();
    double guess = _PETRTM.getBP0InRun(0).x;
    _errorBPnd->setText(analyzeString(truth,guess));
    // k2
    truth = _simulator.getk2();
    guess = _PETRTM.getk2InRun(0).x;
    _errork2->setText(analyzeString(truth,guess));

    // R1
    truth = _simulator.getR1();
    guess = _PETRTM.getR1InRun(0).x;
    _errorR1->setText(analyzeString(truth,guess));
    FUNC_INFO << "** R1" << truth << guess;

    // tau2Ref
    truth = _simulator.getTau2Ref();
    guess = _PETRTM.getTau2RefInRun(0);
    _errorTau2Ref->setText(analyzeString(truth,guess));
    if ( _fitk2RefRadioButton->isChecked() )
        _tau2RefAnalysis->setText(QString::number(guess,'g',2));

    // DV
    truth = _simulator.getDV();
    guess = _PETRTM.getDVInRun(0);
    _errorDV->setText(analyzeString(truth,guess));
    if ( _fitDVRadioButton->isChecked() )
        _DVAnalysis->setText(QString::number(guess,'g',2));

    // challenge
    if ( _fitChallenge->isChecked() )
    {
        truth = _simulator.getChallengeMagPercent();  // delta_BPnd abs
        double deltaBPAbs;
        guess = getChallengeMagPercentFromAnalysis(deltaBPAbs);
        QString valueStringPercent;  valueStringPercent.setNum(guess,'g',2);
        QString valueStringAbs;      valueStringAbs.setNum(deltaBPAbs,'g',2);
        QString diffString;
        if ( guess == 0. )
            diffString = "----";
        else
            diffString.setNum(guess-truth,'g',2);
        if ( analyzeRealData() )
            _errorChallenge->setText(valueStringPercent + " %");
        else
            _errorChallenge->setText(valueStringPercent + " %, abs = " + valueStringAbs + " diff = " + diffString + " %");
    }

    // k4
    if ( _fitk4RadioButton->isChecked() )
    {
        truth = _simulator.getTau4();  // delta_BPnd abs
        guess = _PETRTM.getTau4InRun(0);
        QString valueString;  valueString.setNum(guess,'g',2);
        if ( guess == 0. )
            _errork4->setText("----");
        else
            _errork4->setText(analyzeString(truth,guess));
        _tau4Analysis->setText(QString::number(guess,'g',2));
    }

    // sigma2
    double sigma;
    if ( _PETRTM.isForwardModel() )
        sigma = qSqrt(_PETRTM.getSimulationSigma2(0));
    else
        sigma = qSqrt(_PETRTM.getSigma2());
    QString valueString;  valueString.setNum(sigma,'g',3);
    _sigma->setText(valueString);

    if ( analyzeRealData() )
    {
//        if ( changedAnalysisRegion )
//        {
            double BPnd = _PETRTM.getBP0InRun(0).x;
            double R1   = _PETRTM.getR1InRun(0).x;
            _simulator.setBP0(BPnd);  _simulator.setR1(R1);
            QString text; _BPnd->setText(text.setNum(BPnd));  _R1->setText(text.setNum(R1));
            _simulator.generateTargetTAC();
//        }
        addSimulationTarget();
        updateTargetGraph();
        _plotTarget->conclude(0,true);
        _plotTarget->plotDataAndFit(true);
    }

}
double SimWindow::getChallengeMagPercentFromAnalysis(double &deltaBPAbs)
{
    FUNC_ENTER;
    _PETRTM.updateConditions();
    int nConditions = _PETRTM.getNumberConditions();
    if ( nConditions == 2 )
    {
        _PETRTM.setCurrentCondition(1);
        _PETRTM.evaluateCurrentCondition();
        double BP0 = _PETRTM.getBP0InRun(0).x;
        deltaBPAbs = _PETRTM.getBPndInCurrentCondition();
        double percentChange = 100. * deltaBPAbs / BP0;
        return percentChange;
    }
    else
        return 0.;
}
QString SimWindow::analyzeString(double truth, double guess)
{
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
    _nThreads = indexInBox + 1;
    QString numberString;
    _nSamplesBPnd->setText((numberString.setNum(_nThreads*_numberSimulationsPerThread)));
}

void SimWindow::changedNumberBins()
{
    FUNC_ENTER;
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
    FUNC_ENTER;
    int iBin = indexPlusOne-1;
    double duration = _simulator.getDurationPerBinSec(iBin);
    QString text; _binDuration->setText(text.setNum(duration));
}
void SimWindow::changedBinDuration()
{
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
void SimWindow::changedInfusionDelay()
{
    FUNC_ENTER;
    QString stringEntered = _KBolDelay->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setKBolDelay(value);  scaleSimulationToDataAverage();
        updateAllGraphs();
    }
    else
        _KBolDelay->setText(stringEntered.setNum(_simulator.getKBolDelay()));
}
void SimWindow::changedTau2Ref()
{
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
    QString stringEntered = _BPnd->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setBP0(value);
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _BPnd->setText(stringEntered.setNum(_simulator.getBP0()));
}
void SimWindow::changedR1()
{
    FUNC_ENTER;
    QString stringEntered = _R1->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setR1(value);
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _R1->setText(stringEntered.setNum(_simulator.getR1()));
}
void SimWindow::changedTau4()
{
    FUNC_ENTER;
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
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _tau4->setText(stringEntered.setNum(_simulator.getTau4()));
}
void SimWindow::changedDV()
{
    FUNC_ENTER;
    QString stringEntered = _DV->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setDV(value);
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _DV->setText(stringEntered.setNum(_simulator.getDV()));
}
void SimWindow::changedChallengeTime()
{
    FUNC_ENTER;
    QString stringEntered = _challengeTime->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeTime(value);
        int indexChallenge = _PETRTM.getEventIndex('c');
        _PETRTM.setChallengeOnset(indexChallenge,0,_simulator.getChallengeTime());
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _challengeTime->setText(stringEntered.setNum(_simulator.getChallengeTime()));
}
void SimWindow::changedChallengeMag()
{
    FUNC_ENTER;
    QString stringEntered = _challengeMagPercent->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeMag(value);
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _challengeMagPercent->setText(stringEntered.setNum(_simulator.getChallengeMagPercent()));
    QChar delta = QChar(0x0394);
    _challengeMagLabel->setText(QString("%1BPnd in % (abs = %2)").arg(delta).arg(_simulator.getChallengeMagAbs()));
}
void SimWindow::changedNoiseTar()
{
    FUNC_ENTER;
    QString utilityString = _noiseTar->text();
    bool ok;
    double value = utilityString.toDouble(&ok);
    if ( ok )
    {
        _simulator.setNoiseTar(value);
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _noiseTar->setText(utilityString.setNum(_simulator.getNoiseTar()));
    bool noisy = _simulator.getNoiseRef() != 0. || _simulator.getNoiseTar() != 0.;
    setThreadVisibility(noisy);
}
void SimWindow::changedIgnoreString()
{
    FUNC_ENTER;
    QString stringEntered = _ignoreString->text();
    _PETRTM.setIgnoredPoints(0,true,stringEntered);
    if ( _PETRTM.isSRTM() ) updateAllGraphs();
}
void SimWindow::changedTau2RefAnalysis()
{
    FUNC_ENTER;
    QString stringEntered = _tau2RefAnalysis->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _PETRTM.setTau2Ref(0,value);
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _tau2RefAnalysis->setText(stringEntered.setNum(_PETRTM.getTau2RefInRun(0)));
}
void SimWindow::changedTau4Analysis()
{
    FUNC_ENTER;
    QString stringEntered = _tau4Analysis->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _PETRTM.setTau4(0,value);
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _tau4Analysis->setText(stringEntered.setNum(_PETRTM.getTau4InRun(0)));
}
void SimWindow::changedDVAnalysis()
{
    FUNC_ENTER;
    QString stringEntered = _DVAnalysis->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _PETRTM.setDV(0,value);
        if ( _PETRTM.isSRTM() ) updateAllGraphs();
    }
    else
        _DVAnalysis->setText(stringEntered.setNum(_PETRTM.getDVInRun(0)));
}
void SimWindow::changedParSearchStart()
{
    QString stringEntered = _parSearchStart->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok && _fitk4RadioButton->isChecked() )
        _PETRTM.setSearchParTau4_X(value);
    else if ( ok && _fitDVRadioButton->isChecked() )
        _PETRTM.setSearchParDV_X(value);
    else
        _PETRTM.setSearchParTau2Ref_X(value);
}
void SimWindow::changedParSearchEnd()
{
    QString stringEntered = _parSearchEnd->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok && _fitk4RadioButton->isChecked() )
        _PETRTM.setSearchParTau4_Y(value);
    else if ( ok && _fitDVRadioButton->isChecked() )
        _PETRTM.setSearchParDV_Y(value);
    else
        _PETRTM.setSearchParTau2Ref_Y(value);
}
void SimWindow::changedParSearchStep()
{
    QString stringEntered = _parSearchStep->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok && _fitk4RadioButton->isChecked() )
        _PETRTM.setSearchParTau4_Z(value);
    else if ( ok && _fitDVRadioButton->isChecked() )
        _PETRTM.setSearchParDV_Z(value);
    else
        _PETRTM.setSearchParTau2Ref_Z(value);
}
void SimWindow::changedCheckBoxFitk4(bool state)
{
    FUNC_ENTER;
    _errork4Label->setVisible(state);
    _errork4->setVisible(state);
    _checkBoxk4ErrGraph->setChecked(state);
    _checkBoxk4ErrGraph->setVisible(state);
    _checkBoxDVErrGraph->setChecked(!state);
    _checkBoxDVErrGraph->setVisible(!state);

    dPoint3D parSearch = _PETRTM.getSearchParTau4();
    QString number;
    _parSearchStart->setText(number.setNum(parSearch.x));
    _parSearchEnd->setText(number.setNum(parSearch.y));
    _parSearchStep->setText(number.setNum(parSearch.z));
    _parSearchStart->setVisible(state);
    _parSearchEnd->setVisible(state);
    _parSearchStep->setVisible(state);

    //    clearTimeCurves();
    _PETRTM.setFitk4(state);
    _PETRTM.setPrepared(false);
    if ( _PETRTM.isSRTM() ) updateAllGraphs();
}
void SimWindow::changedCheckBoxFitDV(bool state)
{
    FUNC_ENTER;
    _errorDVLabel->setVisible(state);
    _errorDV->setVisible(state);
    _checkBoxDVErrGraph->setChecked(state);
    _checkBoxDVErrGraph->setVisible(state);
    _checkBoxk4ErrGraph->setChecked(!state);
    _checkBoxk4ErrGraph->setVisible(!state);

    dPoint3D parSearch = _PETRTM.getSearchParDV();
    QString number;
    _parSearchStart->setText(number.setNum(parSearch.x));
    _parSearchEnd->setText(number.setNum(parSearch.y));
    _parSearchStep->setText(number.setNum(parSearch.z));
    _parSearchStart->setVisible(state);
    _parSearchEnd->setVisible(state);
    _parSearchStep->setVisible(state);

//    clearTimeCurves();
    _PETRTM.setFitDV(state);
    _PETRTM.setPrepared(false);
    if ( _PETRTM.isSRTM() ) updateAllGraphs();
}
void SimWindow::changedCheckBoxFitk2Ref(bool state)
{
    FUNC_ENTER;
    _errorTau2RefLabel->setVisible(state);
    _errorTau2Ref->setVisible(state);

    dPoint3D parSearch = _PETRTM.getSearchParTau2Ref();
    QString number;
    _parSearchStart->setText(number.setNum(parSearch.x));
    _parSearchEnd->setText(number.setNum(parSearch.y));
    _parSearchStep->setText(number.setNum(parSearch.z));
    _parSearchStart->setVisible(state);
    _parSearchEnd->setVisible(state);
    _parSearchStep->setVisible(state);

//    clearTimeCurves();
    _PETRTM.setFitk2Ref(state);
    _PETRTM.setPrepared(false);
    if ( _PETRTM.isSRTM() ) updateAllGraphs();
}
void SimWindow::changedCheckBoxChallenge(bool state)
{
    FUNC_ENTER;
    _errorChallengeLabel->setVisible(state);
    _errorChallenge->setVisible(state);
    _checkBoxChallErrVsBPnd->setChecked(state);
    _checkBoxChallErrVsBPnd->setVisible(state);
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
    if ( _PETRTM.isSRTM() ) updateAllGraphs();
}

void SimWindow::changedBPndLow()
{
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
    QString stringEntered = _BPndStep->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
        _BPndStepValue = value;
    else
        _BPndStep->setText(stringEntered.setNum(_BPndStepValue));
}
void SimWindow::changedTau4High()
{
    FUNC_ENTER;
    QString stringEntered = _tau4High->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
        _tau4HighValue = value;
    else
        _tau4High->setText(stringEntered.setNum(_tau4HighValue));
}
void SimWindow::changedTau4Step()
{
    FUNC_ENTER;
    QString stringEntered = _tau4Step->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
        _tau4StepValue = value;
    else
        _tau4Step->setText(stringEntered.setNum(_tau4StepValue));
}
void SimWindow::changedNumberSimulationsBPnd()
{
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
    _plotErrBPndVsBPnd->getPlotSurface()->setVisible(_checkBoxBPndErrVsBPnd->isChecked());
    _plotErrChallVsBPnd->getPlotSurface()->setVisible(_checkBoxChallErrVsBPnd->isChecked());
    _plotTau2RefVsBPnd->getPlotSurface()->setVisible(_checkBoxTau2RefGraph->isChecked());
    _plotErrk4VsBPnd->getPlotSurface()->setVisible(_checkBoxk4ErrGraph->isChecked());
    _plotErrDVVsBPnd->getPlotSurface()->setVisible(_checkBoxDVErrGraph->isChecked());
}
void SimWindow::changedVersusTau4Graphs()
{
    FUNC_ENTER;
    _plotErrBPndVsTau4->getPlotSurface()->setVisible(_checkBoxBPndErrVsTau4->isChecked());
    _plotErrChallVsTau4->getPlotSurface()->setVisible(_checkBoxChallErrVsTau4->isChecked());
    _plotAICVsTau4->getPlotSurface()->setVisible(_checkBoxAICVsTau4->isChecked());
}
void SimWindow::clearBPndCurves()
{
    FUNC_ENTER;
    QCPRange xRange;
    xRange.lower = _BPndLowValue  - _BPndStepValue;
    xRange.upper = _BPndHighValue + _BPndStepValue;

    _plotErrBPndVsBPnd->init();
    _plotErrBPndVsBPnd->setLabelXAxis("BPnd");
//    _plotErrBPndVsBPnd->setLabelXAxis("");
    _plotErrBPndVsBPnd->setLabelYAxis("BPnd % err");
    _plotErrBPndVsBPnd->plotDataAndFit(true);
    _plotErrBPndVsBPnd->setXRange(xRange);

    _plotErrChallVsBPnd->init();
    _plotErrChallVsBPnd->setLabelXAxis("BPnd");
    _plotErrChallVsBPnd->setLabelYAxis(QString("%1BPnd abs err").arg(QChar(0x0394)));
    _plotErrChallVsBPnd->plotDataAndFit(true);
    _plotErrChallVsBPnd->setXRange(xRange);

    _plotTau2RefVsBPnd->init();
    _plotTau2RefVsBPnd->setLabelXAxis("BPnd");
    _plotTau2RefVsBPnd->setLabelYAxis("1/k2'");
    _plotTau2RefVsBPnd->plotDataAndFit(true);
    _plotTau2RefVsBPnd->setXRange(xRange);

    _plotErrk4VsBPnd->init();
    _plotErrk4VsBPnd->setLabelXAxis("BPnd");
    _plotErrk4VsBPnd->setLabelYAxis("1/k4");
    _plotErrk4VsBPnd->plotDataAndFit(true);
    _plotErrk4VsBPnd->setXRange(xRange);

    _plotErrDVVsBPnd->init();
    _plotErrDVVsBPnd->setLabelXAxis("BPnd");
    _plotErrDVVsBPnd->setLabelYAxis("DV");
    _plotErrDVVsBPnd->plotDataAndFit(true);
    _plotErrDVVsBPnd->setXRange(xRange);
}

void SimWindow::clearTau4Curves()
{
    FUNC_ENTER;
    QCPRange xRange = {0.,_tau4HighValue};

    _plotErrBPndVsTau4->init();
    _plotErrBPndVsTau4->setLabelXAxis("tau4 (min)");
//    _plotErrBPndVsBPnd->setLabelXAxis("");
    _plotErrBPndVsTau4->setLabelYAxis("BPnd % err");
    _plotErrBPndVsTau4->plotDataAndFit(true);
    _plotErrBPndVsTau4->setXRange(xRange);

    _plotErrChallVsTau4->init();
    _plotErrChallVsTau4->setLabelXAxis("tau4 (min)");
    _plotErrChallVsTau4->setLabelYAxis(QString("%1BPnd abs err").arg(QChar(0x0394)));
    _plotErrChallVsTau4->plotDataAndFit(true);
    _plotErrChallVsTau4->setXRange(xRange);

    _plotAICVsTau4->init();
    _plotAICVsTau4->setLabelXAxis("tau4 (min)");
    _plotAICVsTau4->setLabelYAxis("AIC");
    _plotAICVsTau4->plotDataAndFit(true);
    _plotAICVsTau4->setXRange(xRange);
}

void SimWindow::calculateTau4Curves()
{  // sweep tau4 truth or tau4 analysis??
    FUNC_ENTER;
    bool noisy = _simulator.getNoiseRef() != 0. || _simulator.getNoiseTar() != 0.;
//    if ( noisy || _PETRTM.isForwardFitk4() )
    if ( noisy )
    {
        calculateTau4CurvesInThreads();
        return;
    }

    double saveTau4    = _PETRTM.getTau4InRun(0); // tau4 will be changed, so save the value for later restoration
    double saveTau2Ref = _PETRTM.getTau2RefFRTMInRun(0);
    QString numberString;
    int nTime = _simulator.getNumberTimeBinsCoarse();
    dMatrix refRegion;      refRegion.resize(1);    refRegion[0].resize(nTime);
    dMatrix tissueVector;   tissueVector.resize(1); tissueVector[0].resize(nTime);
    dMatrix fitVector;      fitVector.resize(1);    fitVector[0].resize(nTime);

    dVector xVector, errBPnd, errChall, AIC;
    double sigma2Sum=0.;

    double BP0 = _simulator.getBP0();

    for (int jt=0; jt<nTime; jt++)
        tissueVector[0][jt] = _simulator.getCtCoarse(jt);
    _PETRTM.setTissueRegion(tissueVector);

    for (double tau4=0; tau4<=_tau4HighValue; tau4 += _tau4StepValue)
    {
        xVector.append(tau4);
        _PETRTM.setTau4(0,tau4);
//        _simulator.setTau4(tau4);  // set BP0 and run a series of samples

        if ( _PETRTM.isRTM2() )
        {
            double bestSigma2 = 1.e20;  double bestTau2Ref = 0.;
            for (double tau2Ref=0; tau2Ref<=2.*saveTau2Ref; tau2Ref += 0.05)
            {
                _PETRTM.setTau2RefFRTM(0,tau2Ref);
                _PETRTM.prepare();
                _PETRTM.fitData(tissueVector,fitVector);
                double sigma2 = _PETRTM.getSimulationSigma2(0);
                if ( sigma2 < bestSigma2 )
                {
                    bestSigma2 = sigma2;
                    bestTau2Ref = tau2Ref;
                }
            }
            _PETRTM.setTau2RefFRTM(0,bestTau2Ref);
        }

        _PETRTM.prepare();
        _PETRTM.fitData(tissueVector,fitVector);

        // get simulator AIC from _PETRM ("getAIC" returns GLM AIC, which is not the right one)
        AIC.append(_PETRTM.getSimulationAIC(0));

        if ( _PETRTM.isForwardModel() )
            sigma2Sum += _PETRTM.getSimulationSigma2(0);
        else
            sigma2Sum += _PETRTM.getSigma2();

        // update the BP error
        double guess = _PETRTM.getBP0InRun(0).x;
        errBPnd.append(percentageError(guess,BP0));
//        errBPnd.append(guess-BP0);
        // update the challenge error
        double truth = _simulator.getChallengeMagPercent();
        double deltaBPAbs;
        guess = getChallengeMagPercentFromAnalysis(deltaBPAbs);
        errChall.append(guess - truth);

        /*
        // update the 1/k4 error
        if ( _PETRTM.isForwardFitk4() )
        {
            truth = _simulator.getTau4();
            guess = _PETRTM.getTau4InRun(0);
            FUNC_INFO << "** tau4 **" << truth << guess;
            errTau4.append(percentageError(guess,truth));
        }
        */

    } // loop over tau4 values
    _PETRTM.setTau4(0,saveTau4);
    _PETRTM.setTau2RefFRTM(0,saveTau2Ref);

    sigma2Sum /= static_cast<double>(xVector.size());
    _sigma2Label->setText(QString("sigma2 = %1").arg(sigma2Sum));

    // restore the value of tau4
    _tau4->setText(numberString.setNum(saveTau4));  changedBPND();

    QColor colors[10] = {Qt::black, Qt::red, Qt::blue, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};
    int nCurves = _plotErrBPndVsTau4->getNumberCurves();
    int iColor  = nCurves%10;

    _plotErrBPndVsTau4->addCurve(0,"error_BPnd");
    _plotErrBPndVsTau4->setPointSize(5);
    _plotErrBPndVsTau4->setColor(colors[iColor]);
    _plotErrBPndVsTau4->setData(xVector,errBPnd);
    _plotErrBPndVsTau4->conclude(0,true);
    _plotErrBPndVsTau4->plotDataAndFit(true);

    _plotErrChallVsTau4->addCurve(0,"error_Challenge");
    _plotErrChallVsTau4->setPointSize(5);
    _plotErrChallVsTau4->setColor(colors[iColor]);
    _plotErrChallVsTau4->setData(xVector,errChall);
    _plotErrChallVsTau4->conclude(0,true);
    _plotErrChallVsTau4->plotDataAndFit(true);

    _plotAICVsTau4->addCurve(0,"AIC");
    _plotAICVsTau4->setPointSize(5);
    _plotAICVsTau4->setColor(colors[iColor]);
    _plotAICVsTau4->setData(xVector,AIC);
    _plotAICVsTau4->conclude(0,true);
    _plotAICVsTau4->plotDataAndFit(true);
}

void SimWindow::calculateBPndCurves()
{
    FUNC_ENTER;
    bool noisy = _simulator.getNoiseRef() != 0. || _simulator.getNoiseTar() != 0.;
//    if ( noisy || _PETRTM.isForwardFitk4() )
    if ( noisy )
    {
        calculateBPndCurvesInThreads();
        return;
    }

    double saveBPnd = _simulator.getBP0(); // BPnd will be changed, so save the value for later restoration
    QString numberString;
    int nTime = _simulator.getNumberTimeBinsCoarse();
    dMatrix refRegion;      refRegion.resize(1);    refRegion[0].resize(nTime);
    dMatrix tissueVector;   tissueVector.resize(1); tissueVector[0].resize(nTime);
    dMatrix fitVector;      fitVector.resize(1);    fitVector[0].resize(nTime);

    dVector xVector, errBPnd, errChall, tau2Ref, errTau4, errDV, errBPndRTM2, errChallRTM2;
    double sigma2Sum=0.;

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
        _PETRTM.setTissueRegion(tissueVector);

        _PETRTM.prepare();
        _PETRTM.fitData(tissueVector,fitVector);

        if ( _PETRTM.isForwardModel() )
            sigma2Sum += _PETRTM.getSimulationSigma2(0);
        else
            sigma2Sum += _PETRTM.getSigma2();

        // update the BP error
        double guess = _PETRTM.getBP0InRun(0).x;
        errBPnd.append(percentageError(guess,BP0));
//        errBPnd.append(guess-BP0);
        // update the challenge error
        double truth = _simulator.getChallengeMagPercent();
        double deltaBPAbs;
        guess = getChallengeMagPercentFromAnalysis(deltaBPAbs);
        errChall.append(guess - truth);
//        errChall.append(percentageError(guess,truth));
        // update the 1/k4 error
        if ( _PETRTM.getFitk4() )
        {
            truth = _simulator.getTau4();
            guess = _PETRTM.getTau4InRun(0);
            FUNC_INFO << "** tau4 **" << truth << guess;
            errTau4.append(percentageError(guess,truth));
        }
        else if ( _PETRTM.getFitDV() )
        {
            truth = _simulator.getDV();
            guess = _PETRTM.getDVInRun(0);
            FUNC_INFO << "** DV **" << truth << guess;
            errDV.append(percentageError(guess,truth));
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
            guess = _PETRTM.getBP0InRun(0).x;
            errBPndRTM2.append(percentageError(guess,truth));
            truth = _simulator.getChallengeMagPercent();
            double deltaBPAbs;
            guess = getChallengeMagPercentFromAnalysis(deltaBPAbs);
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
    int nCurves = _plotErrBPndVsBPnd->getNumberCurves();
    int iColor  = nCurves%10;

    _plotErrBPndVsBPnd->addCurve(0,"error_BPnd");
    _plotErrBPndVsBPnd->setPointSize(5);
    _plotErrBPndVsBPnd->setColor(colors[iColor]);

    _plotErrChallVsBPnd->addCurve(0,"error_Challenge");
    _plotErrChallVsBPnd->setPointSize(5);
    _plotErrChallVsBPnd->setColor(colors[iColor]);

    _plotTau2RefVsBPnd->addCurve(0,"1/k2'");
    _plotTau2RefVsBPnd->setPointSize(5);
    _plotTau2RefVsBPnd->setColor(colors[iColor]);
    if ( _PETRTM.isRTM2() )
        _plotTau2RefVsBPnd->setPointStyle(QCPScatterStyle::ssCross);

    _plotErrk4VsBPnd->addCurve(0,"1/k4");
    _plotErrk4VsBPnd->setPointSize(5);
    _plotErrk4VsBPnd->setColor(colors[iColor]);

    _plotErrDVVsBPnd->addCurve(0,"DV");
    _plotErrDVVsBPnd->setPointSize(5);
    _plotErrDVVsBPnd->setColor(colors[iColor]);

    _plotErrBPndVsBPnd->setData(xVector,errBPnd);
    _plotErrChallVsBPnd->setData(xVector,errChall);
    bool calculateTau2Ref = !_PETRTM.isRTM2() || _checkBoxTau2RefGraph->isChecked();
    if ( calculateTau2Ref )
        _plotTau2RefVsBPnd->setData(xVector,tau2Ref);
    if ( _PETRTM.getFitk4() )
        _plotErrk4VsBPnd->setData(xVector,errTau4);
    else if ( _PETRTM.getFitDV() )
        _plotErrDVVsBPnd->setData(xVector,errDV);

    if ( _PETRTM.isRTM2() && _checkBoxTau2RefGraph->isChecked() )
    { // add a new curves to the two error plots
        _plotErrBPndVsBPnd->addCurve(0,"error_BPnd");
        _plotErrBPndVsBPnd->setPointSize(5);
        _plotErrBPndVsBPnd->setPointStyle(QCPScatterStyle::ssCross);
        _plotErrBPndVsBPnd->setColor(colors[iColor]);
        _plotErrBPndVsBPnd->setData(xVector,errBPndRTM2);
        _plotErrChallVsBPnd->addCurve(0,"error_Challenge");
        _plotErrChallVsBPnd->setPointSize(5);
        _plotErrChallVsBPnd->setPointStyle(QCPScatterStyle::ssCross);
        _plotErrChallVsBPnd->setColor(colors[iColor]);
        _plotErrChallVsBPnd->setData(xVector,errChallRTM2);
    }

    _plotErrBPndVsBPnd->conclude(0,true);
    _plotErrChallVsBPnd->conclude(0,true);
    _plotTau2RefVsBPnd->conclude(0,true);
    _plotErrk4VsBPnd->conclude(0,true);
    _plotErrDVVsBPnd->conclude(0,true);

    _plotErrBPndVsBPnd->plotDataAndFit(true);
    _plotErrChallVsBPnd->plotDataAndFit(true);
    _plotTau2RefVsBPnd->plotDataAndFit(true);
    _plotErrk4VsBPnd->plotDataAndFit(true);
    _plotErrDVVsBPnd->plotDataAndFit(true);
}

double SimWindow::bestTau2RefForRTM2()
{ // search for root
    FUNC_ENTER;
    double tau2Ref = _PETRTM.getTau2RefInRun(0);
    double trueBP0 = _simulator.getBP0();
    double estBP0  = _PETRTM.getBP0InRun(0).x;
    double error   = estBP0 - trueBP0;
    QString numberString;
    if ( error > 0. )
    {
        while ( error > 0. )
        {
            tau2Ref -= 0.05;
            _tau2RefAnalysis->setText(numberString.setNum(tau2Ref));  changedTau2RefAnalysis();
            estBP0 = _PETRTM.getBP0InRun(0).x;
            error = estBP0 - trueBP0;
        }
    }
    else
    {
        while ( error < 0. )
        {
            tau2Ref += 0.05;
            _tau2RefAnalysis->setText(numberString.setNum(tau2Ref));  changedTau2RefAnalysis();
            estBP0 = _PETRTM.getBP0InRun(0).x;
            error = estBP0 - trueBP0;
        }
    }
    return tau2Ref;
}

void SimWindow::clearTimeCurves()
{
    FUNC_ENTER;
    QCPRange xRange;
    xRange.lower = _timeLowValue;
    xRange.upper = _timeHighValue;

    _plotErrBPndOrChallVsTime->init();
    _plotErrBPndOrChallVsTime->setLabelXAxis("time");
    if ( _fitChallenge->isChecked() )
        _plotErrBPndOrChallVsTime->setLabelYAxis("Challenge err (abs)");
    else
        _plotErrBPndOrChallVsTime->setLabelYAxis("BPnd % err");
    _plotErrBPndOrChallVsTime->plotDataAndFit(true);
    _plotErrBPndOrChallVsTime->setXRange(xRange);
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
    FUNC_ENTER;
    QColor colors[10] = {Qt::black, Qt::red, Qt::blue, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};
    int nCurves = _plotErrBPndOrChallVsTime->getNumberCurves();
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
                double truth = _simulator.getChallengeMagPercent();
                double deltaBPAbs;
                double guess = getChallengeMagPercentFromAnalysis(deltaBPAbs);
                errVector.append(guess - truth);
            }
            else
            { // calculate the percent error in BPnd
                _ignoreString->setText(QString("%1-n").arg((jBin)));
                _PETRTM.setIgnoredPoints(0,true,_ignoreString->text());
                analyzeTAC();

                double truth = _simulator.getBP0();
                double guess = _PETRTM.getBP0InRun(0).x;
                errVector.append(percentageError(guess,truth));

                _ignoreString->setText("");
                _PETRTM.setIgnoredPoints(0,true,_ignoreString->text());
                analyzeTAC();
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

    _plotErrBPndOrChallVsTime->addCurve(0,"error_BPnd");
    _plotErrBPndOrChallVsTime->setPointSize(5);
    _plotErrBPndOrChallVsTime->setColor(colors[iColor]);
    _plotErrBPndOrChallVsTime->setData(xVector,errVector);
    _plotErrBPndOrChallVsTime->conclude(0,true);
    _plotErrBPndOrChallVsTime->plotDataAndFit(true);
}

void SimWindow::calculateTau4CurvesInThreads()
{
    FUNC_ENTER;
    updateLieDetectorProgress(10*_nThreads);

    clearErrorMatrices(false);

    /////////////////////////////////////////////////////////
    // Create the average volume.
    qRegisterMetaType<dVector>("dMatrix");
    QVector<lieDetectorTau4 *> simSegment;
    simSegment.resize(_nThreads);
    for (int jThread=0; jThread<_nThreads; jThread++)
    {
        simSegment[jThread] = new lieDetectorTau4(_numberSimulationsPerThread, _tau4Vector, _simulator, _PETRTM);
        connect(simSegment[jThread], SIGNAL(progressLieDetector(int)), this, SLOT(updateLieDetectorProgress(int)));
        connect(simSegment[jThread], SIGNAL(finishedLieDetector(dMatrix,dMatrix,dMatrix,dMatrix,double)),
                this,SLOT(finishedLieDetectorTau4OneThread(dMatrix,dMatrix,dMatrix,dMatrix,double)));
        QThreadPool::globalInstance()->start(simSegment[jThread]);
    }
}

void SimWindow::clearErrorMatrices(bool sweepBP)
{
    _BP0Vector.clear();
    _tau4Vector.clear();

    _errBPndMatrix.clear();
    _errChallMatrix.clear();
    _errTau4Matrix.clear();
    _errDVMatrix.clear();
    _tau2RefMatrix.clear();
    _AICMatrix.clear();

    if ( sweepBP )
    {
        for (double BP0=_BPndLowValue; BP0<=_BPndHighValue; BP0 += _BPndStepValue)
            _BP0Vector.append(BP0);
        int nBP0Values = _BP0Vector.size();
        _errBPndMatrix.resize(nBP0Values);
        _errChallMatrix.resize(nBP0Values);
        _errTau4Matrix.resize(nBP0Values);
        _errDVMatrix.resize(nBP0Values);
        _tau2RefMatrix.resize(nBP0Values);
        _AICMatrix.resize(nBP0Values);
    }
    else
    {
        for (double tau4=0; tau4<=_tau4HighValue; tau4 += _tau4StepValue)
            _tau4Vector.append(tau4);
        int nTau4Values = _tau4Vector.size();
        _errBPndMatrix.resize(nTau4Values);
        _errChallMatrix.resize(nTau4Values);
        _errTau4Matrix.resize(nTau4Values);
        _errDVMatrix.resize(nTau4Values);
        _tau2RefMatrix.resize(nTau4Values);
        _AICMatrix.resize(nTau4Values);
    }
}

void SimWindow::calculateBPndCurvesInThreads()
{
    FUNC_ENTER;
    updateLieDetectorProgress(10*_nThreads);

    clearErrorMatrices(true);

    /////////////////////////////////////////////////////////
    // Create the average volume.
    qRegisterMetaType<dVector>("dMatrix");
    QVector<lieDetectorBPnd *> simSegment;
    simSegment.resize(_nThreads);
    for (int jThread=0; jThread<_nThreads; jThread++)
    {
        simSegment[jThread] = new lieDetectorBPnd(_numberSimulationsPerThread, _BP0Vector, _simulator, _PETRTM);
        connect(simSegment[jThread], SIGNAL(progresslieDetector(int)), this, SLOT(updateLieDetectorProgress(int)));
        connect(simSegment[jThread], SIGNAL(finishedlieDetector(dMatrix,dMatrix,dMatrix,dMatrix,dMatrix,double)),
                this,SLOT(finishedLieDetectorBPndOneThread(dMatrix,dMatrix,dMatrix,dMatrix,dMatrix,double)));
        QThreadPool::globalInstance()->start(simSegment[jThread]);
    }
}

void SimWindow::updateLieDetectorProgress(int iProgress)
{
    FUNC_ENTER;
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

void SimWindow::finishedLieDetectorBPndOneThread(dMatrix errBPnd, dMatrix errChall,
                                                 dMatrix tau2Ref, dMatrix errTau4, dMatrix errDV,
                                                 double sigma2)
{
    qInfo() << "finishedLieDetectorBPndOneThread enter";
    FUNC_ENTER;
    static int nThreadsFinished=0;
    static double sigma2Sum=0.;

    if ( errBPnd.size() != _BP0Vector.size() )
        qFatal("Error: mismatched vector sizes in SimWindow::finishedLieDetectorBPndOneThread");

    int nBP0Values = _BP0Vector.size();
    int nSamples   = errBPnd[0].size();

    QMutex mutex;
    mutex.lock();
    for (int jValue=0; jValue<nBP0Values; jValue++)
    {
        FUNC_INFO << "** jValue errtau4" << jValue << errTau4[jValue][0];
        for ( int jSample=0; jSample<nSamples; jSample++)
        {
            _errBPndMatrix[jValue].append(errBPnd[jValue][jSample]);
            _errChallMatrix[jValue].append(errChall[jValue][jSample]);
            _tau2RefMatrix[jValue].append(tau2Ref[jValue][jSample]);
            _errTau4Matrix[jValue].append(errTau4[jValue][jSample]);
            _errDVMatrix[jValue].append(errDV[jValue][jSample]);
        }
    }
    mutex.unlock();

    nThreadsFinished++;
    sigma2Sum += sigma2;
    if ( nThreadsFinished == _nThreads )
    {
        qInfo() << "all threads finished";
        _progressBar->reset();
        _progressBar->setMaximum(1);  // indicates that the bar is available
        sigma2Sum /= static_cast<double>(_nThreads);
        _sigma2Label->setText(QString("sigma2 = %1").arg(sigma2Sum));
        nThreadsFinished = 0;  // reset for next time
        sigma2Sum = 0.;
        finishedLieDetectorBPndAllThreads();
    }
    qInfo() << "finishedLieDetectorBPndOneThread exit";
}

void SimWindow::finishedLieDetectorTau4OneThread(dMatrix errBPnd, dMatrix errChall, dMatrix AIC, double sigma2)
{
    FUNC_ENTER;
    static int nThreadsFinished=0;
    static double sigma2Sum=0.;

    if ( errBPnd.size() != _tau4Vector.size() )
        qFatal("Error: mismatched vector sizes in SimWindow::finishedLieDetectorTau4OneThread");

    int nTau4Values = _tau4Vector.size();
    int nSamples    = errBPnd[0].size();

    QMutex mutex;
    mutex.lock();
    for (int jValue=0; jValue<nTau4Values; jValue++)
    {
        for ( int jSample=0; jSample<nSamples; jSample++)
        {
            _errBPndMatrix[jValue].append(errBPnd[jValue][jSample]);
            _errChallMatrix[jValue].append(errChall[jValue][jSample]);
            _AICMatrix[jValue].append(AIC[jValue][jSample]);
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
        finishedLieDetectorTau4AllThreads();
    }
}

double SimWindow::calculateMean(dVector vec)
{
    FUNC_ENTER;
    double mean = 0.;
    for (int j=0; j<vec.size(); j++)
        mean += vec[j];
    mean /= static_cast<double>(vec.size());
    return mean;
}
double SimWindow::calculateStDev(double mean, dVector vec)
{
    FUNC_ENTER;
    double stdev = 0.;
    for (int j=0; j<vec.size(); j++)
        stdev += SQR(vec[j] - mean);
    stdev /= static_cast<double>(vec.size()-1);
    stdev = qSqrt(stdev);
    return stdev;
}

void SimWindow::finishedLieDetectorBPndAllThreads()
{
    qInfo() << "finishedLieDetectorBPndAllThreads enter";
    FUNC_ENTER;
    int nBP0Values = _BP0Vector.size();
    int nSamples   = _nThreads * _numberSimulationsPerThread;

    qInfo() << "nSamples" << nSamples;

    dVector errBP, errChall, tau2Ref, errTau4, errDV, errBPSTDEV, errChallSTDEV, tau2RefSTDEV, errTau4STDEV, errDVSTDEV;
    // determine the means
    for (int jValue=0; jValue<nBP0Values; jValue++)
    {
        FUNC_INFO << "** errTau4Matrix" << _errTau4Matrix[jValue];

        qInfo() << "calc mean using" << _errBPndMatrix[jValue].size() << "samples for value" << jValue;
        errBP.append(calculateMean(_errBPndMatrix[jValue]));
        errChall.append(calculateMean(_errChallMatrix[jValue]));
        tau2Ref.append(calculateMean(_tau2RefMatrix[jValue]));
        errTau4.append(calculateMean(_errTau4Matrix[jValue]));
        errDV.append(calculateMean(_errDVMatrix[jValue]));

        errBPSTDEV.append(calculateStDev(errBP.last(),       _errBPndMatrix[jValue]));
        errChallSTDEV.append(calculateStDev(errChall.last(), _errChallMatrix[jValue]));
        tau2RefSTDEV.append(calculateStDev(tau2Ref.last(),   _tau2RefMatrix[jValue]));
        errTau4STDEV.append(calculateStDev(errTau4.last(),   _errTau4Matrix[jValue]));
        errDVSTDEV.append(calculateStDev(errDV.last(),       _errDVMatrix[jValue]));
    }

    // restore the value of BPnd in simulator, plus regenerate graphs
    changedBPND();  // this will use the value saved in the text field of _BPnd, which should have changed during the simulation

    QColor colors[10] = {Qt::black, Qt::red, Qt::blue, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};
    int nCurves = _plotErrBPndVsBPnd->getNumberCurves();
    int iColor  = nCurves%10;

    _plotErrBPndVsBPnd->addCurve(0,"error_BPnd");
    _plotErrBPndVsBPnd->setPointSize(5);
    _plotErrBPndVsBPnd->setColor(colors[iColor]);
    _plotErrBPndVsBPnd->setErrorBars(1);

    _plotErrChallVsBPnd->addCurve(0,"error_Challenge");
    _plotErrChallVsBPnd->setPointSize(5);
    _plotErrChallVsBPnd->setColor(colors[iColor]);
    _plotErrChallVsBPnd->setErrorBars(1);

    _plotTau2RefVsBPnd->addCurve(0,"Tau2Ref");
    _plotTau2RefVsBPnd->setPointSize(5);
    _plotTau2RefVsBPnd->setColor(colors[iColor]);
    _plotTau2RefVsBPnd->setErrorBars(1);
    if ( _PETRTM.isRTM2() )
        _plotTau2RefVsBPnd->setPointStyle(QCPScatterStyle::ssCross);

    _plotErrk4VsBPnd->addCurve(0,"1/k4");
    _plotErrk4VsBPnd->setPointSize(5);
    _plotErrk4VsBPnd->setColor(colors[iColor]);
    _plotErrk4VsBPnd->setErrorBars(1);

    _plotErrDVVsBPnd->addCurve(0,"DV");
    _plotErrDVVsBPnd->setPointSize(5);
    _plotErrDVVsBPnd->setColor(colors[iColor]);
    _plotErrDVVsBPnd->setErrorBars(1);

    dVector xVector;
    for (double BP0=_BPndLowValue; BP0<=_BPndHighValue; BP0 += _BPndStepValue)
        xVector.append(BP0);

    FUNC_INFO << "** set data" << _BP0Vector << errTau4;

    _plotErrBPndVsBPnd->setData(_BP0Vector,errBP,errBPSTDEV);
    _plotErrChallVsBPnd->setData(_BP0Vector,errChall,errChallSTDEV);
    if ( !_PETRTM.isRTM2() )
        _plotTau2RefVsBPnd->setData(_BP0Vector,tau2Ref,tau2RefSTDEV);
    if ( _PETRTM.getFitk4() )
        _plotErrk4VsBPnd->setData(_BP0Vector,errTau4,errTau4STDEV);
    else if ( _PETRTM.getFitDV() )
        _plotErrDVVsBPnd->setData(_BP0Vector,errDV,errDVSTDEV);

    _plotErrBPndVsBPnd->conclude(0,true);
    _plotErrChallVsBPnd->conclude(0,true);
    _plotTau2RefVsBPnd->conclude(0,true);
    _plotErrk4VsBPnd->conclude(0,true);
    _plotErrDVVsBPnd->conclude(0,true);

    _plotErrBPndVsBPnd->plotDataAndFit(true);
    _plotErrChallVsBPnd->plotDataAndFit(true);
    _plotTau2RefVsBPnd->plotDataAndFit(true);
    _plotErrk4VsBPnd->plotDataAndFit(true);
    _plotErrDVVsBPnd->plotDataAndFit(true);

    _progressBar->reset();
    qInfo() << "finishedLieDetectorBPndAllThreads exit";
}

void SimWindow::finishedLieDetectorTau4AllThreads()
{
    FUNC_ENTER;
    int nTau4Values = _tau4Vector.size();
    int nSamples    = _nThreads * _numberSimulationsPerThread;

    FUNC_INFO << "nSamples" << nSamples;

    dVector errBP, errChall, AIC, errBPSTDEV, errChallSTDEV, AICSTDEV;
    // determine the means
    for (int jValue=0; jValue<nTau4Values; jValue++)
    {
        FUNC_INFO << "** errTau4Matrix" << _errTau4Matrix[jValue];

        errBP.append(calculateMean(_errBPndMatrix[jValue]));
        errChall.append(calculateMean(_errChallMatrix[jValue]));
        AIC.append(calculateMean(_AICMatrix[jValue]));

        errBPSTDEV.append(calculateStDev(errBP[jValue],       _errBPndMatrix[jValue]));
        errChallSTDEV.append(calculateStDev(errChall[jValue], _errChallMatrix[jValue]));
        AICSTDEV.append(calculateStDev(AIC[jValue],           _AICMatrix[jValue]));
    }

    // restore the value of BPnd in simulator, plus regenerate graphs
    changedBPND();  // this will use the value saved in the text field of _BPnd, which should have changed during the simulation

    QColor colors[10] = {Qt::black, Qt::red, Qt::blue, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};
    int nCurves = _plotErrBPndVsBPnd->getNumberCurves();
    int iColor  = nCurves%10;

    _plotErrBPndVsTau4->addCurve(0,"error_BPnd");
    _plotErrBPndVsTau4->setPointSize(5);
    _plotErrBPndVsTau4->setColor(colors[iColor]);
    _plotErrBPndVsTau4->setErrorBars(1);

    _plotErrChallVsTau4->addCurve(0,"error_Challenge");
    _plotErrChallVsTau4->setPointSize(5);
    _plotErrChallVsTau4->setColor(colors[iColor]);
    _plotErrChallVsTau4->setErrorBars(1);

    _plotAICVsTau4->addCurve(0,"AIC");
    _plotAICVsTau4->setPointSize(5);
    _plotAICVsTau4->setColor(colors[iColor]);
    _plotAICVsTau4->setErrorBars(1);

    dVector xVector;
    for (double tau4=0; tau4<=_tau4HighValue; tau4 += _tau4StepValue)
        xVector.append(tau4);

    _plotErrBPndVsTau4->setData(_tau4Vector,errBP,errBPSTDEV);
    _plotErrChallVsTau4->setData(_tau4Vector,errChall,errChallSTDEV);
    _plotAICVsTau4->setData(_tau4Vector,AIC,AICSTDEV);

    _plotErrBPndVsTau4->conclude(0,true);
    _plotErrChallVsTau4->conclude(0,true);
    _plotAICVsTau4->conclude(0,true);

    _plotErrBPndVsTau4->plotDataAndFit(true);
    _plotErrChallVsTau4->plotDataAndFit(true);
    _plotAICVsTau4->plotDataAndFit(true);

    _progressBar->reset();
}

void SimWindow::calculateForwardTau4andTau2Ref()
{  // move this function into PETRTM in order to use it with simulations
    FUNC_ENTER;
    if ( !_PETRTM.isRTM2() || !_PETRTM.isForwardModel() ) return;

    _globalFit.tau4Vector.clear();  _globalFit.tau2RefVector.clear();
    _globalFit.DVVector.clear();    _globalFit.AICVector1D.clear();

    // k4
    dPoint3D searchParTau = _PETRTM.getSearchParTau4();
    double tau4Min  = searchParTau.x;
    double tau4Max  = searchParTau.y;
    double tau4Step = searchParTau.z;
    double tau4 = tau4Min;
    while ( tau4 < tau4Max )
    {
        _globalFit.tau4Vector.append(tau4);
        tau4 += tau4Step;
    }

    // k2'
    dPoint3D searchParTau2Ref = _PETRTM.getSearchParTau2Ref();
    double tau2RefMin  = searchParTau2Ref.x;
    double tau2RefMax  = searchParTau2Ref.y;
    double tau2RefStep = searchParTau2Ref.z;
    double tau2Ref = tau2RefMin;
    while ( tau2Ref < tau2RefMax )
    {
        _globalFit.tau2RefVector.append(tau2Ref);
        tau2Ref += tau2RefStep;
    }

    // DV
    double DV = _PETRTM.getDVInRun(0);
    _globalFit.DVVector.append(DV);

    dVector ROI = _simulator.getCtCoarse();

    int nTau4    = _globalFit.tau4Vector.size();
    int nTau2Ref = _globalFit.tau2RefVector.size();
    int nDV      = _globalFit.DVVector.size();

    _globalFit.AICMatrix.resize(nTau4);
    for (int jTau4=0; jTau4<nTau4; jTau4++)
    {
        _globalFit.AICMatrix[jTau4].resize(nTau2Ref);
        for (int jTau2Ref=0; jTau2Ref<nTau2Ref; jTau2Ref++)
            _globalFit.AICMatrix[jTau4][jTau2Ref].fill(1.,nDV);
    }

    QVector<globalk4Fitter *> fitter;
    fitter.resize(nTau4);
    QThreadPool::globalInstance()->setMaxThreadCount(_nThreads);
    globalk4FitterFinished(0);
    qRegisterMetaType<dPoint4D>("dPoint4D");
    for (int jTau4=0; jTau4<nTau4; jTau4++)
    {
        FUNC_INFO << "jTau4" << _globalFit.tau4Vector[jTau4];
        fitter[jTau4] = new globalk4Fitter(_globalFit.tau4Vector[jTau4], _globalFit.tau2RefVector, _globalFit.DVVector,
                                           ROI, _PETRTM);
        connect(fitter[jTau4], SIGNAL(update(dPoint4D)), this, SLOT(globalk4FitterUpdate(dPoint4D)));
        connect(fitter[jTau4], SIGNAL(finished(int)),    this, SLOT(globalk4FitterFinished(int)));
        QThreadPool::globalInstance()->start(fitter[jTau4]);
    }
    FUNC_EXIT;
}

void SimWindow::globalk4FitterUpdate(dPoint4D tau4Tau2RefDVAIC)
{   // tau4Tau2RefDVAIC: (x,y,z,t) = (tau4, tau2Ref, DV, AIC)
    FUNC_ENTER;

    // loop over all points and put the AIC for the (tau4,tau2Ref) point in the correct position in matrix _globalFit.AICMatrix
    // find the position in the AICMatrix
    int iTau4=0;
    for (int jTau4=0; jTau4<_globalFit.tau4Vector.size(); jTau4++)
        if ( _globalFit.tau4Vector[jTau4] == tau4Tau2RefDVAIC.x )
        {
            iTau4 = jTau4; break;
        }
    if ( iTau4 >= _globalFit.tau4Vector.size() ) qFatal("Programming error 1 in timePage::globalk4FitterUpdate");

    int iTau2Ref=0;
    for (int jTau2Ref=0; jTau2Ref<_globalFit.tau2RefVector.size(); jTau2Ref++)
        if ( _globalFit.tau2RefVector[jTau2Ref] == tau4Tau2RefDVAIC.y )
        {
            iTau2Ref = jTau2Ref; break;
        }
    if ( iTau2Ref >= _globalFit.tau2RefVector.size() ) qFatal("Programming error 2 in timePage::globalk4FitterUpdate");

    int iDV=0;
    for (int jDV=0; jDV<_globalFit.DVVector.size(); jDV++)
        if ( _globalFit.DVVector[jDV] == tau4Tau2RefDVAIC.z )
        {
            iDV = jDV; break;
        }
    if ( iDV >= _globalFit.DVVector.size() ) qFatal("Programming error 3 in timePage::globalk4FitterUpdate");

    double AIC = tau4Tau2RefDVAIC.t;
    if ( qIsNaN(AIC) ) qInfo() << "nan for triplet"
                               << _globalFit.tau4Vector[iTau4] << _globalFit.tau2RefVector[iTau2Ref] << _globalFit.DVVector[iTau2Ref];
    if ( qIsInf(AIC) ) qInfo() << "inf for pair"
                               << _globalFit.tau4Vector[iTau4] << _globalFit.tau2RefVector[iTau2Ref] << _globalFit.DVVector[iTau2Ref];
    if ( qIsNaN(AIC) || qIsInf(AIC) ) AIC = 0;
//    qInfo() << "update" << iTau4 << iTau2Ref << "with" << AIC;
    // Put the value into the global AIC matrix
    _globalFit.AICMatrix[iTau4][iTau2Ref][iDV] = AIC;
    // qInfo() << "update" << iTau4 << iTau2Ref << iDV << "with" << AIC;
}

void SimWindow::globalk4FitterFinished(int mode)
{
    static int nThreadsFinished = 0;
    FUNC_ENTER << mode;
    if ( mode <= 0 )
    {
        nThreadsFinished = 0;
        _progressBar->setMaximum(_globalFit.tau4Vector.size());
        _runSimulationAndAnalysis->setEnabled(false);
        FUNC_INFO << "initialize progress to" << _globalFit.tau4Vector.size();
    }
    else
    {
        nThreadsFinished++;
        _progressBar->setValue(qMin(nThreadsFinished,_progressBar->maximum()));
        FUNC_INFO << "progress bar value" << _progressBar->value() << "of" << _progressBar->maximum();
        if ( nThreadsFinished == _globalFit.tau4Vector.size() )
        { // if all threads are finished, clean up and emit and all-finished signal
            nThreadsFinished = 0;  // reset for next time
            _runSimulationAndAnalysis->setEnabled(true);

            int nTau4    = _globalFit.tau4Vector.size();
            int nTau2Ref = _globalFit.tau2RefVector.size();
            int nDV      = _globalFit.DVVector.size();
            // find the best AIC per k4 and also sum over the 2d matrix
            _globalFit.AICVector1D.fill(1.e40,nTau4);
            double AICMax = -1.e10;     double bestTau4=0.;
            double AICMin = 1.e100;     double bestTau2Ref=0.;
            for (int jTau4=0; jTau4<nTau4; jTau4++)
            {
                for (int jTau2Ref=0; jTau2Ref<nTau2Ref; jTau2Ref++)
                {
                    for (int jDV=0; jDV<nDV; jDV++)
                    {
                        double AIC = _globalFit.AICMatrix[jTau4][jTau2Ref][jDV];
                        // qInfo() << "AICMatrix" << jTau4 << jTau2Ref << jDV << AIC;
                        if ( AIC != 0. )
                        {
                            _globalFit.AICVector1D[jTau4] = qMin(_globalFit.AICVector1D[jTau4],AIC);
                            if ( AIC < AICMin )
                            {
                                AICMin = AIC;
                                bestTau4    = _globalFit.tau4Vector[jTau4];
                                bestTau2Ref = _globalFit.tau2RefVector[jTau2Ref];
                            }
                        }
                        AICMax = qMax(AICMax,AIC);
//                        qInfo() << "\n" << _globalFit.tau4Vector[jTau4] << _globalFit.tau2RefVector[jTau2Ref] << _globalFit.DVVector[jDV];
//                        qInfo() << AIC << _globalFit.AICVector1D[jTau4];
                    }
                }
                // qInfo() << "jTau4, AIC1d" << jTau4 << _globalFit.AICVector1D[jTau4];
            }
            for (int jTau4=0; jTau4<nTau4; jTau4++)
            {
                for (int jTau2Ref=0; jTau2Ref<nTau2Ref; jTau2Ref++)
                {
                    for (int jDV=0; jDV<nDV; jDV++)
                        if ( _globalFit.AICVector1D[jTau4] == 0. ) _globalFit.AICVector1D[jTau4] = AICMax;
                }
            }

            // qInfo() << "AICVector1D" << _globalFit.AICVector1D;

            // make a 2D image of AIC, and determine best Tau4 and Tau2Ref
            ImageIO AICImage;
            dPoint3D res, origin;  res.z = 1.;  origin.x = origin.y = origin.z = 0.;
            res.x = _globalFit.tau4Vector[1]    - _globalFit.tau4Vector[0];
            res.y = _globalFit.tau2RefVector[1] - _globalFit.tau2RefVector[0];
            AICImage.setDimensions(nTau4, nTau2Ref, 1, 1);
            AICImage.setResolution(res);
            AICImage.setOrigin(origin);
            double sigma2Tau4    = SQR(4.);
            double sigma2Tau2Ref = SQR(2.);
            double bestTau4Weighted=0.;  double bestTau2RefWeighted=0.;   double bestDV=0.;  double probSum=0.;
            for (int jTau4=0; jTau4<nTau4; jTau4++)
            {
                double tau4 = _globalFit.tau4Vector[jTau4];
                for (int jTau2Ref=0; jTau2Ref<nTau2Ref; jTau2Ref++)
                {
                    double tau2Ref = _globalFit.tau2RefVector[jTau2Ref];
                    for (int jDV=0; jDV<nDV; jDV++)
                    {
                        double DV = _globalFit.DVVector[jDV];
                        double AIC  = _globalFit.AICMatrix[jTau4][jTau2Ref][jDV];
                        double prob = qExp(AICMin-AIC)
                                *     qExp(-SQR(tau4-bestTau4)/sigma2Tau4)      // try to avoid 2-peak problem
                                *     qExp(-SQR(tau4-bestTau4)/sigma2Tau2Ref);  // try to avoid 2-peak problem
                        if ( AIC == 0. ) prob = 0.;   // convergence failure
                        bestTau4Weighted    += prob * tau4;
                        bestTau2RefWeighted += prob * tau2Ref;
                        // if ( prob > 0.01 ) qInfo() << "prob, DV, AIC, tau4, tau2Ref" << prob << DV << AIC << tau4 << tau2Ref;
                        bestDV      += prob * DV;
                        probSum += prob;
                        AICImage.setStackValue(jTau4, jTau2Ref, 0, 0, AIC);
                    }
                }
            }
            bestTau4Weighted    /= probSum;
            bestTau2RefWeighted /= probSum;
            bestDV /= probSum;
            FUNC_INFO << "best set" << bestTau4Weighted << bestTau2RefWeighted << bestDV;
            _tau4Analysis->setText(QString::number(bestTau4Weighted,'g',4));
            _tau2RefAnalysis->setText(QString::number(bestTau2RefWeighted));
            _PETRTM.setTau4(0,bestTau4Weighted);
            _PETRTM.setTau2Ref(0,bestTau2RefWeighted);
            _PETRTM.setDV(0,bestDV);

            if ( nDV > 1 )
                qInfo() << "Over 3D grid search, best (1/k4, 1/k2', DV) = (" << bestTau4Weighted << "," << bestTau2RefWeighted << "," << bestDV << ")";
            else
                qInfo() << "Over 2D grid search, best (1/k4, 1/k2') = (" << bestTau4Weighted << "," << bestTau2RefWeighted << ")";

            // write the image
            AICImage.setFileName("aic");
            AICImage.writeNiftiFile(".",true);


            _progressBar->reset();
            _progressBar->setMaximum(1);  // indicates that the bar is available

            updateAllGraphs();
            showAICPlot();
        } // all done (all threads)
    }
    FUNC_EXIT;
}


globalk4Fitter::globalk4Fitter(double tau4, dVector tau2RefVector, dVector DVVector, dVector ROI, const PETRTM petRTM)
{
    FUNC_ENTER;
    _tau4    = tau4;
    _tau2RefVector = tau2RefVector;
    _DVVector = DVVector;
    _ROI     = ROI;
    _petRTM  = petRTM;
    FUNC_EXIT;
}

void globalk4Fitter::run()
{
    FUNC_ENTER << _petRTM.isForwardModel() << _petRTM.isRTM2();
    if ( !_petRTM.isForwardModel() || !_petRTM.isRTM2() ) return;

    _petRTM.setTau4(0,_tau4);
    for (int jTau2Ref=0; jTau2Ref<_tau2RefVector.size(); jTau2Ref++)
    {
        _petRTM.setTau2Ref(0, _tau2RefVector[jTau2Ref]);
        for (int jDV=0; jDV<_DVVector.size(); jDV++)
        {
            _petRTM.setDV(0, _DVVector[jDV]);
            double globalAIC = 0.;
            dMatrix tissueVector;   tissueVector.resize(1);
            tissueVector[0] = _ROI;
            _petRTM.setTissueRegion(tissueVector);
            _petRTM.prepare();
            _petRTM.fitByForwardModel(true);
            globalAIC += _petRTM._simulator[0].getAIC();
            FUNC_INFO << "globalAIC" << globalAIC << "for" << _tau4 << _tau2RefVector[jTau2Ref];
            dPoint4D tau4Tau2RefDVAIC = {_tau4, _tau2RefVector[jTau2Ref], _DVVector[jDV], globalAIC};
            // qInfo() << "tau4Tau2RefDVAIC" << tau4Tau2RefDVAIC.x << tau4Tau2RefDVAIC.y << tau4Tau2RefDVAIC.z << tau4Tau2RefDVAIC.t;
            emit update(tau4Tau2RefDVAIC);
        }
    }
    emit finished(1);
    FUNC_EXIT;
}
