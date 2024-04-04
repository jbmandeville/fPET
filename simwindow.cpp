#include <QtWidgets>
#include <QDebug>
#include <QVector>
#include <QObject>
#include <QFileDialog>
#include <QString>
#include "simwindow.h"

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
    createPlasmaPage();
    createTargetPage();
    _tabTimeSpace->addTab(_setupPage, tr("Plasma"));
    _tabTimeSpace->addTab(_targetPage, tr("Target Region"));

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
            disconnect(_dataTargetRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(changedDataTargetRegion()));
            _dataTargetRegion->clear();
            for (int jColumn=0; jColumn<_dataColumnNames.count(); jColumn++)
            {
                _dataTargetRegion->addItem(_dataColumnNames.at(jColumn));
            }
            // set the default target region to the first region after the bin info
            int iColumnFirstRegion = qMin(_dataColumnNames.count()-1,qMax(iColumnBinSize,iColumnBinTime)+1);
            _dataTargetRegion->setCurrentIndex(iColumnFirstRegion);
            connect(_dataTargetRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(changedDataTargetRegion()));

            _ROIFileName->setToolTip(fileName);
            _ROIFileName->setText(utilString::getFileNameWithoutDirectory(fileName));
            // Resample the simulation to match the binning of the real data
            _simulator.setDurationBins(_dtBinsSec);
            QString text;
            _numberTimeBins->setText(text.setNum(_dtBinsSec.size())); changedNumberBins();
            // update the bin duration
            int iBin = _binIndex->value() - 1;
            _binDuration->setText(text.setNum(_simulator.getDurationPerBinSec(iBin)));  // this will update graphs
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

void SimWindow::createPlasmaPage()
{
    FUNC_ENTER;
    _setupPage = new QWidget();
    _plotPlasma  = new plotData(0);
    _plotTarget1 = new plotData(1);
    _plotTarget1->getPlotSurface()->setVisible(false);

    //////// The setup page status bar
    _plasmaStatusBar = new QStatusBar;  // must be global so it doesn't go out of scope
    _plasmaStatusBar->setStyleSheet("color:blue");
    _plotPlasma->setQCStatusBar(_plasmaStatusBar);

    auto *setupPlotLayout = new QVBoxLayout();
    setupPlotLayout->addWidget(_plotTarget1->getPlotSurface());
    setupPlotLayout->addWidget(_plotPlasma->getPlotSurface());
    setupPlotLayout->addWidget(_plasmaStatusBar);
    setupPlotLayout->setStretch(0,20);
    setupPlotLayout->setStretch(1,20);
    setupPlotLayout->setStretch(2,1);
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

    auto *setupLayout = new QGridLayout();
    setupLayout->addWidget(numberTimeBinsLabel,0,0);
    setupLayout->addWidget(_numberTimeBins,0,1);
    setupLayout->addWidget(binIndexLabel,1,0);
    setupLayout->addWidget(_binIndex,1,1);
    setupLayout->addWidget(binDurationLabel,2,0);
    setupLayout->addWidget(_binDuration,2,1);
    setupLayout->addWidget(subSampleLabel,3,0);
    setupLayout->addWidget(_subSample,3,1);
    setupGroupBox->setLayout(setupLayout);
    setupLayout->setSpacing(0);
    connect(_numberTimeBins, SIGNAL(editingFinished()), this, SLOT(changedNumberBins()));
    connect(_binDuration, SIGNAL(editingFinished()), this, SLOT(changedBinDuration()));
    connect(_subSample,   SIGNAL(editingFinished()), this, SLOT(changedSubSample()));
    connect(_binIndex, SIGNAL(valueChanged(int)), this, SLOT(changedBinIndex(int)));

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
    plasmaInGroupBox->setLayout(plasmaInLayout);
    plasmaInLayout->setSpacing(0);
    connect(_bolusMag, SIGNAL(editingFinished()), this, SLOT(changedBolusMag()));
    connect(_tauDecay, SIGNAL(editingFinished()), this, SLOT(changedTauBolus()));
    connect(_KBol, SIGNAL(editingFinished()), this, SLOT(changedInfusion()));
    connect(_KBolDelay, SIGNAL(editingFinished()), this, SLOT(changedInfusionDelay()));
    connect(_fastTau, SIGNAL(editingFinished()), this, SLOT(changedFastElimination()));
    connect(_slowTau, SIGNAL(editingFinished()), this, SLOT(changedSlowElimination()));
    connect(_fastFraction, SIGNAL(editingFinished()), this, SLOT(changedFastEliminationFraction()));

    QPixmap pixmapOpenAlignment(":/My-Icons/openFile.png");
    QIcon openIcon(pixmapOpenAlignment);
    _readROIFile = new QPushButton(openIcon,"open",_setupPage);
    _ROIFileName = new QLabel("No real data",_setupPage);
    auto *realDataGroupBox = new QGroupBox("Real data (imported ROIs from table)");

    auto *realDataLayout = new QGridLayout();
    realDataLayout->addWidget(_readROIFile,0,0);
    realDataLayout->addWidget(_ROIFileName,0,1);

    auto *realDataAndStartingLayout = new QVBoxLayout();
    realDataAndStartingLayout->addLayout(realDataLayout);
    realDataGroupBox->setLayout(realDataAndStartingLayout);

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(setupGroupBox);
    rightLayout->addWidget(plasmaInGroupBox);
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

    auto *checkBoxShowTarget  = new QCheckBox("Show Target",_setupPage);
    connect(checkBoxShowTarget,  SIGNAL(clicked(bool)), this, SLOT(showTarget(bool)));
    auto *showWidget = new QWidget();
    auto *showHBoxLayout = new QHBoxLayout();
    showHBoxLayout->addWidget(checkBoxShowTarget);
    showWidget->setLayout(showHBoxLayout);
    showWidget->setStyleSheet("color:Darkred");

    QToolBar *graphToolBar = new QToolBar("graph tool bar");
    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addWidget(showWidget);
    graphToolBar->addAction(rescaleXYAction);

    fullLayout->setMenuBar(graphToolBar);

    // Make toolbar connections
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotPlasma, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotPlasma, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotPlasma, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotPlasma, SLOT(setSelectPoints(bool)));

    connect(dragXAction,     SIGNAL(triggered(bool)), _plotTarget1, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotTarget1, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotTarget1, SLOT(autoScale(bool)));

//    connect(crossCursorAct,  SIGNAL(triggered(bool)), _plotPlasma, SLOT(setSelectPoints()));

    connect(_readROIFile, SIGNAL(pressed()), this, SLOT(getTableDataFile()));

//    _plotPlasma->getPlotSurface()->setVisible(true);
//    _plasmaStatusBar->setVisible(true);

    _setupPage->setLayout(fullLayout);
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

void SimWindow::changedWeightType(int indexInBox)
{
    FUNC_ENTER;
    updateAllGraphs();
}

void SimWindow::changedModelType(int indexInBox)
{
    FUNC_ENTER;
}
void SimWindow::changedBaselineTerms(int indexInBox)
{
    FUNC_ENTER;
    int nBase = indexInBox + 2;
    _fMRIGLM.setBaselineTerms(nBase);
    updateAllGraphs();
}

void SimWindow::createTargetPage()
{
    FUNC_ENTER;
    _targetPage = new QWidget();

    _TRStatusBar = new QStatusBar;  // must be global so it doesn't go out of scope
    _TRStatusBar->setStyleSheet("color:blue");

    _plotBasis  = new plotData(2);
    _plotTarget = new plotData(3);
    _plotTarget->setQCStatusBar(_TRStatusBar);
    _plotBasis->getPlotSurface()->setVisible(false);

    auto *plotTargetLayout = new QVBoxLayout();
    plotTargetLayout->addWidget(_plotTarget->getPlotSurface());
    plotTargetLayout->addWidget(_plotBasis->getPlotSurface());
    plotTargetLayout->addWidget(_TRStatusBar);
    plotTargetLayout->setStretch(0,20);
    plotTargetLayout->setStretch(1,20);
    plotTargetLayout->setStretch(2,1);

    QString numberString;
    int editTextSize=80;

    // target region: input
    _targetSimulationGroupBox = new QGroupBox("Target Simulation: truth");
    QLabel *K1Label = new QLabel("K1");
    QLabel *k2Label = new QLabel("k2");
    QLabel *k3Label = new QLabel("k3");
    QLabel *challengeOnsetTimeLabel1  = new QLabel("Challenge #1 onset");
    QLabel *challengeOffsetTimeLabel1 = new QLabel("Challenge #1 offset");
    QLabel *challengeOnsetTimeLabel2  = new QLabel("Challenge #2 onset");
    QLabel *challengeOffsetTimeLabel2 = new QLabel("Challenge #2 offset");
    QChar delta = QChar(0x0394);
    QLabel *plasmaFraclabel    = new QLabel("Plasma %");
    QLabel *noiseLabel = new QLabel("Noise");
    auto *challengeMagLabel1 = new QLabel(QString("challenge mag (%1k3 in %)").arg(delta));
    auto *challengeMagLabel2 = new QLabel(QString("challenge mag (%1k3 in %)").arg(delta));
    _K1            = new QLineEdit();
    _k2            = new QLineEdit();
    _k3            = new QLineEdit();
    _challengeOnsetTime1  = new QLineEdit();
    _challengeOffsetTime1 = new QLineEdit();
    _challengeOnsetTime2  = new QLineEdit();
    _challengeOffsetTime2 = new QLineEdit();
    _challengeMagPercent1 = new QLineEdit();
    _challengeMagPercent2 = new QLineEdit();
    _plasmaFracTar = new QLineEdit();
    _noiseTar      = new QLineEdit();
    _K1->setText(numberString.setNum(_simulator.getK1()));
    _k2->setText(numberString.setNum(_simulator.getk2()));
    _k3->setText(numberString.setNum(_simulator.getk3()));
    _challengeMagPercent1->setText(numberString.setNum(_simulator.getChallengeMagPercent1()));
    _challengeMagPercent2->setText(numberString.setNum(_simulator.getChallengeMagPercent2()));
    _challengeOnsetTime1->setText(numberString.setNum(_simulator.getChallengeOnsetTime1()));
    _challengeOffsetTime1->setText(numberString.setNum(_simulator.getChallengeOffsetTime1()));
    _challengeOnsetTime2->setText(numberString.setNum(_simulator.getChallengeOnsetTime2()));
    _challengeOffsetTime2->setText(numberString.setNum(_simulator.getChallengeOffsetTime2()));
    _noiseTar->setText(numberString.setNum(_simulator.getNoiseTar()));
    _plasmaFracTar->setText(numberString.setNum(_simulator.getPlasmaPercentTar()));
    _K1->setFixedWidth(editTextSize);
    _k2->setFixedWidth(editTextSize);
    _k3->setFixedWidth(editTextSize);
    _challengeOnsetTime1->setFixedWidth(editTextSize);
    _challengeOffsetTime1->setFixedWidth(editTextSize);
    _challengeOnsetTime2->setFixedWidth(editTextSize);
    _challengeOffsetTime2->setFixedWidth(editTextSize);
    _challengeMagPercent1->setFixedWidth(editTextSize);
    _challengeMagPercent2->setFixedWidth(editTextSize);
    _noiseTar->setFixedWidth(editTextSize);
    _plasmaFracTar->setFixedWidth(editTextSize);

    auto *targetSimulationLayout = new QGridLayout();
    targetSimulationLayout->addWidget(K1Label,0,0);
    targetSimulationLayout->addWidget(_K1,0,1);
    targetSimulationLayout->addWidget(k2Label,1,0);
    targetSimulationLayout->addWidget(_k2,1,1);
    targetSimulationLayout->addWidget(k3Label,2,0);
    targetSimulationLayout->addWidget(_k3,2,1);

    targetSimulationLayout->addWidget(challengeOnsetTimeLabel1,3,0);
    targetSimulationLayout->addWidget(_challengeOnsetTime1,3,1);
    targetSimulationLayout->addWidget(challengeOffsetTimeLabel1,4,0);
    targetSimulationLayout->addWidget(_challengeOffsetTime1,4,1);
    targetSimulationLayout->addWidget(challengeMagLabel1,5,0);
    targetSimulationLayout->addWidget(_challengeMagPercent1,5,1);

    targetSimulationLayout->addWidget(challengeOnsetTimeLabel2,6,0);
    targetSimulationLayout->addWidget(_challengeOnsetTime2,6,1);
    targetSimulationLayout->addWidget(challengeOffsetTimeLabel2,7,0);
    targetSimulationLayout->addWidget(_challengeOffsetTime2,7,1);
    targetSimulationLayout->addWidget(challengeMagLabel2,8,0);
    targetSimulationLayout->addWidget(_challengeMagPercent2,8,1);

    targetSimulationLayout->addWidget(plasmaFraclabel,9,0);
    targetSimulationLayout->addWidget(_plasmaFracTar,9,1);
    targetSimulationLayout->addWidget(noiseLabel,10,0);
    targetSimulationLayout->addWidget(_noiseTar,10,1);

    _targetSimulationGroupBox->setLayout(targetSimulationLayout);
    targetSimulationLayout->setSpacing(0);
    connect(_K1,          SIGNAL(editingFinished()), this, SLOT(changedK1()));
    connect(_k2,          SIGNAL(editingFinished()), this, SLOT(changedk2()));
    connect(_k3,          SIGNAL(editingFinished()), this, SLOT(changedk3()));
    connect(_challengeOnsetTime1, SIGNAL(editingFinished()), this, SLOT(changedChallengeOnsetTime1()));
    connect(_challengeOffsetTime1,SIGNAL(editingFinished()), this, SLOT(changedChallengeOffsetTime1()));
    connect(_challengeOnsetTime2, SIGNAL(editingFinished()), this, SLOT(changedChallengeOnsetTime2()));
    connect(_challengeOffsetTime2,SIGNAL(editingFinished()), this, SLOT(changedChallengeOffsetTime2()));
    connect(_challengeMagPercent1,  SIGNAL(editingFinished()), this, SLOT(changedChallengeMag1()));
    connect(_challengeMagPercent2,  SIGNAL(editingFinished()), this, SLOT(changedChallengeMag2()));
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
    QLabel *modelLabel    = new QLabel("analysis type");
    QLabel *baselineLabel = new QLabel("baseline terms");
    QLabel *weightLabel   = new QLabel("Weighting scheme");
    QLabel *ignoreLabel   = new QLabel("Ignore points");
    QLabel *tauLabel      = new QLabel("tau (response)");
    _modelType = new QComboBox();
    _modelType->addItem("GLM");
    _modelType->addItem("forward");
    _modelType->setCurrentIndex(0);
    _baselineTerms = new QComboBox();
    _baselineTerms->addItem("2");
    _baselineTerms->addItem("3");
    _baselineTerms->addItem("4");
    _baselineTerms->addItem("5");
    _baselineTerms->setCurrentIndex(1);
    _weightType  = new QComboBox();
    _weightType->addItem("Uniform weights");
    _weightType->addItem("11C-Noiseless");
    _weightType->addItem("11C");
    _weightType->addItem("Custom");
    _weightType->setToolTip("Set a weighting scheme for WLS");
    _weightType->setEnabled(false);
    _tauAnalysis  = new QLineEdit();
    _tauAnalysis->setText("0");
    _ignoreString = new QLineEdit();
    _ignoreString->setText("");
    connect(_modelType,    SIGNAL(currentIndexChanged(int)), this, SLOT(changedModelType(int)));
    connect(_baselineTerms,SIGNAL(currentIndexChanged(int)), this, SLOT(changedBaselineTerms(int)));
    connect(_weightType,   SIGNAL(currentIndexChanged(int)), this, SLOT(changedWeightType(int)));
    auto *fitLabel = new QLabel("fit challenge?");
    _fitChallenge = new QCheckBox("fit");
    _fitChallenge->setChecked(true);
    _removeBaseline = new QCheckBox("baseline?");
    _removeBaseline->setChecked(false);
    _signalRaw = new QRadioButton("accumulation");
    _signalDerivative = new QRadioButton("derivative");
    _signalRaw->setChecked(true);

    auto *showLabel = new QLabel("show:");

    auto *analysisLayout = new QGridLayout();
    analysisLayout->addWidget(modelLabel,0,0);
    analysisLayout->addWidget(_modelType,0,1);

    analysisLayout->addWidget(baselineLabel,1,0);
    analysisLayout->addWidget(_baselineTerms,1,1);

    analysisLayout->addWidget(tauLabel,2,0);
    analysisLayout->addWidget(_tauAnalysis,2,1);

    analysisLayout->addWidget(weightLabel,3,0);
    analysisLayout->addWidget(_weightType,3,1);

    analysisLayout->addWidget(ignoreLabel,4,0);
    analysisLayout->addWidget(_ignoreString,4,1);

    analysisLayout->addWidget(fitLabel,5,0);
    analysisLayout->addWidget(_fitChallenge,5,1);

    analysisLayout->addWidget(showLabel,6,0);

    analysisLayout->addWidget(_signalRaw,7,0);
    analysisLayout->addWidget(_signalDerivative,7,1);

    analysisLayout->addWidget(_removeBaseline,8,1);

    _analysisGroupBox->setLayout(analysisLayout);
    analysisLayout->setSpacing(0);
    connect(_tauAnalysis,    SIGNAL(editingFinished()), this, SLOT(changedTauAnalysis()));
    connect(_ignoreString,   SIGNAL(editingFinished()), this, SLOT(changedIgnoreString()));
    connect(_fitChallenge,   SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxChallenge(bool)));
    connect(_removeBaseline, SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxBaseline(bool)));
    connect(_signalDerivative, SIGNAL(toggled(bool)),   this, SLOT(changedDerivativeRadioButton(bool)));
    // target region: errors
    auto *errorGroupBox  = new QGroupBox("Target Region: Percentage Errors");
    QLabel *errork3Label = new QLabel("k3 ");
    QLabel *errork2Label = new QLabel("k2   ");
    _errorK1Label        = new QLabel("K1   ");
    _errorChallengeLabel = new QLabel("Challenge");
    _sigmaLabel          = new QLabel("sigma");
    _errork3             = new QLabel();
    _errork2             = new QLabel();
    _errorK1             = new QLabel();
    _errorChallenge      = new QLabel();
    _sigma               = new QLabel();
    _errorChallengeLabel->setVisible(false);
    _errorChallenge->setVisible(false);
    auto *errorLayout = new QGridLayout();
    errorLayout->addWidget(errork3Label,0,0);
    errorLayout->addWidget(_errork3,0,1);
    errorLayout->addWidget(errork2Label,1,0);
    errorLayout->addWidget(_errork2,1,1);
    errorLayout->addWidget(_errorK1Label,2,0);
    errorLayout->addWidget(_errorK1,2,1);
    errorLayout->addWidget(_errorChallengeLabel,5,0);
    errorLayout->addWidget(_errorChallenge,5,1);
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

    QLabel *analyzeLabel = new QLabel("analyze:",_targetPage);
    _analyzeSimulation = new QRadioButton("Simulation(s)",_targetPage);
    _analyzeRealData = new QRadioButton("ROI Data",_targetPage);
    QButtonGroup *analyzeGroup = new QButtonGroup(_targetPage);
    analyzeGroup->addButton(_analyzeSimulation);
    analyzeGroup->addButton(_analyzeRealData);
    _analyzeSimulation->setChecked(true);

    auto *checkBoxShowBasis  = new QCheckBox("Show Basis",_targetPage);
    auto *showWidget = new QWidget();
    auto *showHBoxLayout = new QHBoxLayout();
    showHBoxLayout->addWidget(checkBoxShowBasis);
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
    connect(checkBoxShowBasis,  SIGNAL(clicked(bool)), this, SLOT(showBasis(bool)));
    connect(_analyzeSimulation, SIGNAL(clicked(bool)), this, SLOT(clickedAnalyzeStimulation(bool)));
    connect(_analyzeRealData,   SIGNAL(clicked(bool)), this, SLOT(clickedAnalyzeRealData(bool)));

    _targetPage->setLayout(fullLayout);
}

void SimWindow::clickedAnalyzeStimulation(bool state)
{
    FUNC_ENTER;
    updateAllGraphs();
}
void SimWindow::clickedAnalyzeRealData(bool state)
{
    FUNC_ENTER;
    updateAllGraphs();
}

void SimWindow::showBasis(bool state)
{
    FUNC_ENTER;
    _plotBasis->getPlotSurface()->setVisible(state);
}
void SimWindow::showTarget(bool state)
{
    FUNC_ENTER;
    _plotTarget1->getPlotSurface()->setVisible(state);
}

void SimWindow::updatePlasmaGraph()
{
    FUNC_ENTER;
    // update the plot
    _plotPlasma->init();
    _plotPlasma->setLegendOn(true);

    _plotPlasma->addCurve(0,"plasma");
    int nTime = _simulator.getNumberTimeBinsCoarse();
    FUNC_INFO << "nTime" << nTime;
    dVector xTime; xTime.resize(nTime);
    dVector yTAC;  yTAC.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = _simulator.getTimeCoarse(jt);
        yTAC[jt]  = _simulator.getCpCoarse(jt);
    }
    FUNC_INFO << xTime;
    FUNC_INFO << yTAC;
    _plotPlasma->setData(xTime,yTAC);
    FUNC_INFO << "here 1";

    _plotPlasma->setPointSize(2);
    FUNC_INFO << "here 2";
    _plotPlasma->conclude(0,true);
    FUNC_INFO << "here 3";
    _plotPlasma->plotDataAndFit(true);
    FUNC_EXIT;
}

void SimWindow::updateTargetGraph(dMatrix yData, dMatrix yFit)
{
    FUNC_ENTER;
    // update the plot
    _plotTarget->init();
    _plotTarget->setLegendOn(true);
    _plotTarget->setLegendPostition(1);
    _plotTarget1->init();
    _plotTarget1->setLegendOn(true);
    _plotTarget1->setLegendPostition(1);

    // add data
    dVector xTime = _simulator.getTimeCourse();
    dVector sim = yData[0];
    bVector ignoreRescale = _fMRIGLM.getIgnoredPointsInRun(0);
    _plotTarget->addCurve(0,"simulation");
    _plotTarget->setColor(Qt::black);
    _plotTarget->setPointSize(3);
    _plotTarget1->addCurve(0,"simulation");
    _plotTarget1->setColor(Qt::black);
    _plotTarget1->setPointSize(3);
    _plotTarget->setData(xTime,sim,ignoreRescale);
    _plotTarget1->setData(xTime,sim,ignoreRescale);

    // add fit
    dVector fit = yFit[0];
    _plotTarget->addCurve(0,"fit");
    _plotTarget->setLineThickness(2);
    _plotTarget->setColor(Qt::red);
    _plotTarget->setData(xTime,fit,ignoreRescale);
    _plotTarget1->addCurve(0,"fit");
    _plotTarget1->setLineThickness(2);
    _plotTarget1->setColor(Qt::red);
    _plotTarget1->setData(xTime,fit,ignoreRescale);

    // add exclusions
    _plotTarget->addCurve(0,"exclusions");
    _plotTarget->setColor(Qt::black);
    _plotTarget->setLineThickness(0);
    _plotTarget->setPointSize(20);
    _plotTarget->setPointStyle(QCPScatterStyle::ssStar);
    _plotTarget->setExport(false);
    _plotTarget1->addCurve(0,"exclusions");
    _plotTarget1->setColor(Qt::black);
    _plotTarget1->setLineThickness(0);
    _plotTarget1->setPointSize(20);
    _plotTarget1->setPointStyle(QCPScatterStyle::ssStar);
    _plotTarget1->setExport(false);
    // set time vectors
    dVector xDataIgnore, yDataIgnore;
    bVector ignorePoint;
    int nTime = xTime.size();
    for (int jt=0; jt<nTime; jt++) // loop over ALL points
    {
        if ( !(_fMRIGLM.getWeight(jt) != 0.) )
        {
            xDataIgnore.append(xTime[jt]);
            yDataIgnore.append(sim[jt]);
            ignorePoint.append(true);
        }
    }

    _plotTarget->setData(xDataIgnore, yDataIgnore, ignorePoint);
    _plotTarget->conclude(0,true);
    _plotTarget->plotDataAndFit(true);
    _plotTarget1->setData(xDataIgnore, yDataIgnore, ignorePoint);
    _plotTarget1->conclude(0,true);
    _plotTarget1->plotDataAndFit(true);

    FUNC_EXIT;
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
    for (int jt=0; jt<_fMRIGLM.getNumberTimePoints(); jt++)
    {
        xData.append(jt);
        yData.append(_fMRIGLM.getWeight(jt));
    }
    _plotBasis->setData(xData, yData);

    int nBasis = _fMRIGLM.getNumberCoefficients();
    for ( int jBasis=0; jBasis<nBasis; jBasis++)
    {
        if ( _fMRIGLM.getEventShape(jBasis) == Shape_Square )
        {
            _plotBasis->addCurve(0,"challenge");
            _plotBasis->setColor(Qt::red);
        }
        else
        {
            _plotBasis->addCurve(0,"drift");
            _plotBasis->setColor(Qt::black);
        }
        // add data
        yData.clear();
        for (int jt=0; jt<_fMRIGLM.getNumberTimePoints(); jt++)
            yData.append(_fMRIGLM.getBasisPoint(jBasis,jt));
        _plotBasis->setData(xData, yData);
    }

    _plotBasis->conclude(0,true);
    _plotBasis->plotDataAndFit(true);
}

void SimWindow::updateAllGraphs()
{
    FUNC_ENTER;
    // run the simulation
    _simulator.run();

    updatePlasmaGraph();
    analyzeTAC();
    updateBasisGraph();  // update basis graph AFTER target graph, because basis functions use target curve
    double AIC=0.;
    _analysisGroupBox->setTitle(QString("Target Region: analysis, AIC = %1").arg(AIC));
}
void SimWindow::analyzeTAC()
{
    FUNC_ENTER;
    int nTime = _simulator.getNumberTimeBinsCoarse();
    // define analysis
    if ( _fMRIGLM.getNumberRuns() != 1 )
    {
        _fMRIGLM.setNumberRuns(1);
        _fMRIGLM.setTimePointsInRun(0,nTime);
        _fMRIGLM.setTimeOrigin(0,0.);
        _fMRIGLM.setTimeStep(0,1.);
        _fMRIGLM.setBaselineTerms(3);
        int indexEvent1 = _fMRIGLM.getEventIndex('1');
        int indexEvent2 = _fMRIGLM.getEventIndex('2');
        _fMRIGLM.setStimulusRun(indexEvent1,0,0);
        _fMRIGLM.setStimulusRun(indexEvent2,0,0);
        _fMRIGLM.setEventShape(indexEvent1,Shape_Square);
        _fMRIGLM.setEventShape(indexEvent2,Shape_Square);
        _fMRIGLM.setStimulusTau(indexEvent1,0.);
        _fMRIGLM.setStimulusTau(indexEvent2,0.);
        if ( _fitChallenge->isChecked() )
        {
            _fMRIGLM.setStimulusMagnitude(indexEvent1,0,1.);
            _fMRIGLM.setStimulusMagnitude(indexEvent2,0,1.);
        }
        else
        {
            _fMRIGLM.setStimulusMagnitude(indexEvent1,0,0.);
            _fMRIGLM.setStimulusMagnitude(indexEvent2,0,1.);
        }
        bool ok;
        double onsetTime  = _challengeOnsetTime1->text().toDouble(&ok);
        double offsetTime = _challengeOffsetTime1->text().toDouble(&ok);
        _fMRIGLM.setStimulusOnset(indexEvent1,0,onsetTime);
        _fMRIGLM.setStimulusOffset(indexEvent1,0,offsetTime);
        onsetTime  = _challengeOnsetTime2->text().toDouble(&ok);
        offsetTime = _challengeOffsetTime2->text().toDouble(&ok);
        _fMRIGLM.setStimulusOnset(indexEvent2,0,onsetTime);
        _fMRIGLM.setStimulusOffset(indexEvent2,0,offsetTime);
    }

    dMatrix dataMatrix; dataMatrix.resize(1);
    if ( analyzeRealData() )
    {
        nTime = _dataTable[0].size();
        int indexTarget = _dataTargetRegion->currentIndex();
        dataMatrix[0] = _dataTable[indexTarget];
    }
    else
    {
        dataMatrix[0].resize(nTime);
        for (int jt=0; jt<nTime; jt++)
            dataMatrix[0][jt] = _simulator.getCtCoarse(jt);
    }
    if ( _signalDerivative->isChecked() )
    {
        differentiateSignal(dataMatrix);
        _fMRIGLM.setfPET(false);
    }
    else
        _fMRIGLM.setfPET(true);

    int indexEvent1 = _fMRIGLM.getEventIndex('1');
    int indexEvent2 = _fMRIGLM.getEventIndex('2');
    if ( _fitChallenge->isChecked() )
    {
        FUNC_INFO << "fit challenge";
        _fMRIGLM.setStimulusMagnitude(indexEvent1,0,1.);
        _fMRIGLM.setStimulusMagnitude(indexEvent2,0,1.);
        _fMRIGLM.prepare();
        _fMRIGLM.defineConditions("1");
    }
    else
    {
        FUNC_INFO << "do not fit challenge";
        _fMRIGLM.setStimulusMagnitude(indexEvent1,0,0.);
        _fMRIGLM.setStimulusMagnitude(indexEvent2,0,0.);
        _fMRIGLM.prepare();
    }
    dMatrix fitMatrix;     fitMatrix.resize(1);     fitMatrix[0].resize(nTime);
    FUNC_INFO << "challenge run" << _fMRIGLM.getStimulusRun(indexEvent1,0);
    _fMRIGLM.fitData(dataMatrix,fitMatrix);
    if ( _removeBaseline->isChecked() )
    {
        _fMRIGLM.convertAbsToPercent(dataMatrix,fitMatrix);
        _fMRIGLM.fitData(dataMatrix, fitMatrix);
    }
    _fMRIGLM.filterEvents(_removeBaseline->isChecked(), dataMatrix, fitMatrix);
    updateTargetGraph(dataMatrix, fitMatrix);
}

void SimWindow::differentiateSignal(dMatrix &signal)
{
    int nRuns = signal.size();
    for (int jRun=0; jRun<nRuns; jRun++)
    {
        dVector runVectorOriginal = signal[jRun];
        int nTime = runVectorOriginal.size();
        dVector diffVector; diffVector.fill(0.,nTime);
        for (int jt=1; jt<nTime; jt++)
        {
            double dt = _simulator.getTimeCoarse(jt) - _simulator.getTimeCoarse(jt-1);
            diffVector[jt] = (runVectorOriginal[jt] - runVectorOriginal[jt-1]) / dt;
        }
        signal[jRun] = diffVector;
    }
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

///////////////////////////////////////
// Slots
///////////////////////////////////////
void SimWindow::changedNumberThreads(int indexInBox)
{
    FUNC_ENTER;
    _nThreads = indexInBox + 1;
    QString numberString;
}

void SimWindow::changedNumberBins()
{
    FUNC_ENTER;
    bool ok;
    int numberBins = _numberTimeBins->text().toInt(&ok);
    if ( ok && numberBins != _simulator.getNumberBins() )
    {
        _simulator.setNumberBins(numberBins);
        if ( _fMRIGLM.getNumberTimePointsInRun(0) != numberBins ) _fMRIGLM.setTimePointsInRun(0,numberBins);
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
        _simulator.setDurationBin(iBin,duration);  // convert sec to min
        _fMRIGLM.setTimeStep(0,duration);
        updateAllGraphs();
    }
    else
    {
        QString text;
        _binDuration->setText(text.setNum(_simulator.getDurationPerBinSec(iBin)));
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
        _simulator.setTauBolus(value);
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
        _simulator.setKBol(value);
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
        _simulator.setKBolDelay(value);
        updateAllGraphs();
    }
    else
        _KBolDelay->setText(stringEntered.setNum(_simulator.getKBolDelay()));
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
void SimWindow::changedFastElimination()
{
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
void SimWindow::changedk2()
{
    FUNC_ENTER;
    QString stringEntered = _k2->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setk2(value);
        updateAllGraphs();
    }
    else
        _k3->setText(stringEntered.setNum(_simulator.getk3()));
}
void SimWindow::changedk3()
{
    FUNC_ENTER;
    QString stringEntered = _k3->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setk3(value);
        updateAllGraphs();
    }
    else
        _k3->setText(stringEntered.setNum(_simulator.getk3()));
}
void SimWindow::changedK1()
{
    FUNC_ENTER;
    QString stringEntered = _K1->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setK1(value);
        updateAllGraphs();
    }
    else
        _K1->setText(stringEntered.setNum(_simulator.getK1()));
}
void SimWindow::changedChallengeOnsetTime1()
{
    FUNC_ENTER;
    QString stringEntered = _challengeOnsetTime1->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeOnsetTime1(value);
        int indexEvent1 = _fMRIGLM.getEventIndex('1');
        _fMRIGLM.setStimulusOnset(indexEvent1,0,value);
        updateAllGraphs();
    }
    else
        _challengeOnsetTime1->setText(stringEntered.setNum(_simulator.getChallengeOnsetTime1()));
}
void SimWindow::changedChallengeOnsetTime2()
{
    FUNC_ENTER;
    QString stringEntered = _challengeOnsetTime2->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeOnsetTime2(value);
        int indexEvent2 = _fMRIGLM.getEventIndex('2');
        _fMRIGLM.setStimulusOnset(indexEvent2,0,value);
        updateAllGraphs();
    }
    else
        _challengeOnsetTime2->setText(stringEntered.setNum(_simulator.getChallengeOnsetTime2()));
}
void SimWindow::changedChallengeOffsetTime1()
{
    FUNC_ENTER;
    QString stringEntered = _challengeOffsetTime1->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeOffsetTime1(value);
        int indexEvent1 = _fMRIGLM.getEventIndex('1');
        _fMRIGLM.setStimulusOffset(indexEvent1,0,value);
        updateAllGraphs();
    }
    else
        _challengeOffsetTime1->setText(stringEntered.setNum(_simulator.getChallengeOffsetTime1()));
}
void SimWindow::changedChallengeOffsetTime2()
{
    FUNC_ENTER;
    QString stringEntered = _challengeOffsetTime2->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeOffsetTime2(value);
        int indexEvent2 = _fMRIGLM.getEventIndex('2');
        _fMRIGLM.setStimulusOffset(indexEvent2,0,value);
        updateAllGraphs();
    }
    else
        _challengeOffsetTime1->setText(stringEntered.setNum(_simulator.getChallengeOffsetTime2()));
}
void SimWindow::changedChallengeMag1()
{
    FUNC_ENTER;
    QString stringEntered = _challengeMagPercent1->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeMag1(value);
        updateAllGraphs();
    }
    else
        _challengeMagPercent1->setText(stringEntered.setNum(_simulator.getChallengeMagPercent1()));
}
void SimWindow::changedChallengeMag2()
{
    FUNC_ENTER;
    QString stringEntered = _challengeMagPercent2->text();
    bool ok;
    double value = stringEntered.toDouble(&ok);
    if ( ok )
    {
        _simulator.setChallengeMag2(value);
        updateAllGraphs();
    }
    else
        _challengeMagPercent2->setText(stringEntered.setNum(_simulator.getChallengeMagPercent2()));
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
        updateAllGraphs();
    }
    else
        _noiseTar->setText(utilityString.setNum(_simulator.getNoiseTar()));
    bool noisy = _simulator.getNoiseTar() != 0.;
    setThreadVisibility(noisy);
}
void SimWindow::changedTauAnalysis()
{
    FUNC_ENTER;
    QString utilityString = _tauAnalysis->text();
    bool ok;
    double value = utilityString.toDouble(&ok);
    int indexEvent1 = _fMRIGLM.getEventIndex('1');
    int indexEvent2 = _fMRIGLM.getEventIndex('2');
    if ( ok )
    {
        _fMRIGLM.setStimulusTau(indexEvent1,value);
        _fMRIGLM.setStimulusTau(indexEvent2,value);
        updateAllGraphs();
    }
    else
    {
        _tauAnalysis->setText(utilityString.setNum(_fMRIGLM.getStimulusTau(indexEvent1)));
        _tauAnalysis->setText(utilityString.setNum(_fMRIGLM.getStimulusTau(indexEvent2)));
    }
}
void SimWindow::changedIgnoreString()
{
    FUNC_ENTER;
    QString stringEntered = _ignoreString->text();
    _fMRIGLM.setIgnoredPoints(0,true,stringEntered);
    updateAllGraphs();
}
void SimWindow::changedCheckBoxChallenge(bool state)
{
    FUNC_ENTER << state;
    _errorChallengeLabel->setVisible(state);
    _errorChallenge->setVisible(state);

    int indexEvent1 = _fMRIGLM.getEventIndex('1');
    int indexEvent2 = _fMRIGLM.getEventIndex('2');
    if ( state )
    {
        _fMRIGLM.setStimulusMagnitude(indexEvent1,0,1.);
        _fMRIGLM.setStimulusMagnitude(indexEvent2,0,1.);
        _fMRIGLM.prepare();
        _fMRIGLM.defineConditions("1");
    }
    else
    {
        _fMRIGLM.setStimulusMagnitude(indexEvent1,0,0.);
        _fMRIGLM.setStimulusMagnitude(indexEvent2,0,0.);
        _fMRIGLM.prepare();
    }
    _fMRIGLM.setPrepared(false);
    updateAllGraphs();
}
void SimWindow::changedCheckBoxBaseline(bool state)
{
    FUNC_ENTER << state;
    updateAllGraphs();
}
void SimWindow::changedDerivativeRadioButton(bool state)
{
    updateAllGraphs();
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
