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
    qDebug() << "SimWindow::SimWindow enter";
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
    _tabTimeSpace->addTab(_setupPage, tr("setup"));
    _tabTimeSpace->addTab(_targetPage, tr("target"));
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
    setThreadVisibility(false);
    updateAllGraphs();

    qDebug() << "SimWindow::SimWindow exit";
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
        qDebug() << "readTableFile";
        QString errorString = utilIO::readTableFile(fileName, _dataColumnNames, _dataTable);
        qDebug() << "_dataColumnNames" << _dataColumnNames;
        qDebug() << "_dataTable size" << _dataTable.size();
        qDebug() << "_dataTable column0" << _dataTable[0];
        qDebug() << "index put" << _dataColumnNames.indexOf(QRegExp("put"), Qt::CaseInsensitive);
        qDebug() << "index put" << _dataColumnNames.indexOf(QRegExp("putamen", Qt::CaseInsensitive));
        qDebug() << "index put" << _dataColumnNames.indexOf(QRegExp("put*", Qt::CaseInsensitive));
        qDebug() << "index put" << _dataColumnNames.indexOf(QRegExp("put*", Qt::CaseInsensitive, QRegExp::Wildcard));
        qDebug() << "index put" << _dataColumnNames.indexOf(QRegExp("*put", Qt::CaseInsensitive, QRegExp::Wildcard));
        if ( !errorString.isEmpty() )
        {
            QMessageBox msgBox;
            msgBox.setText(errorString);
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.exec();
        }
        else
        {
            qDebug() << "getTableDataFile 1";
            QStringList validBinSizeName = {"dt","deltaTime","delta-time","delta_time",
                                            "binSize","bin-size","bin_size",
                                            "binTime","bin-time","bin_time",
                                            "frameDuration","frame-duration","frame_duration",
                                            "frameSize","frame-size","frame_size"};
            QStringList validBinTimeName = {"t","time","time-point","time_point",
                                            "frameTime","frame-time","frame_time"};
            int iColumnBinSize=-1;  int iColumnBinTime=-1;
            if ( !defineTimeBinsFromBinSize(validBinSizeName, iColumnBinSize) )
            {
                qDebug() << "getTableDataFile 2";
                if ( !defineTimeBinsFromTimePoints(validBinTimeName, iColumnBinTime) )
                {
                    qDebug() << "getTableDataFile 3";
                    qDebug() << "getTableDataFile 4";
                    QString binList = validBinSizeName.join(", ");
                    QString binTime = validBinTimeName.join(", ");
                    QMessageBox msgBox;
                    errorString = QString("The table must include either frame durations (1) or time points (2).\n(1)%1\n(2)%2").
                            arg(binList).arg(binTime);
                    msgBox.setText(errorString);
                    msgBox.setIcon(QMessageBox::Critical);
                    msgBox.exec();
                    return;
                }
            }
            // update the GUI to refect the new information
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
            qDebug() << "iColumns" << iColumnBinSize << iColumnBinTime << iColumnFirstRegion;
            _dataTargetRegion->setCurrentIndex(iColumnFirstRegion);
            enablePlasmaMatching(true);
            _ROIFileName->setToolTip(fileName);
            _ROIFileName->setText(utilString::getFileNameWithoutDirectory(fileName));
            // Resample the simulation to match the binning of the real data
            QString text;
            _numberTimeBins->setText(text.setNum(_dtBins.size())); changedNumberBins();
            qDebug() << "SimWindow::getTableDataFile setDurationBins" << _dtBins;
            _simulator.setDurationBins(_dtBins);
            updateAllGraphs();
        }
    }
}

bool SimWindow::defineTimeBinsFromBinSize(QStringList validBinName, int &iColumn)
{
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
    _dtBins.clear();  _timeBins.clear();
    double time=_dataTable[iColumn][0]/60./2.;  // center of 1st bin in min
    for ( int jt=0; jt<_dataTable[iColumn].size(); jt++)
    {
        _dtBins.append(_dataTable[iColumn][jt]/60.); // seconds to minutes
        _timeBins.append(time);
        time += _dtBins.last();
    }
    qDebug() << "SimWindow::defineTimeBinsFromBinSize iColumn" << iColumn << _dataColumnNames.at(iColumn);
    return true;
}

bool SimWindow::defineTimeBinsFromTimePoints(QStringList validBinName, int &iColumn)
{
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
    _dtBins.clear();  _timeBins.clear();
    for ( int jt=0; jt<_dataTable[iColumn].size(); jt++)
    {
        double time = _dataTable[iColumn][jt];
        _timeBins.append(time);
        if ( jt == 0 )
            _dtBins.append(2.*time); // seconds to minutes
        else
        {
            double lowerEdge = _dataTable[iColumn][jt-1] + _dtBins.last()/2.;
            _dtBins.append(2.*(time-lowerEdge));
        }
    }
    qDebug() << "SimWindow::defineTimeBinsFromTimePoints iColumn" << iColumn << _dataColumnNames.at(iColumn);
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
    _RRStatusBar = new QStatusBar;  // must be global so it doesn't go out of scope
    _RRStatusBar->setStyleSheet("color:blue");
    _RRPlot->setQCStatusBar(_RRStatusBar);
//    _RRPlot->setMainPage(this);

    _setupPlotLayout = new QVBoxLayout();
    _setupPlotLayout->addWidget(_plasmaPlot->getPlotSurface());
    _setupPlotLayout->addWidget(_RRPlot->getPlotSurface());
    _setupPlotLayout->addWidget(_RRStatusBar);
    _setupPlotLayout->setStretch(0,1);     // basis plot
    _setupPlotLayout->setStretch(1,1);     // time plot
    showPlasmaRR();
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
    _binDuration->setText(text.setNum(_simulator.getDurationPerBin(0)*60.));
    _subSample->setText(numberString.setNum(_simulator.getSamplesPerBin(0)));

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
    connect(_binIndex, SIGNAL(valueChanged(int)), this, SLOT(changedBinIndex(int)));
    connect(_binDuration, SIGNAL(editingFinished()), this, SLOT(changedBinDuration()));
    connect(_subSample,   SIGNAL(editingFinished()), this, SLOT(changedSubSample()));

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
    auto *realDataLayout = new QGridLayout();
    realDataLayout->addWidget(_readROIFile,0,0);
    realDataLayout->addWidget(_ROIFileName,0,1);
    realDataLayout->addWidget(dataLabel,1,0);
    realDataLayout->addWidget(_dataRefRegion,1,1);
    realDataGroupBox->setLayout(realDataLayout);

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
    auto *radioShowPlasmaRR  = new QRadioButton("Plasma/RR",_setupPage);
    auto *radioShowPlasma    = new QRadioButton("Plasma",_setupPage);
    auto *radioShowRR        = new QRadioButton("RR",_setupPage);
    radioShowPlasmaRR->setChecked(true);

    QToolBar *graphToolBar = new QToolBar("graph tool bar");
    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addAction(rescaleXYAction);
//    graphToolBar->addAction(crossCursorAct);
    graphToolBar->addWidget(separator1);
    graphToolBar->addWidget(showLabel);
    graphToolBar->addWidget(radioShowPlasmaRR);
    graphToolBar->addWidget(radioShowPlasma);
    graphToolBar->addWidget(radioShowRR);
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

    connect(radioShowPlasmaRR, SIGNAL(clicked(bool)), this, SLOT(showPlasmaRR()));
    connect(radioShowPlasma,   SIGNAL(clicked(bool)), this, SLOT(showPlasma()));
    connect(radioShowRR,       SIGNAL(clicked(bool)), this, SLOT(showRR()));

    connect(_tau2Ref,  SIGNAL(editingFinished()), this, SLOT(changedTau2Ref()));
    connect(_tau1Ref,  SIGNAL(editingFinished()), this, SLOT(changedTau1Ref()));
    connect(_plasmaFracRef, SIGNAL(editingFinished()), this, SLOT(changedPlasmaFracRef()));
    connect(_noiseRef, SIGNAL(editingFinished()), this, SLOT(changedNoiseRef()));

    connect(_readROIFile, SIGNAL(pressed()), this, SLOT(getTableDataFile()));
    connect(_dataRefRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(updateAllGraphs()));

    _setupPage->setLayout(fullLayout);
}

void SimWindow::enableComboBoxItem(QComboBox *comboBox, int itemNumber, bool enable)
{
    QStandardItemModel* model = qobject_cast<QStandardItemModel*>(comboBox->model());
    QStandardItem* item= model->item(itemNumber);
    if ( enable )
        item->setFlags(item->flags() | Qt::ItemIsEnabled);   // enable
    else
        item->setFlags(item->flags() & ~Qt::ItemIsEnabled);  // disable
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
    _PETRTM.setPrepared(false);
    updateAllGraphs();

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
    _tau4AnalysisLabel->setVisible(_PETRTM.isFRTM());
    _tau4Analysis->setVisible(_PETRTM.isFRTM());
}

void SimWindow::createTargetPage()
{
    qDebug() << "SimWindow::createTargetPage enter";
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
    connect(_dataTargetRegion, SIGNAL(currentIndexChanged(int)), this, SLOT(updateAllGraphs()));

    // target region: analysis
    auto *analysisGroupBox = new QGroupBox("Target Region: analysis");
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
    _fitk4 = new QCheckBox("fit k4?");
    _fitk4->setChecked(false);
//    _fitk4->setEnabled(false);  // xxx tmp
    _fitChallenge = new QCheckBox("fit challenge?");
    _fitChallenge->setChecked(false);

    auto *analysisLayout = new QGridLayout();
    analysisLayout->addWidget(modelLabel,0,0);
    analysisLayout->addWidget(_modelType,0,1);
    analysisLayout->addWidget(_fitChallenge,1,0);
    analysisLayout->addWidget(_fitk4,1,1);
    analysisLayout->addWidget(weightLabel,2,0);
    analysisLayout->addWidget(_weightType,2,1);
    analysisLayout->addWidget(ignoreLabel,3,0);
    analysisLayout->addWidget(_ignoreString,3,1);
    analysisLayout->addWidget(_tau2RefAnalysisLabel,4,0);
    analysisLayout->addWidget(_tau2RefAnalysis,4,1);
    analysisLayout->addWidget(_tau4AnalysisLabel,5,0);
    analysisLayout->addWidget(_tau4Analysis,5,1);
    analysisGroupBox->setLayout(analysisLayout);
    analysisLayout->setSpacing(0);
    connect(_ignoreString,     SIGNAL(editingFinished()), this, SLOT(changedIgnoreString()));
    connect(_tau2RefAnalysis,  SIGNAL(editingFinished()), this, SLOT(changedTau2RefAnalysis()));
    connect(_tau4Analysis,     SIGNAL(editingFinished()), this, SLOT(changedTau4Analysis()));
    connect(_fitChallenge,     SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxChallenge(bool)));
    connect(_fitk4,            SIGNAL(toggled(bool)),     this, SLOT(changedCheckBoxFitk4(bool)));
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
    rightLayout->addWidget(analysisGroupBox);
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
    QButtonGroup *showGroup = new QButtonGroup(_targetPage);
    showGroup->addButton(radioShowBasisTarget);
    showGroup->addButton(radioShowBasis);
    showGroup->addButton(radioShowTarget);
    radioShowBasisTarget->setChecked(true);

    QLabel *analyzeLabel = new QLabel("analyze:",_targetPage);
    _analyzeSimulation = new QRadioButton("Simulation(s)",_targetPage);
    _analyzeRealData = new QRadioButton("ROI Data",_targetPage);
    QButtonGroup *analyzeGroup = new QButtonGroup(_targetPage);
    analyzeGroup->addButton(_analyzeSimulation);
    analyzeGroup->addButton(_analyzeRealData);
    _analyzeSimulation->setChecked(true);

    QToolBar *graphToolBar = new QToolBar("graph tool bar");
    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addAction(rescaleXYAction);
    graphToolBar->addWidget(separator1);
    graphToolBar->addWidget(showLabel);
    graphToolBar->addWidget(radioShowBasisTarget);
    graphToolBar->addWidget(radioShowBasis);
    graphToolBar->addWidget(radioShowTarget);
    graphToolBar->addWidget(separator2);
    graphToolBar->addWidget(analyzeLabel);
    graphToolBar->addWidget(_analyzeSimulation);
    graphToolBar->addWidget(_analyzeRealData);
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

    _targetPage->setLayout(fullLayout);
    qDebug() << "SimWindow::createTargetPage exit";
}

void SimWindow::createSweepBPndPage()
{
    qDebug() << "SimWindow::createSweepBPndPage enter";
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
    qDebug() << "SimWindow::createSweepBPndPage 1";
    auto *plotLayout = new QVBoxLayout();
    plotLayout->addWidget(_errBPndPlot->getPlotSurface());
    plotLayout->addWidget(_errChallPlot->getPlotSurface());
    plotLayout->addWidget(_tau2RefPlot->getPlotSurface());
    plotLayout->addWidget(_errk4Plot->getPlotSurface());
    QString numberString;
    int editTextSize=80;
    qDebug() << "SimWindow::createSweepBPndPage 2";

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
    qDebug() << "SimWindow::createSweepBPndPage 3";

    clearBPndCurves();
    qDebug() << "SimWindow::createSweepBPndPage 4";
    changedVersusBPndGraphs();

    _sweepBPndPage->setLayout(fullLayout);
    qDebug() << "SimWindow::createSweepBPndPage exit";
}

void SimWindow::createVersusTimePage()
{
    qDebug() << "SimWindow::createVersusTimePage enter";
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
    qDebug() << "SimWindow::createTargetPage exit";

}

void SimWindow::showPlasmaRR()
{
    _plasmaPlot->getPlotSurface()->setVisible(true);
    _RRPlot->getPlotSurface()->setVisible(true);
}
void SimWindow::showPlasma()
{
    _plasmaPlot->getPlotSurface()->setVisible(true);
    _RRPlot->getPlotSurface()->setVisible(false);
}
void SimWindow::showRR()
{
    _plasmaPlot->getPlotSurface()->setVisible(false);
    _RRPlot->getPlotSurface()->setVisible(true);
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
    int nTime = _simulator.getNumberTimeBinsFine();
    dVector xTime; xTime.resize(nTime);
    dVector yTAC;  yTAC.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = _simulator.getTimeFine(jt);
        yTAC[jt]  = _simulator.getCp(jt);
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

    if ( realDataAvailable() )
    {
        qDebug() << "read data not available";
        addDataCurveRR();
        addSimulationCurveRR();
    }
    else
    {
        qDebug() << "read data available";
        addSimulationCurveRR();
        addDataCurveRR();
    }

    _RRPlot->conclude(0,true);
    _RRPlot->plotDataAndFit(true);

    qDebug() << "SimWindow::updateReferenceGraph exit";
}

void SimWindow::addSimulationCurveRR()
{
    if ( realDataAvailable() )
        _RRPlot->addCurve(0,"RR: simulation");
    else
        _RRPlot->addCurve(0,"RR");
    _RRPlot->setColor(Qt::red);
    _RRPlot->setPointSize(5);
    int nTime = _simulator.getNumberTimeBinsCoarse();
    dVector xTime;   xTime.resize(nTime);
    dVector yTAC;    yTAC.resize(nTime);
    dMatrix refRegion;  refRegion.resize(1);  refRegion[0].resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = _simulator.getTimeCoarse(jt);
        yTAC[jt]  = refRegion[0][jt] = _simulator.getCrDown(jt);
    }
    qDebug() << "SimWindow::addSimulationCurveRR bins" << xTime;
//    qDebug() << "SimWindow::addSimulationCurveRR RR"   << yTAC;
    _RRPlot->setData(xTime,yTAC);
}
void SimWindow::addDataCurveRR()
{
    if ( realDataAvailable() )
    {
        int indexInBox = _dataRefRegion->currentIndex();
        if ( indexInBox >= _dataTable.size() )
            qFatal("Fatal Error: the combo-box index exceeds the table index.");
        _RRPlot->addCurve(0,"RR: ROI data");
        _RRPlot->setColor(Qt::gray);
        _RRPlot->setPointSize(5);
        _RRPlot->setData(_timeBins,_dataTable[indexInBox]);
        qDebug() << "SimWindow::addDataCurveRR bins" << _timeBins;
//        qDebug() << "SimWindow::addDataCurveRR RR" << _dataTable[indexInBox];
        // enable buttons
        _analyzeSimulation->setEnabled(true);
        _analyzeRealData->setEnabled(true);
        _targetDataGroupBox->setVisible(true);
    }
    else
    {
        _analyzeSimulation->setEnabled(false);
        _analyzeRealData->setEnabled(false);
        _targetDataGroupBox->setVisible(false);
    }
}
void SimWindow::addSimulationCurveTarget()
{
    if ( realDataAvailable() )
        _targetPlot->addCurve(0,"target: simulation");
    else
        _targetPlot->addCurve(0,"target");
    _targetPlot->setColor(Qt::red);
    _targetPlot->setPointSize(5);
    int nTime = _simulator.getNumberTimeBinsCoarse();
    dVector xTime;   xTime.resize(nTime);
    dVector yTAC;    yTAC.resize(nTime);
    dMatrix tissueVector;  tissueVector.resize(1);  tissueVector[0].resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        xTime[jt] = _simulator.getTimeCoarse(jt);
        yTAC[jt]  = tissueVector[0][jt] = _simulator.getCtDown(jt);
    }
    _PETRTM.setTissueVector(tissueVector);
    qDebug() << "SimWindow::addSimulationCurveRR bins" << xTime;
//    qDebug() << "SimWindow::addSimulationCurveRR RR"   << yTAC;
    _targetPlot->setData(xTime,yTAC);

    dMatrix fitVector;     fitVector.resize(1);     fitVector[0].resize(nTime);
    defineRTMModel();
    // update the RTM model
    _PETRTM.definePETConditions("a c"); // don't define R1, which is not valid for RTM2
    _PETRTM.prepare();
    _PETRTM.fitData(tissueVector,fitVector);

}
void SimWindow::addDataCurveTarget()
{
    if ( realDataAvailable() )
    {
        qDebug() << "SimWindow::addDataCurveTarget 1";
        _targetPlot->addCurve(0,"target: ROI data");
        _targetPlot->setColor(Qt::gray);
        _targetPlot->setPointSize(5);
        qDebug() << "SimWindow::addDataCurveTarget 2";
        int indexInBox = _dataTargetRegion->currentIndex();
        if ( indexInBox >= _dataTable.size() )
            qFatal("Fatal Error: the combo-box index exceeds the table index.");
        else if ( indexInBox < 0 ) return;
        qDebug() << "SimWindow::addDataCurveTarget 3";
        qDebug() << "SimWindow::addDataCurveTarget 3.1" << _timeBins;
        qDebug() << "SimWindow::addDataCurveTarget 3.2" << indexInBox;
        qDebug() << "SimWindow::addDataCurveTarget 3.3" << _dataTable[indexInBox];
        _targetPlot->setData(_timeBins,_dataTable[indexInBox]);
        qDebug() << "SimWindow::addDataCurveTarget 3.5";
        dMatrix tissueVector; tissueVector.resize(1); tissueVector[0] = _dataTable[indexInBox];
        qDebug() << "SimWindow::addDataCurveTarget 4";
        _PETRTM.setTissueVector(tissueVector);
        // enable buttons
        int nTime = _dataTable[0].size();
        dMatrix fitVector;     fitVector.resize(1);     fitVector[0].resize(nTime);
        defineRTMModel();
        // update the RTM model
        _PETRTM.definePETConditions("a c"); // don't define R1, which is not valid for RTM2
        _PETRTM.prepare();
        _PETRTM.fitData(tissueVector,fitVector);
    }
}

void SimWindow::updateTargetGraph()
{
    qDebug() << "SimWindow::updateTargetGraph enter";
    // run the simulation
    _simulator.generateTargetTAC();

    // update the plot
    _targetPlot->init();
    _targetPlot->setLegendOn(true);

    if ( realDataAvailable() )
    {
        qDebug() << "real data available";
        addDataCurveTarget();
        addSimulationCurveTarget();
    }
    else
    {
        qDebug() << "read data not available";
        addSimulationCurveTarget();
        addDataCurveTarget();
    }

    analyzeSimulatedTAC();

    int nTime = _simulator.getNumberTimeBinsCoarse();
    dVector xTime;   xTime.resize(nTime);
    dVector yTAC;    yTAC.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
        xTime[jt] = _simulator.getTimeCoarse(jt);

    // add fit
    _targetPlot->addCurve(0,"fit");
    _targetPlot->setLineThickness(2);
    _targetPlot->setColor(Qt::blue);
    for (int jt=0; jt<nTime; jt++)
        yTAC[jt] = _PETRTM.getFit(jt);
    _targetPlot->setData(xTime,yTAC);

    // add RR
    _targetPlot->addCurve(0,"RR");
    _targetPlot->setColor(Qt::gray);
    for (int jt=0; jt<nTime; jt++)
        yTAC[jt]  = _simulator.getCrDown(jt);
    _targetPlot->setData(xTime,yTAC);

    // add FRTM convolution
    if ( _PETRTM.isFRTM() )
    {
        if ( _PETRTM.isFRTMNew() )
            _targetPlot->addCurve(0,"BP0*convolution");
        else
            _targetPlot->addCurve(0,"convolution");
        _targetPlot->setColor(Qt::magenta);
        for (int jt=0; jt<nTime; jt++)
            yTAC[jt]  = _simulator.getBP0() * _PETRTM.getFRTMConvolution(0,jt);
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
    qDebug() << "SimWindow::updateTargetGraph exit";
}

void SimWindow::defineRTMModel()
{
    qDebug() << "SimWindow::defineRTMModel enter";
    static int firstTime=true;
    if ( analyzeRealData() )
    {
        // matrices: [nRuns][nTime]
        dMatrix timeBins;  timeBins.resize(1);
        dMatrix refRegion; refRegion.resize(1);
        timeBins[0]  = _dataTable[0];         // _dataTable[nColums]
        refRegion[0] = _RRPlot->getYData(0);  // getYData(iCurve)
        _PETRTM.setReferenceRegion(timeBins,refRegion);
    }
    else
    {
        qDebug() << "SimWindow::defineRTMModel 1";
        int nTime = _simulator.getNumberTimeBinsCoarse();
        if ( firstTime || nTime != _PETRTM.getNumberTimePoints() )
        {
            qDebug() << "SimWindow::defineRTMModel 2";
            dMatrix timeBins;  timeBins.resize(1);
            timeBins[0] = _simulator.getTimeBinVector();
            dMatrix refRegion; refRegion.resize(1);
            refRegion[0].resize(nTime);
            for (int jt=0; jt<nTime; jt++)
                refRegion[0][jt] = _simulator.getCrDown(jt);
            qDebug() << "timeBins =" << timeBins;
            qDebug() << "refRegion =" << refRegion;
            _PETRTM.setReferenceRegion(timeBins,refRegion);
        }
        qDebug() << "SimWindow::defineRTMModel 3";

        // Assign these IDs:
        _PETRTM.setR1EventID(0,'R');
        _PETRTM.setk2EventID(0,'k');
        _PETRTM.setk2aEventID(0,'a');
        int indexChallenge = _PETRTM.getEventIndex('c');
        _PETRTM.setChallengeShape(indexChallenge,Challenge_Sigmoid);
        _PETRTM.setChallengeTau(indexChallenge,0.1);
        qDebug() << "SimWindow::defineRTMModel 4";
        if ( firstTime )
        {
            qDebug() << "SimWindow::defineRTMModel 5" << indexChallenge;
            _PETRTM.setChallengeOnset(indexChallenge,0,_simulator.getChallengeTime());
            qDebug() << "SimWindow::defineRTMModel 5.1";
            _PETRTM.setTau2RefSRTM(0,_simulator.getTau2Ref());
            qDebug() << "SimWindow::defineRTMModel 5.2";
            _PETRTM.setTau2RefFRTM(0,_simulator.getTau2Ref());
            qDebug() << "SimWindow::defineRTMModel 5.3";
            _PETRTM.setTau4(0,_simulator.getTau4());
        }
        qDebug() << "SimWindow::defineRTMModel 6";
        firstTime = false;
    }
    qDebug() << "SimWindow::defineRTMModel exit";
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
    qDebug() << "SimWindow::updateBasisGraph exit";
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
            qDebug() << "compare dk2a" << truth << guess;
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
    qDebug() << "k2Ref =" << k2Ref.x << "±" << k2Ref.y << "with relative error" << k2Ref.y / k2Ref.x * 100.;

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
        _errorChallenge->setText(valueString + " %, diff = " + diffString + " %");
    }

    // k4
    if ( _fitk4->isChecked() )
    {
        truth = _simulator.getTau4();  // delta_BPnd abs
        guess = _PETRTM.getTau4InRun(0);
        valueString.setNum(guess,'g',2);
        if ( guess == 0. )
            _errork4->setText("----");
        else
            _errork4->setText(analyzeString(truth,guess));
    }

    /*
    int nCoeff = _PETRTM.getNumberCoefficients();
    qDebug() << "***" << nCoeff << "coefficients ***";
    for (int jCoeff=0; jCoeff<nCoeff; jCoeff++)
        qDebug() << "coefficient" << jCoeff << "=" << _PETRTM.getBeta(jCoeff);
        */

    // sigma2
    double sigma = qSqrt(_PETRTM.getSigma2());
    valueString.setNum(sigma,'g',1);
    _sigma->setText(valueString);

    qDebug() << "SimWindow::analyzeSimulatedTAC exit";
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
        if ( _PETRTM.isFRTMNew() )
            percentChange *= -1.;
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
    if ( ok )
    {
        _simulator.setNumberBins(numberBins);
        _PETRTM.setTimeBins(0,_simulator.getTimeBinVector());
        updateAllGraphs();
    }
    else
    {
        QString text;
        _binDuration->setText(text.setNum(_simulator.getDuration()));
    }
}
void SimWindow::changedBinDuration()
{
    bool ok;
    double duration = _binDuration->text().toDouble(&ok);
    if ( ok )
    {
        int iBin = _binIndex->text().toInt(&ok);
        _simulator.setDurationBin(iBin,duration/60.);  // convert sec to min
        _PETRTM.setTimeBins(0,_simulator.getTimeBinVector());
        updateAllGraphs();
    }
    else
    {
        QString text;
        _binDuration->setText(text.setNum(_simulator.getDuration()));
    }
}
void SimWindow::changedSubSample()
{
    QString stringEntered = _subSample->text();
    bool ok;
    int lSubSample = stringEntered.toInt(&ok);
    int iBin = _binIndex->text().toInt(&ok);
    if ( ok )
    {
        _simulator.setSamplesPerBin(iBin,lSubSample);
        _PETRTM.setTimeBins(0,_simulator.getTimeBinVector());
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

void SimWindow::changedPlasmaFracRef()
{
    QString utilityString = _plasmaFracRef->text();
    bool ok;
    double value = utilityString.toDouble(&ok);
    if ( ok )
    {
        _simulator.setPlasmaPercentRef(value);
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
void SimWindow::changedCheckBoxFitk4(bool state)
{
    qDebug() << "SimWindow::changedCheckBoxFitk4" << state;
    _errork4Label->setVisible(state);
    _errork4->setVisible(state);
    _checkBoxk4ErrGraph->setChecked(state);
    _checkBoxk4ErrGraph->setVisible(state);
//    clearTimeCurves();
    _PETRTM.setFitk4State(state);
    _PETRTM.setPrepared(false);
    qDebug() << "SimWindow::changedCheckBoxFitk4 update";
    updateAllGraphs();
    qDebug() << "SimWindow::changedCheckBoxFitk4 exit";
}
void SimWindow::changedCheckBoxChallenge(bool state)
{
    qDebug() << "SimWindow::changedCheckBoxChallenge" << state;
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
    qDebug() << "SimWindow::changedVersusBPndGraphs enter";
    _errBPndPlot->getPlotSurface()->setVisible(_checkBoxBPndErrGraph->isChecked());
    qDebug() << "SimWindow::changedVersusBPndGraphs 1";
    _errChallPlot->getPlotSurface()->setVisible(_checkBoxChallErrGraph->isChecked());
    qDebug() << "SimWindow::changedVersusBPndGraphs 2";
    _tau2RefPlot->getPlotSurface()->setVisible(_checkBoxTau2RefGraph->isChecked());
    qDebug() << "SimWindow::changedVersusBPndGraphs 3";
    _errk4Plot->getPlotSurface()->setVisible(_checkBoxk4ErrGraph->isChecked());
    qDebug() << "SimWindow::changedVersusBPndGraphs exit";
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
    qDebug() << "SimWindow::calculateBPndCurves enter" << noisy << _simulator.getNoiseRef() << _simulator.getNoiseTar();
    if ( noisy )
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
    for (double BP0=_BPndLowValue; BP0<=_BPndHighValue; BP0 += _BPndStepValue)
    {
        qDebug() << "BP0 =" << BP0;
        xVector.append(BP0);
        _simulator.setBP0(BP0);  // set BP0 and run a series of samples

        _simulator.run();  // run the simulations with randomized noise
        // Perform the analysis
        for (int jt=0; jt<nTime; jt++)
        {
            refRegion[0][jt]    = _simulator.getCrDown(jt);
            tissueVector[0][jt] = _simulator.getCtDown(jt);
        }
        _PETRTM.setReferenceRegion(refRegion);
        _PETRTM.setTissueVector(tissueVector);
        _PETRTM.prepare();
        _PETRTM.fitData(tissueVector,fitVector);

        // update the BP error
        double guess = _PETRTM.getBP0InRun(0);
        errBPnd.append(percentageError(guess,BP0));
        // update the challenge error
        double truth = _simulator.getChallengeMag();
        guess = getChallengeMagFromAnalysis();
        errChall.append(guess - truth);
        // update the 1/k4 error
        if ( _PETRTM.getFitk4State() )
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
        _errk4Plot->setData(xVector,errTau4);

    if ( _PETRTM.isRTM2() && _checkBoxTau2RefGraph->isChecked() )
    { // add a new curves to the two error plots
        _errBPndPlot->addCurve(0,"error_BPnd");
        _errBPndPlot->setPointSize(5);
        _errBPndPlot->setPointStyle(QCPScatterStyle::ssCross);
        _errBPndPlot->setColor(colors[iColor]);
        qDebug() << "setData4" << errBPndRTM2.size();
        _errBPndPlot->setData(xVector,errBPndRTM2);
        _errChallPlot->addCurve(0,"error_Challenge");
        _errChallPlot->setPointSize(5);
        _errChallPlot->setPointStyle(QCPScatterStyle::ssCross);
        _errChallPlot->setColor(colors[iColor]);
        qDebug() << "setData5" << errChallRTM2.size();
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
    qDebug() << "SimWindow::bestTau2RefForRTM2 enter";
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
void SimWindow::calculateTimeCurves()
{
    QMessageBox msgBox;
    msgBox.setText("Function calculateTimeCurves() currently disabled");
    msgBox.setIcon(QMessageBox::Critical);
    msgBox.exec();
}
/*
void SimWindow::calculateTimeCurves()
{
    QColor colors[10] = {Qt::black, Qt::red, Qt::blue, Qt::green,
                         Qt::darkCyan, Qt::darkYellow, Qt::darkMagenta, Qt::darkRed, Qt::darkBlue, Qt::darkGreen};
    int nCurves = _errBPndOrChallVsTimePlot->getNumberCurves();
    int iColor  = nCurves%10;

    _errBPndOrChallVsTimePlot->addCurve(0,"error_BPnd");
    _errBPndOrChallVsTimePlot->setPointSize(5);
    _errBPndOrChallVsTimePlot->setColor(colors[iColor]);

    dVector xVector, errVector;
    double dt = _simulator.getBinResolution();
    double saveTimeDuration  = _simulator.getDuration();
    double saveChallengeTime = _simulator.getChallengeTime();
    for (double time=_timeLowValue; time<=_timeHighValue; time += dt)
    {
        xVector.append(time);
        QString numberString;  numberString.setNum(time);
        if ( _fitChallenge->isChecked() )
        { // calculate the error in the challenge magnitude
            _challengeTime->setText(numberString);  changedChallengeTime();
            double truth = _simulator.getChallengeMag();
            double guess = getChallengeMagFromAnalysis();
            errVector.append(guess - truth);
        }
        else
        { // calculate the percent error in BPnd
//            _binDuration->setText(numberString);  changedTimeDuration();
            double truth = _simulator.getBP0();
            double guess = _PETRTM.getBP0InRun(0);
            errVector.append(percentageError(guess,truth));
            qDebug() << "change time to" << numberString << "with result" << guess;
        }
    }
    // restore the time duration and challenge time
    QString numberString;
    _binDuration->setText(numberString.setNum(saveTimeDuration));
    _challengeTime->setText(numberString.setNum(saveChallengeTime));
    changedChallengeTime();  // this will update the simulation and also the challenge time in analysis

    _errBPndOrChallVsTimePlot->setData(xVector,errVector);
    _errBPndOrChallVsTimePlot->conclude(0,true);
    _errBPndOrChallVsTimePlot->plotDataAndFit(true);
}
*/
void SimWindow::calculateBPndCurvesInThreads()
{
    qDebug() << "SimWindow::calculateBPndCurvesInThreads enter";
    updateLieDetectorProgress(10*_nThreads);
    qDebug() << "SimWindow::calculateBPndCurvesInThreads 1";

    _BP0Vector.clear();
    for (double BP0=_BPndLowValue; BP0<=_BPndHighValue; BP0 += _BPndStepValue)
        _BP0Vector.append(BP0);
    _errBPndMatrix.clear(); _errChallMatrix.clear(); _tau2RefMatrix.clear(); _errTau4Matrix.clear();
    qDebug() << "SimWindow::calculateBPndCurvesInThreads 2";

    /////////////////////////////////////////////////////////
    // Create the average volume.
    qRegisterMetaType<dVector>("dMatrix");
    QVector<lieDetector *> simSegment;
    simSegment.resize(_nThreads);
    for (int jThread=0; jThread<_nThreads; jThread++)
    {
        qDebug() << "SimWindow::calculateBPndCurvesInThreads jThread" << jThread;
        simSegment[jThread] = new lieDetector(_numberSimulationsPerThread, _BP0Vector, _simulator, _PETRTM);
        connect(simSegment[jThread], SIGNAL(progressLieDetector(int)), this, SLOT(updateLieDetectorProgress(int)));
        connect(simSegment[jThread], SIGNAL(finishedLieDetector(dMatrix,dMatrix,dMatrix,dMatrix)),this,SLOT(finishedLieDetectorOneThread(dMatrix,dMatrix,dMatrix,dMatrix)));
        QThreadPool::globalInstance()->start(simSegment[jThread]);
    }
    qDebug() << "SimWindow::calculateBPndCurvesInThreads exit";
}

void SimWindow::updateLieDetectorProgress(int iProgress)
{
    static int progressCounter=0;
    // qDebug() << "SimWindow::updateLieDetectorProgress enter" << progressCounter << "of" << _progressBar->maximum();
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

void SimWindow::finishedLieDetectorOneThread(dMatrix errBPnd, dMatrix errChall, dMatrix tau2Ref, dMatrix errTau4)
{
    qDebug() << "SimWindow::finishedLieDetectorOneThread enter";
    static int nThreadsFinished=0;

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
    qDebug() << "SimWindow::finishedLieDetectorOneThread 3";
    for (int jBP=0; jBP<nBP0Values; jBP++)
    {
        qDebug() << "sizes1[" << jBP << "] =" << errBPnd[jBP].size() << errChall[jBP].size() << tau2Ref[jBP].size() << errTau4[jBP].size();
        qDebug() << "sizes2[" << jBP << "] =" << _errBPndMatrix.size() << _errChallMatrix.size() << _tau2RefMatrix.size() << _errTau4Matrix.size();
        for ( int jSample=0; jSample<nSamples; jSample++)
        {
            _errBPndMatrix[jBP].append(errBPnd[jBP][jSample]);
            _errChallMatrix[jBP].append(errChall[jBP][jSample]);
            _tau2RefMatrix[jBP].append(tau2Ref[jBP][jSample]);
            _errTau4Matrix[jBP].append(errTau4[jBP][jSample]);
        }
    }
    qDebug() << "SimWindow::finishedLieDetectorOneThread 4";
    mutex.unlock();

    nThreadsFinished++;
    if ( nThreadsFinished == _nThreads )
    {
        _progressBar->reset();
        _progressBar->setMaximum(1);  // indicates that the bar is available
        nThreadsFinished = 0;  // reset for next time
        finishedLieDetectorAllThreads();
    }
    qDebug() << "SimWindow::finishedLieDetectorOneThread exit";
}

double SimWindow::calculateMean(dVector vec)
{
    qDebug() << "SimWindow::calculateMean enter" << vec.size();
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
    qDebug() << "SimWindow::finishedLieDetectorAllThreads enter" << nBP0Values << nSamples;

    dVector errBP, errChall, tau2Ref, errTau4, errBPSEM, errChallSEM, tau2RefSEM, errTau4SEM;
    // determine the means
    qDebug() << "SimWindow::finishedLieDetectorAllThreads 1";
    for (int jBP=0; jBP<nBP0Values; jBP++)
    {
        qDebug() << "append to errBP the value" << calculateMean(_errBPndMatrix[jBP]);
        errBP.append(calculateMean(_errBPndMatrix[jBP]));
        errChall.append(calculateMean(_errChallMatrix[jBP]));
        tau2Ref.append(calculateMean(_tau2RefMatrix[jBP]));
        errTau4.append(calculateMean(_errTau4Matrix[jBP]));

        errBPSEM.append(calculateStDev(errBP[jBP],    _errBPndMatrix[jBP])  / qSqrt(nSamples));
        errChallSEM.append(calculateStDev(errChall[jBP], _errChallMatrix[jBP]) / qSqrt(nSamples));
        tau2RefSEM.append(calculateStDev(tau2Ref[jBP],  _tau2RefMatrix[jBP])  / qSqrt(nSamples));
        errTau4SEM.append(calculateStDev(errTau4[jBP], _errTau4Matrix[jBP]) / qSqrt(nSamples));
    }
    qDebug() << "SimWindow::finishedLieDetectorAllThreads 2" << errBP.size() << errChall.size() << tau2Ref.size();


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
    qDebug() << "SimWindow::finishedLieDetectorAllThreads enter";
}
