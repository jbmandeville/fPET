#include <QtWidgets>
#include <QDebug>
#include <QVector>
#include <QObject>
#include <QWidget>
#include <QFileDialog>
#include <QString>
#include <QDateTime>
#include "simwindow.h"

///////////////////////////////////////
/// market average: 6% after inflation
/// market st dev : 24% after inflation (Gaussian)
///

SimWindow::SimWindow()
{
    FUNC_ENTER;
    _tabTimeSpace = new QTabWidget();
    createResultsPage();
    createPlanPage();
    _tabTimeSpace->addTab(_planPage, tr("plan"));
    _tabTimeSpace->addTab(_resultsPage, tr("results"));

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
    connect(aboutAppAct,     &QAction::triggered, this, &SimWindow::aboutApp);

    QAction *quitAction = openMenu->addAction(tr("&Quit"));
    // short-cuts and tooltips
    quitAction->setShortcut(Qt::ControlModifier + Qt::Key_Q);
    connect(quitAction, &QAction::triggered, this, &SimWindow::exitApp);

    QSize defaultWindowSize;
    QRect rec = QApplication::desktop()->screenGeometry();
//    QRect rec = QApplication::desktop()->screenGeometry();
    defaultWindowSize.setWidth(rec.width()*3/4);
    defaultWindowSize.setHeight(rec.height()*3/4);
    resize(defaultWindowSize);

    QStringList columnNames;
    QString errorString = utilIO::readTimeTableFile(":/My-Text/returns.txt",columnNames,_historicalReturnsTable);
    if ( !errorString.isEmpty() )
    {
        qInfo() << errorString;
        exitApp();
    }
    else
    {
//        int nTime = _historicalReturnsTable.size();
//        qInfo() << "nTime =" << nTime;
        /*
        int nColumns = _historicalReturnsTable[0].size();
        qInfo() << "nColumns =" << nColumns;
        for (int jt=0; jt<nTime; jt++)
            qInfo() << _historicalReturnsTable[jt][0] << _historicalReturnsTable[jt][1];
        */
    }

    // read tax brackets
    errorString = utilIO::readTimeTableFile(":/My-Text/brackets.txt",columnNames,_taxBrackets);
    if ( !errorString.isEmpty() )
    {
        qInfo() << errorString;
        exitApp();
    }
    else
    {
//        int nBrackets = _taxBrackets.size();
//        qInfo() << "nBrackets =" << nBrackets;
        /*
        int nColumns = _taxBrackets[0].size();
        qInfo() << "nColumns =" << nColumns;
        for (int jt=0; jt<nBrackets; jt++)
            qInfo() << _taxBrackets[jt][0] << _taxBrackets[jt][1];
*/
    }

    updateAnnuityWithdrawalRates();
    updateGraphsAndPanels();
}

QGroupBox *SimWindow::createJoeGroupBox()
{
    FUNC_ENTER;
    QString numberString;

    auto *retireAgeLabel    = new QLabel("retirement age");
    auto *ssAgeLabel        = new QLabel("SS starting age");
    auto *currentMoneyLabel = new QLabel("market amount");
    auto *FAAmountLabel     = new QLabel("FA amount (#1)");
    auto *RMDAgeLabel       = new QLabel("'RMD' age");
    auto *ageLimitLabel     = new QLabel("Age limit (Joe)");
    _retireAgeJoeString     = new QLineEdit(numberString.setNum(_retireAgeJoe));
    _ssAgeJoeString         = new QLineEdit(numberString.setNum(_ssAgeJoe));
    _marketMoneyJoeString   = new QLineEdit(numberString.setNum(_marketMoneyJoe));
    _annuityFAStartMoneyString = new QLineEdit(numberString.setNum(_annuityFAStartMoney));
    _ageRMDJoeString        = new QLineEdit(numberString.setNum(_ageRMDJoe));

    _ageLimitString = new QLineEdit(numberString.setNum(_ageLimit));
    connect(_ageLimitString, SIGNAL(editingFinished()), this, SLOT(changedAgeLimit()));

    _retireAgeJoeString->setFixedWidth(_editTextWidth);
    _ssAgeJoeString->setFixedWidth(_editTextWidth);
    _marketMoneyJoeString->setFixedWidth(_editTextWidth);
    _annuityFAStartMoneyString->setFixedWidth(_editTextWidth);
    _ageRMDJoeString->setFixedWidth(_editTextWidth);
    _ageLimitString->setFixedWidth(_editTextWidth);

    auto *mainLayout = new QGridLayout();
    mainLayout->setSpacing(0);
    int iRow=0;
    mainLayout->addWidget(retireAgeLabel,iRow,0);
    mainLayout->addWidget(_retireAgeJoeString,iRow++,1);
    mainLayout->addWidget(ssAgeLabel,iRow,0);
    mainLayout->addWidget(_ssAgeJoeString,iRow++,1);
    mainLayout->addWidget(currentMoneyLabel,iRow,0);
    mainLayout->addWidget(_marketMoneyJoeString,iRow++,1);
    mainLayout->addWidget(FAAmountLabel,iRow,0);
    mainLayout->addWidget(_annuityFAStartMoneyString,iRow++,1);
    mainLayout->addWidget(RMDAgeLabel,iRow,0);
    mainLayout->addWidget(_ageRMDJoeString,iRow++,1);
    mainLayout->addWidget(ageLimitLabel,iRow,0);
    mainLayout->addWidget(_ageLimitString,iRow++,1);

    QString title = QString("Joe: age %1").arg(_currentAgeJoe);
    auto *groupBox = new QGroupBox(title);
    groupBox->setLayout(mainLayout);
    connect(_retireAgeJoeString, SIGNAL(editingFinished()),         this, SLOT(changedRetireAgeJoe()));
    connect(_ssAgeJoeString, SIGNAL(editingFinished()),             this, SLOT(changedSSAgeJoe()));
    connect(_marketMoneyJoeString, SIGNAL(editingFinished()),       this, SLOT(changedCurrentMoneyJoe()));
    connect(_annuityFAStartMoneyString, SIGNAL(editingFinished()),  this, SLOT(changedFAStartingAmount()));
    connect(_ageRMDJoeString, SIGNAL(editingFinished()),            this, SLOT(changedRMDAgeJoe()));

    return groupBox;
}

QGroupBox *SimWindow::createEmiriGroupBox()
{
    FUNC_ENTER;
    QString numberString;

    auto *retireAgeLabel    = new QLabel("retirement age");
    auto *ssAgeLabel        = new QLabel("SS starting age");
    auto *currentMoneyLabel = new QLabel("market amount");
    auto *monthlyMoneyLabel = new QLabel("monthly addition");
    _retireAgeEmiriString       = new QLineEdit(numberString.setNum(_retireAgeEmiri));
    _ssAgeEmiriString           = new QLineEdit(numberString.setNum(_ssAgeEmiri));
    _marketMoneyEmiriString     = new QLineEdit(numberString.setNum(_marketMoneyEmiri));
    _monthlyAdditionEmiriString = new QLineEdit(numberString.setNum(_monthlyAdditionEmiri));

    _retireAgeEmiriString->setFixedWidth(_editTextWidth);
    _ssAgeEmiriString->setFixedWidth(_editTextWidth);
    _marketMoneyEmiriString->setFixedWidth(_editTextWidth);
    _monthlyAdditionEmiriString->setFixedWidth(_editTextWidth);

    auto *mainLayout = new QGridLayout();
    mainLayout->setSpacing(0);
    int iRow=0;
    mainLayout->addWidget(retireAgeLabel,iRow,0);
    mainLayout->addWidget(_retireAgeEmiriString,iRow++,1);
    mainLayout->addWidget(ssAgeLabel,iRow,0);
    mainLayout->addWidget(_ssAgeEmiriString,iRow++,1);
    mainLayout->addWidget(currentMoneyLabel,iRow,0);
    mainLayout->addWidget(_marketMoneyEmiriString,iRow++,1);
    mainLayout->addWidget(monthlyMoneyLabel,iRow,0);
    mainLayout->addWidget(_monthlyAdditionEmiriString,iRow++,1);

    QString title = QString("Emiri: age %1").arg(_currentAgeEmiri);
    auto *groupBox = new QGroupBox(title);
    groupBox->setLayout(mainLayout);

    connect(_retireAgeEmiriString,  SIGNAL(editingFinished()), this, SLOT(changedRetireAgeEmiri()));
    connect(_ssAgeEmiriString,      SIGNAL(editingFinished()), this, SLOT(changedSSAgeEmiri()));
    connect(_marketMoneyEmiriString, SIGNAL(editingFinished()), this, SLOT(changedCurrentMoneyEmiri()));
    connect(_monthlyAdditionEmiriString, SIGNAL(editingFinished()), this, SLOT(changedMonthlyAdditionEmiri()));

    return groupBox;
}

QGroupBox *SimWindow::createRatesSSGroupBox()
{
    FUNC_ENTER;
    QString numberString;

    auto *ssJoeLabel     = new QLabel("Joe");
    auto *ssLilikaLabel  = new QLabel("Lilika");
    auto *ssEmiriLabel     = new QLabel("Emiri");
    auto *ssSpousalLabel    = new QLabel("Spousal");

    _ssJoeLabel             = new QLabel(numberString.setNum(_ssJoe));
    _ssLilikaLabel          = new QLabel(numberString.setNum(_ssLilika));

    _ssEmiriLabel         = new QLabel(numberString.setNum(_ssEmiri));
    _ssSpousalLabel       = new QLabel(numberString.setNum(_ssSpousal));

    auto *ssCrisisYear = new QLabel("SS crisis year");
    auto *ssCrisisDrop = new QLabel("SS percent drop");
    _ssCrisisYearString  = new QLineEdit(numberString.setNum(_ssCrisisYear));
    _ssCrisisDropString  = new QLineEdit(numberString.setNum(_ssCrisisDrop));

    _ssCrisisYearString->setFixedWidth(_editTextWidth);
    _ssCrisisDropString->setFixedWidth(_editTextWidth);

    auto *mainLayout = new QGridLayout();
    int iRow=0;
    mainLayout->addWidget(ssJoeLabel,iRow,0);
    mainLayout->addWidget(_ssJoeLabel,iRow++,1);

    mainLayout->addWidget(ssLilikaLabel,iRow,0);
    mainLayout->addWidget(_ssLilikaLabel,iRow++,1);

    mainLayout->addWidget(ssEmiriLabel,iRow,0);
    mainLayout->addWidget(_ssEmiriLabel,iRow++,1);
    mainLayout->addWidget(ssSpousalLabel,iRow,0);
    mainLayout->addWidget(_ssSpousalLabel,iRow++,1);

    mainLayout->addWidget(ssCrisisYear,iRow,0);
    mainLayout->addWidget(_ssCrisisYearString,iRow++,1);
    mainLayout->addWidget(ssCrisisDrop,iRow,0);
    mainLayout->addWidget(_ssCrisisDropString,iRow++,1);
    mainLayout->setSpacing(0);
    _ssJoeLabel->setAlignment(Qt::AlignCenter);
    _ssLilikaLabel->setAlignment(Qt::AlignCenter);
    _ssEmiriLabel->setAlignment(Qt::AlignCenter);
    _ssSpousalLabel->setAlignment(Qt::AlignCenter);

    auto *groupBox = new QGroupBox("Social Security $/year");
    groupBox->setLayout(mainLayout);

    connect(_ssCrisisYearString, SIGNAL(editingFinished()), this, SLOT(changedSSCrisisYear()));
    connect(_ssCrisisDropString, SIGNAL(editingFinished()), this, SLOT(changedSSCrisisPercent()));

    return groupBox;
}

QGroupBox *SimWindow::createHealthCareGroupBox()
{
    FUNC_ENTER;
    QString numberString;

    auto *preRetireLabel     = new QLabel("Pre-retirement (family)");
    auto *privateAdultLabel  = new QLabel("private: adult");
    auto *privateChildLabel  = new QLabel("private: child");
    auto *extraMedicareLabel = new QLabel("Medicare supplement");

    _healthCarePreRetireString  = new QLineEdit(numberString.setNum(_healthCarePreRetire));
    _healthCarePrivateAdultString  = new QLineEdit(numberString.setNum(_healthCarePrivateAdult));
    _healthCarePrivateChildString  = new QLineEdit(numberString.setNum(_healthCarePrivateChild));
    _extraMedicareString  = new QLineEdit(numberString.setNum(_extraMedicare));

    auto *mainLayout = new QGridLayout();
    int iRow=0;
    mainLayout->addWidget(preRetireLabel,iRow,0);
    mainLayout->addWidget(_healthCarePreRetireString,iRow++,1);

    mainLayout->addWidget(privateAdultLabel,iRow,0);
    mainLayout->addWidget(_healthCarePrivateAdultString,iRow++,1);

    mainLayout->addWidget(privateChildLabel,iRow,0);
    mainLayout->addWidget(_healthCarePrivateChildString,iRow++,1);

    mainLayout->addWidget(extraMedicareLabel,iRow,0);
    mainLayout->addWidget(_extraMedicareString,iRow++,1);

    mainLayout->setSpacing(0);
    _ssJoeLabel->setAlignment(Qt::AlignCenter);
    _ssLilikaLabel->setAlignment(Qt::AlignCenter);
    _ssEmiriLabel->setAlignment(Qt::AlignCenter);
    _ssSpousalLabel->setAlignment(Qt::AlignCenter);

    auto *groupBox = new QGroupBox("Health Care costs ($k)");
    groupBox->setLayout(mainLayout);

    connect(_healthCarePreRetireString, SIGNAL(editingFinished()), this, SLOT(changedHCPreRetire()));
    connect(_healthCarePrivateAdultString, SIGNAL(editingFinished()), this, SLOT(changedHCPrivateAdult()));
    connect(_healthCarePrivateChildString, SIGNAL(editingFinished()), this, SLOT(changedHCPrivateChild()));
    connect(_extraMedicareString, SIGNAL(editingFinished()), this, SLOT(changedHCExtraMedicare()));

    return groupBox;
}

QGroupBox *SimWindow::createRetirementStartingGroupBox()
{
    FUNC_ENTER;

    auto *moneyJoeLabel    = new QLabel("$ Joe @ retire");
    auto *moneyEmiriLabel  = new QLabel("$ Emiri @ retire");

    QString numberString;
    _moneyEndWorkJoeLabel   = new QLabel("");
    _moneyEndWorkEmiriLabel = new QLabel("");

    auto *remainderAfterFALabel = new QLabel("$ after FA");
    auto *remainderAfterVALabel = new QLabel("$ after VA");
    _marketStartRetireJoeLabel    = new QLabel("");
    _marketRemainderStartJoeLabel = new QLabel("");

    auto *inflationLabel = new QLabel("% inflation");
    _inflationString = new QLineEdit(numberString.setNum(_inflation));
    _inflationString->setFixedWidth(_editTextWidth);
    connect(_inflationString,SIGNAL(editingFinished()),  this, SLOT(changedInflation()));

    auto *gridLayout = new QGridLayout();
    gridLayout->setSpacing(0);
    int iRow=0;
    gridLayout->addWidget(moneyJoeLabel,iRow,0);
    gridLayout->addWidget(_moneyEndWorkJoeLabel,iRow++,1);
    gridLayout->addWidget(moneyEmiriLabel,iRow,0);
    gridLayout->addWidget(_moneyEndWorkEmiriLabel,iRow++,1);

    gridLayout->addWidget(remainderAfterFALabel,iRow,0);
    gridLayout->addWidget(_marketStartRetireJoeLabel,iRow++,1);

    gridLayout->addWidget(remainderAfterVALabel,iRow,0);
    gridLayout->addWidget(_marketRemainderStartJoeLabel,iRow++,1);

    gridLayout->addWidget(inflationLabel,iRow,0);
    gridLayout->addWidget(_inflationString,iRow++,1);


    auto *groupBox = new QGroupBox("Outcomes:");
    groupBox->setLayout(gridLayout);

    return groupBox;
}

void SimWindow::createResultsPage()
{
    FUNC_ENTER;
    _resultsPage  = new QWidget();

    //////// The setup page status bar
    _statusBarJoePre = new QStatusBar;  // must be global so it doesn't go out of scope
    _statusBarJoePre->setStyleSheet("color:blue");
    _statusBarEmiriPre = new QStatusBar;  // must be global so it doesn't go out of scope
    _statusBarEmiriPre->setStyleSheet("color:blue");


    auto *layoutPlotJoe = new QVBoxLayout();
    layoutPlotJoe->addWidget(_statusBarJoePre);

    auto *plotLayout = new QHBoxLayout();
    plotLayout->addLayout(layoutPlotJoe);

    QHBoxLayout *fullLayout = new QHBoxLayout();
    fullLayout->addLayout(plotLayout);
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

    QToolBar *graphToolBar = new QToolBar("graph tool bar");
    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addAction(rescaleXYAction);

    fullLayout->setMenuBar(graphToolBar);

    // Make toolbar connections

    _resultsPage->setLayout(fullLayout);
}

void SimWindow::createPlanPage()
{
    FUNC_ENTER;
    _planPage = new QWidget();

    //////// The setup page status bar
    _statusBarIncome = new QStatusBar;  // must be global so it doesn't go out of scope
    _statusBarIncome->setStyleSheet("color:blue");
    _plotIncome   = new plotData(0);
    _plotIncome->setLabelXAxis("age Joe");
    _plotIncome->setLabelYAxis("income ($k)");
    _plotIncome->setQCStatusBar(_statusBarIncome);
    auto *layoutPlotIncome = new QVBoxLayout();
    layoutPlotIncome->addWidget(_plotIncome->getPlotSurface());
    layoutPlotIncome->addWidget(_statusBarIncome);
    layoutPlotIncome->setStretch(0,20);
    layoutPlotIncome->setStretch(1,1);

    _statusBarSpending = new QStatusBar;  // must be global so it doesn't go out of scope
    _statusBarSpending->setStyleSheet("color:blue");
    _statusBarSavings = new QStatusBar;  // must be global so it doesn't go out of scope
    _statusBarSavings->setStyleSheet("color:blue");
    _plotSpending   = new plotData(0);
    _plotSpending->setLabelXAxis("age Joe");
    _plotSpending->setLabelYAxis("spending ($k)");
    _plotSpending->setQCStatusBar(_statusBarSpending);
    _plotSavings   = new plotData(0);
    _plotSavings->setLabelXAxis("age Joe");
    _plotSavings->setLabelYAxis("savings ($k)");
    _plotSavings->setQCStatusBar(_statusBarSavings);
    auto *layoutPlotSpendingAndSavings = new QVBoxLayout();
    layoutPlotSpendingAndSavings->addWidget(_plotSpending->getPlotSurface());
    layoutPlotSpendingAndSavings->addWidget(_statusBarSpending);
    layoutPlotSpendingAndSavings->addWidget(_plotSavings->getPlotSurface());
    layoutPlotSpendingAndSavings->addWidget(_statusBarSavings);
    layoutPlotSpendingAndSavings->setStretch(0,20);
    layoutPlotSpendingAndSavings->setStretch(1,1);
    layoutPlotSpendingAndSavings->setStretch(2,20);
    layoutPlotSpendingAndSavings->setStretch(3,1);
    _plotSpendingAndSavingsWidget = new QWidget();
    _plotSpendingAndSavingsWidget->setLayout(layoutPlotSpendingAndSavings);

    _statusBarWR = new QStatusBar;  // must be global so it doesn't go out of scope
    _statusBarWR->setStyleSheet("color:blue");
    _plotWR   = new plotData(0);
    _plotWR->setLabelXAxis("age Joe");
    _plotWR->setLabelYAxis("withdrawal rate");
    _plotWR->setQCStatusBar(_statusBarWR);
    auto *layoutPlotWR = new QVBoxLayout();
    layoutPlotWR->addWidget(_plotWR->getPlotSurface());
    layoutPlotWR->addWidget(_statusBarWR);
    layoutPlotWR->setStretch(0,20);
    layoutPlotWR->setStretch(1,1);

    _statusBarMarketReturnDist = new QStatusBar;  // must be global so it doesn't go out of scope
    _statusBarMarketReturnDist->setStyleSheet("color:blue");
    _plotMarketReturnDist  = new plotData(0);
    _plotMarketReturnDist->setLabelXAxis("percent return");
    _plotMarketReturnDist->setLabelYAxis("probability");
    _plotMarketReturnDist->setQCStatusBar(_statusBarMarketReturnDist);
    auto *layoutPlotMRD = new QVBoxLayout();
    layoutPlotMRD->addWidget(_plotMarketReturnDist->getPlotSurface());
    layoutPlotMRD->addWidget(_statusBarMarketReturnDist);
    layoutPlotMRD->setStretch(0,20);
    layoutPlotMRD->setStretch(1,1);
    _MRDWidget = new QWidget();
    _MRDWidget->setLayout(layoutPlotMRD);

    _statusBarMarketReturns = new QStatusBar;  // must be global so it doesn't go out of scope
    _statusBarMarketReturns->setStyleSheet("color:blue");
    _plotMarketReturns  = new plotData(0);
    _plotMarketReturns->setLabelXAxis("age Joe");
    _plotMarketReturns->setLabelYAxis("return");
    _plotMarketReturns->setQCStatusBar(_statusBarMarketReturns);
    auto *layoutPlotMR = new QVBoxLayout();
    layoutPlotMR->addWidget(_plotMarketReturns->getPlotSurface());
    layoutPlotMR->addWidget(_statusBarMarketReturns);
    layoutPlotMR->setStretch(0,20);
    layoutPlotMR->setStretch(1,1);
    _MRWidget = new QWidget();
    _MRWidget->setLayout(layoutPlotMR);

    auto *layoutPlot2 = new QVBoxLayout();
    layoutPlot2->addWidget(_MRDWidget);
    layoutPlot2->addWidget(_MRWidget);
    layoutPlot2->addLayout(layoutPlotWR);
    _plotWRandMarketReturnsWidget = new QWidget();
    _plotWRandMarketReturnsWidget->setLayout(layoutPlot2);
    _MRDWidget->setVisible(true);
    _MRWidget->setVisible(false);

    auto *plotLayout = new QHBoxLayout();
    plotLayout->addLayout(layoutPlotIncome);
    plotLayout->addWidget(_plotSpendingAndSavingsWidget);
    plotLayout->addWidget(_plotWRandMarketReturnsWidget);

    auto *leftLayout = new QVBoxLayout();
    leftLayout->addWidget(createJoeGroupBox());
    leftLayout->addWidget(createEmiriGroupBox());
    leftLayout->addWidget(createRatesSSGroupBox());
    leftLayout->addWidget(createHealthCareGroupBox());
    leftLayout->addWidget(createRetirementStartingGroupBox());
    updateSS();

    auto *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(createSpendingGroupBox());
    rightLayout->addWidget(createApproachJoeGroupBox());
    rightLayout->addWidget(createMarketRatesAndReturnsGroupBox());
    rightLayout->addWidget(createAnnuityGroupBox());
    rightLayout->addWidget(createAvIncomeGroupBox());

    QHBoxLayout *fullLayout = new QHBoxLayout();
    fullLayout->addLayout(leftLayout);
    fullLayout->addLayout(plotLayout);
    fullLayout->addLayout(rightLayout);
    fullLayout->setStretch(0,1);
    fullLayout->setStretch(1,100);
    fullLayout->setStretch(2,1);

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

    QFrame* separator1 = new QFrame(_planPage);
    separator1->setFrameShape(QFrame::VLine);
    separator1->setLineWidth(3);
    separator1->setFrameShadow(QFrame::Raised);
    QFrame* separator2 = new QFrame(_planPage);
    separator2->setFrameShape(QFrame::VLine);
    separator2->setLineWidth(3);
    separator2->setFrameShadow(QFrame::Raised);
    QFrame* separator3 = new QFrame(_planPage);
    separator3->setFrameShape(QFrame::VLine);
    separator3->setLineWidth(3);
    separator3->setFrameShadow(QFrame::Raised);

    auto *showOption = new QLabel("show:");
    _radioShowAll = new QRadioButton("all");
    _radioShowIncome = new QRadioButton("income only");
    auto *showGroup = new QButtonGroup(_planPage);
    showGroup->addButton(_radioShowAll);
    showGroup->addButton(_radioShowIncome);
    _radioShowAll->setChecked(true);
    connect(_radioShowAll,   SIGNAL(clicked(bool)),  this, SLOT(changedShowPlots()));
    connect(_radioShowIncome,SIGNAL(clicked(bool)),  this, SLOT(changedShowPlots()));

    auto *incomeType = new QLabel("final income:");
    _radioTotalIncome = new QRadioButton("before tax $");
    _radioAfterTaxHealthMortgage = new QRadioButton("after tax $");
    _radioDiscretionaryIncome = new QRadioButton("discretionary $");
    auto *incomeRedGroup = new QButtonGroup(_planPage);
    incomeRedGroup->addButton(_radioTotalIncome);
    incomeRedGroup->addButton(_radioAfterTaxHealthMortgage);
    incomeRedGroup->addButton(_radioDiscretionaryIncome);
    _radioTotalIncome->setChecked(true);
    connect(_radioTotalIncome,SIGNAL(clicked(bool)), this, SLOT(updateGraphsAndPanels()));
    connect(_radioAfterTaxHealthMortgage,SIGNAL(clicked(bool)),   this, SLOT(updateGraphsAndPanels()));
    connect(_radioDiscretionaryIncome,SIGNAL(clicked(bool)),   this, SLOT(updateGraphsAndPanels()));

    auto *incomeShow = new QLabel("show income:");
    _radioShowAllIncome = new QRadioButton("all");
    _radioShowFixedVariableIncome = new QRadioButton("fixed/variable");
    _radioShowIncomeCategories = new QRadioButton("categories");
    auto *incomeTypeGroup = new QButtonGroup(_planPage);
    incomeTypeGroup->addButton(_radioShowAllIncome);
    incomeTypeGroup->addButton(_radioShowFixedVariableIncome);
    incomeTypeGroup->addButton(_radioShowIncomeCategories);
    _radioShowAllIncome->setChecked(true);
    connect(_radioShowAllIncome,SIGNAL(clicked(bool)),    this, SLOT(updateGraphsAndPanels()));
    connect(_radioShowFixedVariableIncome,SIGNAL(clicked(bool)),  this, SLOT(updateGraphsAndPanels()));
    connect(_radioShowIncomeCategories,SIGNAL(clicked(bool)), this, SLOT(updateGraphsAndPanels()));

    auto *incomeShowLayout = new QHBoxLayout();
    incomeShowLayout->addWidget(incomeShow);
    incomeShowLayout->addWidget(_radioShowAllIncome);
    incomeShowLayout->addWidget(_radioShowFixedVariableIncome);
    incomeShowLayout->addWidget(_radioShowIncomeCategories);
    _incomeShowWidget = new QWidget();
    _incomeShowWidget->setLayout(incomeShowLayout);

    auto *incomeTypeLayout = new QHBoxLayout();
    incomeTypeLayout->addWidget(incomeType);
    incomeTypeLayout->addWidget(_radioTotalIncome);
    incomeTypeLayout->addWidget(_radioAfterTaxHealthMortgage);
    incomeTypeLayout->addWidget(_radioDiscretionaryIncome);
    _incomeTypeWidget = new QWidget();
    _incomeTypeWidget->setLayout(incomeTypeLayout);

    QToolBar *graphToolBar = new QToolBar("graph tool bar");
    graphToolBar->addAction(dragXAction);
    graphToolBar->addAction(dragYAction);
    graphToolBar->addAction(rescaleXYAction);
    graphToolBar->addWidget(separator1);
    graphToolBar->addWidget(showOption);
    graphToolBar->addWidget(_radioShowAll);
    graphToolBar->addWidget(_radioShowIncome);
    graphToolBar->addWidget(separator2);
    graphToolBar->addWidget(_incomeTypeWidget);
    graphToolBar->addWidget(separator3);
    graphToolBar->addWidget(_incomeShowWidget);

    fullLayout->setMenuBar(graphToolBar);

    QCPRange xRange, yRange;
    xRange.lower = 60;   xRange.upper = 85;
    yRange.lower = 0.;          yRange.upper = 300.;
    _plotIncome->setAxisRanges(xRange, yRange);

    xRange.lower = 60.;     xRange.upper = _ageLimit;
    yRange.lower = 0.;      yRange.upper = 5000.;
    _plotSavings->setAxisRanges(xRange, yRange);

    yRange.lower = -10.;      yRange.upper = 300.;
    _plotSpending->setAxisRanges(xRange, yRange);

    xRange.lower = 60.;     xRange.upper = 85.;
    yRange.lower = 0.;      yRange.upper = 7.;
    _plotWR->setAxisRanges(xRange, yRange);

    yRange.lower = 0.;      yRange.upper = 100.;
    _plotSpending->setAxisRanges(xRange, yRange);

    // Make toolbar connections
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotIncome,  SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotIncome,  SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotIncome,  SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotIncome,  SLOT(setSelectPoints(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotSpending, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotSpending, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotSpending, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotSpending, SLOT(setSelectPoints(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotSavings, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotSavings, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotSavings, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotSavings, SLOT(setSelectPoints(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotWR,      SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotWR,      SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotWR,      SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotWR,      SLOT(setSelectPoints(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotMarketReturnDist, SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotMarketReturnDist, SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotMarketReturnDist, SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotMarketReturnDist, SLOT(setSelectPoints(bool)));
    connect(dragXAction,     SIGNAL(triggered(bool)), _plotMarketReturns,    SLOT(setXZoom()));
    connect(dragYAction,     SIGNAL(triggered(bool)), _plotMarketReturns,    SLOT(setYZoom()));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotMarketReturns,    SLOT(autoScale(bool)));
    connect(rescaleXYAction, SIGNAL(triggered(bool)), _plotMarketReturns,    SLOT(setSelectPoints(bool)));

    _planPage->setLayout(fullLayout);
}

QGroupBox *SimWindow::createSpendingGroupBox()
{
    FUNC_ENTER;

    QString numberString;

    auto *year0Label = new QLabel("$/year at start");
    _spendingYear0String = new QLineEdit(numberString.setNum(_spendingYear0));
    _spendingYear0String->setFixedWidth(_editTextWidth);
    _spendingYear0String->setToolTip("Spending year 0 including mortgage");
    connect(_spendingYear0String,  SIGNAL(editingFinished()), this, SLOT(changedSpendingYear0()));

    auto *percentIncrease = new QLabel("%/year increase");
    _spendingPerYearIncreasePercentString = new QLineEdit(numberString.setNum(_spendingPerYearIncreasePercent));
    _spendingPerYearIncreasePercentString->setFixedWidth(_editTextWidth);
    connect(_spendingPerYearIncreasePercentString,  SIGNAL(editingFinished()), this, SLOT(changedSpendingSlope()));

    auto *maxSpendLabel = new QLabel("max discretionary");
    _maxDiscretionarySpendingString = new QLineEdit(numberString.setNum(_maxDiscretionarySpending));
    _maxDiscretionarySpendingString->setFixedWidth(_editTextWidth);
    _maxDiscretionarySpendingString->setToolTip("Spending year 0 including mortgage");
    connect(_maxDiscretionarySpendingString,  SIGNAL(editingFinished()), this, SLOT(changedMaxSpending()));

    auto *maxBrokerageLabel = new QLabel("max brokerage");
    _maxBrokerageString = new QLineEdit(numberString.setNum(_maxBrokerage));
    _maxBrokerageString->setFixedWidth(_editTextWidth);
    _maxBrokerageString->setToolTip("Spending year 0 including mortgage");
    connect(_maxBrokerageString,  SIGNAL(editingFinished()), this, SLOT(changedMaxBrokerage()));

    auto *gridLayout = new QGridLayout();
    gridLayout->setSpacing(3);
    int iRow=0;
    gridLayout->addWidget(year0Label,iRow,0);
    gridLayout->addWidget(_spendingYear0String,iRow++,1);

    gridLayout->addWidget(percentIncrease,iRow,0);
    gridLayout->addWidget(_spendingPerYearIncreasePercentString,iRow++,1);

    gridLayout->addWidget(maxSpendLabel,iRow,0);
    gridLayout->addWidget(_maxDiscretionarySpendingString,iRow++,1);

    gridLayout->addWidget(maxBrokerageLabel,iRow,0);
    gridLayout->addWidget(_maxBrokerageString,iRow++,1);

    auto *groupBox = new QGroupBox("Spending after tax");
    groupBox->setLayout(gridLayout);
    return groupBox;
}

QGroupBox *SimWindow::createApproachJoeGroupBox()
{
    FUNC_ENTER;

    QString numberString;

    auto *radioPercentWithdraw = new QRadioButton("% + RMDs");
    _radioTargetIncome         = new QRadioButton("$ target");
    _radioGuardRails           = new QRadioButton("guard-rails");
    connect(radioPercentWithdraw,SIGNAL(clicked(bool)), this, SLOT(changedTargetIncomeRadioButton()));
    connect(_radioTargetIncome,SIGNAL(clicked(bool)), this, SLOT(changedTargetIncomeRadioButton()));
    connect(_radioGuardRails,SIGNAL(clicked(bool)), this, SLOT(changedTargetIncomeRadioButton()));

    auto *buttonGroup = new QButtonGroup();
    buttonGroup->addButton(radioPercentWithdraw);
    buttonGroup->addButton(_radioTargetIncome);
    buttonGroup->addButton(_radioGuardRails);
    radioPercentWithdraw->setChecked(true);

    _incomeFloorTargetLabel = new QLabel("floor");
    _incomeFloorTargetString  = new QLineEdit(numberString.setNum(_incomeFloorTarget));
    _incomeFloorTargetString->setFixedWidth(_editTextWidth);
    connect(_incomeFloorTargetString, SIGNAL(editingFinished()), this, SLOT(changedTargetIncome()));

    _incomeFloorCheckBox = new QCheckBox("enforce floor");
    _incomeFloorCheckBox->setVisible(false);
    _incomeFloorCheckBox->setChecked(true);
    connect(_incomeFloorCheckBox,  SIGNAL(toggled(bool)), this, SLOT(updateGraphsAndPanels()));

    _marketWithDrawLabel = new QLabel("withdrawal rate");
    _marketWithDrawLabel->setToolTip("based upon Joe's $$ only");

    _marketWRString = new QLineEdit(numberString.setNum(_marketWR));
    _marketWRString->setFixedWidth(_editTextWidth);
    _marketWRString->setToolTip("based upon Joe's $$ only");
    connect(_marketWRString, SIGNAL(editingFinished()), this, SLOT(changedMarketWR()));

    auto *bufferLabel = new QLabel("ca$h buffer");
    _bufferString = new QLineEdit(numberString.setNum(_bufferDesired));
    connect(_bufferString,  SIGNAL(editingFinished()), this, SLOT(changedCashBuffer()));

    auto *spiaDetailsLabel  = new QLabel("amount & age");

    auto *radioNoSpia = new QRadioButton("no SPIA");
    _radioSpia5   = new QRadioButton("5-year");
    _radioSpia10  = new QRadioButton("10-year");
    connect(radioNoSpia,  SIGNAL(clicked(bool)), this, SLOT(updateGraphsAndPanels()));
    connect(_radioSpia5,  SIGNAL(clicked(bool)), this, SLOT(setSpia5()));
    connect(_radioSpia10, SIGNAL(clicked(bool)), this, SLOT(setSpia10()));

    auto *spiaGroup = new QButtonGroup();
    spiaGroup->addButton(radioNoSpia);
    spiaGroup->addButton(_radioSpia5);
    spiaGroup->addButton(_radioSpia10);
    radioNoSpia->setChecked(true);

    _spiaAmountString = new QLineEdit(numberString.setNum(_spiaAmount));
    _spiaAgeString    = new QLineEdit(numberString.setNum(_spiaAge));
    _spiaAmountString->setFixedWidth(_editTextWidth);
    _spiaAgeString->setFixedWidth(_editTextWidth);
    connect(_spiaAmountString,SIGNAL(editingFinished()), this, SLOT(changedSpiaAmount()));
    connect(_spiaAgeString,SIGNAL(editingFinished()),    this, SLOT(changedSpiaAge()));

    auto *gridLayout = new QGridLayout();
    int iRow=0;
    gridLayout->addWidget(radioPercentWithdraw,iRow,0);
    gridLayout->addWidget(_radioTargetIncome,iRow,1);
    gridLayout->addWidget(_radioGuardRails,iRow++,2);

    gridLayout->addWidget(_incomeFloorTargetLabel,iRow,0);
    gridLayout->addWidget(_incomeFloorTargetString,iRow++,1);

    gridLayout->addWidget(_marketWithDrawLabel,iRow,0);
    gridLayout->addWidget(_marketWRString,iRow,1);
    gridLayout->addWidget(_incomeFloorCheckBox,iRow++,2);

    gridLayout->addWidget(bufferLabel,iRow,0);
    gridLayout->addWidget(_bufferString,iRow++,1);

    gridLayout->addWidget(radioNoSpia,iRow,0);
    gridLayout->addWidget(_radioSpia5,iRow,1);
    gridLayout->addWidget(_radioSpia10,iRow++,2);
    gridLayout->addWidget(spiaDetailsLabel,iRow,0);
    gridLayout->addWidget(_spiaAmountString,iRow,1);
    gridLayout->addWidget(_spiaAgeString,iRow++,2);

    gridLayout->setSpacing(3);

    _lastSavings = new QLabel("last savings =");
    _successRate = new QLabel("success rate =");

    auto *outcomeLayout = new QVBoxLayout();
    outcomeLayout->addWidget(_lastSavings);
    outcomeLayout->addWidget(_successRate);
    _lastSavings->setVisible(false);
    _successRate->setVisible(false);

    auto *mainLayout = new QVBoxLayout();
    mainLayout->addLayout(gridLayout);
    mainLayout->addLayout(outcomeLayout);

    setIncomeTargetFloorVisibility();

    auto *groupBox = new QGroupBox("withdrawal approach for JOE:");
    groupBox->setLayout(mainLayout);
    return groupBox;
}

QGroupBox *SimWindow::createMarketRatesAndReturnsGroupBox()
{
    FUNC_ENTER;

    QString numberString;

    auto *radioFixedReturn = new QRadioButton("fixed return");
    _radioRand1    = new QRadioButton("rand 1x");
    _radioRand     = new QRadioButton("rand");
    connect(radioFixedReturn,SIGNAL(clicked(bool)), this, SLOT(changedRandomization()));
    connect(_radioRand1,     SIGNAL(clicked(bool)), this, SLOT(changedRandomization()));
    connect(_radioRand,      SIGNAL(clicked(bool)), this, SLOT(changedRandomization()));

    auto *buttonG1 = new QButtonGroup();
    buttonG1->addButton(radioFixedReturn);
    buttonG1->addButton(_radioRand1);
    buttonG1->addButton(_radioRand);
    radioFixedReturn->setChecked(true);

    _radioHistoricalMarket = new QRadioButton("historical");
    _radioHistoricalMarket->setVisible(false);
    _radioHistoricalMarket->setToolTip("Use historical rates after inflation");
    connect(_radioHistoricalMarket,SIGNAL(clicked(bool)), this, SLOT(changedMarketRandomization()));

    _radioRandomMarket = new QRadioButton("random");
    _radioRandomMarket->setVisible(false);
    connect(_radioRandomMarket,SIGNAL(clicked(bool)), this, SLOT(changedMarketRandomization()));

    auto *buttonG2 = new QButtonGroup();
    buttonG2->addButton(_radioHistoricalMarket);
    buttonG2->addButton(_radioRandomMarket);
    _radioHistoricalMarket->setChecked(true);

    _lockMarketReturns = new QCheckBox("lock");
    _lockMarketReturns->setVisible(false);
    connect(_lockMarketReturns,SIGNAL(clicked(bool)), this, SLOT(updateGraphsAndPanels()));

    _marketReturnsRateLabel   = new QLabel("average return");
    _marketReturnsRateString  = new QLineEdit(numberString.setNum(_marketReturnsRate));
    _marketReturnsRateString->setFixedWidth(_editTextWidth);
    _marketReturnsRateString->setToolTip("Average market return AFTER INFLATION");
    connect(_marketReturnsRateString, SIGNAL(editingFinished()), this, SLOT(changedMarketRate()));

    _percentBondsLabel = new QLabel("% bond");
    _percentBondsString = new QLineEdit(numberString.setNum(_percentBonds));
    _percentBondsString->setFixedWidth(_editTextWidth);
    connect(_percentBondsString, SIGNAL(editingFinished()), this, SLOT(changedBondPercentage()));

    auto *marketFeesLabel = new QLabel("fees");
    _marketFeesString  = new QLineEdit(numberString.setNum(_marketFees));
    _marketFeesString->setFixedWidth(_editTextWidth);
    connect(_marketFeesString, SIGNAL(editingFinished()), this, SLOT(changedMarketFees()));

    auto *marketPerYearLabel    = new QLabel("av/year");
    _marketIncomePerYearAvLabel = new QLabel("");

    auto *gridLayout = new QGridLayout();
    int iRow=0;
    gridLayout->addWidget(radioFixedReturn,iRow,0);
    gridLayout->addWidget(_radioRand1,iRow,1);
    gridLayout->addWidget(_radioRand,iRow++,2);

    gridLayout->addWidget(_radioHistoricalMarket,iRow,0);
    gridLayout->addWidget(_radioRandomMarket,iRow,1);
    gridLayout->addWidget(_lockMarketReturns,iRow++,2);

    gridLayout->addWidget(_marketReturnsRateLabel,iRow,0);
    gridLayout->addWidget(_marketReturnsRateString,iRow++,1);

    gridLayout->addWidget(_percentBondsLabel,iRow,0);
    gridLayout->addWidget(_percentBondsString,iRow++,1);

    gridLayout->addWidget(marketFeesLabel,iRow,0);
    gridLayout->addWidget(_marketFeesString,iRow++,1);

    gridLayout->addWidget(marketPerYearLabel,iRow,0);
    gridLayout->addWidget(_marketIncomePerYearAvLabel,iRow++,1);
    gridLayout->setSpacing(0);
    auto *groupBox = new QGroupBox("market");
    groupBox->setLayout(gridLayout);

    return groupBox;
}

QGroupBox *SimWindow::createAnnuityGroupBox()
{
    FUNC_ENTER;

    QString numberString;

    auto *annuityStartAgeLabel   = new QLabel("age start #1/#2");

    _annuityFAStartAgeString     = new QLineEdit(numberString.setNum(_annuityFAStartAge));
    _annuityFAStartAgeString->setFixedWidth(_editTextWidth);
    _annuityFAStartAgeString->setToolTip("age to start Joe's FA");
    connect(_annuityFAStartAgeString, SIGNAL(editingFinished()), this, SLOT(changedFAAnnuityAge()));

    _annuity2StartAgeString     = new QLineEdit(numberString.setNum(_annuity2StartAge));
    _annuity2StartAgeString->setFixedWidth(_editTextWidth);
    _annuity2StartAgeString->setToolTip("age to start Joe's VA");
    connect(_annuity2StartAgeString, SIGNAL(editingFinished()), this, SLOT(changedAnnuity2Age()));

    _annuity2IsFACheckBox = new QCheckBox("annuity #2 = FA");
    connect(_annuity2IsFACheckBox,SIGNAL(clicked(bool)), this, SLOT(updateGraphsAndPanels()));

    auto *annuitizedAmountLabel  = new QLabel("annuity $");
    auto *annuityFeesLabel       = new QLabel("fees");
    auto *annuityWRLabel         = new QLabel("start rates");
    //    auto *annuityHurdleLabel   = new QLabel("VA hurdle");
    auto *annuitizedPercentLabel = new QLabel("annuity %");
    auto *annuityPerYearLabel    = new QLabel("starting & av/year");
    _annuitize2AmountLabel       = new QLabel("");
    _annuitizeFAAmountLabel      = new QLabel("");
    _annuityIncomePerYearAvLabel = new QLabel("");
    _annuityFAIncomePerYearAvLabel = new QLabel("");

    _annuity2FeesString = new QLineEdit(numberString.setNum(_annuity2Fees));
    _annuity2FeesString->setFixedWidth(_editTextWidth);
    connect(_annuity2FeesString, SIGNAL(editingFinished()), this, SLOT(changedAnnuity2Fees()));

    _annuityFAFeesString = new QLineEdit(numberString.setNum(_annuityFAFees));
    _annuityFAFeesString->setFixedWidth(_editTextWidth);
    connect(_annuity2FeesString, SIGNAL(editingFinished()), this, SLOT(changedAnnuityFAFees()));

    _annuity2WRStartingString = new QLineEdit(numberString.setNum(_annuity2WRStarting));
    _annuity2WRStartingString->setFixedWidth(_editTextWidth);
    connect(_annuity2WRStartingString, SIGNAL(editingFinished()), this, SLOT(changedAnnuity2WRStarting()));

    _annuityFAWRStartingString = new QLineEdit(numberString.setNum(_annuityFAWRStarting));
    _annuityFAWRStartingString->setFixedWidth(_editTextWidth);
    connect(_annuityFAWRStartingString, SIGNAL(editingFinished()), this, SLOT(changedAnnuityFAWRStarting()));

    _annuity2PercentageString = new QLineEdit(numberString.setNum(_annuity2Percentage));
    _annuity2PercentageString->setFixedWidth(_editTextWidth);
    _annuity2PercentageString->setToolTip("percentage of Joe's market+buffer $$");
    connect(_annuity2PercentageString, SIGNAL(editingFinished()), this, SLOT(changedAnnuitized2Percent()));

    _annuity2HurdleRateString = new QLineEdit(numberString.setNum(_annuity2HurdleRate));
    _annuity2HurdleRateString->setFixedWidth(_editTextWidth);
    _annuity2HurdleRateString->setEnabled(false);
    connect(_annuity2HurdleRateString, SIGNAL(editingFinished()), this, SLOT(changedAnnuityHurdleRate()));

    auto *demiseLabel = new QLabel("demise/reduction");
    _firstDemiseAgeString = new QLineEdit(numberString.setNum(_firstDemiseAge));
    _firstDemiseAgeString->setFixedWidth(_editTextWidth);
    connect(_firstDemiseAgeString, SIGNAL(editingFinished()), this, SLOT(changedFirstDemiseAge()));

    _reducedBenefitPercentString = new QLineEdit(numberString.setNum(_reducedBenefitPercent));
    _reducedBenefitPercentString->setFixedWidth(_editTextWidth);
    connect(_reducedBenefitPercentString, SIGNAL(editingFinished()), this, SLOT(changedReducedBenefitPercent()));

    auto *gridLayout = new QGridLayout();
    int iRow=0;
    gridLayout->addWidget(annuityStartAgeLabel,iRow,0);
    gridLayout->addWidget(_annuityFAStartAgeString,iRow,1);
    gridLayout->addWidget(_annuity2StartAgeString,iRow++,2);

    gridLayout->addWidget(_annuity2IsFACheckBox,iRow,0);
    gridLayout->addWidget(annuitizedPercentLabel,iRow,1);
    gridLayout->addWidget(_annuity2PercentageString,iRow++,2);

    gridLayout->addWidget(annuitizedAmountLabel,iRow,0);
    gridLayout->addWidget(_annuitizeFAAmountLabel,iRow,1);
    gridLayout->addWidget(_annuitize2AmountLabel,iRow++,2);

    gridLayout->addWidget(annuityFeesLabel,iRow,0);
    gridLayout->addWidget(_annuityFAFeesString,iRow,1);
    gridLayout->addWidget(_annuity2FeesString,iRow++,2);

    gridLayout->addWidget(annuityWRLabel,iRow,0);
    gridLayout->addWidget(_annuityFAWRStartingString,iRow,1);
    gridLayout->addWidget(_annuity2WRStartingString,iRow++,2);

    //    gridLayout->addWidget(annuityHurdleLabel,iRow,0);
    //    gridLayout->addWidget(_annuity2HurdleRateString,iRow++,1);

    gridLayout->addWidget(demiseLabel,iRow,0);
    gridLayout->addWidget(_firstDemiseAgeString,iRow,1);
    gridLayout->addWidget(_reducedBenefitPercentString,iRow++,2);

    gridLayout->addWidget(annuityPerYearLabel,iRow,0);
    gridLayout->addWidget(_annuityIncomePerYearAvLabel,iRow++,1);
    //    gridLayout->setSpacing(0);

    auto *groupBox = new QGroupBox("Lifetime annuity: VA and FA");
    groupBox->setLayout(gridLayout);

    return groupBox;
}

QGroupBox *SimWindow::createAvIncomeGroupBox()
{
    FUNC_ENTER;

    auto *label1 = new QLabel("$$ age 65-74");
    auto *label2 = new QLabel("$$ age 75-84");
    auto *label3 = new QLabel("$$ age 85-94");
    _incomeAverage65_74 = new QLabel("");
    _incomeAverage75_84 = new QLabel("");
    _incomeAverage85_94 = new QLabel("");

    auto *gridLayout = new QGridLayout();
    gridLayout->addWidget(label1,0,0);
    gridLayout->addWidget(_incomeAverage65_74,0,1);
    gridLayout->addWidget(label2,1,0);
    gridLayout->addWidget(_incomeAverage75_84,1,1);
    gridLayout->addWidget(label3,2,0);
    gridLayout->addWidget(_incomeAverage85_94,2,1);
    gridLayout->setSpacing(0);

    _averageIncomeBox = new QGroupBox("average income");
    _averageIncomeBox->setLayout(gridLayout);
    return _averageIncomeBox;
}

void SimWindow::plotMarketReturns()
{
    FUNC_ENTER;
    int iRange = 60;
    double mean  = _marketReturnsRate;
    double sigma = _stockReturnsSigma;
    double sigma2 = SQR(sigma);

    dVector returnsVectorPara;
    dVector probVectorPara;
    double sum = 0.;
    for (int j=-iRange; j<=iRange; j++)
    {
        double percent = j;
        double prob = 0.;
        if ( _radioRand->isChecked() )
            prob = 1./sigma/2./PI * qExp(-SQR(mean-percent)/2./sigma2);
        else if ( qAbs(percent-_marketReturnsRate) <= 0.5 )
            prob = 1.;
        returnsVectorPara.append(percent);
        probVectorPara.append(prob);
        sum += prob;
    }
    for (int j=0; j<probVectorPara.size(); j++)
        probVectorPara[j] /= sum;

    _plotMarketReturnDist->init();
    _plotMarketReturnDist->addCurve(0,"returns over inflation");
    FUNC_INFO << 1;
    _plotMarketReturnDist->setData(returnsVectorPara, probVectorPara);
    _plotMarketReturnDist->setLineThickness(2);
    _plotMarketReturnDist->setColor(Qt::red);

    if ( _radioRand1->isChecked() )
    {
        _plotMarketReturns->init();
        int iStartYear = _yearStartHistoricalMarket + 1928;
        int iEndYear   = iStartYear + (_ageLimit - _ageJoeVector.first());
        if ( iEndYear > 2024 ) iEndYear = 1928 + (iEndYear - 2024);
        QString startingYear = "random";
        if ( _radioHistoricalMarket->isChecked() )
            startingYear = QString("%1 - %2").arg(iStartYear).arg(iEndYear);
        _plotMarketReturns->addCurve(0,startingYear);
        FUNC_INFO << 3;
        _plotMarketReturns->setData(_ageJoeVector, _stockBondReturnsVector);
 //       _plotMarketReturns->setPointSize(2);
        _plotMarketReturns->setLineThickness(2);
        _plotMarketReturns->conclude(0,true);
        _plotMarketReturns->plotDataAndFit(true);
        _plotMarketReturns->showLegend();
    }
    else if ( _radioRand->isChecked() )
    {
        _plotMarketReturnDist->addCurve(0,"returns over inflation");
        FUNC_INFO << 2;
        _plotMarketReturnDist->setData(returnsVectorPara, _marketReturnsProbability);
        _plotMarketReturnDist->setLineThickness(2);
        _plotMarketReturnDist->setColor(Qt::black);
    }

    _plotMarketReturnDist->conclude(0,true);
    _plotMarketReturnDist->plotDataAndFit(true);
    FUNC_EXIT;
}

void SimWindow::plotResults()
{
    FUNC_ENTER;
    plotWithdrawalRate();
    plotIncome();
    plotSpending();
    plotSavings();
    FUNC_EXIT;
}

void SimWindow::plotWithdrawalRate()
{
    FUNC_ENTER;
    _plotWR->init();
    _plotWR->addCurve(0,"withdrawal rate");
    _plotWR->setData(_ageJoeVector, _withdrawalRateVector);
    _plotWR->setLineThickness(2);
    _plotWR->conclude(0,true);
    _plotWR->plotDataAndFit(true);
    FUNC_EXIT;
}
void SimWindow::plotSpending()
{
    // mortgage savings
    _plotSpending->init();
    addPlot(_plotSpending, "extra health care", _healthCare, Qt::green, true, false);
    addPlot(_plotSpending, "mortgage", _mortgageVector,   Qt::darkGray, true, false);
    addPlot(_plotSpending, "discretionary", _discretionarySpending, Qt::red, true, false);
    addPlot(_plotSpending, "required", _spendingVectorRequired, Qt::black, true, false);
//    addPlot(_plotSpending, "baseline", _spendingVectorBaseline, Qt::red, true, true);
    addPlot(_plotSpending, "total", _spendingVectorTotal, Qt::blue, true, false);
    _plotSpending->conclude(0,true);
    _plotSpending->plotDataAndFit(true);
}
void SimWindow::plotSavings()
{
    FUNC_ENTER << _savingsJoeMarketVector.size() << _savingsEmiriVector.size() << _savingsTotal.size();
    _plotSavings->init();
    addPlot(_plotSavings, "market: Joe",    _savingsJoeMarketVector,Qt::blue,      true, true);
    addPlot(_plotSavings, "market: Emiri",  _savingsEmiriVector,    Qt::darkGreen, true, true);
    addPlot(_plotSavings, "market: total",  _savingsTotal,          Qt::red,       true, false);
    addPlot(_plotSavings, "FA: Joe",        _savingsJoeFAVector,    Qt::darkCyan,  true, false);
    addPlot(_plotSavings, "buffer: Joe",    _bufferVector,          Qt::magenta,   true, false);
    addPlot(_plotSavings, "brokerage",      _savingsBrokerage,      Qt::darkYellow,true, false);
    _plotSavings->conclude(0,true);
    _plotSavings->plotDataAndFit(true);
    FUNC_EXIT;
}
void SimWindow::plotRandomSavings()
{
    _plotSavings->init();
    addPlot(_plotSavings, "bottom 10% savings",  _bottom10SavingsVector,    Qt::blue, false, false);
    addPlot(_plotSavings, "bottom quart savings",_bottomQuartSavingsVector, Qt::red,  false, false);
    addPlot(_plotSavings, "top quart savings",   _topQuartSavingsVector,    Qt::red,  false, true);
    addPlot(_plotSavings, "median savings",    _medianSavingsVector,        Qt::red,  true,  false);

    addPlot(_plotSavings, "bottom quart buffer",_bottomQuarterBufferVector, Qt::magenta,  false, true);
    addPlot(_plotSavings, "top quart buffer",   _topQuarterBufferVector,    Qt::magenta,  false,  true);

    _plotSavings->conclude(0,true);
    _plotSavings->plotDataAndFit(true);
}
void SimWindow::plotRandomSpending()
{
    bool lengendOn = _plotIncome->isLegendOn();
    _plotSpending->init();
    _plotSpending->setLegendOn(lengendOn);
    addPlot(_plotSpending, "bottom 10% spending",  _bottom10SpendingVector,    Qt::blue, false, false);
    addPlot(_plotSpending, "bottom quart spending",_bottomQuartSpendingVector, Qt::red,  false, false);
    addPlot(_plotSpending, "top quart spending",   _topQuartSpendingVector,    Qt::red,  false, true);
    addPlot(_plotSpending, "median spending",      _medianSpendingVector,      Qt::red,  true, false);
    //    plotMilestones(_plotIncome);
    _plotSpending->conclude(0,true);
    _plotSpending->plotDataAndFit(true);
}
void SimWindow::plotRandomIncome()
{
    bool lengendOn = _plotIncome->isLegendOn();
    _plotIncome->init();
    _plotIncome->setLegendOn(lengendOn);
    addPlot(_plotIncome, "bottom 10% income",  _bottom10IncomeVector,    Qt::blue, false, false);
    addPlot(_plotIncome, "bottom quart income",_bottomQuartIncomeVector, Qt::red,  false, false);
    addPlot(_plotIncome, "top quart income",   _topQuartIncomeVector,    Qt::red,  false, true);
    addPlot(_plotIncome, "median income",      _medianIncomeVector,      Qt::red,  true, false);
//    plotMilestones(_plotIncome);
    _plotIncome->conclude(0,true);
    _plotIncome->plotDataAndFit(true);
}
void SimWindow::plotRandomWR()
{
    _plotWR->init();
    addPlot(_plotWR, "bottom quart WR",_bottomQuartWRVector, Qt::red, false, false);
    addPlot(_plotWR, "median WR",      _medianWRVector,      Qt::red, true,  false);
    addPlot(_plotWR, "top quart WR",   _topQuartWRVector,    Qt::red, false, false);
    _plotWR->conclude(0,true);
    _plotWR->plotDataAndFit(true);
}

void SimWindow::plotMilestones(plotData *plot)
{
    dVector vertical = {0.,300.};
    // Emiri Medicare: 72
    // RMD Emiri: 82
    dVector milestones = {72., 82.};
    sVector legends = {"Emiri Medicare", "Emiri RMD"};
    // Lilika starts college: 73
    // Mortgage paid: 74
    for (int j=0; j<milestones.size(); j++)
    {
        plot->addCurve(0,legends.at(j));
        dVector milestone = {milestones.at(j), milestones.at(j)};
        plot->setData(milestone, vertical);
        plot->setColor(Qt::gray);
        plot->setLineThickness(4);
    }
}

void SimWindow::plotIncome()
{
    FUNC_ENTER;
    bool lengendOn = _plotIncome->isLegendOn();
    _plotIncome->init();
    _plotIncome->setLegendOn(lengendOn);

    addPlot(_plotIncome, "required spending", _spendingVectorRequired, Qt::black, true, false);
    // logicals = (thick, dashed)
    // thick: false = income; true = withdrawal
    // dashed: true = variable; false = fixed
    if ( _radioShowAllIncome->isChecked() )
    {
        addPlot(_plotIncome, "salary Joe",     _salaryJoeVector, Qt::blue, false, false);
        addPlot(_plotIncome, "salary Emiri",   _salaryEmiriVector, Qt::darkGreen, false, false);

        addPlot(_plotIncome, "SS Spouse", _ssSpousalChildInCareVector, Qt::green, true, false);
        addPlot(_plotIncome, "SS Lilika",      _ssLilikaVector, Qt::magenta, true, false);
        addPlot(_plotIncome, "SS Emiri",       _ssEmiriVector,     Qt::darkGreen, true, false);
        addPlot(_plotIncome, "SS Joe",         _ssJoeVector,     Qt::blue, true, false);

        if ( !_radioTotalIncome->isChecked() )
        {
            // convert taxes to tax savings
            dVector taxSavings;  taxSavings.resize(_taxes.size());
            for (int j=0; j<_taxes.size(); j++) taxSavings[j] = _taxes.at(0) - _taxes.at(j);
            addPlot(_plotIncome, "tax savings", taxSavings,   Qt::darkRed, true, true);

        }

        if ( useSPIA() )
            addPlot(_plotIncome, "SPIA Gap", _spiaVector, Qt::yellow, true, false);

        addPlot(_plotIncome, "Annuity Joe",  _annuityTotalVector, Qt::cyan, true, true);

        addPlot(_plotIncome, "Market Joe",   _withdrawJoeMarketVector,   Qt::blue, true, true);
        addPlot(_plotIncome, "Market Emiri", _withdrawEmiriMarketVector, Qt::green, true, true);

        dVector whichIncome = whichFinalIncome();
        addPlot(_plotIncome, "total",  whichIncome,            Qt::red,     true, false);

    }
    else if ( _radioShowFixedVariableIncome->isChecked() )
    {
        addPlot(_plotIncome, "fixed",  _totalIncomeFixed,       Qt::blue, true, false);
        addPlot(_plotIncome, "variable",  _totalIncomeVariable, Qt::blue, true, true);

        dVector whichIncome = whichFinalIncome();
        addPlot(_plotIncome, "total",  whichIncome,  Qt::red,     true, false);
    }
    else if ( _radioShowIncomeCategories->isChecked()  )
    {
        dVector salaryVector, ssVector, annuityVector, marketVector;
        for (int j=0; j<_ageJoeVector.size(); j++)
        {
            salaryVector.append(_salaryJoeVector.at(j) + _salaryEmiriVector.at(j));
            ssVector.append(_ssJoeVector.at(j) + _ssEmiriVector.at(j) + _ssLilikaVector.at(j) + _ssSpousalChildInCareVector.at(j));
            annuityVector.append(_annuityTotalVector.at(j) + _spiaVector.at(j));
            marketVector.append(_withdrawJoeMarketVector.at(j) + _withdrawEmiriMarketVector.at(j));
        }
        addPlot(_plotIncome, "salaries",  salaryVector, Qt::darkGray, true, false);
        addPlot(_plotIncome, "SS",        ssVector,     Qt::darkGreen, true, false);
        addPlot(_plotIncome, "annuities", annuityVector,Qt::cyan, true, false);
        addPlot(_plotIncome, "market",    marketVector, Qt::blue, true, true);

        dVector whichIncome = whichFinalIncome();
        addPlot(_plotIncome, "total",  whichIncome,  Qt::red,     true, false);
    }

//    plotMilestones(_plotIncome);
    _plotIncome->conclude(0,true);
    _plotIncome->plotDataAndFit(true);
    FUNC_EXIT;
}

void SimWindow::addPlot(plotData *plot, QString legend, dVector data, QColor color, bool thick, bool dashed )
{
    if ( data.size() != _ageJoeVector.size() )
    {
        qInfo() << legend << data.size() << _ageJoeVector.size();
        qFatal("Error: vectors have different dimensions in addPlot");
    }
    else if ( !allZero(data) )
    {
        plot->addCurve(0,legend);
        plot->setData(_ageJoeVector, data);
        plot->setColor(color);
//        plot->setPointSize(2);
        if ( thick )  plot->setLineThickness(2);
        if ( dashed ) plot->setLineDashed(true);
    }
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
void SimWindow::exitApp()
{
    FUNC_ENTER;
    QCoreApplication::exit(0);
}

void SimWindow::changedQLineEdit(QLineEdit *widget, int &value)
{
    bool ok;
    int newValue = widget->text().toInt(&ok);
    if ( ok )
        value = newValue;
    else
    {
        QString numberString;
        widget->setText(numberString.setNum(value));
    }
    updateGraphsAndPanels();
}
void SimWindow::changedQLineEdit(QLineEdit *widget, double &value)
{
    bool ok;
    double newValue = widget->text().toDouble(&ok);
    if ( ok )
        value = newValue;
    else
    {
        QString numberString;
        widget->setText(numberString.setNum(value));
    }
    updateGraphsAndPanels();
}

void SimWindow::changedShowPlots()
{
    _plotSpendingAndSavingsWidget->setVisible(_radioShowAll->isChecked());
    _plotWRandMarketReturnsWidget->setVisible(_radioShowAll->isChecked());
}

void SimWindow::changedTargetIncomeRadioButton()
{
    setIncomeTargetFloorVisibility();
    _marketWithDrawLabel->setVisible(!_radioTargetIncome->isChecked());
    _marketWRString->setVisible(!_radioTargetIncome->isChecked());
    updateGraphsAndPanels();
}
void SimWindow::changedMarketRandomization()
{
    bool showMarketReturnAverage = _radioRandomMarket->isChecked() || !useRandomization();
    _marketReturnsRateLabel->setVisible (showMarketReturnAverage);
    _marketReturnsRateString->setVisible(showMarketReturnAverage);
    updateGraphsAndPanels();
}
void SimWindow::changedRandomization()
{
    _radioRandomMarket->setVisible(useRandomization());
    _radioHistoricalMarket->setVisible(useRandomization());
    _lockMarketReturns->setVisible(useRandomization());
    _lastSavings->setVisible(useRandomization());
    _successRate->setVisible(useRandomization());
    _incomeFloorCheckBox->setVisible(useRandomization() && !_radioTargetIncome->isChecked());
    _percentBondsLabel->setVisible(useRandomization());
    _percentBondsString->setVisible(useRandomization());
    setIncomeTargetFloorVisibility();

    bool showMarketReturnAverage = _radioRandomMarket->isChecked() || !useRandomization();
    _marketReturnsRateLabel->setVisible (showMarketReturnAverage);
    _marketReturnsRateString->setVisible(showMarketReturnAverage);

    bool showReturns = _radioRand1->isChecked();
    _MRWidget->setVisible(showReturns);
    _MRDWidget->setVisible(!showReturns);

    _incomeShowWidget->setEnabled(!_radioRand->isChecked());

    updateGraphsAndPanels();
}

void SimWindow::setIncomeTargetFloorVisibility()
{
    FUNC_ENTER << _radioTargetIncome << _incomeFloorTargetLabel;
    bool useTarget = _radioTargetIncome->isChecked();
    _incomeFloorCheckBox->setVisible(!useTarget);

    if ( useTarget )
        _incomeFloorTargetLabel->setText("target");
    else
        _incomeFloorTargetLabel->setText("floor");
}
