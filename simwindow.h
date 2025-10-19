#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <QMainWindow>
#include <QWidget>
#include "plot.h"

class SimWindow : public QMainWindow
{
    Q_OBJECT
public:
    SimWindow();
private:
    dVector _ageJoeVector;
    dVector _ageEmiriVector;
    dVector _ageLilikaVector;

    double _ssBaseJoe62=28;
    double _ssBaseEmiri62=17;
    int _ageRMD=75;
    int _ageRMDJoe=75;  // adjustable when Joe goes on RMD schedule (to reduce for Lilika's college purposes)
    // RMD rate = 100/years, where years = slope * (age-75) + offset;
    double _RMDSlope = -0.7945;
    double _RMDOffset = 24.236;
    double _annuityVABase67=6.3;
    double _annuityFABase67=7.2;
    double _annuityWRDeltaPerYear=0.15;
    QLineEdit *_ageRMDJoeString;

    // money vectors: fixed
    dVector _salaryJoeVector, _salaryEmiriVector;
    dVector _ssJoeVector, _ssEmiriVector, _ssLilikaVector, _ssSpousalChildInCareVector, _mortgageVector, _spiaVector, _annuityFAJoeVector;
    dVector _taxes;
    dVector _healthCare;
    dVector _bufferVector;
    dVector _spendingVectorBaseline, _spendingVectorRequired, _spendingVectorTotal;
    // money vectors: variable
    dVector _annuity2JoeVector;
    dVector _withdrawJoeMarketVector,  _withdrawEmiriMarketVector;
    // money vectors: summary
    dVector _annuityTotalVector, _annuityFixedVector, _annuityVariableVector; // VA + FA and 0.6VA + FA
    // market savings vectors
    dVector _savingsJoeMarketVector, _savingsJoeFAVector, _savingsEmiriVector, _savingsTotal, _savingsBrokerage;
    dVector _withdrawalRateVector, _stockReturnsVector, _bondReturnsVector, _stockBondReturnsVector;
    dVector _bottom10SavingsVector, _bottomQuartSavingsVector, _medianSavingsVector, _topQuartSavingsVector;
    dVector _bottom10WRVector, _bottomQuartWRVector, _medianWRVector, _topQuartWRVector;
    dVector _bottom10BufferVector, _bottomQuarterBufferVector, _medianBufferVector, _topQuarterBufferVector;
    // incomes:
    dVector _totalIncomeFixed, _totalIncomeVariable, _totalIncome, _incomeAfterTax, _discretionaryIncome, _discretionarySpending;
    dVector _bottom10IncomeVector, _bottomQuartIncomeVector, _medianIncomeVector, _topQuartIncomeVector;
    dVector _bottom10SpendingVector, _bottomQuartSpendingVector, _medianSpendingVector, _topQuartSpendingVector;

    // random market return
    double _stockReturnsFWHM = 41.;
    double _bondReturnsFWHM  = 20.;
    double _bondReturnsMean  = 3.25;
    double _stockReturnsSigma = _stockReturnsFWHM / 2.355;
    double _bondReturnsSigma  = _bondReturnsFWHM / 2.355;
    dVector _marketReturnsProbability;
    double _yearlyStockReturn;      // appears to be all equities and not allow bonds
    double _yearlyBondReturn;
    double _yearlyStockBondReturn;  // could be stocks/bonds
    int _yearStartHistoricalMarket=0;
    int _yearIndexHistoricalMarket=0;
    dMatrix _historicalReturnsTable;

    // taxes
    double _standardDeduction=30.; // $30k
    double _extraDeductionMA=8.8;  // $8.8k
    double _SSTaxableFed=0.85;     // fraction taxable
    double _stateTaxPercent=5.;    // %, but not on SS
    double _FICATaxPercent=7.65;   // %, before retirement
    dMatrix _taxBrackets;          // [percent][cutoff]

    // health care
    double _healthCarePreRetire=12.;   // Emiri: 90k - 12k ira - 6k reimbursable benefits = 72k but gets only 60k, so 12k health care
    double _healthCarePrivateAdult=15.; // $k per year until age 65
    double _healthCarePrivateChild=5.;   // $k per year until age 24
    double _extraMedicare=5.;      // $3k/year/person medicare premiums + 2k out of pocket
    QLineEdit *_healthCarePreRetireString;
    QLineEdit *_healthCarePrivateAdultString;
    QLineEdit *_healthCarePrivateChildString;
    QLineEdit *_extraMedicareString;

    QStatusBar *_statusBar;
    QStatusBar *_statusBarJoePre;
    QStatusBar *_statusBarEmiriPre;
    QProgressBar *_progressBar;

    QStatusBar *_statusBarIncome;
    QStatusBar *_statusBarSpending;
    QStatusBar *_statusBarSavings;
    QStatusBar *_statusBarWR;
    QStatusBar *_statusBarMarketReturnDist;
    QStatusBar *_statusBarMarketReturns;

    int _ageStart=60;
    int _ageLimit=100;
    QLineEdit *_ageLimitString;

    // joe
    QDate _BDJoe = QDate(1962, 5, 5);
    int _currentAgeJoe = _BDJoe.daysTo(QDate::currentDate())/365;
    int _retireAgeJoe=64;
    int _ssAgeJoe=64;
    int _annuityFAStartMoney=700;
    int _marketMoneyJoe=1360;  // subtract 150 for HELOC
    double _ssJoe=41;
    int _salaryJoe62=110;
    int _salaryJoe=60;
    QLineEdit *_retireAgeJoeString;
    QLineEdit *_ssAgeJoeString;
    QLineEdit *_annuityFAStartMoneyString;
    QLineEdit *_marketMoneyJoeString;
    QLabel *_ssJoeLabel;

    // annuity: VA (although could be FA)
    int _annuity2StartAge=65;
    double _annuity2Fees=0.5;
    double _annuity2WRStarting=_annuityVABase67;
    double _annuity2HurdleRate=3.9;
    double _annuity2Percentage=0.;
    double _annuity2StartMoney=0.;
    double _annuity2YearlyIncomeStart=0.;

    // annuity: FA
    int _annuityFAStartAge=65;
    double _annuityFAFees=0.03;
    double _annuityFAWRStarting=_annuityFABase67;
    double _annuityFAAccumulationAmount=0.;
    double _annuityFAYearlyIncomeStart=0.;

    // Emiri
    QDate _BDEmiri = QDate(1969, 6, 27);
    int _currentAgeEmiri = _BDEmiri.daysTo(QDate::currentDate())/365;
    int _retireAgeEmiri=58;
    int _ssAgeEmiri=64;
    int _monthlyAdditionEmiri=1000;
    int _marketMoneyEmiri=380;
    double _ssEmiri=15.;
    double _ssSpousal;
    int _salaryEmiri=80;
    QLineEdit *_retireAgeEmiriString;
    QLineEdit *_ssAgeEmiriString;
    QLineEdit *_monthlyAdditionEmiriString;
    QLineEdit *_marketMoneyEmiriString;
    QLabel *_ssEmiriLabel;
    QLabel *_ssSpousalLabel;

    // Lilika
    QDate _BDLilika = QDate(2017, 6, 3);
    int _currentAgeLilika = _BDLilika.daysTo(QDate::currentDate())/365;
    int _collegeAgeLilika=18;
    double _ssLilika=21.;
    QLabel *_ssLilikaLabel;

    // SS crisis
    int _ssCrisisYear = 2033;
    int _ssCrisisDrop = 0;   // 23% alternative
    int _ssCrisisAgeJoe = _currentAgeJoe + _ssCrisisYear - QDate::currentDate().year();
    QLineEdit *_ssCrisisYearString;
    QLineEdit *_ssCrisisDropString;

    // cash buffer
    double _bufferDesired=0.;
    QLineEdit *_bufferString;
    double _bufferAmount=0.;
    // spending
    QLineEdit *_spendingYear0String;
    QLineEdit *_spendingPerYearIncreasePercentString;
    QLineEdit *_maxDiscretionarySpendingString;
    QLineEdit *_maxBrokerageString;
    double _spendingYear0 = 110.;  // spending including mortgage = 12 * (5 CC + 3 house + 1 extra) = 108
    double _spendingPerYearIncreasePercent = 0.;
    double _mortgage=25;
    int _mortgageAgeJoe=74; // 2036
    double _maxDiscretionarySpending=200.;
    double _maxBrokerage=200.;

    // rates & fees
    double _marketReturnsRate=8.;
    double _marketFees=0.1;
    double _marketWR=6.;
    double _marketIncomePerYearAv;
    int _percentBonds = 0.;
    QLineEdit *_marketReturnsRateString;
    QLineEdit *_marketFeesString;
    QLineEdit *_marketWRString;
    QLabel *_percentBondsLabel;
    QLineEdit *_percentBondsString;
    QLabel *_marketIncomePerYearAvLabel;
    QLabel *_marketWithDrawLabel;
    QLabel *_marketReturnsRateLabel;

    double _reducedBenefitPercent=67.;
    double _annuityIncomePerYearAv;
    int _firstDemiseAge=93;

    // FA (annuity #1)
    QLineEdit *_annuityFAStartAgeString;
    QLineEdit *_annuityFAFeesString;
    QLineEdit *_annuityFAWRStartingString;
    QLabel *_annuityFAIncomePerYearAvLabel;

    // annuity #2
    QLineEdit *_annuity2StartAgeString;
    QLineEdit *_annuity2FeesString;
    QLineEdit *_annuity2WRStartingString;
    QLineEdit *_annuity2HurdleRateString;
    QLineEdit *_annuity2PercentageString;
    QLabel *_annuityIncomePerYearAvLabel;
    QCheckBox *_annuity2IsFACheckBox;

    // either VA or FA
    QLineEdit *_reducedBenefitPercentString;
    QLineEdit *_firstDemiseAgeString;

    // SPIA
    double _spiaAmount=300;
    int _spiaAge=65;
    double _spiaWR=13.1;
    double _inflation=3.;
    QRadioButton *_radioSpia10;
    QRadioButton *_radioSpia5;
    QLineEdit *_spiaAmountString;
    QLineEdit *_spiaAgeString;
    QLineEdit *_inflationString;

    int _editTextWidth=60;
    QWidget *_incomeTypeWidget;
    QRadioButton *_radioShowAll;
    QRadioButton *_radioShowIncome;

    QWidget *_incomeShowWidget;
    QRadioButton *_radioTotalIncome;
    QRadioButton *_radioAfterTaxHealthMortgage;
    QRadioButton *_radioDiscretionaryIncome;
    QRadioButton *_radioShowAllIncome;
    QRadioButton *_radioShowFixedVariableIncome;
    QRadioButton *_radioShowIncomeCategories;

    QRadioButton *_radioTargetIncome;
    QRadioButton *_radioGuardRails;
    QRadioButton *_radioRand;
    QRadioButton *_radioRand1;
    QRadioButton *_radioHistoricalMarket;
    QRadioButton *_radioRandomMarket;
    QCheckBox *_lockMarketReturns;
    QLabel *_incomeFloorTargetLabel;
    int _incomeFloorTarget=150;
    QLineEdit *_incomeFloorTargetString;
    QCheckBox *_incomeFloorCheckBox;
    QLabel *_lastSavings;
    QLabel *_successRate;

    QTabWidget *_tabTimeSpace;
    QWidget *_planPage;
    QWidget *_resultsPage;

    // outcomes
    QLabel *_moneyEndWorkJoeLabel;
    QLabel *_moneyEndWorkEmiriLabel;
    QLabel *_moneyEndWorkTotalLabel;

    // outcomes
    QLabel *_marketStartRetireJoeLabel;
    QLabel *_annuitize2AmountLabel;
    QLabel *_annuitizeFAAmountLabel;
    QLabel *_marketRemainderStartJoeLabel;
    QLabel *_incomeAverage65_74;
    QLabel *_incomeAverage75_84;
    QLabel *_incomeAverage85_94;
    QGroupBox *_averageIncomeBox;

    // plots

    plotData *_plotMarketReturnDist;
    plotData *_plotMarketReturns;
    plotData *_plotWR;
    plotData *_plotIncome;
    plotData *_plotSpending;
    plotData *_plotSavings;
    QWidget *_plotSpendingAndSavingsWidget;
    QWidget *_plotWRandMarketReturnsWidget;
    QWidget *_MRDWidget;
    QWidget *_MRWidget;

    void changedQLineEdit(QLineEdit *widget, int &value);
    void changedQLineEdit(QLineEdit *widget, double &value);
    double getRandomStockReturn();
    double getRandomBondReturn();
    inline void randomizeHistoricalStartingPoint() {_yearStartHistoricalMarket = _yearIndexHistoricalMarket = _historicalReturnsTable.size() * RAND_0_1;}
    dPoint2D getHistoricalMarketReturn();
    double GaussianRandomizer(double mean, double sigma, double cutoff);

    void createResultsPage();
    QGroupBox *createJoeGroupBox();
    QGroupBox *createEmiriGroupBox();
    QGroupBox *createRatesSSGroupBox();
    QGroupBox *createHealthCareGroupBox();
    QGroupBox *createRetirementStartingGroupBox();
    QGroupBox *createSpendingGroupBox();
    QGroupBox *createApproachJoeGroupBox();
    QGroupBox *createMarketRatesAndReturnsGroupBox();
    QGroupBox *createAnnuityGroupBox();
    QGroupBox *createAvIncomeGroupBox();

    void createPlanPage();

    void clearAllVectors();
    void updateYearlyMarketReturns(int iAge);
    void runRetirement();
    void randomRetirementAndSummarize();
    double payTaxes(int whichAge, double totalIncome);
    double getFederalTax(double income);
    double payHealthCare(int whichAge);

    void plotMarketReturns();
    void createTimeVectors();
    void createFixedMoneyVectors();
    void createSummaryInfo();
    void plotResults();
    void plotWithdrawalRate();
    void plotIncome();
    void plotSpending();
    void plotSavings();
    void plotRandomSavings();
    void plotRandomSpending();
    void plotRandomIncome();
    void plotRandomWR();
    void createRandomSummaryInfo();
    void extractQuartiles(dMatrix &inputVector, dVector &bottom10, dVector &bottomQuart, dVector &medianVector, dVector &topQuart);
    double findTimeAverageFromVector(dVector ageVector, dVector moneyVector, double age1, double age2);

    void createSSVectors();
    dVector createStaticMoneyVector(dVector ageVector, int ageStart, int ageEnd, double amount);
    dVector createStaticMoneyVector(dVector ageVector, int ageStart, int ageEnd, double amount, double inflation);
    dVector createStaticMoneyVector(dVector ageVector, int ageStart, int ageEnd, int ageEligible, double ageShift, double amount, double newAmount,
                                    double inflationPercent);
    dVector createStaticMoneyVector(dVector ageVector, int ageStart, int ageEnd, int ageEligible,
                                    double ageShift, double amount, double newAmount);
    dVector addMoneyVectors(dVector moneyVector1, dVector moneyVector2);

    dVector createTimeVector(int startAge, int endAge);
    dVector createMoneyVectorPre(dVector ageVector, double startMoney, double monthlyAddition, int endAddtion,
                                 double rate, double fees);
    double withdrawMarketIncomeWithGuardRails(double age, double ageRMD, double nonMarketIncome, int ageIncomeStarts, double targetWRPreRMD, double &guardRailDollars, double &savings);
    double withdrawMarketIncome(double age, double ageRMD, double salaryLikeIncome, int ageIncomeStarts, double withdrawalRatePreRMD, double &qualifiedSavings);
    double withdrawMarketIncomeFloor(double age, double ageRMD, double salaryLikeIncome, double incomeTarget,
                                     int ageIncomeStarts, double &qualifiedSavings);
    void addPlot(plotData *plot, QString legend, dVector data, QColor color, bool thick, bool dashed );
    void plotMilestones(plotData *plot);
    bool allZero(dVector data);

    void updateSS();
    void updateAnnuities();
    void updateGraphs();
    void setIncomeTargetFloorVisibility();

    inline bool useRandomization() {return _radioRand1->isChecked() || _radioRand->isChecked();}
    inline bool useSPIA() {return _radioSpia5->isChecked() || _radioSpia10->isChecked();}

    inline double getRequiredTotalIncome()
    {
        double requiredTotalIncome=0.;  // default
        if ( _incomeFloorCheckBox->isChecked() ) requiredTotalIncome = _incomeFloorTarget;
        return requiredTotalIncome;
    }

    inline double getSpendingPerYear(double age)
    {
        double ratioPerYear = 1. + _spendingPerYearIncreasePercent/100.;
        double fractionalIncrease = pow(ratioPerYear,age-_ageStart);
        return _spendingYear0 * fractionalIncrease;
    }

    inline void updateAnnuityWithdrawalRates()
    {
        double delta = (_annuity2StartAge - 67) * _annuityWRDeltaPerYear;
        if ( _annuity2IsFACheckBox->isChecked() )
        {
            _annuity2WRStarting = _annuityFABase67+delta;
            _annuity2Fees = 0.03;
        }
        else
        {
            _annuity2WRStarting = _annuityVABase67+delta;
            _annuity2Fees = 0.5;
        }
        QString number;
        _annuity2WRStartingString->setText(number.setNum(_annuity2WRStarting,'g',3));
        _annuity2FeesString->setText(number.setNum(_annuity2Fees,'g',3));

        delta = (_annuityFAStartAge - 67) * _annuityWRDeltaPerYear;
        _annuityFAWRStarting = _annuityFABase67+delta;
        _annuityFAWRStartingString->setText(number.setNum(_annuityFAWRStarting,'g',3));
    }
    inline double mixedMarketReturn(double stockReturn, double bondReturn)
    {
        double marketReturn = (1.-_percentBonds/100.) * stockReturn +  _percentBonds/100. * bondReturn;
        return marketReturn;
    }

    // getters
    inline dVector whichFinalIncome()
    {
        dVector whichIncome;
        if ( _radioTotalIncome->isChecked() )
            whichIncome = _totalIncome;
        else if ( _radioAfterTaxHealthMortgage->isChecked() )
            whichIncome = _incomeAfterTax;
        else
            whichIncome = _discretionaryIncome;
        return whichIncome;
    }

private slots:
    void exitApp();
    void aboutApp();
    inline void updateGraphsAndPanels() {updatePanels(); updateGraphs();}
    void changedShowPlots();

    inline void updatePanels()
    {
        updateAnnuityWithdrawalRates();
        updateSS();
        updateAnnuities();
    }

    inline void changedRetireAgeJoe()
    {
        changedQLineEdit(_retireAgeJoeString, _retireAgeJoe);
        QString numberString;
        _ssAgeJoeString->setText(numberString.setNum(qMax(_ssAgeJoe,62)));
        changedQLineEdit(_ssAgeJoeString, _ssAgeJoe);
        _spiaAgeString->setText(numberString.setNum(qMax(62,_retireAgeJoe)));
        changedQLineEdit(_spiaAgeString, _spiaAge);
        updateAnnuityWithdrawalRates();
        updateGraphsAndPanels();
    }
    inline void changedSSAgeJoe()         {changedQLineEdit(_ssAgeJoeString, _ssAgeJoe);}
    inline void changedCurrentMoneyJoe()  {changedQLineEdit(_marketMoneyJoeString, _marketMoneyJoe);}
    inline void changedFAStartingAmount() {changedQLineEdit(_annuityFAStartMoneyString, _annuityFAStartMoney);}
    inline void changedRMDAgeJoe()        {changedQLineEdit(_ageRMDJoeString, _ageRMDJoe);}

    inline void changedRetireAgeEmiri()
    {
        changedQLineEdit(_retireAgeEmiriString, _retireAgeEmiri);
        QString numberString;
        _ssAgeEmiriString->setText(numberString.setNum(qMax(64,_retireAgeEmiri)));
        changedQLineEdit(_ssAgeEmiriString, _ssAgeEmiri);
    }
    inline void changedSSAgeEmiri()             {changedQLineEdit(_ssAgeEmiriString, _ssAgeEmiri);}
    inline void changedCurrentMoneyEmiri()      {changedQLineEdit(_marketMoneyEmiriString, _marketMoneyEmiri);}
    inline void changedMonthlyAdditionEmiri()   {changedQLineEdit(_monthlyAdditionEmiriString, _monthlyAdditionEmiri);}

    inline void changedAgeLimit()               {changedQLineEdit(_ageLimitString, _ageLimit);}
    inline void changedSpiaAmount()             {changedQLineEdit(_spiaAmountString, _spiaAmount);}
    inline void changedSpiaAge()                {changedQLineEdit(_spiaAgeString, _spiaAge);}
    inline void changedAnnuity2Age()
    {
        changedQLineEdit(_annuity2StartAgeString, _annuity2StartAge);
        _retireAgeJoe = qMin(_retireAgeJoe,_annuity2StartAge);
        QString numberString;
        _retireAgeJoeString->setText(numberString.setNum(_retireAgeJoe));
        updateAnnuityWithdrawalRates();
        updateGraphsAndPanels();
    }
    inline void changedFAAnnuityAge()
    {
        changedQLineEdit(_annuityFAStartAgeString, _annuityFAStartAge);
        _retireAgeJoe = qMin(_retireAgeJoe,_annuityFAStartAge);
        QString numberString;
        _retireAgeJoeString->setText(numberString.setNum(_retireAgeJoe));
        updateAnnuityWithdrawalRates();
        updateGraphsAndPanels();
    }
    inline void changedInflation()       {changedQLineEdit(_inflationString, _inflation);}

    inline void changedSSCrisisYear()    {changedQLineEdit(_ssCrisisYearString, _ssCrisisYear);}
    inline void changedSSCrisisPercent() {changedQLineEdit(_ssCrisisDropString, _ssCrisisDrop);}

    inline void changedHCPreRetire()     {changedQLineEdit(_healthCarePreRetireString, _healthCarePreRetire);}
    inline void changedHCPrivateAdult()  {changedQLineEdit(_healthCarePrivateAdultString, _healthCarePrivateAdult);}
    inline void changedHCPrivateChild()  {changedQLineEdit(_healthCarePrivateChildString, _healthCarePrivateChild);}
    inline void changedHCExtraMedicare() {changedQLineEdit(_extraMedicareString, _extraMedicare);}

    inline void changedTargetIncome()    {FUNC_ENTER; changedQLineEdit(_incomeFloorTargetString, _incomeFloorTarget);}
    inline void changedCashBuffer()      {changedQLineEdit(_bufferString, _bufferDesired);}

    inline void changedSpendingYear0() {changedQLineEdit(_spendingYear0String, _spendingYear0);}
    inline void changedSpendingSlope() {changedQLineEdit(_spendingPerYearIncreasePercentString, _spendingPerYearIncreasePercent);}
    inline void changedMaxSpending() {changedQLineEdit(_maxDiscretionarySpendingString, _maxDiscretionarySpending);}
    inline void changedMaxBrokerage() {changedQLineEdit(_maxBrokerageString, _maxBrokerage);}
    inline void changedMarketRate()    {changedQLineEdit(_marketReturnsRateString, _marketReturnsRate);}
    inline void changedMarketFees()    {changedQLineEdit(_marketFeesString, _marketFees);}
    inline void changedMarketWR()              {changedQLineEdit(_marketWRString, _marketWR);}
    inline void changedBondPercentage()        {changedQLineEdit(_percentBondsString, _percentBonds);}
    inline void changedAnnuity2Fees()          {changedQLineEdit(_annuity2FeesString, _annuity2Fees);}
    inline void changedAnnuityFAFees()         {changedQLineEdit(_annuityFAFeesString, _annuityFAFees);}
    inline void changedAnnuity2WRStarting()    {changedQLineEdit(_annuity2WRStartingString, _annuity2WRStarting);}
    inline void changedAnnuityFAWRStarting()   {changedQLineEdit(_annuityFAWRStartingString, _annuityFAWRStarting);}
    inline void changedAnnuityHurdleRate()     {changedQLineEdit(_annuity2HurdleRateString, _annuity2HurdleRate);}
    inline void changedAnnuitized2Percent()    {changedQLineEdit(_annuity2PercentageString, _annuity2Percentage);}
    inline void changedFirstDemiseAge()        {changedQLineEdit(_firstDemiseAgeString, _firstDemiseAge);}
    inline void changedReducedBenefitPercent() {changedQLineEdit(_reducedBenefitPercentString, _reducedBenefitPercent);}

    inline void setSpia5()
    {
        _spiaWR = 22.3;
        updateGraphsAndPanels();
    }
    inline void setSpia10()
    {
        _spiaWR = 13.1;
        updateGraphsAndPanels();
    }

    void changedTargetIncomeRadioButton();
    void changedMarketRandomization();
    void changedRandomization();
public slots:
};
#endif // SIMWINDOW_H

