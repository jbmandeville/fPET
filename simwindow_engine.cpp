#include <QtWidgets>
#include <QDebug>
#include <QVector>
#include <QObject>
#include <QWidget>
#include <QFileDialog>
#include <QString>
#include <QDateTime>
#include "simwindow.h"

void SimWindow::runRetirement()
{
    FUNC_ENTER;

    QString numberString;

    if ( !_lockMarketReturns->isChecked() ) randomizeHistoricalStartingPoint();
    createTimeVectors();
    // Create income vector that have no dependence upon current savings (salary, social security, mortgage)
    createFixedMoneyVectors();

    double marketMoneyJoe   = _marketMoneyJoe;
    double marketMoneyEmiri = _marketMoneyEmiri;
    double yearlyFAIncome   = 0.;
    double yearlyAnnuity2Income = 0.;
    double healthCareYear0  = 0.;
    double savingsBrokerage = 0.;

    // start by moving some money to a buffer
    _bufferAmount   = _bufferDesired;
    marketMoneyJoe -= _bufferAmount;

    // define guard rails, which depend upon an initial withdrawal rate converted to dollars
    double guardRailDollarsJoe   = _marketWR/100. * marketMoneyJoe;
    double guardRailDollarsEmiri = _marketWR/100. * marketMoneyEmiri;

    _annuity2StartMoney = 0.;
    _annuityFAAccumulationAmount = _annuityFAStartMoney;

    clearAllVectors();

    for (int jAge=0; jAge<_ageJoeVector.size(); jAge++)
    {
        // year 0 is based on pre-retirement
        double ageJoe   = _ageJoeVector.at(jAge);
        double ageEmiri = _ageEmiriVector.at(jAge);
        double yearlyMarketIncomeJoe   = 0.;
        double yearlyMarketIncomeEmiri = 0.;

        // next year: market return etc
        // calculate the market return over the coming year
        updateYearlyMarketReturns(jAge);

        if ( ageJoe == _retireAgeJoe )
        {
            int iMarketJoe = marketMoneyJoe + _bufferAmount;  int iAnnuity = _annuityFAAccumulationAmount;
            int totalJoe = iMarketJoe + iAnnuity;
            int percentMarket = iMarketJoe * 100 / totalJoe;
            _marketStartRetireJoeLabel->setText(QString("%1 = %2 %").arg(numberString.setNum(marketMoneyJoe+_bufferAmount)).arg(numberString.setNum(percentMarket)));
            _moneyEndWorkJoeLabel->setText(numberString.setNum(totalJoe) + " , " + numberString.setNum(iMarketJoe) + "/" + numberString.setNum(iAnnuity));
        }
        if ( ageEmiri == _retireAgeEmiri )
        {
            int iMarketEmiri = marketMoneyEmiri;
            _moneyEndWorkEmiriLabel->setText(numberString.setNum(iMarketEmiri));
        }

        // annuity #1 is always use FA, but age is variable
        if ( ageJoe == _annuityFAStartAge )
        {
            // fixed annuity
            yearlyFAIncome = _annuityFAAccumulationAmount * _annuityFAWRStarting / 100.;
            _annuityFAYearlyIncomeStart = yearlyFAIncome;
            if ( _radioRand->isChecked() )
                _annuitizeFAAmountLabel->setText("NA");
            else
            {
                QString numberString;
                _annuitizeFAAmountLabel->setText(numberString.setNum(static_cast<int>(_annuityFAAccumulationAmount)));
            }
            _annuityFAAccumulationAmount = 0.;
        }

        // annuity #2
        if ( ageJoe == _annuity2StartAge )
        {
            _annuity2StartMoney = _annuity2Percentage/100. * (marketMoneyJoe+_bufferAmount);  // assume buffer is in qualified account
            marketMoneyJoe -= _annuity2StartMoney;
            yearlyAnnuity2Income = _annuity2StartMoney * _annuity2WRStarting / 100.;
            _annuity2YearlyIncomeStart = yearlyAnnuity2Income;
            if ( _radioRand->isChecked() )
                _annuitize2AmountLabel->setText("NA");
            else
            {
                QString numberString;
                _annuitize2AmountLabel->setText(numberString.setNum(static_cast<int>(_annuity2StartMoney)));
                _marketRemainderStartJoeLabel->setText(numberString.setNum(marketMoneyJoe+_bufferAmount));
            }
        }

        if ( useSPIA() && ageJoe == _spiaAge ) marketMoneyJoe -= _spiaAmount;

        // Joe's annuity income
        if ( ageJoe == _firstDemiseAge )
        {
            yearlyFAIncome *= _reducedBenefitPercent/100.;
            yearlyAnnuity2Income *= _reducedBenefitPercent/100.;
        }

        double totalSalary = _salaryJoeVector[jAge] + _salaryEmiriVector[jAge];
        double totalSS     = _ssJoeVector[jAge]     + _ssEmiriVector[jAge] + _ssLilikaVector[jAge] + _ssSpousalChildInCareVector[jAge];
        double totalAnnuity = yearlyAnnuity2Income + yearlyFAIncome;
        double salaryLikeIncome = totalSalary + totalSS + totalAnnuity + _spiaVector[jAge];
        if ( _radioTargetIncome->isChecked() )
        { // income floor
            yearlyMarketIncomeJoe = withdrawMarketIncomeFloor(ageJoe, _ageRMDJoe, salaryLikeIncome, _incomeFloorTarget,
                                                              _retireAgeJoe, marketMoneyJoe);
            yearlyMarketIncomeEmiri = withdrawMarketIncomeFloor(ageEmiri, _ageRMD, salaryLikeIncome + yearlyMarketIncomeJoe,
                                                                _incomeFloorTarget, _retireAgeEmiri, marketMoneyEmiri);
        } // income floor
        else if ( _radioGuardRails->isChecked() )
        { // follow a guardrails approach
            // Joe
            yearlyMarketIncomeJoe = withdrawMarketIncomeWithGuardRails(ageJoe, _ageRMDJoe, salaryLikeIncome, _retireAgeJoe, _marketWR,
                                                                       guardRailDollarsJoe ,marketMoneyJoe);
            // Emiri
            yearlyMarketIncomeEmiri = withdrawMarketIncomeWithGuardRails(ageEmiri, _ageRMD, salaryLikeIncome + yearlyMarketIncomeJoe,
                                                                         _retireAgeEmiri, 0, guardRailDollarsEmiri ,marketMoneyEmiri);
        }
        else
        { // follow a percentage withdrawal rule
            // Joe
            double marketWRJoePreRMD   = _marketWR;
            yearlyMarketIncomeJoe   = withdrawMarketIncome(ageJoe, _ageRMDJoe, salaryLikeIncome, _retireAgeJoe, marketWRJoePreRMD, marketMoneyJoe);

            double marketWREmiriPreRMD = 0.;
            yearlyMarketIncomeEmiri = withdrawMarketIncome(ageEmiri, _ageRMD, salaryLikeIncome + yearlyMarketIncomeJoe,
                                                           _retireAgeEmiri, marketWREmiriPreRMD, marketMoneyEmiri);
        }

        double fixedAnnuity=0.;  double variableAnnuity=0.;
        double fixedFactor = 0.7;
        double varFactor = 1.-fixedFactor;
        fixedAnnuity     = fixedFactor * yearlyAnnuity2Income + yearlyFAIncome;
        variableAnnuity  = varFactor * yearlyAnnuity2Income;

        double totalIncomeFixed = totalSalary + totalSS + fixedAnnuity + _spiaVector[jAge];
        double totalIncomeVariable = yearlyMarketIncomeJoe + yearlyMarketIncomeEmiri + variableAnnuity;

        double totalIncome  = totalIncomeFixed + totalIncomeVariable;
        double totalSavings = marketMoneyJoe + marketMoneyEmiri + _annuityFAAccumulationAmount;

        // Assume everything is ok (enough money) and pay taxes
        double taxes = payTaxes(jAge, totalIncome);  // uses _totalIncome, _ssVectors
        double health = payHealthCare(jAge);
        if ( jAge == 0 ) healthCareYear0 = health;
        health -= healthCareYear0;
        // convert health care to "extra health care" relative to pre-retirement
        double incomeAfterTax = totalIncome - taxes;
        double mortageReduction = _mortgage - _mortgageVector[jAge];

        double spendingRequired = getSpendingPerYear(ageJoe) + health - mortageReduction;
        double additionalAmountNeeded = spendingRequired-incomeAfterTax;
        double taxEfficiency = 0.75;
        if ( additionalAmountNeeded > 0. && _bufferAmount > 0. )
        {
            // withdraw from buffer
            double withdrawPreTax = additionalAmountNeeded/taxEfficiency;
            double incomeAdditional = taxEfficiency * qMin(withdrawPreTax,_bufferAmount);
            withdrawPreTax = incomeAdditional/taxEfficiency;  // could be altered by low buffer amount
            _bufferAmount -= withdrawPreTax;  // this much lost inside IRA
            additionalAmountNeeded -= incomeAdditional;  // this much received outside IRA
            yearlyMarketIncomeJoe += incomeAdditional;
            incomeAfterTax += incomeAdditional;
        }
        if ( additionalAmountNeeded > 0. && marketMoneyJoe > 0. )
        {
            // withdraw from qualified market funds
            double withdrawPreTax = additionalAmountNeeded/taxEfficiency;
            double incomeAdditional = taxEfficiency * qMin(withdrawPreTax,marketMoneyJoe);
            withdrawPreTax = incomeAdditional/taxEfficiency;  // could be altered by low market amount
            marketMoneyJoe -= withdrawPreTax;  // this much lost inside IRA
            additionalAmountNeeded -= incomeAdditional;  // this much received outside IRA
            yearlyMarketIncomeJoe += incomeAdditional;
            incomeAfterTax += incomeAdditional;
        }
        if ( additionalAmountNeeded > 0. && marketMoneyEmiri > 0. )
        {
            // withdraw from qualified market funds
            double withdrawPreTax = additionalAmountNeeded/taxEfficiency;
            double incomeAdditional = taxEfficiency * qMin(withdrawPreTax,marketMoneyEmiri);
            withdrawPreTax = incomeAdditional/taxEfficiency;  // could be altered by low market amount
            marketMoneyEmiri -= withdrawPreTax;  // this much lost inside IRA
            additionalAmountNeeded -= incomeAdditional;  // this much received outside IRA
            incomeAfterTax += incomeAdditional;
            yearlyMarketIncomeEmiri += incomeAdditional;
        }
        // By taking additional market income above, we may have changed totalIncomeVariable, marketMoneyJoe, marketMoneyEmiri
        totalIncomeVariable = yearlyMarketIncomeJoe + yearlyMarketIncomeEmiri + variableAnnuity;
        totalIncome  = totalIncomeFixed + totalIncomeVariable;

        double discretionaryIncomeAvailable = incomeAfterTax - spendingRequired;
        double discretionarySpending = discretionaryIncomeAvailable;
        double brokerageIncome=0.;
        if ( discretionaryIncomeAvailable < 0. )
        {
            if ( savingsBrokerage > 0. )
            {
                double moneyNeeded = qMax(-discretionaryIncomeAvailable,savingsBrokerage*_marketWR/100.);
                brokerageIncome = qMin(savingsBrokerage,moneyNeeded/0.85); // assume 15% Long-term capital gain
                savingsBrokerage -= brokerageIncome;
                discretionarySpending += brokerageIncome;
            }
        }
        else
        { // money may be available to deposit into brokerage
            discretionarySpending = qMin(discretionaryIncomeAvailable,_maxDiscretionarySpending);
            double depositIntoBrokerage = qMax(0.,discretionaryIncomeAvailable - discretionarySpending);
            if ( savingsBrokerage < _maxBrokerage )
                savingsBrokerage += depositIntoBrokerage;  // deposit into brokerage
            else
                discretionarySpending += depositIntoBrokerage;  // spend it on something
        }
        totalSavings = marketMoneyJoe + marketMoneyEmiri + _annuityFAAccumulationAmount + savingsBrokerage;
        double totalWR = 0.;
        double totalMarketIncome = yearlyMarketIncomeJoe + yearlyMarketIncomeEmiri + brokerageIncome;
        if ( totalSavings > 0. )
            totalWR = 100. * totalMarketIncome / totalSavings;

        _taxes.append(taxes);
        _healthCare.append(health);

        _incomeAfterTax.append(incomeAfterTax);
        _totalIncomeFixed.append(totalIncomeFixed);
        _totalIncomeVariable.append(totalIncomeVariable+brokerageIncome);
        _totalIncome.append(totalIncome+brokerageIncome);
        _discretionaryIncome.append(discretionaryIncomeAvailable);
        _discretionarySpending.append(discretionarySpending);

        _spendingVectorBaseline.append(getSpendingPerYear(ageJoe));
        _spendingVectorRequired.append(spendingRequired);
        _spendingVectorTotal.append(spendingRequired+discretionarySpending);

        _withdrawJoeMarketVector.append(yearlyMarketIncomeJoe);
        _withdrawEmiriMarketVector.append(yearlyMarketIncomeEmiri);
        _withdrawalRateVector.append(totalWR);

        _annuityFAJoeVector.append(yearlyFAIncome);
        _annuity2JoeVector.append(yearlyAnnuity2Income);
        _annuityTotalVector.append(totalAnnuity);
        _annuityFixedVector.append(fixedAnnuity);
        _annuityVariableVector.append(variableAnnuity);

        _savingsJoeMarketVector.append(marketMoneyJoe);
        _savingsEmiriVector.append(marketMoneyEmiri);
        _savingsTotal.append(totalSavings);
        _bufferVector.append(_bufferAmount);
        _savingsBrokerage.append(savingsBrokerage);

        if ( ageJoe >= _currentAgeJoe )
        {
            // update for next year
            marketMoneyJoe   *= 1. + (_yearlyStockBondReturn-_marketFees)/100.;
            savingsBrokerage *= 1. + (_yearlyStockBondReturn-_marketFees)/100.;
            if ( ageJoe < _annuityFAStartAge ) _annuityFAAccumulationAmount *= 1.02;  // assume 2% over inflation
            marketMoneyEmiri *= 1. + (_yearlyStockBondReturn-_marketFees)/100.;
            if ( !_annuity2IsFACheckBox->isChecked() )
            {
                double hurdleRatio = (1.+_yearlyStockReturn/100.-_annuity2Fees/100.) / (1.+_annuity2HurdleRate/100.);
                yearlyAnnuity2Income *= hurdleRatio;          // assume VA and all market rates are already corrected for inflation
            }
            else
                yearlyAnnuity2Income *= 1.-_inflation/100. - _annuity2Fees/100.;   // correct FA down for inflation and fees
            yearlyFAIncome *= 1.-_inflation/100. - _annuityFAFees/100.;   // correct FA down for inflation and fees
        }
        _savingsJoeFAVector.append(_annuityFAAccumulationAmount);

        double averageLast2Years = _stockBondReturnsVector.at(jAge);
        if ( jAge > 0 ) averageLast2Years = (_stockBondReturnsVector.at(jAge-1) + _stockBondReturnsVector.at(jAge))/2.;
        // replenish buffer only for good market and low non-IRA savings
        if ( averageLast2Years > 0. && ageJoe < _ageRMD )
        //        if ( _yearlyStockBondReturn > 0. )// && ageJoe < _ageRMD )
        {
            double fillBufferAmount = qMin(_bufferDesired - _bufferAmount, marketMoneyJoe);
            _bufferAmount  += fillBufferAmount;
            marketMoneyJoe -= fillBufferAmount;

            fillBufferAmount = qMin(_bufferDesired - _bufferAmount, marketMoneyEmiri);
            _bufferAmount  += fillBufferAmount;
            marketMoneyEmiri -= fillBufferAmount;
        }
    } // year

    FUNC_EXIT;
}

void SimWindow::randomRetirementAndSummarize()
{
    int nTrials = 10000;
    dMatrix savingsTotalPerTrial;  savingsTotalPerTrial.resize(nTrials);
    dMatrix incomeTotalPerTrial;   incomeTotalPerTrial.resize(nTrials);
    dMatrix WRTotalPerTrial;       WRTotalPerTrial.resize(nTrials);
    dMatrix bufferPerTrial;        bufferPerTrial.resize(nTrials);
    dMatrix discretionary;         discretionary.resize(nTrials);

    double marketIncomePerYearAv=0.;
    double annuityIncomePerYearAv=0.;

    int nFailures = 0;
    int iRange = 60;

    _marketReturnsProbability.fill(0.,2*iRange+1);
    for (int jTrial=0; jTrial<nTrials; jTrial++)
    {
        runRetirement();
        createSummaryInfo();
        incomeTotalPerTrial[jTrial]  = whichFinalIncome();
        savingsTotalPerTrial[jTrial] = _savingsTotal;
        discretionary[jTrial] = _discretionarySpending;
        WRTotalPerTrial[jTrial] = _withdrawalRateVector;
        bufferPerTrial[jTrial] = _bufferVector;
        marketIncomePerYearAv += _marketIncomePerYearAv;
        annuityIncomePerYearAv+= _annuityIncomePerYearAv;
        // What is a failure? Could choose running out of savings, but instead choose salary floor
        bool tooLowIncome = false;
        for (int jAge=0; jAge<_ageJoeVector.size(); jAge++)
            if ( _ageJoeVector.at(jAge) >= _retireAgeJoe ) tooLowIncome |= _discretionarySpending[jAge] < 0.;
        if ( tooLowIncome ) nFailures++;
        //            if ( _savingsTotal.last() <= 0 )  // should never happen for % withdrawal
        //                nFailures++;
    } // jTrial
    double failureProb = static_cast<double>(nFailures) / static_cast<double>(nTrials) * 100.;
    double successProb = 100. - failureProb;
    QString numberString;
    numberString.setNum(successProb,'g',3);
    _successRate->setText(QString("success rate = %1 %").arg(numberString));

    marketIncomePerYearAv  /= static_cast<double>(nTrials);
    annuityIncomePerYearAv /= static_cast<double>(nTrials);
    annuityIncomePerYearAv = qRound(annuityIncomePerYearAv);
    _marketIncomePerYearAvLabel->setText(numberString.setNum(marketIncomePerYearAv,'g',3));
    _annuityIncomePerYearAvLabel->setText(numberString.setNum(annuityIncomePerYearAv,'g',3));

    extractQuartiles(savingsTotalPerTrial,   _bottom10SavingsVector, _bottomQuartSavingsVector,  _medianSavingsVector,   _topQuartSavingsVector);
    extractQuartiles(incomeTotalPerTrial,    _bottom10IncomeVector,  _bottomQuartIncomeVector,   _medianIncomeVector,    _topQuartIncomeVector);
    extractQuartiles(discretionary,          _bottom10SpendingVector,_bottomQuartSpendingVector, _medianSpendingVector,  _topQuartSpendingVector);
    extractQuartiles(WRTotalPerTrial,        _bottom10WRVector,      _bottomQuartWRVector,       _medianWRVector,        _topQuartWRVector);
    extractQuartiles(bufferPerTrial,         _bottom10BufferVector,  _bottomQuarterBufferVector, _medianBufferVector,    _topQuarterBufferVector);

    dMatrix bufferVsTime;          bufferVsTime.resize(_ageJoeVector.size());
    dMatrix discretionaryVsTime;   discretionaryVsTime.resize(_ageJoeVector.size());

    numberString = QString("last savings = %1 / %2 / %3 / %4 ")
                       .arg(static_cast<int>(_bottom10SavingsVector.last()))
                       .arg(static_cast<int>(_bottomQuartSavingsVector.last()))
                       .arg(static_cast<int>(_medianSavingsVector.last()))
                       .arg(static_cast<int>(_topQuartSavingsVector.last()));
    _lastSavings->setText(numberString);

    double sum = 0.;
    for (int j=0; j<_marketReturnsProbability.size(); j++)
        sum += _marketReturnsProbability.at(j);
    for (int j=0; j<_marketReturnsProbability.size(); j++)
        _marketReturnsProbability[j] /= sum;
    plotRandomSpending();
    plotRandomSavings();
    plotRandomIncome();
    plotRandomWR();
    createRandomSummaryInfo();
}

void SimWindow::clearAllVectors()
{
    _withdrawJoeMarketVector.clear();
    _withdrawEmiriMarketVector.clear();

    _annuityFAJoeVector.clear();
    _annuity2JoeVector.clear();

    _savingsJoeMarketVector.clear();
    _savingsJoeFAVector.clear();
    _savingsEmiriVector.clear();
    _savingsTotal.clear();
    _savingsBrokerage.clear();

    _withdrawalRateVector.clear();

    _bufferVector.clear();

    _spendingVectorBaseline.clear();
    _spendingVectorRequired.clear();
    _spendingVectorTotal.clear();

    _totalIncomeFixed.clear();
    _totalIncomeVariable.clear();
    _totalIncome.clear();

    _incomeAfterTax.clear();
    _taxes.clear();
    _healthCare.clear();
    _discretionaryIncome.clear();
    _discretionarySpending.clear();

    _annuityTotalVector.clear();
    _annuityFixedVector.clear();
    _annuityVariableVector.clear();
    if ( !_lockMarketReturns->isChecked() )
    {
        _stockReturnsVector.clear();
        _bondReturnsVector.clear();
        _stockBondReturnsVector.clear();
    }
}

void SimWindow::extractQuartiles(dMatrix &inputVector, dVector &bottom10, dVector &bottomQuart, dVector &medianVector, dVector &topQuart)
{
    // inputVector[nTrials][nTime]
    // output vectors: sorted by trials

    bottom10.clear();
    bottomQuart.clear();
    medianVector.clear();
    topQuart.clear();

    int nTrials = inputVector.size();
    int nTime = _ageJoeVector.size();
    // input vector is [nTrials][nTime]; we need [nTime][nTrials]
    dMatrix swapInputVector;  swapInputVector.resize(nTime);

    for (int jTime=0; jTime<nTime; jTime++)
    {
        swapInputVector[jTime].resize(nTrials);
        for (int jTrial=0; jTrial<nTrials; jTrial++)
            swapInputVector[jTime][jTrial] = inputVector[jTrial][jTime];
        // savings
        utilMath::topDownMergeSort(swapInputVector[jTime]);
        int nTrials10   = nTrials * 0.10;
        int nTrialsLow  = nTrials * 0.25;
        int nTrialsMid  = nTrials * 0.50;
        int nTrialsHigh = nTrials * 0.75;
        double low10     = swapInputVector[jTime][nTrials10];
        double lowQuart  = swapInputVector[jTime][nTrialsLow];
        double median    = swapInputVector[jTime][nTrialsMid];
        double highQuart = swapInputVector[jTime][nTrialsHigh];
        bottom10.append(low10);
        bottomQuart.append(lowQuart);
        medianVector.append(median);
        topQuart.append(highQuart);
    }
}

double SimWindow::getRandomStockReturn()
{
    double mean  = _marketReturnsRate;
    double sigma = _stockReturnsSigma;
    double marketReturn = GaussianRandomizer(mean,sigma,60.);
    return marketReturn;
}
double SimWindow::getRandomBondReturn()
{
    double mean  = _bondReturnsMean;
    double sigma = _bondReturnsSigma;
    double bondReturn = GaussianRandomizer(mean,sigma,60.);
    return bondReturn;
}

dPoint2D SimWindow::getHistoricalMarketReturn()  // variable annuity appears not to allow bonds
{ // returns {stocks, bonds} for {variable annuity, market with stock/bond mixture}
    _yearIndexHistoricalMarket++;
    // every nth time or so, reset the starting point
    double oneOverN = 1./1000000.;
    if ( _radioRand->isChecked() && _radioHistoricalMarket->isChecked() && RAND_0_1 < oneOverN ) randomizeHistoricalStartingPoint();
    if ( _yearIndexHistoricalMarket >= _historicalReturnsTable.size() )
        _yearIndexHistoricalMarket = 0;

    dPoint2D marketReturn;
    double stockReturn = _historicalReturnsTable[_yearIndexHistoricalMarket][0];
    double bondReturn  = _historicalReturnsTable[_yearIndexHistoricalMarket][1];
    marketReturn.x = stockReturn;
    marketReturn.y = bondReturn;

    return marketReturn;
}

double SimWindow::GaussianRandomizer(double mean, double sigma, double cutoff)
{
    FUNC_ENTER;
    double yGauss, y, x;
    do
    {
        /* Choose x between +- the cutoff */
        //        x = cutoff * ( -1. + 2. * QRandomGenerator::global()->generateDouble() );
        x = cutoff * ( -1. + 2. * RAND_0_1 );
        /* Compute the Gaussian function of x. */
        double arg = - 0.5 * SQR(x-mean) /SQR(sigma);
        yGauss = qExp( arg );
        /* Choose a random y from 0 to 1. */
        //        y = QRandomGenerator::global()->generateDouble();
        y = RAND_0_1;
    }
    while (y > yGauss);

    return( x );
}

void SimWindow::updateGraphs()
{
    FUNC_ENTER;

    QString numberString;

    if ( !_radioRand->isChecked() )
    { // this is either a fixed return or a single randomization
        runRetirement();
        createSummaryInfo();
        plotResults();
        _annuityIncomePerYearAv = qRound(_annuityIncomePerYearAv);
        QString annuityAmounts = QString("%1 / %2  &  %3").arg(qRound(_annuity2YearlyIncomeStart)).arg(qRound(_annuityFAYearlyIncomeStart)).
                                 arg(numberString.setNum(_annuityIncomePerYearAv,'g',3));
        _annuityIncomePerYearAvLabel->setText(annuityAmounts);
        _marketIncomePerYearAvLabel->setText(numberString.setNum(_marketIncomePerYearAv,'g',3));
    }
    else
    {
        randomRetirementAndSummarize();
        _annuityIncomePerYearAvLabel->setText("NA");
        _marketIncomePerYearAvLabel->setText("NA");
    }

    plotMarketReturns();
}

void SimWindow::updateYearlyMarketReturns(int iAge)
{ // calculate _yearlyStockReturn, _yearlyBondReturn, _yearlyStockBondReturn
    if ( !useRandomization() )
        _yearlyStockReturn = _yearlyBondReturn = _yearlyStockBondReturn = _marketReturnsRate;
    else
    {
        if ( _lockMarketReturns->isChecked() )
        {
            _yearlyStockReturn = _stockReturnsVector.at(iAge);
            _yearlyBondReturn  = _bondReturnsVector.at(iAge);
        }
        else if ( _radioHistoricalMarket->isChecked() )
        {
            dPoint2D returns   = getHistoricalMarketReturn();
            _yearlyStockReturn = returns.x;
            _yearlyBondReturn  = returns.y;
        }
        else
        {
            _yearlyStockReturn = getRandomStockReturn();
            _yearlyBondReturn  = getRandomBondReturn();
        }
        _yearlyStockBondReturn = mixedMarketReturn(_yearlyStockReturn,_yearlyBondReturn);
        int iRange = _marketReturnsProbability.size()/2 - 1;
        int iBin = (_yearlyStockBondReturn + iRange);
        if ( iBin >= 0 && iBin < _marketReturnsProbability.size() )
            _marketReturnsProbability[iBin]++;
    }
    if ( !_lockMarketReturns->isChecked() )
    {
        _stockReturnsVector.append(_yearlyStockReturn);
        _bondReturnsVector.append(_yearlyBondReturn);
        _stockBondReturnsVector.append(_yearlyStockBondReturn);
    }
}

dVector SimWindow::createTimeVector(int startAge, int endAge)
{
    dVector ageVector;
    for (int jAge=startAge; jAge<= endAge; jAge++)
    {
        double dAge = jAge;
        ageVector.append(dAge);
    }
    return ageVector;
}

dVector SimWindow::createMoneyVectorPre(dVector ageVector, double startMoney, double monthlyAddition, int endAddtion,
                                        double rate, double fees)
{
    dVector accumulation;
    double money=startMoney;

    for (int jAge=0; jAge<ageVector.size(); jAge++)
    {
        double age = ageVector.at(jAge);
        accumulation.append(money);
        double addition = 0.;
        if ( age < endAddtion ) addition = monthlyAddition*12 / 1000.;
        double interest = money * (rate-fees)/100.;
        double newMoney = addition + interest;
        money += newMoney;
        FUNC_INFO << "@ age" << ageVector.at(jAge) << "$$ =" << accumulation.at(jAge);
    }
    return accumulation;
}

void SimWindow::createTimeVectors()
{
    FUNC_ENTER << "_currentAgeJoe" << _currentAgeJoe << "ageLimit" << _ageLimit;
    _ageJoeVector   = createTimeVector(_ageStart, _ageLimit);
    _ageEmiriVector = _ageJoeVector;
    _ageLilikaVector= _ageJoeVector;
    for (int jAge=0; jAge<_ageJoeVector.size(); jAge++)
    {
        _ageEmiriVector[jAge]  += _currentAgeEmiri  - _currentAgeJoe;
        _ageLilikaVector[jAge] += _currentAgeLilika - _currentAgeJoe;
    }
}

void SimWindow::createFixedMoneyVectors()
{
    FUNC_ENTER;
    // Add Joe's SS
    FUNC_INFO << "age Emiri" << _ageEmiriVector;
    FUNC_INFO << "age Lilika" << _ageLilikaVector;

    // create static money vectors: salaries and SS
    //                       ageVector, ageStart, ageEnd, ageEligible, amount, ageShift, newAmount
    _salaryJoeVector  = createStaticMoneyVector(_ageJoeVector, _ageJoeVector.first(), _retireAgeJoe, 0, 62, _salaryJoe62, _salaryJoe);
    _salaryEmiriVector= createStaticMoneyVector(_ageEmiriVector, _ageJoeVector.first()-7, _retireAgeEmiri, _salaryEmiri);
    _mortgageVector   = createStaticMoneyVector(_ageJoeVector,_ageStart,_mortgageAgeJoe,_mortgage,_inflation);

    createSSVectors();

    double spiaPerYear = 0.;
    if ( useSPIA() ) spiaPerYear = _spiaAmount * _spiaWR/100.;
    int duration = 5;
    if ( _radioSpia10->isChecked() ) duration = 10;
    //                                        ageVector, ageStart,   ageEnd,       amount, inflation
    _spiaVector = createStaticMoneyVector(_ageJoeVector,_spiaAge,_spiaAge+duration,spiaPerYear,_inflation);
}

void SimWindow::createSSVectors()
{
    int nTime = _ageJoeVector.size();
    _ssJoeVector.resize(nTime);
    _ssEmiriVector.resize(nTime);
    _ssLilikaVector.resize(nTime);
    _ssSpousalChildInCareVector.resize(nTime);
    for (int jAge=0; jAge<nTime; jAge++)
    {
        double ageJoe = _ageJoeVector.at(jAge);
        double ageEmiri = ageJoe - 7;
        double ageLilika = ageJoe - 55;

        double ssJoe=_ssJoe;
        double salaryJoe=0;
        double ssDeductionJoeSalary = 0.;
        if  ( ageJoe < _retireAgeJoe && ageJoe < 67 ) salaryJoe = _salaryJoe;
        if ( ageJoe >= _ssAgeJoe && ageJoe < _firstDemiseAge )
        { // take SS but reduce it according to the salary
            ssDeductionJoeSalary = qMax(0.,(salaryJoe - 23.4)/2.);
            ssJoe = _ssJoe - ssDeductionJoeSalary * _ssJoe/(_ssJoe+_ssLilika);
        }

        // Lilika's SS will be reduced according to Joe's salary but also Emiri's salary
        double ssLilika=0.;
        if ( ageJoe >= _ssAgeJoe && ageLilika < 18 )
            ssLilika = _ssLilika - ssDeductionJoeSalary * _ssLilika/(_ssJoe+_ssLilika);

        // family max = 1.5 * Joe's 67-age benefit = 3 * _ssLilika
        double maxBenefitSpousal = 3. * _ssLilika - ssJoe - ssLilika;
        // Emiri's spousal SS will be reduced according to Joe's salary and Emiri's salary
        double ssSpousal=0.;
        double salaryEmiri=0;
        double ssDeductionSpousalSalary = 0.;
        if  ( ageEmiri < _retireAgeEmiri && ageEmiri < 67 ) salaryEmiri = _salaryEmiri;

        if ( ageJoe >= _ssAgeJoe && ageLilika < 16 )
        {
            // take spousal SS but reduce it according to Emiri's salary
            ssDeductionSpousalSalary = qMax(0.,(salaryEmiri - 23.4)/2.);
            ssSpousal = qMax(0.,maxBenefitSpousal - ssDeductionSpousalSalary);
        }
        // Emiri's ss will be reduced according to her salary
        double ssEmiri=_ssEmiri;
        if ( ageEmiri >= _ssAgeEmiri )
        { // take SS but reduce it according to the salary
            if ( ageEmiri < _firstDemiseAge )
            {
                double ssDeductionEmiriSalary = qMax(0.,(salaryEmiri - 23.4)/2.);
                ssEmiri = _ssEmiri - ssDeductionEmiriSalary;
            }
            else
                ssEmiri = ssJoe;
        }
        else
            ssEmiri = 0.;

        if ( ageJoe > _ssCrisisAgeJoe )
        {
            ssJoe     *= 1. - static_cast<double>(_ssCrisisDrop)/100.;
            ssEmiri   *= 1. - static_cast<double>(_ssCrisisDrop)/100.;
            ssLilika  *= 1. - static_cast<double>(_ssCrisisDrop)/100.;
            ssSpousal *= 1. - static_cast<double>(_ssCrisisDrop)/100.;
        }

        if ( ageJoe < _ssAgeJoe || ageJoe >= _firstDemiseAge )
        {
            if ( ageEmiri >= _ssAgeEmiri ) ssEmiri = ssJoe;
            ssJoe = 0.;
        }

        _ssJoeVector[jAge]                = ssJoe;
        _ssEmiriVector[jAge]              = ssEmiri;
        _ssLilikaVector[jAge]             = ssLilika;
        _ssSpousalChildInCareVector[jAge] = ssSpousal;
    }
}

dVector SimWindow::createStaticMoneyVector(dVector ageVector, int ageStart, int ageEnd, double amount)
{
    int ageEligible=0;
    int ageShift=ageEnd+1;
    double newAmount=0.;
    dVector outputVector = createStaticMoneyVector(ageVector, ageStart, ageEnd, ageEligible, ageShift, amount, newAmount);
    return outputVector;
}

dVector SimWindow::createStaticMoneyVector(dVector ageVector, int ageStart, int ageEnd, double amount, double inflation)
{
    int ageEligible=0;
    int ageShift=ageEnd+1;
    double newAmount=0.;
    dVector outputVector = createStaticMoneyVector(ageVector, ageStart, ageEnd, ageEligible, ageShift, amount, newAmount, inflation);
    return outputVector;
}

dVector SimWindow::createStaticMoneyVector(dVector ageVector, int ageStart, int ageEnd, int ageEligible, double ageShift,
                                           double amount, double newAmount)
{
    double inflation=0.;
    dVector outputVector = createStaticMoneyVector(ageVector, ageStart, ageEnd, ageEligible, ageShift, amount, newAmount, inflation);
    return outputVector;
}

dVector SimWindow::createStaticMoneyVector(dVector ageVector, int ageStart, int ageEnd, int ageEligible, double ageShift, double amount,
                                           double newAmount, double inflationPercent)
{
    dVector moneyVector;
    double inflationAdjustedAmount    = amount;
    double inflationAdjustedNewAmount = newAmount;
    double inflationAdjustment        = 1.-inflationPercent/100.;
    for (int jValue=0; jValue<ageVector.size(); jValue++)
    {
        double age = ageVector.at(jValue);
        if ( age >= ageStart && age >= ageEligible && age < ageEnd )
        {
            if ( age <= ageShift )
            {
                moneyVector.append(inflationAdjustedAmount);
            }
            else
            {
                moneyVector.append(inflationAdjustedNewAmount);
            }
            inflationAdjustedAmount *= inflationAdjustment;
            inflationAdjustedNewAmount *= inflationAdjustment;
        }
        else
            moneyVector.append(0.);
    }
    return moneyVector;
}

double SimWindow::findTimeAverageFromVector(dVector ageVector, dVector moneyVector, double age1, double age2)
{ // average from age1 to age2 (inclusive)
    double average = 0.;
    int nYears=0;
    for (int jAge=0; jAge<ageVector.size(); jAge++)
    {
        double age = ageVector.at(jAge);
        if ( age >= age1 && age <= age2 )
        {
            average += moneyVector[jAge];
            nYears++;
        }
    }
    if ( nYears > 0 ) average /= static_cast<double>(nYears);
    return average;
}

void SimWindow::createSummaryInfo()
{
    FUNC_ENTER;
    int nTime = _ageJoeVector.size();

    // calculate withdrawal rate abnd average income per year
    _marketIncomePerYearAv = 0.;
    _annuityIncomePerYearAv = 0.;
    bool tooLowIncome = false;
    double incomeMin=1.e10;
    for (int jAge=0; jAge<nTime; jAge++)
    {
        double totalMarketIncome  = _withdrawJoeMarketVector[jAge] + _withdrawEmiriMarketVector[jAge];
        double totalAnnuityIncome = _annuity2JoeVector[jAge] + _annuityFAJoeVector[jAge];
        _marketIncomePerYearAv += totalMarketIncome;
        _annuityIncomePerYearAv+= totalAnnuityIncome;
        double income = _totalIncome.at(jAge);
        if ( _ageJoeVector.at(jAge) >= _retireAgeJoe )
            tooLowIncome |= income < _spendingVectorRequired.at(jAge);
        if ( income < incomeMin ) incomeMin = income;
    }
    if ( tooLowIncome )
        _successRate->setText(QString("FAILURE: income %1").arg(incomeMin));
    else
        _successRate->setText("success (> income floor)");
    _marketIncomePerYearAv  /= static_cast<double>(_ageJoeVector.size());
    _annuityIncomePerYearAv /= static_cast<double>(_ageJoeVector.size());

    double average65_74 = findTimeAverageFromVector(_ageJoeVector, whichFinalIncome(), 65., 74.);
    double average75_84 = findTimeAverageFromVector(_ageJoeVector, whichFinalIncome(), 75., 84.);
    double average85_95 = findTimeAverageFromVector(_ageJoeVector, whichFinalIncome(), 85., 94.);
    QString numberString;
    _incomeAverage65_74->setText(numberString.setNum(average65_74,'g',3));
    _incomeAverage75_84->setText(numberString.setNum(average75_84,'g',3));
    _incomeAverage85_94->setText(numberString.setNum(average85_95,'g',3));

    double av = ( average65_74 + average75_84 + average85_95 ) / 3.;
    _averageIncomeBox->setTitle("average income = " + numberString.setNum(av,'g',3));
}

void SimWindow::createRandomSummaryInfo()
{
    FUNC_ENTER;
    QString numberString;

    double age1 = 65.;
    double age2 = 74.;
    double average = findTimeAverageFromVector(_ageJoeVector, _bottom10IncomeVector, age1, age2);
    _incomeAverage65_74->setText(numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _bottomQuartIncomeVector, age1, age2);
    _incomeAverage65_74->setText(_incomeAverage65_74->text() + " | " + numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _medianIncomeVector, age1, age2);
    _incomeAverage65_74->setText(_incomeAverage65_74->text() + " | " + numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _topQuartIncomeVector, age1, age2);
    _incomeAverage65_74->setText(_incomeAverage65_74->text() + " | " + numberString.setNum(average,'g',3));

    age1 = 75.;
    age2 = 84.;
    average = findTimeAverageFromVector(_ageJoeVector, _bottom10IncomeVector, age1, age2);
    _incomeAverage75_84->setText(numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _bottomQuartIncomeVector, age1, age2);
    _incomeAverage75_84->setText(_incomeAverage75_84->text() + " | " + numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _medianIncomeVector, age1, age2);
    _incomeAverage75_84->setText(_incomeAverage75_84->text() + " | " + numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _topQuartIncomeVector, age1, age2);
    _incomeAverage75_84->setText(_incomeAverage75_84->text() + " | " + numberString.setNum(average,'g',3));

    age1 = 85.;
    age2 = 94.;
    average = findTimeAverageFromVector(_ageJoeVector, _bottom10IncomeVector, age1, age2);
    _incomeAverage85_94->setText(numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _bottomQuartIncomeVector, age1, age2);
    _incomeAverage85_94->setText(_incomeAverage85_94->text() + " | " + numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _medianIncomeVector, age1, age2);
    _incomeAverage85_94->setText(_incomeAverage85_94->text() + " | " + numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _topQuartIncomeVector, age1, age2);
    _incomeAverage85_94->setText(_incomeAverage85_94->text() + " | " + numberString.setNum(average,'g',3));

    age1 = _retireAgeJoe;
    age2 = _ageLimit;
    average = findTimeAverageFromVector(_ageJoeVector, _bottom10IncomeVector, age1, age2);
    _averageIncomeBox->setTitle("average income = " + numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _bottomQuartIncomeVector, age1, age2);
    _averageIncomeBox->setTitle(_averageIncomeBox->title()   + " | " + numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _medianIncomeVector, age1, age2);
    _averageIncomeBox->setTitle(_averageIncomeBox->title()   + " | " + numberString.setNum(average,'g',3));
    average = findTimeAverageFromVector(_ageJoeVector, _topQuartIncomeVector, age1, age2);
    _averageIncomeBox->setTitle(_averageIncomeBox->title()   + " | " + numberString.setNum(average,'g',3));
}

bool SimWindow::allZero(dVector data)
{
    bool allZero=true;
    for (int j=0; j<data.size(); j++)
        allZero &= data.at(j) == 0.;
    return allZero;
}

dVector SimWindow::addMoneyVectors(dVector moneyVector1, dVector moneyVector2)
{ // assume a consistent starting point, and add vectors to create one with a length of the shorter of the two
    dVector sumVector;
    int iCounter=0;
    bool endOfVector=false;
    while (!endOfVector)
    {
        sumVector.append(moneyVector1[iCounter]+moneyVector2[iCounter]);
        iCounter++;
        endOfVector = (iCounter== moneyVector1.size()) || (iCounter== moneyVector2.size());
    }
    return sumVector;
}

double SimWindow::withdrawMarketIncomeWithGuardRails(double age, double ageRMD, double nonMarketIncome, int ageIncomeStarts, double targetWRPreRMD,
                                                     double &guardRailDollars, double &savings)
{
    double marketWR = guardRailDollars / savings * 100.;
    double delta = targetWRPreRMD / 10.;
    double lowerRail = targetWRPreRMD - delta;
    double upperRail = targetWRPreRMD + delta;
    if ( marketWR < lowerRail )  // then increase spending
        guardRailDollars *= 1.1;
    else if ( marketWR > upperRail ) // then decrease spending
        guardRailDollars /= 1.1;
    marketWR = guardRailDollars / savings * 100.;
    return withdrawMarketIncome(age, ageRMD, nonMarketIncome, ageIncomeStarts, marketWR, savings);
}

double SimWindow::withdrawMarketIncome(double age, double ageRMD, double salaryLikeIncome, int ageIncomeStarts, double withdrawalRatePreRMD,
                                       double &qualifiedSavings)
{
    FUNC_ENTER;
    double incomeRequired = getRequiredTotalIncome();

    // try to withdraw in this order (except must honor RMD)
    double incomeBuffer=0;
    double incomeQualified=0.;

    double incomeTotal = salaryLikeIncome;
    if ( age >= ageIncomeStarts )
    {
        if ( age >= ageRMD )
        //        if ( age > 0 )
        { // withdraw from qualified market funds first
            // get rid of the buffer (invest it)
            if ( _bufferAmount > 0 )
            {
                qualifiedSavings += _bufferAmount;
                _bufferAmount = 0.;
            }

            double yearsRMD = _RMDSlope * (age - ageRMD) + _RMDOffset;
            double withdrawalRate = 100./yearsRMD;
            incomeQualified = qMin(qualifiedSavings * withdrawalRate/100.,qualifiedSavings);
            qualifiedSavings -= incomeQualified;
            incomeTotal += incomeQualified;
        }
        else // before RMD age
        { // withdraw from buffer first
            double qualifiedIncomeTarget = (qualifiedSavings+_bufferAmount) * withdrawalRatePreRMD/100.;
            double incomeDesired = qMax(incomeRequired,salaryLikeIncome+qualifiedIncomeTarget);
            incomeRequired = qMax(incomeDesired,incomeRequired);

            if ( incomeTotal < incomeDesired )
            {
                // withdraw from qualified buffer if available
                incomeBuffer = qMin(incomeDesired-incomeTotal,_bufferAmount);
                _bufferAmount -= incomeBuffer;
                incomeTotal += incomeBuffer;
            }

            if ( incomeTotal < incomeDesired )
            {
                // withdraw qualified income
                incomeQualified = qMin(incomeDesired-incomeTotal,qualifiedSavings);
                qualifiedSavings -= incomeQualified;
                incomeTotal += incomeQualified;
            }

            if ( incomeTotal < incomeDesired )
            {
                // if income is STILL below the floor, the only choice is to take more from qualified savings
                double extraQualified = qMin(incomeDesired - incomeTotal,qualifiedSavings);
                incomeQualified += extraQualified;
                qualifiedSavings -= extraQualified;
            }
        } // < ageRMD
    } // age >= ageIncomeStarts

    // as a last resort, withdraw more from the qualified non-buffer savings
    if ( incomeTotal < incomeRequired )
    {
        // if income is STILL below the floor, the only choice is to take more from qualified savings
        double extraQualified = qMin(incomeRequired - incomeTotal,qualifiedSavings);
        incomeQualified += extraQualified;
        qualifiedSavings -= extraQualified;
    }

    return incomeQualified + incomeBuffer;
}

double SimWindow::withdrawMarketIncomeFloor(double age, double ageRMD, double salaryLikeIncome, double incomeTarget,
                                            int ageIncomeStarts, double &qualifiedSavings)
{
    double incomeQualified=0.;
    double incomeBuffer=0;

    double incomeTotal = salaryLikeIncome;
    if ( age >= ageIncomeStarts )
    {
        if ( age >= ageRMD )
        {  // withdraw from qualified market funds first
            double yearsRMD = _RMDSlope * (age - ageRMD) + _RMDOffset;
            double withdrawalRate = 100./yearsRMD;
            incomeQualified = qMin(qualifiedSavings * withdrawalRate/100.,qualifiedSavings);
            qualifiedSavings -= incomeQualified;
            incomeTotal += incomeQualified;

            if ( incomeTotal < incomeTarget )
            { // remove $$ from buffer first
                incomeBuffer = qMin(incomeTarget - incomeTotal,_bufferAmount);
                _bufferAmount -= incomeBuffer;
                incomeTotal   += incomeBuffer;
            }
            FUNC_INFO << "@ age" << age << "income" << incomeTotal << incomeBuffer << "qualifiedSavings =" << qualifiedSavings;
        }
        else // before RMD age
        {
            // withdraw from buffer next
            if ( incomeTotal < incomeTarget )
            {
                incomeBuffer = qMin(incomeTarget - incomeTotal,_bufferAmount);
                incomeTotal += incomeBuffer;
                _bufferAmount -= incomeBuffer;
            }
        }
    }

    if ( incomeTotal < incomeTarget )
    { // remove more qualified income
        double extraQualified = qMin(incomeTarget - incomeTotal,qualifiedSavings);
        incomeQualified += extraQualified;
        qualifiedSavings -= extraQualified;
        incomeTotal += extraQualified;
    }

    return incomeQualified + incomeBuffer;
}

void SimWindow::updateSS()
{
    FUNC_ENTER;
    QString numberString;

    double extraYearsJoe = _ssAgeJoe - 62;
    _ssJoe = (10*static_cast<int>(_ssBaseJoe62 * qPow(1.08,extraYearsJoe)))/10.;
    numberString.setNum(_ssJoe,'g',3);
    _ssJoeLabel->setText(numberString);

    double extraYearsEmiri = _ssAgeEmiri - 62;
    _ssEmiri = _ssBaseEmiri62 * qPow(1.08,extraYearsEmiri);
    numberString.setNum(_ssEmiri,'g',3);
    _ssEmiriLabel->setText(numberString);

    double ssJoeAge67 = (10*static_cast<int>(_ssBaseJoe62 * qPow(1.08,5)))/10.;
    _ssLilika  = ssJoeAge67/2.;
    numberString.setNum(_ssLilika,'g',3);
    _ssLilikaLabel->setText(numberString);

    _ssSpousal = (1.5*ssJoeAge67 - _ssJoe)/2.;
    numberString.setNum(_ssSpousal,'g',3);
    _ssSpousalLabel->setText(numberString);
}

void SimWindow::updateAnnuities()
{
    runRetirement();
}

double SimWindow::payTaxes(int whichAge, double totalIncome)
{
    // FICA is no longer withheld, but most of this will be needed for health care premiums
    // 7.6% FICA means ~ $15k for Joe+Emiri
    // This will needed by Emiri prior to age 65

    double taxableIncomeFederal = totalIncome;
    double taxableIncomeState   = taxableIncomeFederal;

    // 15% of SS taxes not taxable at federal level
    taxableIncomeFederal -= (1.-_SSTaxableFed) * _ssJoeVector.at(whichAge);
    taxableIncomeFederal -= (1.-_SSTaxableFed) * _ssEmiriVector.at(whichAge);
    taxableIncomeFederal -= (1.-_SSTaxableFed) * _ssSpousalChildInCareVector.at(whichAge);
    // minor child SS is not taxable (tax based upon child)
    taxableIncomeFederal -= _ssLilikaVector.at(whichAge);

    // no state tax on SS
    taxableIncomeState -= _ssJoeVector.at(whichAge);
    taxableIncomeState -= _ssEmiriVector.at(whichAge);
    taxableIncomeState -= _ssLilikaVector.at(whichAge);
    taxableIncomeState -= _ssSpousalChildInCareVector.at(whichAge);

    // standard deductions
    taxableIncomeFederal-= _standardDeduction;
    taxableIncomeState  -= _standardDeduction;
    taxableIncomeState  -= _extraDeductionMA;

    double taxState = taxableIncomeState * _stateTaxPercent/100.;
    double taxFed   = getFederalTax(taxableIncomeFederal);

    // subtract FICA while working
    double FICA=0;
    if ( _ageJoeVector.at(whichAge) < _retireAgeJoe )
        FICA += _salaryJoeVector.at(whichAge) * _FICATaxPercent/100.;
    if ( _ageEmiriVector.at(whichAge) < _retireAgeEmiri )
        FICA +=_salaryEmiriVector.at(whichAge) * _FICATaxPercent/100.;

    //    QINFO << "tax at age" << _ageJoeVector.at(whichAge);
    //    QINFO << "federal state FICA = " << taxFed << taxState << FICA;

    double taxes = taxFed + taxState + FICA;
    return taxes;
}

double SimWindow::getFederalTax(double income)
{
    //    qInfo() << "income" << income;
    int nBrackets = _taxBrackets.size();
    dVector percentVector, cutoffVector;
    for (int jBracket=0; jBracket<nBrackets; jBracket++)
    {
        percentVector.append(_taxBrackets[jBracket][0]);
        cutoffVector.append(_taxBrackets[jBracket][1]);
    }

    double tax = qMin(income,cutoffVector[0]) * percentVector[0]/100.;
    //    qInfo() << "tax0" << tax;

    if ( income > cutoffVector[0] )
        tax += (qMin(income,cutoffVector[1])-cutoffVector[0]) * percentVector[1] / 100.;
    //    qInfo() << "tax1" << tax;

    if ( income > cutoffVector[1] )
        tax += (qMin(income,cutoffVector[2])-cutoffVector[1]) * percentVector[2] / 100.;
    //    qInfo() << "tax2" << tax;

    if ( income > cutoffVector[2] )
        tax += (qMin(income,cutoffVector[3])-cutoffVector[2]) * percentVector[3] / 100.;
    //    qInfo() << "tax3" << tax;

    if ( income > cutoffVector[3] )
        tax += (qMin(income,cutoffVector[4])-cutoffVector[3]) * percentVector[4] / 100.;
    //    qInfo() << "tax4" << tax;

    if ( income > cutoffVector[4] )
        tax += (qMin(income,cutoffVector[5])-cutoffVector[4]) * percentVector[5] / 100.;
    //    qInfo() << "tax5" << tax;

    if ( income > cutoffVector[5] )
        tax += (qMin(income,cutoffVector[6])-cutoffVector[5]) * percentVector[6] / 100.;
    //    qInfo() << "tax6" << tax;

    return tax;
}

double SimWindow::payHealthCare(int whichAge)
{
    double ageJoe      = _ageJoeVector.at(whichAge);
    double ageEmiri    = _ageEmiriVector.at(whichAge);
    double ageLilika   = _ageLilikaVector.at(whichAge);

    double healthCost = 0.;
    // assume the order is 1) Joe retires, 2) Emiri retires, 3) Lilika reaches age 24
    if ( ageEmiri < _retireAgeEmiri )
    {
        healthCost += _healthCarePreRetire;  // health care for whole family
        if ( ageJoe >= _retireAgeJoe ) healthCost += _extraMedicare;
    }
    else
    {
        if ( ageEmiri < 65 )
            healthCost += _healthCarePrivateAdult;
        else
            healthCost += _extraMedicare; // Emiri
        healthCost += _extraMedicare;     // Joe
        if ( ageLilika < 25 ) healthCost += _healthCarePrivateChild;
    }

    return healthCost;
}
