#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QMessageBox>
#include <qmath.h>
#include <QTextStream>
#include <QDebug>
#include "io.h"
#include "generalglm.h"
#include "petrtm.h"

using namespace utilIO;

dVector PETRTM::getEquilibrationVector(int iFile)  // returns k4 * k2 / k2a(t)
{ // get value and error (into x and y) from 2-parameter RTM model using k2 and k2a
    dVector equilibrationVector;  equilibrationVector.fill(0.,_nTimePerRun[iFile]);
    if ( _modelRTM == RTM_rFRTM2_R1 )
    {
        for (int jt=0; jt<_nTimePerRun[iFile]; jt++)
            equilibrationVector[jt] = (1.+_BPnd[iFile][jt]) / _tau4[iFile];
    }
    else
    {
        // Baseline k2 and k2a
        QChar k2ID  = _k2EventID[iFile];
        QChar k2aID = _k2aEventID[iFile];
        int iCoeffk2 = -1;  int iCoeffk2a = -1;
        for ( int jCoeff=0; jCoeff<getNumberCoefficients(); jCoeff++ )
        {
            if ( _basisID[jCoeff] == k2ID  ) iCoeffk2  = jCoeff;
            if ( _basisID[jCoeff] == k2aID ) iCoeffk2a = jCoeff;
        }
        if ( iCoeffk2 < 0 || iCoeffk2a < 0 ) return equilibrationVector;
        // k2 is constant
        double k2 = getBeta(iCoeffk2);
        // set k2a to be a constant vector in time (for now)
        dVector k2aVector; k2aVector.resize(_nTimePerRun[iFile]);
        for (int jt=0; jt<_nTimePerRun[iFile]; jt++)
            k2aVector[jt] = getBeta(iCoeffk2a);

        // Now challenges
        for ( int jChallenge=0; jChallenge<_maxChallenges; jChallenge++ )
        {
            if ( isGoodChallengeInRun(jChallenge, iFile) )
            {
                QChar challengeID = _challengeEventID[jChallenge];
                int iCoeffChallenge = getEventCoefficient(challengeID);
                if ( iCoeffChallenge < 0 ) continue;  // should never occur
                dVector shape;
                createChallengeShape(iFile, jChallenge, shape);
                for (int jt=0; jt<shape.size(); jt++)
                    k2aVector[jt] += getBeta(iCoeffChallenge) * shape[jt];
            }
        }
        bool allValid = k2 > 0.;
        for (int jt=0; jt<_nTimePerRun[iFile]; jt++)
            allValid &= k2aVector[jt] > 0.;
        if ( allValid )
        {
            for (int jt=0; jt<_nTimePerRun[iFile]; jt++)
                equilibrationVector[jt] = k2 / k2aVector[jt] / _tau4[iFile];
        }
    }
    return equilibrationVector;
}

dVector PETRTM::getBPndVector(int iFile)
{ // get value and error (into x and y)
    dVector BPndVector;  BPndVector.fill(0.,_nTimePerRun[iFile]);
    // Baseline k2 and k2a define baseline BPnd
    dPoint2D BP;        BP.x = BP.y = 0.;
    QChar k2ID  = _k2EventID[iFile];
    QChar k2aID = _k2aEventID[iFile];
    int iCoeffk2 = -1;  int iCoeffk2a = -1;
    for ( int jCoeff=0; jCoeff<getNumberCoefficients(); jCoeff++ )
    {
        if ( _basisID[jCoeff] == k2ID  ) iCoeffk2  = jCoeff;
        if ( _basisID[jCoeff] == k2aID ) iCoeffk2a = jCoeff;
    }
    if ( iCoeffk2 < 0 || iCoeffk2a < 0 ) return BPndVector;

    double k2     = getBeta(iCoeffk2);
    dVector k2aVector; k2aVector.resize(_nTimePerRun[iFile]);
    for (int jt=0; jt<_nTimePerRun[iFile]; jt++)
        k2aVector[jt] = getBeta(iCoeffk2a);

    // Now challenges
    for ( int jChallenge=0; jChallenge<_maxChallenges; jChallenge++ )
    {
        if ( isGoodChallengeInRun(jChallenge, iFile) )
        {
            QChar challengeID = _challengeEventID[jChallenge];
            int iCoeffChallenge = getEventCoefficient(challengeID);
            if ( iCoeffChallenge < 0 ) continue;  // should never occur
            dVector shape;
            createChallengeShape(iFile, jChallenge, shape);
            for (int jt=0; jt<shape.size(); jt++)
                k2aVector[jt] += getBeta(iCoeffChallenge) * shape[jt];
        }
    }
    bool allValid = k2 > 0.;
    for (int jt=0; jt<_nTimePerRun[iFile]; jt++)
        allValid &= k2aVector[jt] > 0.;
    if ( allValid )
    {
        for (int jt=0; jt<_nTimePerRun[iFile]; jt++)
            BPndVector[jt] = k2 / k2aVector[jt] - 1;
    }
    return BPndVector;
}

dPoint2D PETRTM::getBPndVersusTime(int iFile, int iTime)  // VERY INEFFICIENT METHOD: EVERY POINT ITIME CALLS createChallengeShape
{ // get value and error (into x and y)
    // Baseline k2 and k2a define baseline BPnd
    dPoint2D BP;        BP.x = BP.y = 0.;
    QChar k2ID  = _k2EventID[iFile];
    QChar k2aID = _k2aEventID[iFile];
    int iCoeffk2 = -1;  int iCoeffk2a = -1;
    for ( int jCoeff=0; jCoeff<getNumberCoefficients(); jCoeff++ )
    {
        if ( _basisID[jCoeff] == k2ID  ) iCoeffk2  = jCoeff;
        if ( _basisID[jCoeff] == k2aID ) iCoeffk2a = jCoeff;
    }
    if ( iCoeffk2 < 0 || iCoeffk2a < 0 ) return BP;
    double k2     = getBeta(iCoeffk2);
    double k2Err2 = getVar(iCoeffk2);
    double k2a    = getBeta(iCoeffk2a);

    // Now challenges
    double k2aErr2=0.;
    for ( int jChallenge=0; jChallenge<_maxChallenges; jChallenge++ )
    {
        if ( isGoodChallengeInRun(jChallenge, iFile) )
        {
            QChar challengeID = _challengeEventID[jChallenge];
            int iCoeffChallenge = getEventCoefficient(challengeID);
            if ( iCoeffChallenge < 0 ) continue;  // should never occur
            dVector shape;
            createChallengeShape(iFile, jChallenge, shape);
            k2a     += getBeta(iCoeffChallenge) * shape[iTime];
            k2aErr2 += getVar(iCoeffChallenge) * shape[iTime];
        }
    }
    if ( k2 > 0. && k2a > 0. )
    {
        BP.x= k2 / k2a - 1;
        BP.y = BP.x * qSqrt( k2Err2/k2/k2 + k2aErr2/k2a/k2a);
    }
    return BP;
}

dPoint2D PETRTM::getValueAndErrorForCurrentCondition()
{   // Get "k2" for the current condition. k2 is a parameter stored in "beta", whereas conditions are in "effectSize"
    // For a single condition, pointing to a coefficient, return beta. Otherwise return the effect size
    dPoint2D result;
    int iCoeff = getOnlyCoefficientInCondition(getCurrentCondition());
    if ( iCoeff < 0 )
    { // this is not a condition on a single coefficient
        result.x = getEffectSize();
        result.y = getEffectStDev();
    }
    else
    {
        result.x = getBeta(iCoeff);
        result.y = getSEM(iCoeff);
    }
    return result;
}

dPoint2D PETRTM::getBPndInCurrentCondition()
{
    dPoint2D BP;  BP.x = BP.y = 0.;

    qDebug() << "PETRTM::getBPndInCurrentCondition enter";
    // updateConditions();

    int iCondition = getCurrentCondition();
    // This condition should be composed entirely of k2a and challenge events.
    // Compute the BP for each k2a individually, and sum using the contrast matrix
    int nCoeffInCondition = getNumberCoefficientsInCondition(iCondition);
    for (int jCoeffInCondition=0; jCoeffInCondition<nCoeffInCondition; jCoeffInCondition++)
    {
        int iCoeff = _iCoeffInCondition[iCondition][jCoeffInCondition];
        dPoint2D par = averageParameter(_matchingk2InCondition[iCondition][jCoeffInCondition]);
        if ( par.x == 0. ) return BP;
        double k2=par.x;   double k2Err2=par.y;
        double k2a=0.;  double k2aErr2=0.;  double deltak2a = 0.;  double deltak2aErr2 = 0.;

        qDebug() << "PETRTM::getBPndInCurrentCondition type" << iCoeff << _basisShape[iCoeff];
        if ( iCoeff >=0 && _basisShape[iCoeff] == Type_k2a )
        {
            k2a     = getBeta(iCoeff);
            k2aErr2 = getVar(iCoeff);
        }
        else if ( iCoeff >=0 && _basisShape[iCoeff] == Type_challenge )
        {
            deltak2a     = getBeta(iCoeff);
            deltak2aErr2 = getVar(iCoeff);
            dPoint2D par = averageParameter(_matchingk2aInCondition[iCondition][jCoeffInCondition]);
            if ( par.x == 0. ) return BP;
            k2a=par.x;   k2aErr2=par.y;
        }

        if ( k2 !=0. && k2a != 0. )
        {
            if ( _basisShape[iCoeff] == Type_k2a )
            {
                double BPnd = k2/k2a - 1.;
                BP.x += BPnd * _contrastMatrix[iCondition][0][iCoeff];
                BP.y += BPnd * BPnd * ( k2Err2/k2/k2 + k2aErr2/k2a/k2a );
            }
            else // if ( _basisShape[iCoeff] == Type_challenge )
            {
                double totalk2a = k2a + deltak2a;
                double DBPnd = k2/k2a * deltak2a / totalk2a;
                BP.x += DBPnd * _contrastMatrix[iCondition][0][iCoeff];
                BP.y += DBPnd * DBPnd * ( k2Err2/k2/k2 + k2aErr2/k2a/k2a + deltak2aErr2/deltak2a/deltak2a +
                                          (k2aErr2+deltak2aErr2)/totalk2a/totalk2a );

            }
        }
    }
    BP.y = qSqrt(BP.y);
    return BP;
}

dPoint2D PETRTM::averageParameter(iVector iCoeffVector)
{  // x = variance-weighted mean, y = weighted variance
    dPoint2D result; result.x = 0.;  result.y = 0.;
    int nValues = iCoeffVector.size();
    if ( nValues == 0 ) return result;
    else if ( nValues == 1 )
    {
        int iCoeff = iCoeffVector[0];
        if ( iCoeff < 0 ) return result;
        result.x = getBeta(iCoeff);
        result.y = getVar(iCoeff);
        return result;
    }
    else
    {
        double weightSum = 0.;
        for ( int jPar=0; jPar<nValues; jPar++ )
        {
            int iCoeff = iCoeffVector[jPar];
            if ( iCoeff < 0 ) return result;
            double value  = getBeta(iCoeff);
            double weight = 1/getVar(iCoeff);
            weightSum += weight;
            result.x += weight * value;   // mean
        }
        if ( weightSum > 0. )
        {
            result.x /= weightSum;
            result.y = qSqrt(1./weightSum);  // e.g. for var1=var2=x, var = sqrt(x^2+x^2)/2 = var1/sqrt(2)
        }
        return result;
    }
}

double PETRTM::getBP0InRun(int iRun)
{ // make this more efficient
    if ( _k2EventCoefficient[iRun] < 0 || _k2aEventCoefficient[iRun] < 0 ) return 0.;

    double BPnd = 0.;
    if ( _modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1 )
        BPnd = _BPnd[iRun][0]; // BPnd is a fixed parameter
    else
    { // BPnd is a derived parameter
        double k2  = getBeta(_k2EventCoefficient[iRun]);
        double k2a = getBeta(_k2aEventCoefficient[iRun]);
        if ( k2 == 0. || k2a == 0. )
            return 0.;
        else
            BPnd = k2/k2a - 1.;
    }
    return BPnd;
}

dPoint2D PETRTM::getTau2InCurrentCondition()
{ // tau2 = 1 / k2
    dPoint2D tau2;  tau2.x = tau2.y = 0.; // save value and error in x,y
    if ( currentconditionIsk2Type() )
    {
        int iCoeffk2 = getOnlyCoefficientInCondition(getCurrentCondition());
        if ( iCoeffk2 >= 0 )
        {
            double k2    = getBeta(iCoeffk2);
            double k2Err = getSEM(iCoeffk2);
            if ( k2 != 0. )
            {
                tau2.x = 1. / k2;
                tau2.y = qAbs(tau2.x) * k2Err / k2;
            }
        }
    }
    return tau2;
}

dPoint2D PETRTM::getTau2RefInCurrentCondition()
{ // tau2_ref = R1 / k2
    dPoint2D tau2Ref;  tau2Ref.x = tau2Ref.y = 0.; // save value and error in x,y
    int iCondition = getCurrentCondition();
    if ( currentconditionIsR1Type() )
    {
        int nCoeffInCondition = getNumberCoefficientsInCondition(iCondition);
        for (int jCoeffInCondition=0; jCoeffInCondition<nCoeffInCondition; jCoeffInCondition++)
        {
            int iCoeffR1 = _iCoeffInCondition[iCondition][jCoeffInCondition];
            dPoint2D par = averageParameter(_matchingk2InCondition[iCondition][jCoeffInCondition]);
            if ( par.x == 0. ) return tau2Ref;
            double k2=par.x;   double k2Err2=par.y;
            double R1    = getBeta(iCoeffR1);
            double R1Err2 = getVar(iCoeffR1);
            if ( R1 != 0. && k2 != 0. )
            {
                tau2Ref.x = R1 / k2;
                tau2Ref.y = qAbs(tau2Ref.x) * qSqrt( k2Err2/k2/k2 + R1Err2/R1/R1 );
            }
            /*
            if ( iCoeffR1 >=0 && iCoeffk2 >= 0 )
            {
                double R1    = getBeta(iCoeffR1);
                double R1Err2 = getVar(iCoeffR1);
                double k2     = getBeta(iCoeffk2);
                double k2Err2 = getVar(iCoeffk2);
                if ( R1 != 0. && k2 != 0. )
                {
                    tau2Ref.x = R1 / k2;
                    tau2Ref.y = qAbs(tau2Ref.x) * qSqrt( k2Err2/k2/k2 + R1Err2/R1/R1 );
                }
            }
            */
        }
    }
    return tau2Ref;
}

void PETRTM::setTau2Ref(int iRun, double tau2Ref)
{
    qDebug() << "*** changed tau2Ref from" << _tau2RefSRTM[iRun] << "to" << tau2Ref;
    if ( isRTM2() )
    {   // set tau2Ref for the specified run
        if ( isSRTM() )
            _tau2RefSRTM[iRun] = tau2Ref;
        else
            _tau2RefFRTM[iRun] = tau2Ref;
    }
    else
    {
        // set tau2Ref for the specified run AND FOR ALL OTHER RUNS WITH MATCHING R1 and k2
        QChar R1 = getR1EventID(iRun);
        QChar k2 = getk2EventID(iRun);
        for ( int jRun=0; jRun<_nRuns; jRun++ )
        {
            if ( getR1EventID(jRun) == R1 && getk2EventID(jRun) == k2 )
            {
                if ( isSRTM() )
                    _tau2RefSRTM[jRun] = tau2Ref;
                else
                    _tau2RefFRTM[jRun] = tau2Ref;
            }
        }
    }
    setPrepared(false);
};


double PETRTM::getTau2RefInRun(int iRun)
{ // tau2_ref = R1 / k2
    double tau2Ref=0.;
    if ( _modelRTM == RTM_SRTM3 || _modelRTM == RTM_rFRTM3 || _modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1 )
    { // then tau2Ref is a derived parameter
        if ( _R1EventCoefficient[iRun] < 0 || _k2EventCoefficient[iRun] < 0 ) return 0.;
        double R1 = getBeta(_R1EventCoefficient[iRun]);
        double k2 = getBeta(_k2EventCoefficient[iRun]);

        if ( R1 != 0. && k2 != 0. )
            tau2Ref = R1 / k2;
    }
    else // tau2Ref is a fixed parameter
    {
        if ( isSRTM() )
            tau2Ref = _tau2RefSRTM[iRun];
        else
            tau2Ref = _tau2RefFRTM[iRun];
    }
    return tau2Ref;
}

dPoint2D PETRTM::getk2RefInRun(int iRun)
{ // k2_ref = k2 / R1
    dPoint2D k2Ref;  k2Ref.x = k2Ref.y = 0.; // save value and error in x,y
    if ( _modelRTM == RTM_SRTM3 || _modelRTM == RTM_rFRTM3 || _modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1 )
    {
        if ( _R1EventCoefficient[iRun] < 0 || _k2EventCoefficient[iRun] < 0 ) return k2Ref;
        double R1 = getBeta(_R1EventCoefficient[iRun]);
        double k2 = getBeta(_k2EventCoefficient[iRun]);
        double R1err = getSEM(_R1EventCoefficient[iRun]);
        double k2err = getSEM(_k2EventCoefficient[iRun]);
        if ( R1 != 0. && k2 != 0. )
        {
            k2Ref.x = k2 / R1;
            k2Ref.y = k2Ref.x * qSqrt( (R1err/R1)*(R1err/R1) + (k2err/k2)*(k2err/k2) );
        }
    }
    else
    {
        if ( isSRTM() )
            k2Ref.x = 1./_tau2RefSRTM[iRun];
        else
            k2Ref.x = 1./_tau2RefFRTM[iRun];
        k2Ref.y = 0.;
    }
    return k2Ref;
}

dPoint2D PETRTM::getR1InRun(int iRun)
{
    dPoint2D R1;  R1.x = R1.y = 0.; // save value and error in x,y
    if ( _modelRTM == RTM_SRTM3 || _modelRTM == RTM_rFRTM3 || _modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1 )
    {
        if ( _R1EventCoefficient[iRun] < 0 ) return R1;
        R1.x = getBeta(_R1EventCoefficient[iRun]);
        R1.y = getSEM(_R1EventCoefficient[iRun]);
    }
    return R1;
}

dPoint2D PETRTM::getk2InRun(int iRun)
{
    qDebug() << "PETRTM::getk2InRun enter" << iRun;
    dPoint2D k2;  k2.x = k2.y = 0.; // save value and error in x,y
    qDebug() << "PETRTM::getk2InRun size" << _k2EventCoefficient.size();
    qDebug() << "PETRTM::getk2InRun index" << _k2EventCoefficient[iRun] << _k2EventID[iRun];
    k2.x = getBeta(_k2EventCoefficient[iRun]);
    k2.y = getSEM(_k2EventCoefficient[iRun]);
    qDebug() << "PETRTM::getk2InRun exit" << k2.x;
    return k2;
}

dPoint2D PETRTM::getk2aInRun(int iRun)
{
    dPoint2D k2a;  k2a.x = k2a.y = 0.; // save value and error in x,y
    k2a.x = getBeta(_k2aEventCoefficient[iRun]);
    k2a.y = getSEM(_k2aEventCoefficient[iRun]);
    return k2a;
}

dPoint2D PETRTM::getk2RefInCurrentCondition()
{ // tau2_ref = R1 / k2
    dPoint2D k2Ref;  k2Ref.x = k2Ref.y = 0.; // save value and error in x,y
    int iCondition = getCurrentCondition();
    if ( currentconditionIsR1Type() )
    {
        int nCoeffInCondition = getNumberCoefficientsInCondition(iCondition);
        for (int jCoeffInCondition=0; jCoeffInCondition<nCoeffInCondition; jCoeffInCondition++)
        {
            int iCoeffR1 = _iCoeffInCondition[iCondition][jCoeffInCondition];
            dPoint2D par = averageParameter(_matchingk2InCondition[iCondition][jCoeffInCondition]);
            if ( par.x == 0. ) return k2Ref;
            double k2=par.x;   double k2Err2=par.y;
            double R1     = getBeta(iCoeffR1);
            double R1Err2 = getVar(iCoeffR1);
            if ( R1 != 0. && k2 != 0. )
            {
                k2Ref.x = k2 / R1;
                k2Ref.y = qAbs(k2Ref.x) * qSqrt( k2Err2/k2/k2 + R1Err2/R1/R1 );
            }
/*
            int iCoeffk2 = _matchingk2InCondition[iCondition][jCoeffInCondition];
            if ( iCoeffR1 >=0 && iCoeffk2 >= 0 )
            {
                double R1    = getBeta(iCoeffR1);
                double R1Err2 = getVar(iCoeffR1);
                double k2     = getBeta(iCoeffk2);
                double k2Err2 = getVar(iCoeffk2);
                if ( R1 != 0. && k2 != 0. )
                {
                    k2Ref.x = k2 / R1;
                    k2Ref.y = qAbs(k2Ref.x) * qSqrt( k2Err2/k2/k2 + R1Err2/R1/R1 );
                }
            }
*/
        }
    }
    return k2Ref;
}

QString PETRTM::getCurrentConditionTypeString()
{
    QString typeString;
    if ( getCurrentConditionShape() == Type_R1 )
        typeString = "R1";
    else if ( getCurrentConditionShape() == Type_k2 )
        typeString = "k2";
    else if ( getCurrentConditionShape() == Type_k2a )
        typeString = "k2a";
    else if ( getCurrentConditionShape() == Type_dCrdt )
        typeString = "dCr/dt";
    else if ( getCurrentConditionShape() == Type_challenge )
        typeString = "dk2a";
    return typeString;
}

void PETRTM::setTissueVector(bool newTAC, dMatrix tissueRegion)
{
    QMutex mutex;
    mutex.lock();
    _tissRegionRaw = tissueRegion;
    _tissRegion    = tissueRegion;

    if ( !newTAC )
    { // this section is experimental and should not be utilized for now; call via createRunBasisFunction()
        for ( int jRun=0; jRun<_nRuns; jRun++ )
        {
            int iCoeff_dCrdt = _dCrdtEventCoefficient[jRun];
            if ( iCoeff_dCrdt >= 0 )
            {
                for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                    _tissRegion[jRun][jt] -= getRegressor(iCoeff_dCrdt,jt);
            }
        }
    }

    if ( _smoothingScale != 0. )
        fitLoessCurve(_tissRegion);  // xxx
    _tissRegionDeriv = _tissRegion;  // this ensure dimensions are correct
    if ( (_modelRTM == RTM_rFRTM3 || _modelRTM == RTM_rFRTM2 || _modelRTM == RTM_rFRTM2_R1 || _PETWeightingModel == Weights_Custom)
         && _tau4Default != 0. )
        differentiateByRun(_tissRegionDeriv); // rFRTM requires a tissue derivative for the convolution term
    mutex.unlock();
    setPrepared(false);
}

dPoint2D PETRTM::getReferenceRegion(bool useFit, int iRun, int iTime)
{ // pass useFit explicitly, so that this can be called externally to get both raw and fit
    dPoint2D returnValue;
    returnValue.x= getTimeInRun(iRun,iTime);
    if ( useFit )
        returnValue.y= _refRegion[iRun][iTime];     // could be raw or fit
    else
        returnValue.y= _refRegionRaw[iRun][iTime];  // always raw
    return returnValue;
}
dPoint2D PETRTM::getTissueRegion(bool useFit, int iRun, int iTime)
{
    dPoint2D returnValue;
    returnValue.x= getTimeInRun(iRun,iTime);
    if ( useFit )
        returnValue.y= _tissRegion[iRun][iTime];    // could be either raw or fit
    else
        returnValue.y= _tissRegionRaw[iRun][iTime]; // always raw
    return returnValue;
}

dPoint2D PETRTM::getReferenceRegionTimesR1(int iRun, int iTime)
{ // legacy function (xxx jbm)
    dPoint2D returnValue;
    returnValue.x= getTimeInRun(iRun,iTime);

    if ( _modelRTM == RTM_SRTM3 || _modelRTM == RTM_rFRTM3 || _modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1 )
    { // get R1
        QChar eventID = _R1EventID[iRun];
        int iCoeff = -1;
        for ( int jCoeff=0; jCoeff<getNumberCoefficients(); jCoeff++ )
        {
            if ( _basisShape[jCoeff] == Type_R1 && _basisID[jCoeff] == eventID )
                iCoeff = jCoeff;
        }
        if ( iCoeff < 0 )
            qFatal("Fatal error in getReferenceRegionTimesR1: R1 coefficient not found");
        double R1 = getBeta(iCoeff);
        returnValue.y = _refRegion[iRun][iTime] * R1;
    }
    else
    { // R1 = k2 * tau2_ref
        QChar eventID = _k2EventID[iRun];
        int iCoeff = -1;
        for ( int jCoeff=0; jCoeff<getNumberCoefficients(); jCoeff++ )
        {
            if ( _basisShape[jCoeff] == Type_k2 && _basisID[jCoeff] == eventID )
                iCoeff = jCoeff;
        }
        if ( iCoeff < 0 )
            qFatal("Fatal error in getReferenceRegionTimesR1: k2 coefficient not found");
        double k2 = getBeta(iCoeff);
        if ( isSRTM() )
            returnValue.y = _refRegion[iRun][iTime] * k2 * _tau2RefSRTM[iRun];
        else
            returnValue.y = _refRegion[iRun][iTime] * k2 * _tau2RefFRTM[iRun];
    }
    return returnValue;
}

void PETRTM::setReferenceRegionFromTableColumn(int iColumn)
{
    _referenceRegionTableColumn = iColumn;
    if ( iColumn >= _columnNames[0].count() )
        qFatal("Programming error in PETRTM::setReferenceRegionFromTableColumn");
    _refRegionName = "table:" + _columnNames[0].at(iColumn);
    dMatrix referenceRegion;
    referenceRegion.resize(_nRuns);
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        referenceRegion[jRun].fill(0.,_nTimePerRun[jRun]);
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
            referenceRegion[jRun][jt] = _table[jRun][jt][iColumn];
    }
    setReferenceRegion(referenceRegion);
}

void PETRTM::setReferenceRegion(dMatrix timeBins, dMatrix referenceRegionRaw)
{  // set the raw and fit version of the reference region
    qDebug() << "PETRTM::setReferenceRegion enter";
    int nRuns = referenceRegionRaw.size();
    setNumberRuns(nRuns);
    for (int jRun=0; jRun<nRuns; jRun++)
    {
        setTimePointsInRun(jRun,referenceRegionRaw[jRun].size());
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
            _table[jRun][jt][0] = timeBins[jRun][jt];
    }
    updateReferenceRegion();
    setPrepared(false);
    _referenceRegionIsDefined = true;
    qDebug() << "PETRTM::setReferenceRegion exit";
}

void PETRTM::setReferenceRegion(dMatrix referenceRegionRaw)
{  // set the raw and fit version of the reference region
    qDebug() << "PETRTM::setReferenceRegion enter";
    _refRegionRaw = referenceRegionRaw;
    updateReferenceRegion();
    setPrepared(false);
    _referenceRegionIsDefined = true;
    qDebug() << "PETRTM::setReferenceRegion exit";
}

void PETRTM::updateReferenceRegion()
{
    qDebug() << "PETRTM::updateReferenceRegion 1";
    _refRegion = _refRegionRaw;
    if ( _smoothingScale != 0. )
        fitLoessCurve(_refRegion);
    _refRegionIntegral = _refRegion;
    _refRegionDeriv    = _refRegion;
    integrateByRun(_refRegionIntegral);
    differentiateByRun(_refRegionDeriv);
    setPrepared(false);
}

QString PETRTM::createConditions()
{
    int nEvents = countEvents();
    QString conditions = "";
    // Create all single k2a events first
    for ( int jEvent=0; jEvent<nEvents; jEvent ++)
    {
        if ( _basisShape[jEvent] == Type_k2a || _basisShape[jEvent] == Type_challenge )
        {
            conditions.append(_basisID[jEvent]);
            conditions.append(' ');
        }
    }
    // Make all pair-wise subtractions of k2a or challenge events.
    for ( int jEvent=0; jEvent<getNumberCoefficients(); jEvent++ )
    {
        bool BP0_type = _basisShape[jEvent] == Type_k2a || _basisShape[jEvent] == Type_challenge;
        for ( int jEvent1=jEvent+1; jEvent1<getNumberCoefficients(); jEvent1++ )
        {
            bool BP1_type = _basisShape[jEvent1] == Type_k2a || _basisShape[jEvent1] == Type_challenge;
            if ( BP0_type && BP1_type )
            {
                QChar eventID  = _basisID.at(jEvent);
                QChar eventID1 = _basisID.at(jEvent1);
                conditions.append(eventID);
                conditions.append('-');
                conditions.append(eventID1);
                conditions.append(' ');
            }
        }
    }
    // Create all other events
    for ( int jEvent=0; jEvent<nEvents; jEvent ++)
    {
        if ( _basisShape[jEvent] == Type_k2 ||
             _basisShape[jEvent] == Type_R1 ||
             _basisShape[jEvent] == Type_dCrdt )
        {
            conditions.append(_basisID[jEvent]);
            conditions.append(' ');
        }
    }
    definePETConditions(conditions);
    setCurrentCondition(0);
    return conditions;
}

void PETRTM::updateConditions()
{
    countEvents();  // updates coefficient arrays
    QString conditionString = getConditionString();
    definePETConditions(conditionString);
}

void PETRTM::definePETConditions( QString conditionString )
{
    if ( ! (_modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1) )
        defineConditions(conditionString);

    // The rest of this function helps to turn BPnd into a "meta-GLM parameter"
    // by finding k2 and k2a indices that can form BP values together with conditions of type k2a or challenge
    int nConditions = getNumberConditions();
    if ( nConditions == 0 ) return;

    _iCoeffInCondition.resize(nConditions);
    _matchingk2InCondition.resize(nConditions);
    _matchingk2aInCondition.resize(nConditions);
    for (int jCondition=0; jCondition<nConditions; jCondition++)
    {
        int nCoefficientsInCondition = getNumberCoefficientsInCondition(jCondition);
        _iCoeffInCondition[jCondition].resize(nCoefficientsInCondition);
        _matchingk2InCondition[jCondition].resize(nCoefficientsInCondition);
        _matchingk2aInCondition[jCondition].resize(nCoefficientsInCondition);
        for ( int jCoeffInCondition=0; jCoeffInCondition<nCoefficientsInCondition; jCoeffInCondition++ )
        {
            _matchingk2InCondition[jCondition][jCoeffInCondition].resize(0);
            _matchingk2aInCondition[jCondition][jCoeffInCondition].resize(0);
            int iCoeff = _indexCoeffInCondition[jCondition][jCoeffInCondition];
            QChar ID = _basisID[iCoeff];
            _iCoeffInCondition[jCondition][jCoeffInCondition] = iCoeff;
            if ( _basisShape[iCoeff] == Type_k2a || _basisShape[iCoeff] == Type_R1 )
            {
                // Find the run for this event
                for ( int jRun=0; jRun<_nRuns; jRun ++ )
                {
                    if ( _k2aEventID[jRun] == ID || _R1EventID[jRun] == ID )
                    { // k2 and k2a needed for BP (including challenges); k2 need for tau2' with R1
                        _matchingk2InCondition[jCondition][jCoeffInCondition].append(_k2EventCoefficient[jRun]);
                        _matchingk2aInCondition[jCondition][jCoeffInCondition].append(_k2aEventCoefficient[jRun]);
                    }
                }
            }
            else
            { // challenge
                for ( int jChallenge=0; jChallenge<_maxChallenges; jChallenge++)
                {
                    if ( _challengeEventID[jChallenge] == ID )
                    {
                        for ( int jStim=0; jStim<_maxStimuli; jStim++)
                        {
                            if ( isGoodStimulus(jChallenge,jStim) )
                            {
                                int iRun = _challengeRun[jChallenge][jStim];
                                _matchingk2InCondition[jCondition][jCoeffInCondition].append(_k2EventCoefficient[iRun]);
                                _matchingk2aInCondition[jCondition][jCoeffInCondition].append(_k2aEventCoefficient[iRun]);
                            } // good stimulus
                        } // jStim
                    } // matching IDs
                } // jChallenge
            } // is challenge
        } // jCoeffInCondition
    }
    // Potentially reset the current condition
    if ( getCurrentCondition() >= getNumberConditions() ) setCurrentCondition(0);
}
int PETRTM::getNumberChallenges()
{
    int nChallenges=0;
    for ( int jRun=0; jRun<_nRuns; jRun++ )
    {
        for (int jChallenge=0; jChallenge<_maxChallenges; jChallenge++)
        {
            if ( isGoodChallengeInRun(jChallenge,jRun) )
                nChallenges++;
        }
    }
    return nChallenges;
}

int PETRTM::getNumberChallengesInRun(int iRun)
{
    int nChallenges=0;
    for (int jChallenge=0; jChallenge<_maxChallenges; jChallenge++)
    {
        if ( isGoodChallengeInRun(jChallenge,iRun) )
            nChallenges++;
    }
    return nChallenges;
}

QChar PETRTM::getFirstGoodChallenge()
{
    for (int jChallenge=0; jChallenge<_maxChallenges; jChallenge++)
    {
        for ( int jRun=0; jRun<_nRuns; jRun++ )
        {
            if ( isGoodChallengeInRun(jChallenge,jRun) )
            {
                QChar eventID = getEventChar(jChallenge);
                return eventID;
            }
        }
    }
    return -1;
}
int PETRTM::getFirstGoodChallengeIndex()
{
    for ( int jRun=0; jRun<_nRuns; jRun++ )
    {
        for (int jChallenge=0; jChallenge<_maxChallenges; jChallenge++)
        {
            if ( isGoodChallengeInRun(jChallenge,jRun) )
                return jChallenge;
        }
    }
    return -1;
}
int PETRTM::getFirstGoodChallengeIndexInRun( int iRun )
{
    for (int jChallenge=0; jChallenge<_maxChallenges; jChallenge++)
    {
        if ( isGoodChallengeInRun(jChallenge,iRun) )
            return jChallenge;
    }
    return -1;
}

bool PETRTM::isGoodChallenge(QChar ID)
{
    int indexChallenge = getEventIndex(ID);
    return isGoodChallenge(indexChallenge);
}

bool PETRTM::isGoodChallenge( int indexChallenge )
{
    bool goodChallenge = false;
    for ( int jStim=0; jStim<_maxStimuli; jStim++ )
        goodChallenge |= isGoodStimulus(indexChallenge,jStim);
    return goodChallenge;
}

bool PETRTM::isGoodChallengeInRun( QChar ID, int iRun )
{
    int indexChallenge = getEventIndex(ID);
    return isGoodChallengeInRun(indexChallenge, iRun);
}

bool PETRTM::isGoodChallengeInRun( int indexChallenge, int iRun )
{
    bool goodChallenge = false;
    for ( int jStim=0; jStim<_maxStimuli; jStim++ )
        goodChallenge |= ( isGoodStimulus(indexChallenge,jStim) && iRun == _challengeRun[indexChallenge][jStim] );
    return goodChallenge;
}

bool PETRTM::isGoodStimulusInRun( int indexChallenge, int indexStimulus, int iRun)
{
    bool goodStimulusInRun = isGoodStimulus(indexChallenge, indexStimulus) && iRun == _challengeRun[indexChallenge][indexStimulus];
    return goodStimulusInRun;
}

bool PETRTM::isGoodStimulus( int indexChallenge, int indexStimulus )
{
    int iShape = _challengeShape[indexChallenge];
    if ( iShape == Challenge_none ) return false;
    int infoRequired = getChallengeInfoRequired(iShape);
    bool goodStimulus = true;
    int iRun = _challengeRun[indexChallenge][indexStimulus];
    goodStimulus &= iRun >= 0;
    if ( goodStimulus && infoRequired == 1 ) //  including goodStimulus in test prevents iRun=-1 from throwing an error
        goodStimulus &= _challengeOn[indexChallenge][indexStimulus] >= getTimeInRun(iRun,0);
    else if ( goodStimulus && infoRequired == 2 )
    { // stim OFF > ON && ON falls within run timing
        goodStimulus &= _challengeOff[indexChallenge][indexStimulus] > _challengeOn[indexChallenge][indexStimulus];
        goodStimulus &= _challengeOn[indexChallenge][indexStimulus]  >= 0.;
        goodStimulus &= _challengeOn[indexChallenge][indexStimulus]  < getTimeInRun(iRun,_nTimePerRun[iRun]-1);
    }
    return goodStimulus;
}
int PETRTM::getFirstBadStimulus(int indexChallenge)
{ // return the first bad stimulus (first available new stimulus index) for a given event index
    for (int jStimulus=0; jStimulus<_maxStimuli; jStimulus++)
    {
        if ( ! isGoodStimulus(indexChallenge, jStimulus) )
            return jStimulus;
    }
    _maxStimuli++;
    _challengeRun[indexChallenge].append(-1);
    _challengeOn[indexChallenge].append(0.);
    _challengeOff[indexChallenge].append(0.);
    return _maxStimuli-1;
}
int PETRTM::getFirstGoodStimulus(QChar ID)
{ // return the first bad stimulus (first available new stimulus index) for a given event index
    int indexChallenge = getEventIndex(ID);
    return getFirstGoodStimulus(indexChallenge);
}
int PETRTM::getFirstGoodStimulus(int indexChallenge)
{ // return the first bad stimulus (first available new stimulus index) for a given event index
    for (int jStimulus=0; jStimulus<_maxStimuli; jStimulus++)
    {
        if ( isGoodStimulus(indexChallenge, jStimulus) )
            return jStimulus;
    }
    return -1;
}
int PETRTM::getFirstGoodStimulusInRun(int indexChallenge, int iRun)
{ // return the first bad stimulus (first available new stimulus index) for a given event index
    for (int jStimulus=0; jStimulus<_maxStimuli; jStimulus++)
    {
        if ( isGoodStimulusInRun(indexChallenge, jStimulus, iRun) )
            return jStimulus;
    }
    return -1;
}

int PETRTM::getChallengeInfoRequired(int iShape)
{ // 0 = run, 1 = run + onset, 2 = run + onset + offset
    if ( iShape == Challenge_Gamma || iShape == Challenge_Sigmoid )
        return 1; // requires run identifier + onset time
    else if ( iShape == Challenge_Square || iShape == Challenge_RampUp || iShape == Challenge_RampDown )
        return 2; // requires run identifier + onset time + offset time
    else
        return 0;
}

void PETRTM::setR1EventID(int iRun, QChar ID)
{
    if ( iRun < 0 ) return;
    _R1EventID[iRun] = ID;
    setPrepared(false);
    setBasisFunctionsChanged();
};

void PETRTM::setdCrdtEventID(int iRun, QChar ID)
{
    if ( iRun < 0 ) return;
    _dCrdtEventID[iRun] = ID;
    setPrepared(false);
    setBasisFunctionsChanged();
};

void PETRTM::setNumberRuns(int nFiles)
{ // This is called once or if the number of files changes
    _nRuns = nFiles;
    if ( _nRuns == 0 ) return;

    if ( _nRuns != _nTimePerRun.size() )
    { // the following block needs only 1 allocation
        _R1EventID.fill('1',_nRuns);
        _k2EventID.fill('a',_nRuns);
        _k2aEventID.fill('A',_nRuns);
        _dCrdtEventID.fill('v',_nRuns);
        _R1EventCoefficient.resize(_nRuns);
        qDebug() << "resize _k2EventCoefficient";
        _k2EventCoefficient.resize(_nRuns);
        _k2aEventCoefficient.resize(_nRuns);
        _dCrdtEventCoefficient.resize(_nRuns);
        _nTimePerRun.resize(_nRuns);
        _frameFiles.fill("",_nRuns);
        _tau4.fill(0.,_nRuns);
        _tau2RefSRTM.fill(3.5,_nRuns);
        _tau2RefFRTM.fill(1.5,_nRuns);
        _BPnd.resize(_nRuns);
        _weights.resize(_nRuns);
        _ignoreString.resize(_nRuns);
        _table.resize(_nRuns);     // [_nRuns][_nTimeInRun][_nColumns]
        _quadLOESS.resize(_nRuns);     // _nTimeInRun
        _columnNames.resize(_nRuns);

        _refRegion.resize(_nRuns);
        _refRegionRaw.resize(_nRuns);
        _refRegionIntegral.resize(_nRuns);
        _refRegionDeriv.resize(_nRuns);
        _tissRegionRaw.resize(_nRuns);
        _tissRegion.resize(_nRuns);
        _tissRegionDeriv.resize(_nRuns);
        _frtmConvolution.resize(_nRuns);
        _frtmConvolutionRaw.resize(_nRuns);
        _frtmConvolutionRelative.resize(_nRuns);

        _challengeEventID.resize(_maxChallenges);
        for (int jChallenge=0; jChallenge<_maxChallenges; jChallenge++)
        {
            int iAscii;
            if ( jChallenge < 9 )
                iAscii = jChallenge + 49;  // event index 0 = '1' = ascii 49
            else if ( jChallenge < 35)
                iAscii = jChallenge + 88;  // event index 9 = 'a' = ascii 97
            else
                iAscii = jChallenge + 30;  // event index 35 = 'A' = ascii 65
            _challengeEventID[jChallenge] = iAscii;
        }

        _challengeShape.fill(Challenge_none,_maxChallenges);
        _tau.fill(10.,_maxChallenges);
        _alpha.fill(1.,_maxChallenges);

        _challengeRun.resize(_maxChallenges);
        _challengeOn.resize(_maxChallenges);
        _challengeOff.resize(_maxChallenges);
        for ( int jChallenge=0; jChallenge<_maxChallenges; jChallenge++ )
        {
            _challengeRun[jChallenge].fill(-1,_maxStimuli);
            _challengeOn[jChallenge].fill(0.,_maxStimuli);
            _challengeOff[jChallenge].fill(0.,_maxStimuli);
        }

        setWeightingModel(_PETWeightingModel);
    }

    setPrepared(false);
}

int PETRTM::writeReferenceRegion(int iRun, QString fileName, QString RRName)
{
    int nTime = _table[iRun].size();
    if ( nTime != _nTimePerRun[iRun] )
    {
        _warningErrors->append(QString("Error: the number of time points in the table (%1) does not match the number of time points (%2)").
                               arg(nTime).arg(_nTimePerRun[iRun]));
        return(1);
    }
    else if ( !_referenceRegionIsDefined )
    {
        _warningErrors->append(QString("Error: the reference region is not yet defined"));;
        return(1);
    }

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return(1);
    QTextStream out(&file);

    dVector dtVector;  dtVector.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
        dtVector[jt] = _table[iRun][jt][0];
    _table[iRun].clear();  // this is necessary
    _table[iRun].resize(nTime);  // resize table vector
    out << "delta_time " << RRName << "\n";
    for (int jt=0; jt<nTime; jt++)
    {
        double dt = dtVector[jt] * 60.;  // convert to seconds
        if ( jt == 0 ) dt *= 2.;         // reset the first bin
        out << dt << " " << _refRegion[iRun][jt] << "\n";
        _table[iRun][jt].append(dtVector[jt]);
        _table[iRun][jt].append(_refRegion[iRun][jt]);
    }
    file.close();

    // success, so point to the new frame file
    _frameFiles[iRun] = fileName;
    _refRegionName = "table:" + RRName;
    _referenceRegionTableColumn = 1;
    _columnNames[iRun].clear();
    _columnNames[iRun].append("delta_time");
    _columnNames[iRun].append(RRName);

    setPrepared(false);

    return(0);
}

int PETRTM::readTimeBinsFile(int iRun, QString fileName)
{
    if ( iRun >= _frameFiles.size() ) _frameFiles.resize(iRun+1);
    QFileInfo info(fileName);
    _frameFiles[iRun] = info.absoluteFilePath();
//    if ( info.isRelative() )
//        info.makeAbsolute();

    utilIO::readTableFile(iRun, fileName, _columnNames[iRun], _table, *_warningErrors);

    int colon = _refRegionName.lastIndexOf(":");
    bool useTable = colon>0;
    if ( useTable )
    { // double-check the name
        QString tableName = _refRegionName.left(colon);
        useTable = !tableName.compare("table",Qt::CaseInsensitive);
    }
    if ( useTable )
    {
        int length = _refRegionName.length();
        int afterColon = length - colon - 1;
        QString regionName = _refRegionName.right(afterColon);
        int iFound = -1;
        for (int jColumn=0; jColumn<_columnNames[iRun].size(); jColumn++)
        {
            if ( !regionName.compare(_columnNames[iRun].at(jColumn),Qt::CaseInsensitive) )
                iFound = jColumn;
        }
        if ( iFound < 0 )
        {
            _warningErrors->append(QString("Error locating the reference region %1 in table file %2.").arg(regionName).arg(fileName));
            return(1);
        }
        else if ( _referenceRegionTableColumn > 0 && iFound != _referenceRegionTableColumn )
        {
            _warningErrors->append(QString("Error: use consistent column locations for reference regions in tables."));
            return(1);
        }
        else
            _referenceRegionTableColumn = iFound;
    }

    for ( int jt=0; jt<_table[iRun].size(); jt++)
    {
        double dt = _table[iRun][jt][0];
        dt /= 60.;                // seconds to minutes
//        if ( jt == 0 ) dt /= 2.;  // center the first bin
        _table[iRun][jt][0] = dt;
    }

    if ( _smoothingScale != 0. )
        defineLOESSFitting(iRun);

    setPrepared(false);

    return(0);
}

void PETRTM::defineLOESSFitting(int iRun)
{
    // Define local quadratic smoothing for each this run: (one function for each time point)
    _quadLOESS[iRun].resize(_nTimePerRun[iRun]);
    for ( int jt=0; jt<_nTimePerRun[iRun]; jt++)
    {
        double time0 = getTimeInRun(iRun,jt);
        double halfWidth = _smoothingScale;
        int iLowSide=0;  double time = time0;  double width=0.;  double maxWidth=0.;
        while ( jt - iLowSide >= 0 && width < halfWidth )
        {
            time = getTimeInRun(iRun,jt-iLowSide);
            width = qAbs(time-time0);
            if ( width > maxWidth ) maxWidth = width;
            iLowSide++;
        }
        iLowSide--;
        int iHighSide=0;  time = time0;  width = 0.;
        while ( jt + iHighSide < _nTimePerRun[iRun] && width < halfWidth )
        {
            time = getTimeInRun(iRun,jt+iHighSide);
            width = qAbs(time-time0);
            if ( width > maxWidth ) maxWidth = width;
            iHighSide++;
        }
        iHighSide--;
        // Set weights
        dVector weight;  weight.clear();
        int iHalfWidth = qMax(1,qMin(iLowSide,iHighSide));
        for ( int jRel=-iHalfWidth; jRel<=iHalfWidth; jRel++ )
        {
            int iTime = jt + jRel;
            if ( iTime >= 0 && iTime < _nTimePerRun[iRun] )
            {
                time = getTimeInRun(iRun,iTime);
                width = qAbs(time-time0);
                double distance = width / maxWidth;  // maxWidth > 0 always
                if ( distance < 1. )
                    weight.append(qPow(1-qPow(distance,3),3));
            }
        }
        // define (nCoeff, nTime)
        _quadLOESS[iRun][jt].define(qMin(weight.size(),3),weight.size());  // qMax: don't define 3 terms with only 2 points
        _quadLOESS[iRun][jt].setWeights(weight);
    }

}

void PETRTM::saveTimeModelFiles(QString dirName, QStringList dataFileList)
{
    if ( dataFileList.count() != _nRuns )
        qFatal("Programming error: dataFileList does not match _nRuns in PETRTM::saveTimeModelFiles.");
    QDir directory = QDir(dirName);
    QString fileName = dirName + "/timeModel.dat";
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);

    // time-model
    bool allSameTau4=true;
    if ( isFRTM() )
    {
        for (int jRun=0; jRun<_nRuns; jRun++)
            allSameTau4 &= (_tau4[jRun] == _tau4[0]);
        if ( allSameTau4 )
            out << "time-model FRTM " << _tau4[0] << "\n";
        else
            out << "time-model FRTM\n";
    }
    else
        out << "time-model SRTM\n";

    // reference region
    if ( _refRegionName.size() != 0 )
        out << "reference-region " << _refRegionName << "\n";
    if ( _brainRegionName.size() != 0 )
        out << "brain-region " << _brainRegionName << "\n";

    // Smoothing scale
    if ( _smoothingScale != 0. )
        out << "smoothing-scale " << _smoothingScale << "\n";

    // weighting model
    if ( _PETWeightingModel == Weights_11C )
        out << "weights 11C\n";
    else if ( _PETWeightingModel == Weights_11C_Noiseless )
        out << "weights 11C-noiseless\n";
    else if ( _PETWeightingModel == Weights_18F )
        out << "weights 18F\n";
    else if ( _PETWeightingModel == Weights_18F_Noiseless )
        out << "weights 18F-noiseless\n";
    else if ( _PETWeightingModel == Weights_Custom )
        out << "weights Signal\n";

    // conditions
    QString conditionString = getConditionString();
    if ( conditionString.size() != 0 )
        out << "conditions " << conditionString << "\n";

    // list of scans
    out << "scans:\n";
    for ( int jFile=0; jFile<getNumberRuns(); jFile++ )
    {
        QString frameFileName = directory.relativeFilePath(_frameFiles[jFile]);
        out << dataFileList.at(jFile) << " pet" << jFile+1 << ".glm " << frameFileName + "\n";
    }
    file.close();

    for ( int jRun=0; jRun<getNumberRuns(); jRun++ )
    {
        QString runNumber;     runNumber.setNum(jRun+1);
        fileName = dirName + "/pet" + runNumber + ".glm";
        writeGLMFile(jRun,fileName,allSameTau4);
    } // jRun
}

void PETRTM::writeGLMFile(int iRun, QString fileName, bool allSameTau4)
{
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);
    // tau4
    if ( ! allSameTau4 )
        out << "1/k4 " << _tau4[iRun] << "\n\n";
    // tau2Ref
    if ( isRTM2() )
    {
        if ( isFRTM() )
            out << "1/k2' " << _tau2RefFRTM[iRun] << "\n\n";
        else
            out << "1/k2' " << _tau2RefSRTM[iRun] << "\n\n";
    }
    else
    {
        out << "0\n\n";
        out << _R1EventID[iRun] << " R1\n\n";    // R1 with 3-parameter versions
    }
    out << "k2 " << _k2EventID[iRun]  << "\n\n";       // always k2
    out << "k2a " << _k2aEventID[iRun] << "\n\n";      // always k2a
    if ( _dCrdtIncluded )
        out << "dCrdt " << _dCrdtEventID[iRun] << "\n\n";

    // Challenges?
    qInfo() << "# challenges = " << getNumberChallengesInRun(iRun);
    if ( getNumberChallengesInRun(iRun) > 0 )
    {
        for (int jChallenge=0; jChallenge<_maxChallenges; jChallenge++)
        {
            if ( isGoodChallengeInRun(jChallenge,iRun) )
            {
                out << "dk2a " << getEventChar(jChallenge);
                if ( _challengeShape[jChallenge] == Challenge_Constant )
                    out << " constant";
                else if ( _challengeShape[jChallenge] == Challenge_Square )
                    out << " square";
                else if ( _challengeShape[jChallenge] == Challenge_RampUp )
                    out << " ramp-up";
                else if ( _challengeShape[jChallenge] == Challenge_RampDown )
                    out << " ramp-down";
                else if ( _challengeShape[jChallenge] == Challenge_Sigmoid )
                    out << " sigmoid " << _tau[jChallenge];
                else if ( _challengeShape[jChallenge] == Challenge_Gamma )
                {
                    out << " gamma " << _tau[jChallenge];
                    if ( getChallengeAlpha(jChallenge) != 1. )
                        out << " " << _alpha[jChallenge];
                } // challenge shapes
                out << "\n";
                for (int jStimulus=0; jStimulus<_maxStimuli; jStimulus++)
                {
                    if ( isGoodStimulusInRun(jChallenge, jStimulus, iRun) )
                    {
                        if ( _challengeShape[jChallenge] == Challenge_Constant )
                            out << "\n";
                        else if ( _challengeShape[jChallenge] == Challenge_Square ||
                             _challengeShape[jChallenge] == Challenge_RampUp ||
                             _challengeShape[jChallenge] == Challenge_RampDown )
                            out << _challengeOn[jChallenge][jStimulus] << " "
                                << _challengeOff[jChallenge][jStimulus] << "\n";
                        else // gamma or sigmoid
                            out << _challengeOn[jChallenge][jStimulus] << "\n";
                    }
                }
                out << "\n";  // spaces between multiple challenges
            } // good challenge
        } // jChallenge
    } // # challenges > 0

    // ignored points?
    QString ignoredString = getIgnoredString(iRun);
    if ( ignoredString.size() != 0 )
        out << "\nignore " << ignoredString << "\n";
    file.close();
}

void PETRTM::writeGLMFileOldFormat(int iRun, QString fileName)
{
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);
    // tau2Ref
    if ( isRTM2() )
    {
        if ( isFRTM() )
            out << _tau2RefFRTM[iRun] << "\n\n";
        else
            out << _tau2RefSRTM[iRun] << "\n\n";
    }
    else
    {
        out << "0\n\n";
        out << _R1EventID[iRun] << " R1\n\n";    // R1 with 3-parameter versions
    }
    out << _k2EventID[iRun]  << " k2\n\n";       // always k2
    out << _k2aEventID[iRun] << " k2a\n\n";      // always k2a
    if ( _dCrdtIncluded )
        out << _dCrdtEventID[iRun] << " dCrdt\n\n";

    // Challenges?
    qInfo() << "# challenges = " << getNumberChallengesInRun(iRun);
    if ( getNumberChallengesInRun(iRun) > 0 )
    {
        for (int jChallenge=0; jChallenge<_maxChallenges; jChallenge++)
        {
            if ( isGoodChallengeInRun(jChallenge,iRun) )
            {
                out << getEventChar(jChallenge) << " dk2a";
                if ( _challengeShape[jChallenge] == Challenge_Constant )
                    out << " constant";
                else if ( _challengeShape[jChallenge] == Challenge_Square )
                    out << " square";
                else if ( _challengeShape[jChallenge] == Challenge_RampUp )
                    out << " ramp-up";
                else if ( _challengeShape[jChallenge] == Challenge_RampDown )
                    out << " ramp-down";
                else if ( _challengeShape[jChallenge] == Challenge_Sigmoid )
                    out << " sigmoid " << _tau[jChallenge];
                else if ( _challengeShape[jChallenge] == Challenge_Gamma )
                {
                    out << " gamma " << _tau[jChallenge];
                    if ( getChallengeAlpha(jChallenge) != 1. )
                        out << " " << _alpha[jChallenge];
                } // challenge shapes
                out << "\n";
                for (int jStimulus=0; jStimulus<_maxStimuli; jStimulus++)
                {
                    if ( isGoodStimulusInRun(jChallenge, jStimulus, iRun) )
                    {
                        if ( _challengeShape[jChallenge] == Challenge_Constant )
                            out << "\n";
                        else if ( _challengeShape[jChallenge] == Challenge_Square ||
                             _challengeShape[jChallenge] == Challenge_RampUp ||
                             _challengeShape[jChallenge] == Challenge_RampDown )
                            out << _challengeOn[jChallenge][jStimulus] << " "
                                << _challengeOff[jChallenge][jStimulus] << "\n";
                        else // gamma or sigmoid
                            out << _challengeOn[jChallenge][jStimulus] << "\n";
                    }
                }
            } // good challenge
        } // jChallenge
    } // # challenges > 0

    // ignored points?
    QString ignoredString = getIgnoredString(iRun);
    if ( ignoredString.size() != 0 )
        out << "\nignore " << ignoredString << "\n";
    file.close();
}

double PETRTM::getTimeInRun(int iRun, int iTime)
{
    if ( iRun < 0 || iTime < 0 ) return -1.;

    // This first bin is centered at 1/2 the bin width
    double dt0  = _table[iRun][0][0];
    double time = dt0/2.;
    for ( int jt=1; jt<iTime; jt++ )
    {
        double lastDelta = _table[iRun][jt-1][0];
        double thisDelta = _table[iRun][jt][0];
        time += lastDelta/2. + thisDelta/2.;  // half bin width on either side
    }
    return time;
}

void PETRTM::setChallengeShape(int indexChallenge, int iShape)
{
    _challengeShape[indexChallenge] = iShape;
    setPrepared(false);
}
void PETRTM::setChallengeRun(int indexChallenge, int indexStimulus, int iRun)
{
    _challengeRun[indexChallenge][indexStimulus] = iRun;
    setPrepared(false);
}
void PETRTM::setChallengeOnset(int indexChallenge, int indexStimulus, double time)
{
    _challengeOn[indexChallenge][indexStimulus] = time;
    setPrepared(false);
}
void PETRTM::setChallengeOffset(int indexChallenge, int indexStimulus, double time)
{
    _challengeOff[indexChallenge][indexStimulus] = time;
    setPrepared(false);
}
void PETRTM::setChallengeTau(int indexChallenge, double tau)
{
    _tau[indexChallenge] = tau;
    setPrepared(false);
}
void PETRTM::setChallengeAlpha(int indexChallenge, double alpha)
{
    _alpha[indexChallenge] = alpha;
    setPrepared(false);
}

bool PETRTM::getFrameStatus()
{
    bool allFilesRead = true;
    for (int jRun=0; jRun<_nRuns; jRun++)
        allFilesRead &= !_frameFiles[jRun].isEmpty();  // true if no file names are empty
    return allFilesRead;
}

void PETRTM::setTimePointsInRun(int iRun, int nTime)
{
    if ( iRun<0 || iRun >=_nRuns )
    {
        qWarning() << "Error: iRun = " << iRun << " but # runs = " << _nRuns;
        exit(1);
    }
    _nTimePerRun[iRun] = nTime;
    if ( nTime != _table[iRun].size() )
    {
        _table[iRun].resize(nTime);
        _quadLOESS[iRun].resize(nTime);
        for (int jt=0; jt<nTime; jt++)
            _table[iRun][jt].append(1.);  // initialize time bin widths (1st column=frames) to 1.
    }
    if ( _refRegion[iRun].size() != nTime || _tissRegion[iRun].size() != nTime )
    {
        _dCrdtEventCoefficient.fill(-1,_nTimePerRun[iRun]);
        _refRegion[iRun].fill(0.,_nTimePerRun[iRun]);
        _refRegionRaw[iRun].fill(0.,_nTimePerRun[iRun]);
        _refRegionIntegral[iRun].fill(0.,_nTimePerRun[iRun]);
        _refRegionDeriv[iRun].fill(0.,_nTimePerRun[iRun]);
        _tissRegionRaw[iRun].fill(0.,_nTimePerRun[iRun]);
        _tissRegion[iRun].fill(0.,_nTimePerRun[iRun]);
        _frtmConvolution[iRun].fill(0.,_nTimePerRun[iRun]);
        _frtmConvolutionRaw[iRun].fill(0.,_nTimePerRun[iRun]);
        _frtmConvolutionRelative[iRun].fill(0.,_nTimePerRun[iRun]);
        _BPnd[iRun].fill(0.,_nTimePerRun[iRun]);
    }
    setPrepared(false);
}

void PETRTM::prepare()
{
    qDebug() << "PETRTM::prepare enter";
    createAllBasisFunctions();
    qDebug() << "PETRTM::prepare 1";
    calculatePseudoInverse();
    qDebug() << "PETRTM::prepare 2" ;
    if ( ! (_modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1) )
        defineConditions(getConditionString());
    qDebug() << "PETRTM::prepare 3";
    setPrepared(true);
    qDebug() << "PETRTM::prepare exit";
}

void PETRTM::fitData(dMatrix timeSeriesVector, dMatrix &yFit)
{
    QVector<ROI_data> ROIVector;
    int nRuns = timeSeriesVector.size();
    ROIVector.resize(nRuns);
    for ( int jRun=0; jRun<nRuns; jRun++)
    {
        ROIVector[jRun].xTime.clear();
        ROIVector[jRun].ySignal = timeSeriesVector[jRun];
        for ( int jt=0; jt<timeSeriesVector[jRun].size(); jt++)
            ROIVector[jRun].xTime.append(jt+1);
    }
    fitData(ROIVector, yFit);
}


void PETRTM::fitData(QVector<ROI_data> timeSeriesVector, dMatrix &yFit)
{
    if ( _modelRTM == RTM_SRTM3 || _modelRTM == RTM_SRTM2 || _modelRTM == RTM_SRTM2_R1 || _tau4Default == 0. )
        fitDataByGLM(timeSeriesVector, yFit);
    else if ( _modelRTM == RTM_SRTM2x2 )
        fitDataByIterativeAlternating2ParametersModels(timeSeriesVector, yFit);
    else if ( _modelRTM == RTM_rFRTM3 || _modelRTM == RTM_rFRTM2 )
    {
        if ( _calibratingTau4 && _modelRTM == RTM_rFRTM2 )
            fitDataByIterativeAlternating2ParametersModels(timeSeriesVector, yFit);
        else
            fitDataByFRTMBasisFunctions(timeSeriesVector, yFit);
    }
}

void PETRTM::fitDataByGLM(QVector<ROI_data> timeSeriesVector, dMatrix &yFit)
{ // timeSeriesVector follows the ROI_data structure, with 1d vectors for x and y
    if ( !isReady() )
        prepare();
    // Copy the ROI input data into a simple vector
    int nTimeTotal = getNumberTimePoints();
    dVector data;
    data.fill(0.,nTimeTotal);
    int iStartRun = 0;
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
            data[jt+iStartRun] = timeSeriesVector[jRun].ySignal[jt];
        iStartRun += _nTimePerRun[jRun];
    }
    fitWLS(data,true);
    // Attach the fit.
    iStartRun = 0;
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
        {
            int iTimeTotal = jt + iStartRun;
            yFit[jRun][jt] = getFit(iTimeTotal);
        }
        iStartRun += _nTimePerRun[jRun];
    }
}

void PETRTM::fitDataByFRTMBasisFunctions(QVector<ROI_data> timeSeriesVector, dMatrix &yFit)
{
    if ( !isReady() )
        prepare();
    // Copy the ROI input data into a simple vector, and then do 1 fit
    int nTimeTotal = getNumberTimePoints();
    dVector data;   data.fill(0.,nTimeTotal);
    int iStartRun = 0;
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
            data[jt+iStartRun] = timeSeriesVector[jRun].ySignal[jt];
        iStartRun += _nTimePerRun[jRun];
    }
    fitWLS(data,true);

    // Now iterate with modified basis functions
    dVector BP0;  BP0.resize(_nRuns);
    bool goodIteration=true;
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        BP0[jRun] = getBP0InRun(jRun);
        goodIteration &= BP0[jRun] > 0.;
    }
    bool moreIterations=goodIteration;


    dPoint2D deltaBP; deltaBP.x = deltaBP.y = 0.;
    /*
    qInfo() << "\n";
    if ( getNumberConditions() >= 2 )
    {
        setCurrentCondition(2);
        deltaBP = getBPndInCurrentCondition();
    }
    qInfo() << "SRTM BP0 =" << BP0[0] << "deltaBP = " << deltaBP.x << "sigma2 = " << 1.;
    */

    _nIterations=0;
    double lastSigma2 = getSigma2();
    while (goodIteration && moreIterations)
    {
        resetAndCalculateFRTMConvolution(true);
        recreateAllFRTMBasisFunctions();
        fitWLS(data,true);
        double sigma2 = getSigma2();
        double diff = qAbs(sigma2-lastSigma2)/lastSigma2 * 100.;
        lastSigma2 = sigma2;
        moreIterations = diff > 0.05;  // 0.05%
        goodIteration = true;
        for (int jRun=0; jRun<_nRuns; jRun++)
        {
            BP0[jRun] = getBP0InRun(jRun);
            goodIteration &= BP0[jRun] > 0.;
        }
        _nIterations++;
        moreIterations &= _nIterations < 20;
/*
        if ( getNumberConditions() >= 2 )
            deltaBP = getBPndInCurrentCondition();
        qInfo() << "iterate" << _nIterations << "BP0 =" << BP0 << "deltaBP = " << deltaBP.x << deltaBP.y << "sigma2 = " << getSigma2();
*/
    }

    if ( !goodIteration )
        _nIterations *= -1;
    else
    {
        // Attach the fit.
        iStartRun = 0;
        for (int jRun=0; jRun<_nRuns; jRun++)
        {
            for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
            {
                int iTimeTotal = jt + iStartRun;
                yFit[jRun][jt] = getFit(iTimeTotal);
            }
            iStartRun += _nTimePerRun[jRun];
        }
    }
}

void PETRTM::fitDataByIterativeAlternating2ParametersModels(QVector<ROI_data> timeSeriesVector, dMatrix &yFit)
{
    // This should start in the standard rFRTM2 model
    if ( _modelRTM != RTM_SRTM2x2 && _modelRTM != RTM_rFRTM2 )
        exit(1);
    if ( _modelRTM == RTM_SRTM2x2 )
        _modelRTM = RTM_SRTM2;
    else
        _modelRTM = RTM_rFRTM2;
    prepare();

    // Calculate BPnd using k2/k2a formulation starting with _tau2RefSRTM
    if ( isSRTM() )
        fitDataByGLM(timeSeriesVector, yFit);
    else
        fitDataByFRTMBasisFunctions(timeSeriesVector, yFit);

    double sigma2, lastSigma2;
    lastSigma2 = sigma2 = getSigma2();
    int nIterations =0;
    bool moreIterations = true;
    while ( moreIterations )
    {
        // Extract BPnd values from k2/k2a formulation
        for (int jRun=0; jRun<_nRuns; jRun++)
        {
            for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
            {
                dPoint2D BP = getBPndVersusTime(jRun,jt);
                _BPnd[jRun][jt] = BP.x;
            }
        }
        // Switch to the rFRTM2_R1 model to extract _tau2Ref
        switchToAlternateModel_rFRTM2();
        prepare();

        if ( isSRTM() )
            fitDataByGLM(timeSeriesVector, yFit);
        else
            fitDataByFRTMBasisFunctions(timeSeriesVector, yFit);

        for (int jRun=0; jRun<_nRuns; jRun++)
            _tau2RefFRTM[jRun] = getTau2RefInRun(jRun);  // this is a derived parameter from rFRTM2_R1
        // Switch back to "standard" model
        switchToAlternateModel_rFRTM2();
        prepare();
//        resetAndCalculateFRTMConvolution(false);
//        recreateAllFRTMBasisFunctions();
        // Calculate BPnd

        if ( isSRTM() )
            fitDataByGLM(timeSeriesVector, yFit);
        else
            fitDataByFRTMBasisFunctions(timeSeriesVector, yFit);

        sigma2 = getSigma2();
        double diff = qAbs(sigma2-lastSigma2)/lastSigma2 * 100.;
        lastSigma2 = sigma2;
        moreIterations = diff > 0.01;
        moreIterations &= nIterations < 20;
        nIterations++;
    }
    if ( _modelRTM == RTM_SRTM2 )
        _modelRTM = RTM_SRTM2x2;
    else
        _modelRTM = RTM_rFRTM2;
}

void PETRTM::switchToAlternateModel_rFRTM2()
{
    if ( isRTM2_R1() )
    {
        if ( isSRTM() )
            _modelRTM = RTM_SRTM2;
        else
            _modelRTM = RTM_rFRTM2;
    }
    else
    {
        if ( isSRTM() )
            _modelRTM = RTM_SRTM2_R1;
        else
            _modelRTM = RTM_rFRTM2_R1;
    }
}

int PETRTM::getTotalTimeIndex(int iRun, int iTimeInRun)
{
    int iTimeTotal = iTimeInRun;
    for (int jRun=0; jRun<iRun; jRun++)
        iTimeTotal += _nTimePerRun[jRun];
    return iTimeTotal;
}

void PETRTM::readGLMIgnoreBlock(QTextStream *in_stream, int iRun, QString inputString)
{
    bool firstIgnored = true;
    setIgnoredPoints(iRun,firstIgnored,inputString);
    while (!in_stream->atEnd())
    {
        QString line = in_stream->readLine();
        QString unCommented = line.left(line.indexOf("#"));
        if ( unCommented.isEmpty() )
            break;  // end block with empty line
        setIgnoredPoints(iRun,firstIgnored,unCommented);
        firstIgnored = false;
    }
}

void PETRTM::setIgnoredPoints(int iRun, bool resetWeights, QString ignoreString)
{ // the ignored string should have fields like "5-10" separated by commas or spaces
    if ( resetWeights )
    {  // set weights using a specific time model
        int nTimeInRun = nTimeInRun;
        setWeightsInRun(iRun);
    }
    iVector includeVolume; includeVolume.fill(false,_weights[iRun].size());
    int error = utilString::decodeSelectionList(ignoreString, includeVolume);
    if ( error == 0 )
    {  // this section turns weight on/off but does not alter non-binary values
        for (int jt=0; jt<includeVolume.size(); jt++)
        {
            if ( includeVolume[jt] ) _weights[iRun][jt]=0.;
        }
        setPrepared(false);
        _ignoreString[iRun] = ignoreString;
    }
}

QString PETRTM::getIgnoredString(int iRun)
{
    QString ignoreString = "";
    if ( iRun < 0 ) return ignoreString;
    int nTimeInRun=getNumberTimePointsInRun(iRun);
    iVector ignore;  ignore.resize(nTimeInRun);
    for (int jTime=0; jTime<nTimeInRun; jTime++)
    {
        int iTimeTotal = getTotalTimeIndex(iRun,jTime);
        if ( getWeight(iTimeTotal) == 0. )
            ignore[jTime] = 1;
        else
            ignore[jTime] = 0;
    }
    ignoreString = utilString::recodeSelectionList(ignore);
    return ignoreString;
}

void PETRTM::setWeightingModel(int whichWeightingModel)
{
    _PETWeightingModel = whichWeightingModel;

    for ( int jRun=0; jRun<_nRuns; jRun++ )
    {
        setWeightsInRun(jRun);
        setIgnoredPoints(jRun,false,_ignoreString[jRun]);
    }
    setPrepared(false);
}

void PETRTM::setWeightsInRun(int iRun)
{
    _weights[iRun].fill(1.,_nTimePerRun[iRun]);
    if ( _PETWeightingModel != Weights_Uniform )
    {
        double tau;
        if ( _PETWeightingModel == Weights_18F || _PETWeightingModel == Weights_18F_Noiseless )
            tau = 110.;     // 18F, minutes
        else // 11C is the default, which will be used in "custom"
            tau = 20.2334;  // 11C, minutes
        tau *= 1.442695;    // convert from half life to exponential time constant
        double trace=0.;
        for (int jt=0; jt<_nTimePerRun[iRun]; jt++)
        {
            double time = getTimeInRun(iRun,jt);
            double dt = _table[iRun][jt][0];
            _weights[iRun][jt] = dt / qExp(time/tau);
            if ( _PETWeightingModel == Weights_11C || _PETWeightingModel == Weights_18F )
                _weights[iRun][jt] *= _tissRegion[iRun][jt];
            if ( _PETWeightingModel == Weights_Custom )
            {
                double physiologicalWeight = 1. - SQR(_frtmConvolutionRelative[iRun][jt]);
                _weights[iRun][jt] *= physiologicalWeight;
            }
            trace += _weights[iRun][jt];
        }
        // Normalize weights: average value should = 1; this must be true to make F statistic valid.
        // E.g., post-hoc determination of _sigma2 assumes accurate separability of _sigma2 and weights = 1/sigma2
        for (int jt=0; jt<_nTimePerRun[iRun]; jt++)
            _weights[iRun][jt] *= static_cast<double>(_nTimePerRun[iRun])/trace;
    } // 11C or 18F
}

bool PETRTM::isValidID(int iRun, int iType, QChar eventID)
{
    qDebug() << "isValidID enter" << iType << eventID;
    bool valid = true;
    valid &= getEventIndex(eventID) >= 0;
    qDebug() << "isValidID 1" << valid;
    if ( iType == Type_R1 && (_modelRTM == RTM_rFRTM2 || _modelRTM == RTM_SRTM2) )
        valid = false;
    else if ( iType == Type_dCrdt && !_dCrdtIncluded )
        valid = false;
    qDebug() << "isValidID 2" << valid;
    // make sure the reference region is defined for basis functions that need it
    if ( iType == Type_R1 || iType == Type_k2 || iType == Type_dCrdt )
    {
        qDebug() << "isValidID size1" << _refRegionRaw.size() << _nRuns << valid;
        valid &= _refRegionRaw.size() == _nRuns;
        qDebug() << "isValidID size2" << _refRegionRaw[iRun].size() << _nTimePerRun[iRun] << valid;
        if ( valid ) valid &= _refRegionRaw[iRun].size() == _nTimePerRun[iRun];
        qDebug() << "isValidID size3" << valid;
        if ( valid )
        {
            bool nonZero = false;
            for (int jt=0; jt<_nTimePerRun[iRun]; jt++)
                nonZero |= _refRegionRaw[iRun][jt] != 0.;
            valid &= nonZero;
            qDebug() << "isValidID nonZero" << nonZero;
        }
//        if ( iType == Type_dCrdt )
//            valid &= _tau4[iRun] != 0.;
    }
    else if ( iType == Type_k2a && (_modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1) )
        valid = false;
    qDebug() << "isValidID exit" << valid;
    return valid;
};

int PETRTM::countEvents()
{
    qDebug() << "PETRTM::countEvents enter" << getNumberCoefficients();
    // Count events
    _basisID.resize(0);
    _basisShape.resize(0);
    _challengeIndex.resize(0);
    // This method of counting events assumes that IDs are independent across type.
    // e.g., an ID cannot be a k2 and k2a event in different runs
    for ( int jRun=0; jRun<_nRuns; jRun++ )
    {
        qDebug() << "event IDs for run" << jRun << "are" << _k2EventID[jRun] << _k2aEventID[jRun] << _R1EventID[jRun];
        qDebug() << "isValid" <<
                    isValidID(jRun, Type_k2,  _k2EventID[jRun])   <<
                    isValidID(jRun, Type_k2a, _k2aEventID[jRun]) <<
                    isValidID(jRun, Type_R1,  _R1EventID[jRun]);

        qDebug() << "append k2?" << ! _basisID.contains(_k2EventID[jRun]) << isValidID(jRun, Type_k2, _k2EventID[jRun]);

        // Count k2, k2a, and R1 events
        if ( ! _basisID.contains(_k2EventID[jRun]) && isValidID(jRun, Type_k2, _k2EventID[jRun]) )
        {
            _basisID.append(_k2EventID[jRun]);
            _basisShape.append(Type_k2);
            _challengeIndex.append(-1);
            qDebug() << "append basisID" << _k2EventID[jRun];
        }
        if ( ! _basisID.contains(_k2aEventID[jRun]) && isValidID(jRun, Type_k2a, _k2aEventID[jRun]) )
        {
            _basisID.append(_k2aEventID[jRun]);
            _basisShape.append(Type_k2a);
            _challengeIndex.append(-1);
        }
        if ( ! _basisID.contains(_R1EventID[jRun]) && isValidID(jRun, Type_R1, _R1EventID[jRun]) &&
             (_modelRTM == RTM_SRTM3 || _modelRTM == RTM_rFRTM3 || _modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1) )
        { // R1 terms should occur in 3-parameter models
            _basisID.append(_R1EventID[jRun]);
            _basisShape.append(Type_R1);
            _challengeIndex.append(-1);
        }
        if ( ! _basisID.contains(_dCrdtEventID[jRun]) && isValidID(jRun, Type_dCrdt, _dCrdtEventID[jRun]) )
        {
            _basisID.append(_dCrdtEventID[jRun]);
            _basisShape.append(Type_dCrdt);
            _challengeIndex.append(-1);
        }
        // Define mapping of IDs onto coefficients
        _R1EventCoefficient[jRun]  = getEventCoefficient(_R1EventID[jRun]);
        _k2EventCoefficient[jRun]  = getEventCoefficient(_k2EventID[jRun]);
        _k2aEventCoefficient[jRun] = getEventCoefficient(_k2aEventID[jRun]);
        _dCrdtEventCoefficient[jRun] = getEventCoefficient(_dCrdtEventID[jRun]);
        qDebug() << "_k2EventCoefficient for run" << jRun << "=" << _k2EventCoefficient[jRun] << _k2EventID[jRun];
    }

    for ( int jChallenge=0; jChallenge<_maxChallenges; jChallenge++ )
    {
        if ( isGoodChallenge(jChallenge) && ! _basisID.contains(_challengeEventID[jChallenge]) )
        {
            _basisID.append(_challengeEventID[jChallenge]);
            _basisShape.append(Type_challenge);
            _challengeIndex.append(jChallenge);
        }
    }
    int nCoeff = _basisID.size();
    qDebug() << "PETRTM::countEvents exit" << nCoeff;
    return nCoeff;
}

void PETRTM::setRTMModelType(RTMModelTypes model)
{
    _modelRTM = model;
    setPrepared(false);

}

void PETRTM::setRTMModelType(QString model)
{
    qDebug() << "PETRTM::setRTMModelType enter" << model;
    if ( model == "SRTM3" )
        _modelRTM = RTM_SRTM3;
    else if ( model == "SRTM2" )
        _modelRTM = RTM_SRTM2;
    else if ( model == "rFRTM3" )
        _modelRTM = RTM_rFRTM3;
    else if ( model == "rFRTM2" || model == "rFRTM2-1" )
        _modelRTM = RTM_rFRTM2;
    else if ( model == "SRTM2_R1" )
        _modelRTM = RTM_SRTM2_R1;
    else if ( model == "rFRTM2_R1" )
        _modelRTM = RTM_rFRTM2_R1;
    else if ( model == "SRTM2x2" )
        _modelRTM = RTM_SRTM2x2;
    else
    {
        qWarning() << "Error: model not defined: " << model;
        exit(1);
    }

//    prepare();
//    updateConditions();
//    if ( getCurrentCondition() >= getNumberConditions() ) setCurrentCondition(0);
    setPrepared(false);
    qDebug() << "PETRTM::setRTMModelType exit" << _modelRTM;
};

void PETRTM::createAllBasisFunctions()
{
    qDebug() << "PETRTM::createAllBasisFunctions enter";
    int nCoeff = countEvents();
    qInfo() << "*************** create PET-RTM basis functions **************" << nCoeff << _nRuns;
    if ( nCoeff == 0 ) return;

    int nTimeTotal = 0;
    for (int jRun=0; jRun<_nRuns; jRun++)
        nTimeTotal += _nTimePerRun[jRun];

    init(nTimeTotal, nCoeff);

    bool calculateConvolution = false;
    if ( (_modelRTM == RTM_rFRTM2_R1) || (_PETWeightingModel == Weights_Custom) )
        calculateConvolution = true;  // calculate the convolution
    if ( isFRTM() || (_modelRTM == RTM_rFRTM2_R1) || (_PETWeightingModel == Weights_Custom) )
        resetAndCalculateFRTMConvolution(calculateConvolution);

    dVector weights;    weights.resize(nTimeTotal);
    if ( _nRuns != _weights.size() ) _weights.resize(_nRuns);
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        if ( _weights[jRun].size() != _nTimePerRun[jRun] ||
             (_PETWeightingModel == Weights_11C || _PETWeightingModel == Weights_18F || _PETWeightingModel == Weights_Custom) )
        { // update weights if the array size is bad or if signal-dependent weights are being used
            setWeightsInRun(jRun);
            setIgnoredPoints(jRun,false,_ignoreString[jRun]);
        }
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
        {
            int iTimeTotal = getTotalTimeIndex(jRun,jt);
            weights[iTimeTotal] = _weights[jRun][jt];
        }
    }
    setWeights(weights);

    dVector basis;
    for (int jCoeff=0; jCoeff<getNumberCoefficients(); jCoeff++)
    {
        basis.fill(0.,nTimeTotal);  // fill and resize
        QChar eventID = _basisID[jCoeff];
        if ( _basisShape[jCoeff] != Type_challenge )
            createRunBasisFunction(true, eventID, basis);
        else
            createChallengeBasisFunction(jCoeff, basis);
        addOrInsertBasisFunction(jCoeff,basis);
    }
    qDebug() << "PETRTM::createAllBasisFunctions exit";
}

void PETRTM::recreateAllFRTMBasisFunctions()
{
    dVector basis;

    int nTimeTotal = getNumberTimePoints();  // this is retrieved from the underlying GLM, but it shouldn't change during "recreate"
    for (int jCoeff=0; jCoeff<getNumberCoefficients(); jCoeff++)
    {
        basis.fill(0.,nTimeTotal);  // fill and resize
        QChar eventID = _basisID[jCoeff];
        if ( _basisShape[jCoeff] == Type_challenge )
            createChallengeBasisFunction(jCoeff, basis);
        else
            createRunBasisFunction(false, eventID, basis);
        addOrInsertBasisFunction(jCoeff,basis);
//        for (int jt=0; jt<getNumberTimePoints(); jt++)
//            qInfo() << "basis[" << eventID << "][" << jt << "] =" << basis[jt];
    }

    calculatePseudoInverse();
    if ( ! (_modelRTM == RTM_SRTM2_R1 || _modelRTM == RTM_rFRTM2_R1) )
        defineConditions(getConditionString());
    setPrepared(true);
}

void PETRTM::createRunBasisFunction(bool newTAC, QChar eventID, dVector &basis)
{
    qDebug() << "createRunBasisFunction" << isFRTM() << isRTM2() << isSRTM();
    dMatrix basisFunction;   basisFunction.resize(_nRuns);

    // For when not new TACs, reset the tissue vector to subtract the dCR/dt term
    // if ( !newTAC ) setTissueVector(newTAC, _tissRegionRaw);

    // Create integrals
    dMatrix tissueIntegral = _tissRegion;
    integrateByRun(tissueIntegral);
    dMatrix frtmConvolutionIntegral;
    if ( isFRTM() )
    {
        frtmConvolutionIntegral = _frtmConvolution; // [_nRuns][_nTimeInRun]; for use with modified basis functions
        integrateByRun(frtmConvolutionIntegral);
    }

    for ( int jRun=0; jRun<_nRuns; jRun++ )
    {
        basisFunction[jRun].fill(0.,_nTimePerRun[jRun]);

        if ( eventID == _R1EventID[jRun] )  // true for SRTM3, rFRTM3, SRTM2_R1, rFRTM2_R1
        {
            for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                basisFunction[jRun][jt] = static_cast<double>(_refRegion[jRun][jt]);
        }

        else if ( eventID == _dCrdtEventID[jRun] )
        {
            for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                basisFunction[jRun][jt] = static_cast<double>(_refRegionDeriv[jRun][jt]);
        }

        else if ( eventID == _k2EventID[jRun] )
        {
            if ( isFRTM() )
            { // true for rFRTM2, rFRTM3, rFRTM2_R1
                for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                    basisFunction[jRun][jt] = _refRegionIntegral[jRun][jt] - frtmConvolutionIntegral[jRun][jt];
            }
            else // true for SRTM3, SRTM2, SRTM2_R1
            {
                for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                    basisFunction[jRun][jt] = static_cast<double>(_refRegionIntegral[jRun][jt]);
            }

            if ( isRTM2_k2a() )
            { // k2 and k2a model
                bool useSRTMTau2Ref = isSRTM();
                double tau2Ref;
                if ( useSRTMTau2Ref )
                    tau2Ref = _tau2RefSRTM[jRun];
                else
                    tau2Ref = _tau2RefFRTM[jRun];
                for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                    basisFunction[jRun][jt] += tau2Ref * _refRegion[jRun][jt];
            }
            else if ( isRTM2_R1() )
            { // k2 and R1 model
                if ( isFRTM() )
                {
                    for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                        basisFunction[jRun][jt] -= 1./(1.+_BPnd[jRun][0]) * (tissueIntegral[jRun][jt] - frtmConvolutionIntegral[jRun][jt]);
                }
                else
                {
                    for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                        basisFunction[jRun][jt] -= 1./(1.+_BPnd[jRun][0]) * tissueIntegral[jRun][jt];
                }
            }
        }

        else if ( eventID == _k2aEventID[jRun] )
        {
            if ( isFRTM() )
            {
                for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                    basisFunction[jRun][jt] = - ( tissueIntegral[jRun][jt] - frtmConvolutionIntegral[jRun][jt] );
//rfrtm-mod                    basisFunction[jRun][jt] = - tissueIntegral[jRun][jt];
            }
            else
            {
                for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                    basisFunction[jRun][jt] = - tissueIntegral[jRun][jt];
            }

        } // event types

    } // jRun

    // convert matrix "basisFunction" to concatenated vector "basis"
    int iStartRun = 0;
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
            basis[jt+iStartRun] = basisFunction[jRun][jt];
        iStartRun += _nTimePerRun[jRun];
    }
}

void PETRTM::createChallengeBasisFunction(int iCoeff, dVector &basis)
{
    int indexChallenge = _challengeIndex[iCoeff];
    dMatrix basisFunction;   basisFunction.resize(_nRuns);
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        basisFunction[jRun].fill(0.,_nTimePerRun[jRun]);
        if ( isGoodChallengeInRun(indexChallenge,jRun) )
        {
            dVector shape;
            createChallengeShape(jRun, indexChallenge, shape);
            for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                basisFunction[jRun][jt] = _tissRegion[jRun][jt];
//rfrtm-mod
            if ( isFRTM() )
            {
                for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                    basisFunction[jRun][jt] -= _frtmConvolution[jRun][jt];  // Ct - dC/dt X E(BPnd,k4)
            }
//rfrtm-mod
            for ( int jt=0; jt<_nTimePerRun[jRun]; jt++)
                basisFunction[jRun][jt] *= - shape[jt];  // - {Ct - dC/dt X E(BPnd,k4)} * shape(t)
        }
    } // jRun
    // Integrate the basis functions by run
    integrateByRun(basisFunction);

    // convert matrix "basisFunction" to concatenated vector "basis"
    int iStartRun = 0;
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
            basis[jt+iStartRun] = basisFunction[jRun][jt];
        iStartRun += _nTimePerRun[jRun];
    }
}

void PETRTM::resetAndCalculateFRTMConvolution(bool calculateConvolution)
{ // passed matrices (run,time) should be pre-filled with zeros.
    qDebug() << "resetAndCalculateFRTMConvolution enter" << calculateConvolution << _tau4[0];
    // Re-compute the convolution integrals
    for (int jRun=0; jRun<_nRuns; jRun++)
        _frtmConvolutionRaw[jRun].fill(0.,_nTimePerRun[jRun]);
    _frtmConvolution = _frtmConvolutionRaw;
    if ( !calculateConvolution ) return;

    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        dVector equilibrationVector = getEquilibrationVector(jRun);
        if ( _tau4[jRun] != 0. )
        {   // tau4 != 0
            /*
            // rfrtm2-mod
            dPoint2D k2a = getk2aInRun(jRun);
            dPoint2D k2 = getk2InRun(jRun);
            // rfrtm2-mod
            */
            double maxConvolution = 0.;
            for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
            {
                double time = getTimeInRun(jRun,jt);
                for ( int jtPrime=0; jtPrime<=jt; jtPrime++ )
                { // integrate from 0 to current time in run (=iTime_run & time)
                    double timePrime = getTimeInRun(jRun,jtPrime);
                    double dt = _table[jRun][jtPrime][0];
                    double argument = equilibrationVector[jtPrime];
                    double exponential = qExp(-argument * (time-timePrime) );
                    _frtmConvolutionRaw[jRun][jt] += dt * exponential * _tissRegionDeriv[jRun][jtPrime];
                } // jtPrime
                if ( _frtmConvolutionRaw[jRun][jt] > maxConvolution ) maxConvolution = _frtmConvolutionRaw[jRun][jt];
                /*
                // rfrtm2-mod
                double mult = 0.;
                if ( k2.x > 0. && k2a.x > 0. )
                    mult = 1. - k2a.x / k2.x;
                _frtmConvolutionRaw[jRun][jt] *= mult;
                // rfrtm2-mod
                */
            } // jt
            if ( maxConvolution != 0. )
            {
                for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
                    _frtmConvolutionRelative[jRun][jt] = _frtmConvolutionRaw[jRun][jt] / maxConvolution;
            }
        } // tau4 != 0
    } // jRun

    _frtmConvolution = _frtmConvolutionRaw;
    if ( _smoothingScale != 0. )
        // Potentially smooth the convolution to reduce noise
        fitLoessCurve(_frtmConvolution);  // additional smoothing xxx
}

void PETRTM::calculateTau2PrimeForSRTM2ForBiasFreeBPND(int iRun)
{
    dMatrix tissueIntegral          = _tissRegion;
    dMatrix frtmConvolutionIntegral = _frtmConvolution; // [_nRuns][_nTimeInRun]; for use with modified basis functions
    integrateByRun(tissueIntegral);
    integrateByRun(frtmConvolutionIntegral);

    int iChallenge = getFirstGoodChallengeIndexInRun(iRun);
    int iStimulus  = getFirstGoodStimulusInRun(iChallenge,iRun);
    int iTimeLastBaseline=0;
    if ( iChallenge < 0 || iStimulus < 0 )
        iTimeLastBaseline = _nTimePerRun[iRun] - 1;
    else
    {
        double timeChallenge = _challengeOn[iChallenge][iStimulus];
        qDebug() << "timeChallenge =" << timeChallenge;
        for (int jt=0; jt<_nTimePerRun[jt]; jt++)
        {
            if ( getTimeInRun(iRun,jt) < timeChallenge )
                iTimeLastBaseline = jt;
            else
                break;
        }
    }
    qDebug() << "iTimeLastBaseline =" << iTimeLastBaseline;
    double rVal = _refRegion[iRun][iTimeLastBaseline];
    double rInt = _refRegionIntegral[iRun][iTimeLastBaseline];
    double tInt = tissueIntegral[iRun][iTimeLastBaseline];
    double cInt = frtmConvolutionIntegral[iRun][iTimeLastBaseline];
    qDebug() << "tInt" << tInt << "cInt" << cInt;
    double part1 = _tau2RefFRTM[iRun] * tInt / (tInt - cInt);
    qDebug() << "rInt" << rInt << "rVal" << rVal << "_tau2RefFRTM" << _tau2RefFRTM[iRun];
    double part2 = (tInt - rInt) / (tInt - cInt) * cInt / rVal;
    qDebug() << "PETRTM::calculateTau2PrimeForSRTM2ForBiasFreeBPND part1 part2" << part1 << part2;
    qDebug() << "PETRTM::calculateTau2PrimeForSRTM2ForBiasFreeBPND frtm srtm2" << _tau2RefFRTM[iRun] << part1 - part2;
    qDebug() << "PETRTM::calculateTau2PrimeForSRTM2ForBiasFreeBPND part2a part2b" << (tInt - rInt) / (tInt - cInt) << cInt / rVal;
}

void PETRTM::setSmoothingScale(double smoothingScale)
{
    _smoothingScale = smoothingScale;
    if ( _smoothingScale != 0. )
    {
        for ( int jRun=0; jRun<_nRuns; jRun++ )
            defineLOESSFitting(jRun);
    }
    updateReferenceRegion();
}

void PETRTM::fitLoessCurve(dMatrix &runData)
{
    dMatrix originalData = runData;

    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
        {
            int iFullWidth = _quadLOESS[jRun][jt].getNumberTimePoints();
            int iHalfWidth = iFullWidth/2;
            dVector localData;  localData.resize(iFullWidth);
            for ( int jRel=-iHalfWidth; jRel<=iHalfWidth; jRel++ )
            {
                int iTime = jt + jRel;
                if ( iTime >= 0 && iTime <_nTimePerRun[jRun] )
                    localData[jRel+iHalfWidth] = originalData[jRun][iTime];
                else
                    localData[jRel+iHalfWidth] = 0.;  // should also have weight=0.
                _quadLOESS[jRun][jt].fitWLS(localData,true);
                runData[jRun][jt] = _quadLOESS[jRun][jt].getFit(iHalfWidth);
            }
        }
    }
}

void PETRTM::smoothRunData(dMatrix &runData)
{
    // PET frames can vary. Use points that contribute more than 10% of a Gaussian weighting; use at least 1 neighbor on each side
    dMatrix originalData = runData;
    double time, timePrime, sum, weight_sum, weight, fwhm;
    double smoothBinRatio=4.;
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
        {
            time = getTimeInRun(jRun,jt);
            sum = originalData[jRun][jt];
            weight_sum = 1.;  // value of Gauss(0)
            fwhm = smoothBinRatio * _table[jRun][jt][0];
            // Negative side
            bool decentWeight=true;
            int iTime = jt-1;
            while ( decentWeight && iTime >= 0 )
            {
                timePrime = getTimeInRun(jRun,iTime);
                weight = Gauss(time-timePrime,fwhm);
                sum += originalData[jRun][iTime] * weight;
                weight_sum += weight;
                decentWeight = weight > 0.1;
                iTime--;
            }
            decentWeight=true;
            iTime = jt+1;
            while ( decentWeight && iTime < _nTimePerRun[jRun] )
            {
                timePrime = getTimeInRun(jRun,iTime);
                weight = Gauss(time-timePrime,fwhm);
                sum += originalData[jRun][iTime] * weight;
                weight_sum += weight;
                decentWeight = weight > 0.1;
                iTime++;
            }
            runData[jRun][jt] = sum / weight_sum;
        }
    }
}

double PETRTM::Gauss(double x, double fwhm)
{
    double sigma = .42466 * fwhm;
    double value = qExp(-x*x/2./sigma/sigma);
    return value;
}

void PETRTM::createChallengeShape(int iRun, int indexChallenge, dVector &shape)
{
    shape.fill(0,_nTimePerRun[iRun]);
    for ( int jStim=0; jStim<_maxStimuli; jStim++)
    {
        if ( isGoodStimulusInRun(indexChallenge,jStim,iRun) )
        {
            for ( int jt=0; jt<_nTimePerRun[iRun]; jt++)
            {
                double time = getTimeInRun(iRun,jt);
                if ( _challengeShape[indexChallenge] == Challenge_Constant )
                    shape[jt] += 1.;
                else if ( _challengeShape[indexChallenge] == Challenge_Square )
                {
                    if ( time >= _challengeOn[indexChallenge][jStim] && time < _challengeOff[indexChallenge][jStim] )
                        shape[jt] += 1.;
                }
                else if ( _challengeShape[indexChallenge] == Challenge_RampUp )
                {
                    double duration = _challengeOff[indexChallenge][jStim] - _challengeOn[indexChallenge][jStim];
                    if ( time >= _challengeOn[indexChallenge][jStim] && time < _challengeOff[indexChallenge][jStim] )
                        shape[jt] += ( time-_challengeOn[indexChallenge][jStim] )/duration;
                }
                else if ( _challengeShape[indexChallenge] == Challenge_RampDown )
                {
                    double duration = _challengeOff[indexChallenge][jStim] - _challengeOn[indexChallenge][jStim];
                    if ( time >= _challengeOn[indexChallenge][jStim] && time < _challengeOff[indexChallenge][jStim] )
                        shape[jt] += 1. - ( time-_challengeOn[indexChallenge][jStim] )/duration;
                }
                else if ( _challengeShape[indexChallenge] == Challenge_Gamma )
                {
                    double time0 = _challengeOn[indexChallenge][jStim];
                    double tau = _tau[indexChallenge];
                    // double alpha = _alpha[indexChallenge];
                    if ( time >= time0 )
                        shape[jt] += (time-time0)/tau * qExp(1.-(time-time0)/tau);
                }
                else if ( _challengeShape[indexChallenge] == Challenge_Sigmoid )
                {
                    double time0 = _challengeOn[indexChallenge][jStim];
                    double tau = _tau[indexChallenge];
                    // double alpha = _alpha[indexChallenge];
                    if ( time >= time0 )
                        shape[jt] += (time-time0)/tau / qSqrt((1.+ (time-time0)/tau*(time-time0)/tau  ));
                }
            } // jt
        } // isGoodStimulusInRun
    } // jStim
}
void PETRTM::integrateByRun(dMatrix &runMatrix )  // [nRuns][nTimePerRun]
{
    qDebug() << "PETRTM::integrateByRun enter" << _table.size() << _nRuns;
    dMatrix copiedMatrix = runMatrix;
    if ( _table.size() != _nRuns )
    { // this is just to provide a "default" option; bins sizes may differ in time
        _table.resize(_nRuns);
        for ( int jRun=0; jRun<_nRuns; jRun++ )
        {
            qDebug() << "PETRTM::integrateByRun 0" << _nTimePerRun[jRun];
            _table[jRun].resize(_nTimePerRun[jRun]);
            for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
                _table[jRun][jt].append(1.);  // initialize time bin widths (1st column=frames) to 1.
        }
    }
    qDebug() << "PETRTM::integrateByRun 1";
    // Find the integral from t'=0 to t'=t for each point. Reset at the start of each run.
    qDebug() << "PETRTM::integrateByRun 2" << runMatrix.size();
    for ( int jRun=0; jRun<runMatrix.size(); jRun++ )
    {
        qDebug() << "PETRTM::integrateByRun 3" << runMatrix[jRun].size();
        for ( int jt=0; jt<runMatrix[jRun].size(); jt++ )
        {
            double integral = 0.;
            for ( int jtPrime=0; jtPrime<=jt; jtPrime++ )
            {
                double dt = _table[jRun][jtPrime][0];
                integral += copiedMatrix[jRun][jtPrime] * dt;
            }
            runMatrix[jRun][jt] = integral;
        }
    }
    qDebug() << "PETRTM::integrateByRun exit";
}

void PETRTM::differentiateByRun(dMatrix &runMatrix )  // [nRuns][nTimePerRun]
{
    qDebug() << "PETRTM::differentiateByRun enter";
    dMatrix copiedMatrix = runMatrix;
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        for (int jt=0; jt<_nTimePerRun[jRun]; jt++)
        {
            double dt = _table[jRun][jt][0];
            if ( jt != 0. )
                runMatrix[jRun][jt] = (copiedMatrix[jRun][jt] - copiedMatrix[jRun][jt-1]) / dt;
            else
                runMatrix[jRun][jt] = copiedMatrix[jRun][jt] / dt;
        }
    }
    qDebug() << "PETRTM::differentiateByRun exit";
}

void PETRTM::averageStimuli(dMatrix yData, dMatrix yFit)
{
    bool useRatio = true;
    if ( _challengeForStimAv < 0 || _challengeForStimAv >= _maxChallenges) return;
    int length = _nPreForChallAv + _nPostForChallAv + 1;
    _xForChallAv.fill(0.,length);
    _yForChallAv.fill(0.,length);
    _yFitForChallAv.fill(0.,length);
    _ySEMForChallAv.fill(0.,length);
    // find the length of the vector
    int nStimuli = 0;
    for ( int jStim=0; jStim<_maxStimuli; jStim++)
    {
        if ( isGoodStimulus(_challengeForStimAv,jStim) )
        {
            int iRun     = _challengeRun[_challengeForStimAv][jStim];
            if ( jStim != getFirstGoodStimulusInRun(_challengeForStimAv,iRun) ) continue;
            double onSet = _challengeOn[_challengeForStimAv][jStim];
            int iStart = 0;
            for ( int jt=0; jt<_nTimePerRun[iRun]; jt++)
            {
                double time = getTimeInRun(iRun,jt);
                if ( time >= onSet )
                {
                    iStart = jt;
                    break;
                }
            } // jt
            double time0 = getTimeInRun(iRun,iStart);
            int iTime = iStart - _nPreForChallAv;
            for ( int jt=0; jt<_yForChallAv.size(); jt++, iTime++)
            {
                if ( iTime >= 0 && iTime < yData[iRun].size() )
                {
                    double time = getTimeInRun(iRun,iTime);
                    _xForChallAv[jt]    += (time - time0);
                    if ( useRatio )
                    {
                        double value = yData[iRun][iTime] / _refRegion[iRun][iTime] - 1.;
                        _yForChallAv[jt]    += value;
                        double fit   = yFit[iRun][iTime]  / _refRegion[iRun][iTime] - 1.;
                        _yFitForChallAv[jt] += fit;
                    }
                    else
                    {
                        _yForChallAv[jt]    += yData[iRun][iTime];
                        _yFitForChallAv[jt] += yFit[iRun][iTime];
                    }
                }
            }
            nStimuli++;
        } // isGoodStimulus
    }
    for ( int jt=0; jt<_yForChallAv.size(); jt++)
    {
        _xForChallAv[jt]    /= static_cast<double>(nStimuli);
        _yForChallAv[jt]    /= static_cast<double>(nStimuli);
        _yFitForChallAv[jt] /= static_cast<double>(nStimuli);
    }

    // Now compute SEM
    if ( nStimuli <= 1 )
    {
        for ( int jt=0; jt<_yForChallAv.size(); jt++)
            _ySEMForChallAv[jt] = 0.;
    }
    else
    {
        // find the length of the vector
        for ( int jStim=0; jStim<_maxStimuli; jStim++)
        {
            if ( isGoodStimulus(_challengeForStimAv,jStim) )
            {
                int iRun     = _challengeRun[_challengeForStimAv][jStim];
                if ( jStim != getFirstGoodStimulusInRun(_challengeForStimAv,iRun) ) continue;
                double onSet = _challengeOn[_challengeForStimAv][jStim];
                int iStart = 0;
                for ( int jt=0; jt<_nTimePerRun[iRun]; jt++)
                {
                    double time = getTimeInRun(iRun,jt);
                    if ( time >= onSet )
                    {
                        iStart = jt;
                        break;
                    }
                } // jt
                int iTime = iStart - _nPreForChallAv;
                for ( int jt=0; jt<_yForChallAv.size(); jt++, iTime++)
                {
                    if ( iTime >= 0 && iTime < yData[iRun].size() )
                    {
                        if ( useRatio )
                        {
                            double value = yData[iRun][iTime] / _refRegion[iRun][iTime] - 1.;
                            _ySEMForChallAv[jt] += SQR(value - _yForChallAv[jt]);
                        }
                        else
                            _ySEMForChallAv[jt] += SQR(yData[iRun][iTime] - _yForChallAv[jt]);
                    }
                }
            } // isGoodStimulus
        }
        // compute stdev and then sem
        for ( int jt=0; jt<_yForChallAv.size(); jt++)
        {
            _ySEMForChallAv[jt] /= static_cast<double>(nStimuli-1);
            _ySEMForChallAv[jt] = qSqrt(_ySEMForChallAv[jt]);       // stdev
            _ySEMForChallAv[jt] /= static_cast<double>(nStimuli);  // sem
        }
    }
}

void PETRTM::getStimAveragingVectors(dVector &xForChallAv, dVector &yForChallAv, dVector &ySEMForChallAv, dVector &yFitForChallAv)
{
    int length = _yForChallAv.size();
    xForChallAv.resize(length); yForChallAv.resize(length); ySEMForChallAv.resize(length); yFitForChallAv.resize(length);
    xForChallAv = _xForChallAv; yForChallAv = _yForChallAv; ySEMForChallAv = _ySEMForChallAv; yFitForChallAv = _yFitForChallAv;
}
