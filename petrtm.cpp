#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QMessageBox>
#include <qmath.h>
#include <QTextStream>
#include <QDebug>
#include <QMutex>
#include "io.h"
#include "generalglm.h"
#include "petrtm.h"

using namespace utilIO;

int PETRTM::readTimeModel(QString fileName)
{
    // Read the time model file
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return(1);
    QTextStream in_stream(&file);

    QString allConditions = "";
    int iLine=0;
    while (!in_stream.atEnd())
    {
        iLine++;
        QString line = in_stream.readLine();
        QString unCommented = line.left(line.indexOf("#"));
        if ( unCommented.isEmpty() )
            continue;
        QStringList stringList = unCommented.split(QRegExp("[ ]"), QString::SkipEmptyParts);
        QString keyWord = stringList.at(0);
        if ( !keyWord.compare("time-model",Qt::CaseInsensitive) )
        { // set generic model name (SRTM or rFRTM); 2 or 3 parameter variants set later based upon 1st GLM file
            QString modelName = stringList.at(1);
            if ( !modelName.compare("FRTM",Qt::CaseInsensitive)  || !modelName.compare("FRTM3",Qt::CaseInsensitive) ||
                 !modelName.compare("rFRTM",Qt::CaseInsensitive) || !modelName.compare("rFRTM3",Qt::CaseInsensitive) ||
                 !modelName.compare("FRTM2",Qt::CaseInsensitive) || !modelName.compare("rFRTM2",Qt::CaseInsensitive) )
            {
                _modelRTM = RTM_rFRTM3;
            }
            else if ( !modelName.compare("SRTM",Qt::CaseInsensitive) || !modelName.compare("SRTM3",Qt::CaseInsensitive) ||
                      !modelName.compare("SRTM2",Qt::CaseInsensitive) )
            {
                _modelRTM = RTM_SRTM3;
            }
            else
            {
                qInfo() << "Error reading the model name" << modelName;
                return(1);
            }
            if ( isFRTM() )
            { // need to read 1/k4
                if ( stringList.size() > 2 )
                {
                    QString tau4 = stringList.at(2);
                    bool ok;
                    double value;
                    value = tau4.toDouble(&ok);
                    if ( ok )
                    {
                        _tau4Default = value;
                        for ( int jRun=0; jRun<_nRuns; jRun++ )
                            if ( _tau4[jRun] == 0. ) _tau4[jRun] = _tau4Default;
                    }
                }
            }
        }
        else if ( !keyWord.compare("reference-region",Qt::CaseInsensitive) )
        {
            // The reference region might be
            // 1) an explicit file name (should end with .ovl), or
            // 2) a name to be matched in the list of overlays
            _refRegionName = utilString::replaceEnvironmentVariables(stringList.at(1));
        }
        else if ( !keyWord.compare("brain-region",Qt::CaseInsensitive) )
        {
            // The reference region might be
            // 1) an explicit file name (should end with .ovl), or
            // 2) a name to be matched in the list of overlays
            _brainRegionName = utilString::replaceEnvironmentVariables(stringList.at(1));
        }
        else if ( !keyWord.compare("smoothing-scale",Qt::CaseInsensitive) )
        {
            QString smoothingString = stringList.at(1);
            bool ok;
            _smoothingScaleMin = smoothingString.toDouble(&ok);
            if ( !ok )
            {
                qInfo() << "Error reading the value for the smoothing string:" << smoothingString;
                return(1);
            }
        }
        else if ( !keyWord.compare("weights",Qt::CaseInsensitive) )
        {
            QString modelName = stringList.at(1);
            if ( !modelName.compare("11C",Qt::CaseInsensitive)  )
                _PETWeightingModel = Weights_11C;
            else if ( !modelName.compare("18F",Qt::CaseInsensitive)  )
                _PETWeightingModel = Weights_18F;
            else if ( !modelName.compare("11C-noiseless",Qt::CaseInsensitive)  )
                _PETWeightingModel = Weights_11C_Noiseless;
            else if ( !modelName.compare("18F-noiseless",Qt::CaseInsensitive)  )
                _PETWeightingModel = Weights_18F_Noiseless;
            else if ( !modelName.compare("Custom",Qt::CaseInsensitive)  )
                _PETWeightingModel = Weights_noUptake;
            else
                _PETWeightingModel = Weights_Uniform;
        }
        else if ( !keyWord.compare("conditions",Qt::CaseInsensitive) || !keyWord.compare("condition",Qt::CaseInsensitive) )
        {
            for ( int jString=1; jString<stringList.size(); jString++ )
                _conditionList.append(stringList.at(jString));
            allConditions = getConditionString();
        }
        else if ( !keyWord.compare("runs",Qt::CaseInsensitive)  || !keyWord.compare("runs:",Qt::CaseInsensitive) ||
                  !keyWord.compare("scans",Qt::CaseInsensitive) || !keyWord.compare("scans:",Qt::CaseInsensitive) ||
                  !keyWord.compare("files",Qt::CaseInsensitive) || !keyWord.compare("files:",Qt::CaseInsensitive) )
            break;
    } // while !end

    // Now read GLM files and table files
    int iRun=0;
    while (!in_stream.atEnd())
    {
        iLine++;
        QString line = in_stream.readLine();
        QString unCommented = line.left(line.indexOf("#"));
        if ( unCommented.isEmpty() )
            continue;
        QStringList stringList = unCommented.split(QRegExp("[ ]"), QString::SkipEmptyParts);
        QString dataFile = stringList.at(0);
        QString glmFile = stringList.at(1);
        if ( stringList.count() > 2 )
        {
            QString tableFile = stringList.at(2);
            if ( readTimeBinsFile(iRun,tableFile) ) return(1);
        }
        else
        {
            qInfo() << "Error: each file requires a table with time frames (seconds) in the first column";
            qInfo() << "[TAC file name] [GLM file name] [table file name]";
            return(1);
        }
        readGLMFile(iRun,glmFile);
        iRun++;
    }
    file.close();

    // Now that all runs have been read, define this as 2 or 3-parameter RTM
    bool RTM2 = true; // all runs are RTM2 if ends up as true
    bool RTM3 = true; // all runs are RTM3 if ends up as true
    if ( isSRTM() )
    {
        for (int jRun=0; jRun<_nRuns; jRun++)
        {
            RTM2 &= _tau2RefSRTMFixed[jRun] != 0.;
            RTM3 &= _tau2RefSRTMFixed[jRun] == 0.;
        }
    }
    else
    {
        for (int jRun=0; jRun<_nRuns; jRun++)
        {
            RTM2 &= _tau2RefFRTMFixed[jRun] != 0.;
            RTM3 &= _tau2RefFRTMFixed[jRun] == 0.;
        }
    }
    if ( RTM2 && !RTM3 )
    {
        if ( isSRTM() )
        {
            if ( _modelRTM != RTM_SRTM2Fit )
                _modelRTM = RTM_SRTM2;
        }
        else
            _modelRTM = RTM_rFRTM2;
    }
    else if ( RTM3 && !RTM2 )
    { // variable _tau2Ref
        if ( isSRTM() )
            _modelRTM = RTM_SRTM3;
        else
            _modelRTM = RTM_rFRTM3;
    }
    else
    {
        QMessageBox msgBox;  msgBox.setIcon(QMessageBox::Critical);
        QString errorText = QString("This analysis mixed 2-parameter and 3-parameter models!");
        msgBox.setText(errorText);
        msgBox.exec();
        return(1);
    }

    // define conditions at the end, when everything else is established
    prepare();  // this is needed to avoid killing conditions due to undefined basis functions
    _petConditionsFromTimeModelFile = allConditions;
    definePETConditions(allConditions);

    return(0);
}

int PETRTM::readGLMFile(int iRun, QString fileName)
{
    // Delete all existing challenges in this run by setting an invalid run number
    for ( int jChallenge=0; jChallenge<_maxChallenges; jChallenge++ )
    {
        for ( int jStim=0; jStim<_maxStimuli; jStim++ )
        {
            if ( isGoodStimulusInRun(jChallenge, jStim, iRun) )
                _challengeRun[jChallenge][jStim] = -1;
        }
    }

    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qInfo() << "Error opening file " << fileName;
        return(1);
    }
    QTextStream in_stream(&file);

    QString line = in_stream.readLine();
    QString unCommented = line.mid(0, line.indexOf("#"));
    QStringList stringList = unCommented.split(QRegExp("[,\\s]"), QString::SkipEmptyParts);
    int nStrings = stringList.size();
    bool ok = true;

    // In the old format, the first line should be either a fixed value of 1/k2' or 0, with the latter indicating 3-parameter RTM
    if ( nStrings == 0 )
        ok = false;
    else
    {
        QString tau2Ref = stringList.at(0);
        tau2Ref.toDouble(&ok);
    }
    file.close();
    if ( ok )
        return readGLMFileOldFormat(iRun, fileName);
    else
        return readGLMFileNewFormat(iRun, fileName);
}

int PETRTM::readGLMFileOldFormat(int iRun, QString fileName)
{
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qInfo() << "Error opening file " << fileName;
        return(1);
    }
    QTextStream in_stream(&file);

    int iLine=1;
    QString line = in_stream.readLine();
    QString unCommented = line.mid(0, line.indexOf("#"));
    QStringList stringList = unCommented.split(QRegExp("[,\\s]"), QString::SkipEmptyParts);
    int nStrings = stringList.size();
    bool ok = true;

    // 1/k2'
    if ( nStrings == 0 )
        ok = false;
    else
    {
        QString tau2Ref = stringList.at(0);
        if ( isSRTM() )
            _tau2RefSRTMFixed[iRun] = tau2Ref.toDouble(&ok);
        else
            _tau2RefFRTMFixed[iRun] = tau2Ref.toDouble(&ok);
    }
    if ( !ok )
    {
        QString errorText = QString("Error reading the value of 1/k2' for run %1").arg(iRun);
        utilString::errorMessage(errorText);
    }

    setIgnoredPoints(iRun,true,"");
    bool lookForNewEvent = true;
    bool lookForIgnore = false;
    bool firstIgnored = true;
    int iStimulus=0;
    // Start looking for events
    int indexChallenge=0;
    while (!in_stream.atEnd())
    {
        line = in_stream.readLine();  iLine++;
        unCommented = line.left(line.indexOf("#"));
        if ( unCommented.isEmpty() )
        {
            lookForNewEvent = true;
            continue;
        }
        QStringList valueList = unCommented.split(QRegExp("[,\\s]"), QString::SkipEmptyParts);
        int nValues = valueList.size();

        if ( lookForNewEvent )
        { // this should be a event ID and shape
            lookForNewEvent = false;  // code below is for new event, but then it should be false
            QString eventIDString = valueList.at(0);
            int ignore = eventIDString.toInt();
            if ( eventIDString == "ignore" || eventIDString == "Ignore" || ignore == -1 )
            {
                lookForNewEvent = false;
                lookForIgnore = true;
                continue;
            }
            if ( nValues < 2 )
            {
                qInfo() << "Error; there should be at least 2 values on line" << iLine << ": ID and type (k2, k2a, R1, challenge)";
                return(1);
            }
            QChar eventID = eventIDString.at(0);
            QString eventType = valueList.at(1);  // 2nd string is the type: k2, k2a, R1, challenge
            // Set run-specific R1, k2, k2a
            if ( ! eventType.compare("R1",Qt::CaseInsensitive) )
                _R1EventID[iRun] = eventID;
            else if ( ! eventType.compare("k2",Qt::CaseInsensitive) )
                _k2EventID[iRun] = eventID;
            else if ( ! eventType.compare("k2a",Qt::CaseInsensitive) )
                _k2aEventID[iRun] = eventID;
            else if ( ! eventType.compare("dCr/dt",Qt::CaseInsensitive) )
            {
                _dCrdtEventID[iRun] = eventID;
                setInclusionOfdCrdt(true);
            }
            else if ( ! eventType.compare("challenge",Qt::CaseInsensitive) ||
                      ! eventType.compare("dk2a",Qt::CaseInsensitive)      ||
                      ! eventType.compare("k2c",Qt::CaseInsensitive)       ||
                      ! eventType.compare("delta_k2a",Qt::CaseInsensitive) ||
                      ! eventType.compare("delta-k2a",Qt::CaseInsensitive) ||
                      ! eventType.compare("deltak2a",Qt::CaseInsensitive) )
            {
                indexChallenge = getEventIndex(eventID);
                if ( indexChallenge < 0 )
                {
                    QString errorText = QString("Error; event %1 in file %2 is not allowed").arg(eventID).arg(fileName);
                    QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                    return(1);
                }
                if ( indexChallenge < 0 ) continue;
                iStimulus = getFirstBadStimulus(indexChallenge);
                _challengeShape[indexChallenge] = Challenge_Square; // default
                if ( nValues > 2 )
                { // read the challenge shape
                    QString challengeShape = valueList.at(2);
                    if ( ! challengeShape.compare("constant",Qt::CaseInsensitive) )
                    {
                        _challengeShape[indexChallenge] = Challenge_Constant;
                        lookForNewEvent = true;  // constant events need no extra information
                    }
                    else if ( ! challengeShape.compare("square",Qt::CaseInsensitive) )
                        _challengeShape[indexChallenge] = Challenge_Square;
                    else if ( ! challengeShape.compare("ramp-up",Qt::CaseInsensitive) )
                        _challengeShape[indexChallenge] = Challenge_RampUp;
                    else if ( ! challengeShape.compare("ramp-down",Qt::CaseInsensitive) )
                        _challengeShape[indexChallenge] = Challenge_RampDown;
                    else if ( ! challengeShape.compare("sigmoid",Qt::CaseInsensitive) )
                    {
                        _challengeShape[indexChallenge] = Challenge_Sigmoid;
                        if ( nValues < 4 ) // [ID] sigmoid [tau]
                        {
                            QString errorText = QString("Error; specify the sigmoidal 'tau' value on line %1 of file %2").arg(iLine).arg(fileName);
                            QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                            return(1);
                        }
                        QString tau = valueList.at(3);
                        bool ok;
                        _challengeTau[indexChallenge] =tau.toDouble(&ok);
                        if ( !ok )
                        {
                            QString errorText = QString("Error reading the value of tau for the sigmoid challenge on line %1 of file %2").arg(iLine).arg(fileName);
                            QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                            return(1);
                        }
                    }
                    else if ( ! challengeShape.compare("gamma",Qt::CaseInsensitive) )
                    {
                        _challengeShape[indexChallenge] = Challenge_Gamma;
                        if ( nValues < 4 ) // [ID] gamma [tau]
                        {
                            QString errorText = QString("Error; specify the gamma 'tau' value on line %1 of file %2").arg(iLine).arg(fileName);
                            QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                            return(1);
                        }
                        QString tau = valueList.at(3);
                        bool ok;
                        _challengeTau[indexChallenge] =tau.toDouble(&ok);
                        if ( !ok )
                        {
                            QString errorText = QString("Error reading the value of tau for the gamma challenge on line %1 of file %2").arg(iLine).arg(fileName);
                            QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                            return(1);
                        }
                        if ( nValues >= 4 ) // [ID] gamma [tau] [alpha]
                        {
                            QString alpha = valueList.at(3);
                            bool ok;
                            _challengeAlpha[indexChallenge] =alpha.toDouble(&ok);
                            if ( !ok )
                            {
                                QString errorText = QString("Error reading the value of alpha for the gamma challenge on line %1").arg(iLine);
                                QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                                return(1);
                            }
                        }
                    } // shapes
                } // nValues > 2
            } // challenge
            else
            {
                QString errorText = QString("Unrecognized keyword %1 on line %2").arg(eventType).arg(iLine);
                QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                return(1);
            }
        } // lookForNewEvent
        else if ( !lookForIgnore )
        { // this should be a repeated instance of the last event
            // All events should have a start time
            QString eventOn = valueList.at(0);
            _challengeOn[indexChallenge][iStimulus] = eventOn.toDouble(&ok);
            _challengeRun[indexChallenge][iStimulus] = iRun;
            if ( !ok )
            {
                QString errorText = QString("Error reading the ON time for a stimulus from file %1").arg(fileName);
                QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                return(1);
            }
            // Square events must have a 2nd time
            double dNext;
            bool needOnAndOff = _challengeShape[indexChallenge] == Challenge_Square
                    ||          _challengeShape[indexChallenge] == Challenge_RampUp
                    ||          _challengeShape[indexChallenge] == Challenge_RampDown;
            if ( needOnAndOff && nValues < 2 )
            {
                QString errorText = "On and off times both are needed for events of type square or ramp";
                QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                return(1);
            }
            else if ( nValues >= 2 )
            {  // read 2nd value (off or mag)
                QString nextWord = valueList.at(1);
                dNext = nextWord.toDouble(&ok);
                if ( !ok )
                {
                    QString errorText;
                    bool needOn = (_challengeShape[indexChallenge] == Challenge_Gamma)
                            ||    (_challengeShape[indexChallenge] == Challenge_Sigmoid);
                    if ( needOnAndOff )
                        errorText = QString("Error reading the OFF time for event %1 from file %2").arg(getEventChar(indexChallenge)).arg(fileName);
                    else if ( needOn )
                        errorText = QString("Error reading the magnitude from file %1").arg(fileName);
                    QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                    return(1);
                } // !ok
                if ( needOnAndOff )
                    _challengeOff[indexChallenge][iStimulus] = dNext;
            }
            iStimulus = getFirstBadStimulus(indexChallenge);
        }  // stimulus
        else
        { // ignore
            setIgnoredPoints(iRun,firstIgnored,unCommented);
            firstIgnored = false;
        }
    }
    file.close();
    setPrepared(false);
    return(0);
}

int PETRTM::readGLMFileNewFormat(int iRun, QString fileName)
{
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qInfo() << "Error opening file " << fileName;
        return(1);
    }
    QTextStream in_stream(&file);

    setIgnoredPoints(iRun,true,"");
    bool lookForNewEvent = true;
    int iStimulus=0;
    // Start looking for events
    int indexChallenge=0;
    int iLine=0;
    while (!in_stream.atEnd())
    {
        QString line = in_stream.readLine();  iLine++;  // Without the subtraction above, this ++ would put iLine 1 ahead of error messages
        QString unCommented = line.left(line.indexOf("#"));
        if ( unCommented.isEmpty() )
        {
            lookForNewEvent = true;
            continue;
        }
        QStringList valueList = unCommented.split(QRegExp("[,\\s]"), QString::SkipEmptyParts);
        int nValues = valueList.size();

        if ( lookForNewEvent )
        { // this should be a event ID and shape
            if ( nValues < 2 )
            {
                qInfo() << "Error; there should be at least 2 values on line" << iLine << ": ID and type (k2, k2a, R1, challenge)";
                return(1);
            }
            lookForNewEvent = false;  // code below is for new event, but then it should be false

            QString eventType = valueList.at(0);     // 1st string is the type: k2, k2a, R1, challenge, 1/k2', 1/k4
            QString eventIDString = valueList.at(1);  QChar eventID = eventIDString.at(0);
            int ignore = eventIDString.toInt();
            if ( eventType == "ignore" || eventType == "Ignore" || ignore == -1 )
                readGLMIgnoreBlock(&in_stream, iRun, eventIDString);
            // Set run-specific R1, k2, k2a
            else if ( ! eventType.compare("1/k2'",Qt::CaseInsensitive) ||
                      ! eventType.compare("1/k2_ref",Qt::CaseInsensitive) ||
                      ! eventType.compare("tau2'",Qt::CaseInsensitive) ||
                      ! eventType.compare("tau2_ref",Qt::CaseInsensitive) )
            {  // Set a fixed value of tau2_ref
                QString tau2Ref = valueList.at(1);
                bool ok;
                if ( isSRTM() )
                {
                    if ( nValues == 6 )
                    {
                        _modelRTM = RTM_SRTM2Fit;
                        bool ok1, ok2, ok3, ok4, ok5;
                        _tau2RefSRTMFixed[iRun]    = tau2Ref.toDouble(&ok1);
                        _tau2RefSRTMCalOffset = valueList.at(2).toDouble(&ok2);
                        _tau2RefSRTMCal[iRun][0] = valueList.at(3).toDouble(&ok3);
                        _tau2RefSRTMCal[iRun][1] = valueList.at(4).toDouble(&ok4);
                        _tau2RefSRTMCal[iRun][2] = valueList.at(5).toDouble(&ok5);
                        ok = ok1 && ok2 && ok3 && ok4 && ok5;
                    }
                    else
                        _tau2RefSRTMFixed[iRun] = tau2Ref.toDouble(&ok);
                }
                else
                    _tau2RefFRTMFixed[iRun] = tau2Ref.toDouble(&ok);
                if ( !ok )
                {
                    QString errorText = QString("Error reading the value of 1/k2' on line %1 of file %2").arg(iLine).arg(fileName);
                    utilString::errorMessage(errorText);
                }
            }
            else if ( ! eventType.compare("1/k4",Qt::CaseInsensitive) ||
                      ! eventType.compare("tau4",Qt::CaseInsensitive) )
            {  // Set a fixed value of tau4
                QString tau4String = valueList.at(1);
                bool ok;
                double tau4 = tau4String.toDouble(&ok);
                if ( !ok )
                {
                    QString errorText = QString("Error reading the value of 1/k2' on line %1 of file %2").arg(iLine).arg(fileName);
                    utilString::errorMessage(errorText);
                }
                else
                    _tau4[iRun] = _tau4Default = tau4;
            }
            else if ( ! eventType.compare("R1",Qt::CaseInsensitive) )
            {
                _R1EventID[iRun] = eventID;
                _tau2RefSRTMFixed[iRun] = 0.;  // set to 0 to flag SRTM3
            }
            else if ( ! eventType.compare("k2",Qt::CaseInsensitive) )
                _k2EventID[iRun] = eventID;
            else if ( ! eventType.compare("k2a",Qt::CaseInsensitive) )
                _k2aEventID[iRun] = eventID;
            else if ( ! eventType.compare("dCr/dt",Qt::CaseInsensitive) )
            {
                _dCrdtEventID[iRun] = eventID;
                setInclusionOfdCrdt(true);
            }
            else if ( ! eventType.compare("challenge",Qt::CaseInsensitive) ||
                      ! eventType.compare("dk2a",Qt::CaseInsensitive)      ||
                      ! eventType.compare("k2c",Qt::CaseInsensitive)       ||
                      ! eventType.compare("delta_k2a",Qt::CaseInsensitive) ||
                      ! eventType.compare("delta-k2a",Qt::CaseInsensitive) ||
                      ! eventType.compare("deltak2a",Qt::CaseInsensitive) )
            {
                indexChallenge = getEventIndex(eventID);
                if ( indexChallenge < 0 )
                {
                    QString errorText = QString("Error; event %1 in file %2 is not allowed").arg(eventID).arg(fileName);
                    QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                    return(1);
                }
                if ( indexChallenge < 0 ) continue;
                iStimulus = getFirstBadStimulus(indexChallenge);
                _challengeShape[indexChallenge] = Challenge_Square; // default
                if ( nValues > 2 )
                { // read the challenge shape
                    QString challengeShape = valueList.at(2);
                    if ( ! challengeShape.compare("constant",Qt::CaseInsensitive) )
                    {
                        _challengeShape[indexChallenge] = Challenge_Constant;
                        lookForNewEvent = true;  // constant events need no extra information
                    }
                    else if ( ! challengeShape.compare("square",Qt::CaseInsensitive) )
                        _challengeShape[indexChallenge] = Challenge_Square;
                    else if ( ! challengeShape.compare("ramp-up",Qt::CaseInsensitive) )
                        _challengeShape[indexChallenge] = Challenge_RampUp;
                    else if ( ! challengeShape.compare("ramp-down",Qt::CaseInsensitive) )
                        _challengeShape[indexChallenge] = Challenge_RampDown;
                    else if ( ! challengeShape.compare("sigmoid",Qt::CaseInsensitive) )
                    {
                        _challengeShape[indexChallenge] = Challenge_Sigmoid;
                        if ( nValues < 4 ) // [ID] sigmoid [tau]
                        {
                            QString errorText = QString("Error; specify the sigmoidal 'tau' value on line %1 of file %2").arg(iLine).arg(fileName);
                            QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                            return(1);
                        }
                        QString tau = valueList.at(3);
                        bool ok;
                        _challengeTau[indexChallenge] =tau.toDouble(&ok);
                        if ( !ok )
                        {
                            QString errorText = QString("Error reading the value of tau for the sigmoid challenge on line %1 of file %2").arg(iLine).arg(fileName);
                            QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                            return(1);
                        }
                    }
                    else if ( ! challengeShape.compare("gamma",Qt::CaseInsensitive) )
                    {
                        _challengeShape[indexChallenge] = Challenge_Gamma;
                        if ( nValues < 4 ) // [ID] gamma [tau]
                        {
                            QString errorText = QString("Error; specify the gamma 'tau' value on line %1 of file %2").arg(iLine).arg(fileName);
                            QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                            return(1);
                        }
                        QString tau = valueList.at(3);
                        bool ok;
                        _challengeTau[indexChallenge] =tau.toDouble(&ok);
                        if ( !ok )
                        {
                            QString errorText = QString("Error reading the value of tau for the gamma challenge on line %1 of file %2").arg(iLine).arg(fileName);
                            QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                            return(1);
                        }
                        if ( nValues >= 4 ) // [ID] gamma [tau] [alpha]
                        {
                            QString alpha = valueList.at(3);
                            bool ok;
                            _challengeAlpha[indexChallenge] =alpha.toDouble(&ok);
                            if ( !ok )
                            {
                                QString errorText = QString("Error reading the value of alpha for the gamma challenge on line %1").arg(iLine);
                                QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                                return(1);
                            }
                        }
                    } // shapes
                } // nValues > 2
            } // challenge
        } // lookForNewEvent
        else
        { // this should be a repeated instance of the last event
            // All events should have a start time
            QString eventOn = valueList.at(0);  bool ok;
            _challengeOn[indexChallenge][iStimulus] = eventOn.toDouble(&ok);
            _challengeRun[indexChallenge][iStimulus] = iRun;
            if ( !ok )
            {
                QString errorText = QString("Error reading the ON time for a stimulus from file %1").arg(fileName);
                QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                return(1);
            }
            // Square events must have a 2nd time
            double dNext;
            bool needOnAndOff = _challengeShape[indexChallenge] == Challenge_Square
                    ||          _challengeShape[indexChallenge] == Challenge_RampUp
                    ||          _challengeShape[indexChallenge] == Challenge_RampDown;
            if ( needOnAndOff && nValues < 2 )
            {
                QString errorText = "On and off times both are needed for events of type square or ramp";
                QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                return(1);
            }
            else if ( nValues >= 2 )
            {  // read 2nd value (off or mag)
                QString nextWord = valueList.at(1);  bool ok;
                dNext = nextWord.toDouble(&ok);
                if ( !ok )
                {
                    QString errorText;
                    bool needOn = (_challengeShape[indexChallenge] == Challenge_Gamma)
                            ||    (_challengeShape[indexChallenge] == Challenge_Sigmoid);
                    if ( needOnAndOff )
                        errorText = QString("Error reading the OFF time for event %1 from file %2").arg(getEventChar(indexChallenge)).arg(fileName);
                    else if ( needOn )
                        errorText = QString("Error reading the magnitude from file %1").arg(fileName);
                    QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
                    return(1);
                } // !ok
                if ( needOnAndOff )
                    _challengeOff[indexChallenge][iStimulus] = dNext;
            }
            iStimulus = getFirstBadStimulus(indexChallenge);
        }  // stimulus
    }
    file.close();
    setPrepared(false);
    return(0);
}

dVector PETRTM::getEquilibrationVector(int iFile)
{   // Create the vector: k4 * (1 + BPnd) using saved values of k4 and BPnd
    //    1) Fit for BPnd using fixed k4, with BPnd saved from last iteration as _BPndForIterations
    // OR 2) Fit k4 using fixed BPnd (saved from 1st pass in _BPndForIterations)
    int nTimeInRun = _dtBinsSec[iFile].size();
    dVector equilibrationVector;  equilibrationVector.fill(0.,nTimeInRun);
//    qDebug() << "getEquilibrationVector tau4" << _tau4[0];
    for (int jt=0; jt<nTimeInRun; jt++)
        equilibrationVector[jt] = 1./_tau4[iFile] * (1. + _BPndForIterations[iFile] );
    return equilibrationVector;
}

dVector PETRTM::getBPndVector(int iFile)
{ // get value and error (into x and y)
    int nTimeInRun = _dtBinsSec[iFile].size();
    dVector BPndVector;  BPndVector.fill(0.,nTimeInRun);
    // Baseline k2 and k2a define baseline BPnd
    dPoint2D BP;        BP.x = BP.y = 0.;
    // Baseline k2 and k2a
    int iCoeffk2  = _k2EventCoefficient[iFile];
    int iCoeffk2a = _k2aEventCoefficient[iFile];
    if ( iCoeffk2 < 0 || iCoeffk2a < 0 ) return BPndVector;

    double k2     = getBeta(iCoeffk2);
    dVector k2aVector; k2aVector.resize(nTimeInRun);
    for (int jt=0; jt<nTimeInRun; jt++)
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
    for (int jt=0; jt<nTimeInRun; jt++)
        allValid &= k2aVector[jt] > 0.;
    if ( allValid )
    {
        if ( isFRTMNew() )
        {
            for (int jt=0; jt<nTimeInRun; jt++)
//                BPndVector[jt] = _tau4[iFile] * k2aVector[jt]/k2;
                BPndVector[jt] = k2aVector[jt]/k2;

        }
        else
        {
            for (int jt=0; jt<nTimeInRun; jt++)
                BPndVector[jt] = k2 / k2aVector[jt] - 1;
        }
    }
    return BPndVector;
}

dPoint2D PETRTM::getBPndVersusTime(int iFile, int iTime)  // VERY INEFFICIENT METHOD: EVERY POINT ITIME CALLS createChallengeShape
{ // get value and error (into x and y)
    // Baseline k2 and k2a define baseline BPnd
    dPoint2D BP;        BP.x = BP.y = 0.;
    // Baseline k2 and k2a
    int iCoeffk2  = _k2EventCoefficient[iFile];
    int iCoeffk2a = _k2aEventCoefficient[iFile];
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
        if ( isFRTMNew() )
//            BP.x= _tau4[iFile] * k2a/k2;
            BP.x= k2a/k2;
        else
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
//    qDebug() << "PETRTM::getBPndInCurrentCondition enter";
    dPoint2D BP;  BP.x = BP.y = 0.;

    updateConditions();

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
                if ( isFRTMNew() )
                { // BP = "k2a" / k2
                    double BPnd = k2a/k2;
                    BP.x += BPnd * _contrastMatrix[iCondition][0][iCoeff];
                    BP.y += BPnd * BPnd * ( k2Err2/k2/k2 + k2aErr2/k2a/k2a );
                }
                else
                { // BP = k2/k2a - 1
                    double BPnd = k2/k2a - 1.;
                    BP.x += BPnd * _contrastMatrix[iCondition][0][iCoeff];
                    BP.y += BPnd * BPnd * ( k2Err2/k2/k2 + k2aErr2/k2a/k2a );
                }
            }
            else // if ( _basisShape[iCoeff] == Type_challenge )
            {
                // k2a      = k2 * BP  , BP = k2a/k2
                // deltak2a = k2 * dBP / (1+BP+dBP)
                // dBP = deltak2a/k2*(1+k2a/k2) / (1-deltak2a/k2)
                if ( isFRTMNew() ) // "delta_k2a" = delta_k2k3 = k2 * k4 * dBP * BP/(1+BP)
                {
                    double DBPnd = - deltak2a/k2 * (1+k2a/k2) / (1 - deltak2a/k2);
                    BP.x += DBPnd * _contrastMatrix[iCondition][0][iCoeff];
                    BP.y += DBPnd * DBPnd * ( k2Err2/k2/k2 + deltak2aErr2/deltak2a/deltak2a );
                }
                else
                {
                    double totalk2a = k2a + deltak2a;
                    double DBPnd = k2/k2a * deltak2a / totalk2a;
                    BP.x += DBPnd * _contrastMatrix[iCondition][0][iCoeff];
                    BP.y += DBPnd * DBPnd * ( k2Err2/k2/k2 + k2aErr2/k2a/k2a + deltak2aErr2/deltak2a/deltak2a +
                                              (k2aErr2+deltak2aErr2)/totalk2a/totalk2a );
                }
            }
        }
    }
    BP.y = qSqrt(BP.y);
//    qDebug() << "PETRTM::getBPndInCurrentCondition exit" << BP.x << BP.y;
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
//        qDebug() << "PETRTM::averageParameter" << iCoeff << result.x;
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
    double k2  = getBeta(_k2EventCoefficient[iRun]);
    double k2a = getBeta(_k2aEventCoefficient[iRun]);
    if ( k2 == 0. || k2a == 0. )
        return 0.;
    else
    {
        if ( isFRTMNew() )
        {
//            BPnd = _tau4[iRun] * k2a / k2;
            BPnd = k2a / k2;
//            qDebug() << "FRTMNew: k2a =" << k2a << "k2 =" << k2 << "_tau4 =" << _tau4[iRun] << "BPnd =" << BPnd;
        }
        else
        {
            BPnd = k2/k2a - 1.;
//            qDebug() << "FRTMOld: k2a =" << k2a << "k2 =" << k2 << "BPnd =" << BPnd;
        }
    }
    return BPnd;
}

double PETRTM::getTau4InRun(int iRun)
{ // "k2a" = k2 * k4 * BPnd, so 1/k4 = k2/k2a * BPnd
    if ( _k2EventCoefficient[iRun] < 0 || _k2aEventCoefficient[iRun] < 0 ) return 0.;

    double tau4;
    if ( isFRTMFitk4() )
    { // calculate 1/k4 from k2 and "k2a" = k2 * k4
        double k2  = getBeta(_k2EventCoefficient[iRun]);
        double k2a = getBeta(_k2aEventCoefficient[iRun]);
        if ( k2 == 0. || k2a == 0. )
            tau4 = 0.;
        else
            tau4 = k2 / k2a;
//        qDebug() << "PETRTM::getTau4InRun ongoing" << tau4;
    }
    else if (_fitk4UsingFixedBPnd ) // k4 fitting has been completed, so return saved value
    {
        tau4 = _tau4[iRun];
//        qDebug() << "PETRTM::getTau4InRun completed" << tau4;
    }
    else if ( isFRTM() )
        tau4 = _tau4[iRun];
    else
//        tau4 = 0.;
        tau4 = _tau4[iRun];

    return tau4;
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
    if ( isRTM2() )
    {   // set tau2Ref for the specified run
        if ( isSRTM() )
            _tau2RefSRTMFixed[iRun] = tau2Ref;
        else
            _tau2RefFRTMFixed[iRun] = tau2Ref;
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
                    _tau2RefSRTMFixed[jRun] = tau2Ref;
                else
                    _tau2RefFRTMFixed[jRun] = tau2Ref;
            }
        }
    }
    setPrepared(false);
};


double PETRTM::getTau2RefInRun(int iRun)
{ // tau2_ref = R1 / k2
    double tau2Ref=0.;
    if ( isRTM3() )
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
            tau2Ref = _tau2RefSRTMFixed[iRun];
        else
            tau2Ref = _tau2RefFRTMFixed[iRun];
    }
    return tau2Ref;
}

dPoint2D PETRTM::getk2RefInRun(int iRun)
{ // k2_ref = k2 / R1
    dPoint2D k2Ref;  k2Ref.x = k2Ref.y = 0.; // save value and error in x,y
    if ( isRTM3() )
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
            k2Ref.x = 1./_tau2RefSRTMFixed[iRun];
        else
            k2Ref.x = 1./_tau2RefFRTMFixed[iRun];
        k2Ref.y = 0.;
    }
    return k2Ref;
}

dPoint2D PETRTM::getR1InRun(int iRun)
{
    dPoint2D R1;  R1.x = R1.y = 0.; // save value and error in x,y
    if ( isRTM3() )
    {
        if ( _R1EventCoefficient[iRun] < 0 ) return R1;
        R1.x = getBeta(_R1EventCoefficient[iRun]);
        R1.y = getSEM(_R1EventCoefficient[iRun]);
    }
    return R1;
}

dPoint2D PETRTM::getk2InRun(int iRun)
{
    dPoint2D k2;  k2.x = k2.y = 0.; // save value and error in x,y
    k2.x = getBeta(_k2EventCoefficient[iRun]);
    k2.y = getSEM(_k2EventCoefficient[iRun]);
    return k2;
}

dPoint2D PETRTM::getk2aInRun(int iRun)
{
    dPoint2D k2a;  k2a.x = k2a.y = 0.; // save value and error in x,y
    k2a.x = getBeta(_k2aEventCoefficient[iRun]);
    k2a.y = getSEM(_k2aEventCoefficient[iRun]);
    return k2a;
}

dPoint2D PETRTM::getdk2aInRun(QChar challengeID)
{
    dPoint2D dk2a;  dk2a.x = dk2a.y = 0.; // save value and error in x,y
    int iCoeffChallenge = getEventCoefficient(challengeID);
    if ( iCoeffChallenge >= 0 )  // should always be true
    {
        dk2a.x = getBeta(iCoeffChallenge);
        dk2a.y = getSEM(iCoeffChallenge);
    }
    return dk2a;
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

void PETRTM::setTissueVector(dMatrix tissueRegion)
{
    QMutex mutex;
    mutex.lock();
    _tissRegionRaw = tissueRegion;
    _tissRegion    = tissueRegion;

    if ( _smoothingScaleMin != 0. )
        fitLoessCurve(_tissRegion);  // this replaces the tissue region with the fit
    _tissRegionDeriv = _tissRegion;  // this ensure dimensions are correct
    if ( (isFRTM() && _tau4Default != 0.) || _PETWeightingModel == Weights_noUptake )
        differentiateByRun(_tissRegionDeriv); // rFRTM requires a tissue derivative for the convolution term
    mutex.unlock();
    setPrepared(false);
}

dPoint2D PETRTM::getReferenceRegion(bool useFit, int iRun, int iTime)
{ // pass useFit explicitly, so that this can be called externally to get both raw and fit
    dPoint2D returnValue;
//    returnValue.x= getTimeInRun(iRun,iTime);
    returnValue.x = _timeInRun[iRun][iTime];
    if ( useFit )
        returnValue.y = _refRegion[iRun][iTime];     // could be raw or fit
    else
        returnValue.y = _refRegionRaw[iRun][iTime];  // always raw
    return returnValue;
}
dPoint2D PETRTM::getTissueRegion(bool useFit, int iRun, int iTime)
{
    dPoint2D returnValue;
//    returnValue.x = getTimeInRun(iRun,iTime);
    returnValue.x = _timeInRun[iRun][iTime];
    if ( useFit )
        returnValue.y = _tissRegion[iRun][iTime];    // could be either raw or fit
    else
        returnValue.y = _tissRegionRaw[iRun][iTime]; // always raw
    return returnValue;
}

dPoint2D PETRTM::getReferenceRegionTimesR1(int iRun, int iTime)
{ // legacy function (xxx jbm)
    dPoint2D returnValue;
//    returnValue.x = getTimeInRun(iRun,iTime);
    returnValue.x = _timeInRun[iRun][iTime];

    if ( isRTM3() )
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
            returnValue.y = _refRegion[iRun][iTime] * k2 * _tau2RefSRTMFixed[iRun];
        else
            returnValue.y = _refRegion[iRun][iTime] * k2 * _tau2RefFRTMFixed[iRun];
    }
    return returnValue;
}

void PETRTM::defineFrameInterpolation(int iRun)
{
//    qDebug() << "defineFrameInterpolation enter" << iRun << _binSplit.size();
    // First find the minimum bin width
    int nTimeInRun = _dtBinsSec[iRun].size();
//    qDebug() << "defineFrameInterpolation 1" << _binSplit[iRun].size() << nTimeInRun;
    if ( _binSplit[iRun].size() != nTimeInRun )   _binSplit[iRun].resize(nTimeInRun);
//    qDebug() << "defineFrameInterpolation 2" << _dtBinsSec[iRun];
    _minBin[iRun]=1e6;
    for (int jt=0; jt<nTimeInRun; jt++)
    {
        _binSplit[iRun][jt].clear();
         if ( _dtBinsSec[iRun][jt] < _minBin[iRun] ) _minBin[iRun] = _dtBinsSec[iRun][jt];
    }
//    qDebug() << "defineFrameInterpolation minBin" << _minBin[iRun];
    int iBinFine=0;
    for (int jt=0; jt<nTimeInRun; jt++)
    {
        int nSplit = _interpolationFactor * _dtBinsSec[iRun][jt] / _minBin[iRun];
        if ( _dtBinsSec[iRun][jt] % _minBin[iRun] )
        {
            // designate the convolutions (rFRTM) cannot be used by clearing vectors
            _binSplit[iRun].clear();
            return;
        }
        else
        {
            _binSplit[iRun][jt].clear();
            for (int jSplit=0; jSplit<nSplit; jSplit++, iBinFine++)
                _binSplit[iRun][jt].append(iBinFine);
        }
    }
//    for (int jt=0; jt<nTimeInRun; jt++)
//        qDebug() << "binSplit[" << jt << "] =" << _binSplit[iRun][jt];
//    qDebug() << "bdefineFrameInterpolation exit";
}
dVector PETRTM::interpolateTissueVector(int iRun, dVector tissueVector)
{
    int nTime = _dtBinsSec[iRun].size();
    dVector tissueFine;
    for (int jt=0; jt<nTime; jt++)
    {
        int nSplit = _binSplit[iRun][jt].size() ;
        if ( nSplit == 1 )
        {
            tissueFine.append(tissueVector[jt]);
//            int indexFine = _binSplit[iRun][jt][0];  // original value
//            tissueFine[indexFine] = tissueVector[jt];
        }
        else
        {
//            double t1 = getTimeInRun(iRun,jt);
            double t1 = _timeInRun[iRun][jt];
            double y1 = tissueVector[jt];
            double binSize = _dtBinsSec[iRun][jt]/60.;
            double tLowerEdge = t1 - binSize/2.;

            // Lower point at jt-1 unless at endpoint
            double t0, y0;
            if ( jt != 0 )
            {
//                t0 = getTimeInRun(iRun,jt-1);
                t0 = _timeInRun[iRun][jt-1];
                y0 = tissueVector[jt-1];
            }
            else
                t0 = y0 = 0.;

            // Upper point at jt+1 unless at endpoint
            double t2, y2;
            if ( jt != nTime-1 )
            {
//                t2 = getTimeInRun(iRun,jt+1);
                t2 = _timeInRun[iRun][jt+1];
                y2 = tissueVector[jt+1];
            }
            else
            {
                t2 = t1 + binSize/2.;
                y2 = y1;
            }

            double binSizeFine = binSize / static_cast<double>(nSplit);
            for (int jSplit=0; jSplit<nSplit; jSplit++)
            {
                double t = tLowerEdge + binSizeFine * (jSplit + 0.5);
                tissueFine.append(interpolateNewtonDividedDifferencesQuad(t,t0,t1,t2,y0,y1,y2));
            }
        }
    }
//    qDebug() << "input tissue vector" << tissueVector;
//    qDebug() << "fine  tissue vector" << tissueFine;
    return tissueFine;
}
dVector PETRTM::combineFineTissueVector(int iRun, dVector tissueFine)
{
    int nTime = _dtBinsSec[iRun].size();
    dVector tissueVector;  tissueVector.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        int nSplit = _binSplit[iRun][jt].size();
        double sum=0.;
        for (int jSplit=0; jSplit<nSplit; jSplit++)
        {
            int indexFine = _binSplit[iRun][jt][jSplit];
            sum += tissueFine[indexFine];
        }
        tissueVector[jt] = sum / static_cast<double>(nSplit);
    }
    return tissueVector;
}
double PETRTM::interpolateNewtonDividedDifferencesQuad(double x, double x0, double x1, double x2, double f0, double f1, double f2)
{
    double secondDiff = ( (f1-f2)/(x1-x2) - (f0-f1)/(x0-x1) ) / (x2-x1);
    double f = f0 + (x-x0) * (f1-f0)/(x1-x0) + (x-x0)*(x-x1) * secondDiff;
    return f;
}

void PETRTM::setTimeBinsSec(int iRun, iVector timeBinsSec)
{
    if ( iRun < _nRuns && iRun < _refRegion.size() )
    {
        if ( timeBinsSec.size() != _refRegion[iRun].size() )
            qInfo() << "Error: the number of time bins (" << timeBinsSec.size() << ") does not match the reference region (" << _refRegion[iRun].size() << ").";
        else
        {
            _dtBinsSec[iRun] = timeBinsSec;
            defineFrameInterpolation(iRun);
            for (int jt=0; jt<timeBinsSec.size(); jt++)
                _timeInRun[iRun][jt] = getTimeInRun(iRun, jt);
        }
    }
    else
        qFatal("Error: attempting to set time bins for non-existent run.");
    setPrepared(false);
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
        int nTimeInRun = _dtBinsSec[jRun].size();
        referenceRegion[jRun].fill(0.,nTimeInRun);
        for (int jt=0; jt<nTimeInRun; jt++)
            referenceRegion[jRun][jt] = _table[jRun][jt][iColumn]; // set RR from table column
    }
    setReferenceRegion(referenceRegion);
}

void PETRTM::setReferenceRegion(iMatrix timeBinsSec, dMatrix referenceRegionRaw) // [nRuns][nTimePerRun]
{  // timeBins in units of seconds
    int nRuns = referenceRegionRaw.size();
    setNumberRuns(nRuns);
    for (int jRun=0; jRun<nRuns; jRun++)
    {
        setNumberTimePointsInRun(jRun,referenceRegionRaw[jRun].size());
        setTimeBinsSec(jRun,timeBinsSec[jRun]);
        int nTimeInRun   = _dtBinsSec[jRun].size();
        for (int jt=0; jt<nTimeInRun; jt++)
            _table[jRun][jt][0] = static_cast<double>(timeBinsSec[jRun][jt]);  // put the time bins into table column 0
    }
    setReferenceRegion(referenceRegionRaw);
}

void PETRTM::setReferenceRegion(dMatrix referenceRegionRaw)
{  // set the raw and fit version of the reference region
    _refRegionRaw = referenceRegionRaw;
    updateReferenceRegion();
    setPrepared(false);
    _referenceRegionIsDefined = true;
}

void PETRTM::updateReferenceRegion()
{
    _refRegion = _refRegionRaw;
    if ( _smoothingScaleMin != 0. )
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
//        goodStimulus &= _challengeOn[indexChallenge][indexStimulus] >= getTimeInRun(iRun,0);
        goodStimulus &= _challengeOn[indexChallenge][indexStimulus] >= _timeInRun[iRun][0];
    else if ( goodStimulus && infoRequired == 2 )
    { // stim OFF > ON && ON falls within run timing
        goodStimulus &= _challengeOff[indexChallenge][indexStimulus] > _challengeOn[indexChallenge][indexStimulus];
        goodStimulus &= _challengeOn[indexChallenge][indexStimulus]  >= 0.;
        int nTimeInRun = _dtBinsSec[iRun].size();
//        goodStimulus &= _challengeOn[indexChallenge][indexStimulus]  < getTimeInRun(iRun,nTimeInRun-1);
        goodStimulus &= _challengeOn[indexChallenge][indexStimulus]  < _timeInRun[iRun][nTimeInRun-1];
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

    if ( _nRuns != _dtBinsSec.size() )
    { // the following block needs only 1 allocation
        _R1EventID.fill('1',_nRuns);
        _k2EventID.fill('a',_nRuns);
        _k2aEventID.fill('A',_nRuns);
        _dCrdtEventID.fill('v',_nRuns);
        _R1EventCoefficient.resize(_nRuns);
        _k2EventCoefficient.resize(_nRuns);
        _k2aEventCoefficient.resize(_nRuns);
        _dCrdtEventCoefficient.resize(_nRuns);
        _frameFiles.fill("",_nRuns);
        _tau4.fill(0.,_nRuns);
        _tau2RefSRTMFixed.fill(3.5,_nRuns);
        _tau2RefFRTMFixed.fill(1.5,_nRuns);
        _BPndForIterations.fill(0.,_nRuns);
        _tau2RefSRTMCal.resize(_nRuns);
        for (int jRun=0; jRun<_nRuns; jRun++)
        {
            _tau2RefSRTMCal[jRun].resize(3);
            _tau2RefSRTMCal[jRun][0] = _tau2RefSRTMCal[jRun][1] = _tau2RefSRTMCal[jRun][2] = 0.;
        }
        _weights.resize(_nRuns);
        _ignoreString.resize(_nRuns);
        _dtBinsSec.resize(_nRuns);
        _timeInRun.resize(_nRuns);
        _minBin.fill(0,_nRuns);
        _binSplit.resize(_nRuns);
        _table.resize(_nRuns);     // [_nRuns][nTimeInRun][_nColumns]
        _quadLOESS.resize(_nRuns);     // nTimeInRun
        _columnNames.resize(_nRuns);

        _refRegion.resize(_nRuns);
        _refRegionRaw.resize(_nRuns);
        _refRegionIntegral.resize(_nRuns);
        _refRegionDeriv.resize(_nRuns);
        _tissRegionRaw.resize(_nRuns);
        _tissRegion.resize(_nRuns);
        _tissRegionDeriv.resize(_nRuns);
        _frtmConv_dCtdtE.resize(_nRuns);
        _frtmConv_dCtdtERaw.resize(_nRuns);
        _frtmConv_CtE.resize(_nRuns);
        _frtmConvDeWeightUptake.resize(_nRuns);

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
        _challengeTau.fill(10.,_maxChallenges);
        _challengeAlpha.fill(1.,_maxChallenges);

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
    int nTimeInRun = _refRegion[iRun].size();
    if ( nTimeInRun != _dtBinsSec[iRun].size() )
    {
        _warningErrors->append(QString("Error: the number of time points in the table (%1) does not match the number of time points (%2)").
                               arg(nTimeInRun).arg(_dtBinsSec[iRun].size()));
        return(1);
    }
    else if ( !_referenceRegionIsDefined )
    {
        _warningErrors->append(QString("Error: the reference region is not yet defined"));
        return(1);
    }

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return(1);
    QTextStream out(&file);

    // Write the reference region as a table file containing: 0) time bins (convert to sec), 1) reference region
    _table[iRun].clear();  // this is necessary, because the old table file is overwritten in order to import the ref region
    _table[iRun].resize(nTimeInRun);  // resize table vector
    out << "delta_time " << RRName << "\n";
    for (int jt=0; jt<nTimeInRun; jt++)
    {
        out << _dtBinsSec[iRun][jt] << " " << _refRegion[iRun][jt] << "\n";
        _table[iRun][jt].append(_dtBinsSec[iRun][jt]);
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

    QString errorString = utilIO::readTimeTableFile(fileName, _columnNames[iRun], _table[iRun]);
    if ( !errorString.isEmpty() )
    {
        _warningErrors->append(errorString);
        return(1);
    }

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

    int nTime = _table[iRun].size();
    iVector timeBinsSec; timeBinsSec.resize(nTime);
    for ( int jt=0; jt<nTime; jt++)
        timeBinsSec[jt] = static_cast<int>(_table[iRun][jt][0]);
    setTimeBinsSec(iRun,timeBinsSec);

    if ( _smoothingScaleMin != 0. )
        _quadLOESS[iRun].define(_timeInRun[iRun],_smoothingScaleMin);

    setPrepared(false);

    return(0);
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
    if ( _smoothingScaleMin != 0. )
        out << "smoothing-scale " << _smoothingScaleMin << "\n";

    // weighting model
    if ( _PETWeightingModel == Weights_11C )
        out << "weights 11C\n";
    else if ( _PETWeightingModel == Weights_11C_Noiseless )
        out << "weights 11C-noiseless\n";
    else if ( _PETWeightingModel == Weights_18F )
        out << "weights 18F\n";
    else if ( _PETWeightingModel == Weights_18F_Noiseless )
        out << "weights 18F-noiseless\n";
    else if ( _PETWeightingModel == Weights_noUptake )
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
            out << "1/k2' " << _tau2RefFRTMFixed[iRun] << "\n\n";
        else if ( isSRTMReg() )
            out << "1/k2' " << _tau2RefSRTMFixed[iRun] << " " << _tau2RefSRTMCalOffset << " "
                << _tau2RefSRTMCal[iRun][0] << " " << _tau2RefSRTMCal[iRun][1] << " " << _tau2RefSRTMCal[iRun][2] <<  "\n\n";
        else
            out << "1/k2' " << _tau2RefSRTMFixed[iRun] << "\n\n";
    }
    else
        out << "R1 " << _R1EventID[iRun] << "\n\n"; // R1 with 3-parameter versions
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
                    out << " sigmoid " << _challengeTau[jChallenge];
                else if ( _challengeShape[jChallenge] == Challenge_Gamma )
                {
                    out << " gamma " << _challengeTau[jChallenge];
                    if ( getChallengeAlpha(jChallenge) != 1. )
                        out << " " << _challengeAlpha[jChallenge];
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
            out << _tau2RefFRTMFixed[iRun] << "\n\n";
        else
            out << _tau2RefSRTMFixed[iRun] << "\n\n";
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
                    out << " sigmoid " << _challengeTau[jChallenge];
                else if ( _challengeShape[jChallenge] == Challenge_Gamma )
                {
                    out << " gamma " << _challengeTau[jChallenge];
                    if ( getChallengeAlpha(jChallenge) != 1. )
                        out << " " << _challengeAlpha[jChallenge];
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

double PETRTM::getTimeInRun(int iRun, double rBin)
{ // get the time at bin "rBin" included in (0,nTimePerRun); the first bin starts at 0 and ends at 1 with center point 0.5
    if ( iRun < 0 || rBin < 0. ) return -1.;
    int iBin = qFloor(rBin);
    // subtract half the current bin to align definitions of "iBin" and "rBin"
    double deltaTimeInBin = static_cast<double>(_dtBinsSec[iRun][iBin])/60.;
    double timeLastBin = getTimeInRun(iRun,iBin) - deltaTimeInBin/2.;  // switch from bin-centered times
    double deltaBin = rBin - iBin;
    double deltaTime = deltaBin * deltaTimeInBin; // add rest of time bin
    return timeLastBin + deltaTime;
}

double PETRTM::getTimeInRun(int iRun, int iTime)
{ // get the time at time point iTime
    if ( iRun < 0 || iTime < 0 ) return -1.;

    // This first bin is centered at 1/2 the bin width
    double dt0  = static_cast<double>(_dtBinsSec[iRun][0])/60.;
    double time = dt0/2.;
    for ( int jt=1; jt<=iTime; jt++ )
    {
        double lastDelta = static_cast<double>(_dtBinsSec[iRun][jt-1])/60.;
        double thisDelta = static_cast<double>(_dtBinsSec[iRun][jt])/60.;
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
    _challengeTau[indexChallenge] = tau;
    setPrepared(false);
}
void PETRTM::setChallengeAlpha(int indexChallenge, double alpha)
{
    _challengeAlpha[indexChallenge] = alpha;
    setPrepared(false);
}

bool PETRTM::getFrameStatus()
{
    bool allFilesRead = true;
    for (int jRun=0; jRun<_nRuns; jRun++)
        allFilesRead &= !_frameFiles[jRun].isEmpty();  // true if no file names are empty
    return allFilesRead;
}

void PETRTM::setNumberTimePointsInRun(int iRun, int nTime)
{
    if ( iRun<0 || iRun >=_nRuns )
    {
        qWarning() << "Error: iRun = " << iRun << " but # runs = " << _nRuns;
        exit(1);
    }
    if ( nTime != _dtBinsSec[iRun].size() )
    {
        _dtBinsSec[iRun].clear();
        _table[iRun].resize(nTime);
        for (int jt=0; jt<nTime; jt++)
        {
            _dtBinsSec[iRun].append(60);
            _table[iRun][jt].append(60.);  // initialize time bin widths (1st column=frames) to 1.
        }
        _dCrdtEventCoefficient.fill(-1,nTime);
        _refRegion[iRun].fill(0.,nTime);
        _refRegionRaw[iRun].fill(0.,nTime);
        _refRegionIntegral[iRun].fill(0.,nTime);
        _refRegionDeriv[iRun].fill(0.,nTime);
        _tissRegionRaw[iRun].fill(0.,nTime);
        _tissRegion[iRun].fill(0.,nTime);
        _frtmConv_dCtdtE[iRun].fill(0.,nTime);
        _frtmConv_dCtdtERaw[iRun].fill(0.,nTime);
        _frtmConv_CtE[iRun].fill(0.,nTime);
        _frtmConvDeWeightUptake[iRun].fill(0.,nTime);
        _timeInRun[iRun].fill(1.,nTime);
        setTimeBinsSec(iRun,_dtBinsSec[iRun]);
    }
    setPrepared(false);
}

void PETRTM::prepare()
{
    // Perhaps do error checking here first?
//    updateReferenceRegion();

    // For iterative models, add the 1st fit to the "prepare" step
    bool multiPass = isFRTM() || isSRTMReg() || _PETWeightingModel == Weights_noUptake;
    if ( multiPass )
    {
        if ( _PETWeightingModel == Weights_noUptake )
        {
            setWeightingModel(Weights_Uniform);
            createAllBasisFunctions();
            fitData(_tissRegion);
            setWeightingModel(Weights_noUptake);  // now weighting model has BPnd info
        }
        else if ( isSRTMReg() )
        {
            setRTMModelType(RTM_SRTM2);
            createAllBasisFunctions();
            fitData(_tissRegion);
            setRTMModelType(RTM_SRTM2Fit);  // now SRTMReg has BPnd info
        }
        else if ( isFRTM() )
        {
            RTMModelTypes saveModel = _modelRTM;
            if ( isRTM3() )
                setRTMModelType(RTM_SRTM3);
            else
                setRTMModelType(RTM_SRTM2);
            createAllBasisFunctions();
            fitData(_tissRegion);
//            qDebug() << "PETRTM::prepare, BP0 =" << getBP0InRun(0);
            setRTMModelType(saveModel); // now FRTM has BPnd info
        }
    }
    // Now prepare the model using prior information from above if needed
    createAllBasisFunctions();
}

void PETRTM::fitData(QVector<ROI_data> timeSeriesVector, dMatrix &yFit)
{
    // iterative methods
    bool fitByIteration = (isFRTM() && _tau4Default != 0.);
    if (  fitByIteration )
        fitDataByGLMIterationForConsistency(timeSeriesVector, yFit);
    else
        fitDataByGLM(timeSeriesVector, yFit); // 1-step

    // potentially fit k4 as 2nd stage
    if ( _fitk4UsingFixedBPnd )  // k4 fitting has been selected, but the current 1st-stage model could be anything
    {

        RTMModelTypes saveTimeModel = _modelRTM;
        // Change models, prepare new basis functions (with newTAC=true), and fit the model as a starting point
        if ( isRTM3() )
            setRTMModelType(RTM_rFRTM3New);
        else
            setRTMModelType(RTM_rFRTM2New);
        createAllBasisFunctions();
        // Now fit by k4 line scan
        fitDataByGLMPlusLineScan(timeSeriesVector, yFit);

        setRTMModelType(saveTimeModel);
        prepare();
        fitDataByGLM(timeSeriesVector, yFit);
        if (  fitByIteration )
            fitDataByGLMIterationForConsistency(timeSeriesVector, yFit);
    }

}

void PETRTM::fitData(dMatrix timeSeriesVector)
{ // convenience function when yFit can be discarded (e.g., want access to beta values)
    dMatrix yFit = timeSeriesVector;  // resize to correct dimensions
    fitData(timeSeriesVector, yFit);
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
dVector PETRTM::makeVectorFromROIData(QVector<ROI_data> timeSeriesVector)
{
    int nTimeTotal = getNumberTimePoints();
    dVector data;
    data.fill(0.,nTimeTotal);
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        int iStartRun = jRun * _dtBinsSec[jRun].size();
        for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
            data[jt+iStartRun] = timeSeriesVector[jRun].ySignal[jt];
        iStartRun += _dtBinsSec[jRun].size();
    }
    return data;
}

void PETRTM::fitDataByGLM(dMatrix timeSeriesVector)
{ // convenience function when yFit can be discarded (e.g., want access to beta values)
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
    fitDataByGLM(ROIVector);
}
void PETRTM::fitDataByGLM(QVector<ROI_data> timeSeriesVector)
{ // convenience function when yFit can be discarded (e.g., want access to beta values)
    dMatrix yFit;
    int nRuns = timeSeriesVector.size();
    yFit.resize(nRuns);
    for ( int jRun=0; jRun<nRuns; jRun++)
        yFit[jRun].resize(_dtBinsSec[jRun].size());
    fitDataByGLM(timeSeriesVector, yFit);
}
void PETRTM::fitDataByGLM(QVector<ROI_data> timeSeriesVector, dMatrix &yFit)
{ // timeSeriesVector follows the ROI_data structure, with 1d vectors for x and y
    // Copy the ROI input data into a simple vector
    dVector data = makeVectorFromROIData(timeSeriesVector);
    fitWLSForIteration(data);
    // Attach the fit.
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        int iStartRun = jRun * _dtBinsSec[jRun].size();
        for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
        {
            int iTimeTotal = jt + iStartRun;
            yFit[jRun][jt] = getFit(iTimeTotal);
        }
    }
}

void PETRTM::fitWLSForIteration(dVector &data)
{

    // fit by WLS using generalGLM
    fitWLS(data,true);

    // Save values that could be needed for iterative methods
    if ( isFRTMFitk4() )
    {
        for (int jRun=0; jRun<_nRuns; jRun++)
            _tau4[jRun] = getTau4InRun(jRun);
    }
    else
    {
        for (int jRun=0; jRun<_nRuns; jRun++)
            _BPndForIterations[jRun] = getBP0InRun(jRun);
//        qDebug() << "save BPnd for iterations" <<_BPndForIterations[0];
    }

}


void PETRTM::fitDataByGLMIterationForConsistency(QVector<ROI_data> timeSeriesVector, dMatrix &yFit)
{   // Fit by GLM iteration when linear parameters determine a value used in a non-linear part of the expression.
    // E.g., use this for rFRTM as originally published, where BPnd=k2/k2a-1 from one iteration is used to form basis functions for the next iteration.
    // The stopping criteria is consistency: the value of the iterative parameter does not change.
    // Assume that the TAC has already been fit by an initial model at this point
    // Copy the ROI input data into a simple vector
    fitDataByGLM(timeSeriesVector, yFit);

    // Make sure all BP0 values are positive
    bool moreIterations=true;
    for (int jRun=0; jRun<_nRuns; jRun++)
        moreIterations &= getBP0InRun(jRun) > 0.;
    bool goodIteration = moreIterations;

    _nIterations=0;
    int maxIterations = 20;

    // the model should have been fit using fitDataByGLM prior to this function
    double lastSigma2 = getSigma2();

/*
    if ( isFRTMFitk4() )
        qDebug() << "iteration" << _nIterations << ", tau4 =" << _tau4[0];
    else
        qDebug() << "iteration" << _nIterations << ", BP0 =" << getBP0InRun(0);
    qDebug() << "*** starting:" << getNumberCoefficients() << "coefficients ***";
    for (int jCoeff=0; jCoeff<getNumberCoefficients(); jCoeff++)
        qDebug() << "coefficient" << jCoeff << "=" << getBeta(jCoeff);
*/
    moreIterations = false;  // turn off loop: just fit 1ce with new model
    while (goodIteration && moreIterations)
    {
        createAllBasisFunctions();
        fitDataByGLM(timeSeriesVector, yFit); // 1-step

        double sigma2 = getSigma2();
        double diff = qAbs(sigma2-lastSigma2)/lastSigma2 * 100.;
//        qDebug() << "lastSigma2, sigma2, diff" << lastSigma2 << sigma2 << diff;
        lastSigma2 = sigma2;
        moreIterations = diff > 0.01;
//        goodIteration = lastSigma2 > sigma2;
        goodIteration = true;
        for (int jRun=0; jRun<_nRuns; jRun++)
        {
            if ( isFRTMFitk4() ) // make sure all values of 1/k4 > 0
                goodIteration &= _tau4[jRun] > 0.;
            else // make sure all values of BPnd > 0
                goodIteration &= getBP0InRun(jRun) > 0.;
        }
        _nIterations++;
        moreIterations &= _nIterations < maxIterations;
/*
        if ( isFRTMFitk4() )
            qDebug() << "iteration" << _nIterations << ", tau4 =" << _tau4[0];
        else
            qDebug() << "iteration" << _nIterations << ", BP0 =" << getBP0InRun(0);
        qDebug() << "***" << getNumberCoefficients() << "coefficients ***";
        for (int jCoeff=0; jCoeff<getNumberCoefficients(); jCoeff++)
            qDebug() << "coefficient" << jCoeff << "=" << getBeta(jCoeff);
*/
    }

    if ( !goodIteration )
        _nIterations *= -1;
    else
    {
        // Attach the fit.
        for (int jRun=0; jRun<_nRuns; jRun++)
        {
            int iStartRun = jRun * _dtBinsSec[jRun].size();
            for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
            {
                int iTimeTotal = jt + iStartRun;
                yFit[jRun][jt] = getFit(iTimeTotal);
            }
        }
    }
//    qDebug() << "PETRTM::fitDataByGLMIterationForConsistency exit";
}

void PETRTM::fitDataByGLMPlusLineScan(QVector<ROI_data> timeSeriesVector, dMatrix &yFit)
{   // Fit by GLM plus a line scan on an additional parameter.
    // E.g., use this with rFRTM and a fixed value of k4 (or BPnd) plus a line scan on BPnd (or k4)
    // Assume that the TAC has already been fit by an initial model at this point
//    qDebug() << "PETRTM::fitDataByGLMPlusLineScan enter" << isFRTM();
    // Copy the ROI input data into a simple vector
    dVector data = makeVectorFromROIData(timeSeriesVector);

    if ( _fitk4UsingFixedBPnd )
    { // fix BPnd and vary 1/k4
        for ( int jRun=0; jRun<_nRuns; jRun++ )
        {
            _tau4[jRun]      = _tau4Default;
            double increment = _tau4Default / 4.;
            double minimumIncrement = 0.1;
            double sigma2;
            while ( increment > minimumIncrement )
            {
                sigma2 = lineScanTau4(jRun, increment, data);
//                qDebug() << "lineScan1D returns" << _tau4[jRun] << sigma2;
            }
//            qDebug() << "PETRTM::fitDataByGLMPlusLineScan result: " << _tau4[jRun];
        }
    }
    else
    { // fix k4 and vary BPnd
        for ( int jRun=0; jRun<_nRuns; jRun++ )
        {
            _BPndForIterations[jRun] = 0.;
            double increment = 1.;
            double minimumIncrement = 0.1;
            double sigma2;
            while ( increment > minimumIncrement )
            {
//                sigma2 = lineScan1D(jRun, _BPndForIterations, increment, data);
//                qDebug() << "lineScan1D returns" << _tau4[jRun] << sigma2;
            }
//            qDebug() << "PETRTM::fitDataByGLMPlusLineScan result: " << _BPndForIterations[jRun];
        }
    }

    // Attach the fit.
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        int iStartRun = jRun * _dtBinsSec[jRun].size();
        for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
        {
            int iTimeTotal = jt + iStartRun;
            yFit[jRun][jt] = getFit(iTimeTotal);
        }
    }
//    qDebug() << "PETRTM::fitDataByGLMPlusLineScan exit";
}
double PETRTM::updateLineScanFit(dVector data)
{
    createAllBasisFunctions();
    fitWLSForIteration(data);
    return getSigma2();
}

double PETRTM::lineScanTau4( int iRun, double &valueIncrement, dVector data )
{ // return value is sigma2
    double xParabola[3], yParabola[3];

    double valueCenter = _tau4[iRun];
    // Save the initial cost function.
    double sigma2Center = updateLineScanFit(data);
//    qDebug() << "lineScan1D enter" << _tau4[iRun] << valueIncrement << sigma2Center;
    xParabola[1] = valueCenter;
    yParabola[1] = sigma2Center;

    // test the - direction
    double valueLow = valueCenter - valueIncrement;
    _tau4[iRun] = valueLow;  double sigma2Low = updateLineScanFit(data);
//    qDebug() << "lineScan1D low" << _tau4[iRun] << sigma2Low << getSigma2();
    xParabola[0] = valueLow;
    yParabola[0] = sigma2Low;

    // test the + direction
    double valueHigh = valueCenter + valueIncrement;
    _tau4[iRun] = valueHigh; double sigma2High = updateLineScanFit(data);
//    qDebug() << "lineScan1D high" << _tau4[iRun] << sigma2High << getSigma2();
    xParabola[2] = valueHigh;
    yParabola[2] = sigma2High;

    bool minimumIsBracketed = sigma2Center < sigma2Low && sigma2Center < sigma2High;
    bool noChangeInSigma2   = sigma2Center == sigma2Low || sigma2Center == sigma2High;
    bool goInLowDirection   = sigma2Center > sigma2Low;
//    qDebug() << "lineScan1D bools" << minimumIsBracketed << noChangeInSigma2 << goInLowDirection;
    // bool goInHighDirection  = sigma2Center > sigma2High;
    double sigma2;  double delta = 0.;
    if ( minimumIsBracketed )
    { // the increment range contains the maximum of the cost function, so interpolate to get the best estimate
        double xMax, yMax;
        if ( utilMath::ParabolicInterpolation(xParabola, yParabola, xMax, yMax) )
        { // set the final increment to the interpolated value; half the increment range for next time
            sigma2 = yMax;
            delta = xMax - valueCenter;
            valueIncrement /= 2.;
        }
        else
        { // This should never happen.
            sigma2 = sigma2Center;
            valueIncrement = 0.;
        }
    }
    else if ( noChangeInSigma2 )
    { // not enough information somehow
        sigma2 = sigma2Center;
        valueIncrement = 0.;
    }
    else
    {
        double sign;  double tolerance = 1.e-5;
        if ( goInLowDirection )
        { // - direction: keep going below
            sigma2 = sigma2Low;
            _tau4[iRun] = valueLow;
            sign = -1.;
        }
        else // if ( goInHighDirection )
        { // + direction
            sigma2 = sigma2High;
            _tau4[iRun] = valueHigh;
            sign = 1.;
        }
        double sigma2New = sigma2;
        while ( sigma2New <= sigma2 )  // keep going as long as things get better (low sigma2)
        {
            sigma2 = sigma2New;
            _tau4[iRun] += sign * valueIncrement;
            sigma2New = updateLineScanFit(data);
//            qDebug() << "lineScan1D searching" << _tau4[iRun] << sigma2New;
            if ( sigma2New > sigma2 )
            { // The last step did not improve things, so set the delta value
                _tau4[iRun] -= sign * valueIncrement;
                valueIncrement /= 2.;
                sigma2 = updateLineScanFit(data);
//                qDebug() << "lineScan1D exit (reverse)" << sigma2;
                return sigma2;
            }
            else if ( sigma2New/sigma2 - 1. < tolerance )
            { // good enough, so stop
                sigma2 = sigma2New;
                delta = _tau4[iRun] - valueCenter;
                valueIncrement /= 2.;
//                qDebug() << "lineScan1D exit (tolerance)" << sigma2;
                return sigma2;
            }
            else
            {
                // good step but not within tolerance, so interpolate
                if ( sign < 0 )
                { // sign = -1; drop last point and incorporate new point
                    xParabola[2] = xParabola[1];
                    xParabola[1] = xParabola[0];
                    xParabola[0] = _tau4[iRun];
                    yParabola[2] = yParabola[1];
                    yParabola[1] = yParabola[0];
                    yParabola[0] = sigma2New;
                }
                else
                { // sign = +1; drop first point and incorporate new point
                    xParabola[0] = xParabola[1];
                    xParabola[1] = xParabola[2];
                    xParabola[2] = _tau4[iRun];
                    yParabola[0] = yParabola[1];
                    yParabola[1] = yParabola[2];
                    yParabola[2] = sigma2New;
                }
            } // sigma2New <= sigma2
            // reached a value of sigma2New > sigma2; now interpolate on last step
            double xMax, yMax;
            if ( utilMath::ParabolicInterpolation(xParabola, yParabola, xMax, yMax) )
            {
//                qDebug() << "ParabolicInterpolation" << xParabola << yParabola << xMax << yMax;
                sigma2 = yMax;
                delta = xMax - valueCenter;
            }
            else
                sigma2 = sigma2Center;
        }
    }
    _tau4[iRun] = valueCenter + delta;
    sigma2 = updateLineScanFit(data);
//    qDebug() << "lineScan1D exit (bracketed)" << _tau4[iRun] << delta << valueIncrement << sigma2;
    return sigma2;
}

int PETRTM::getTotalTimeIndex(int iRun, int iTimeInRun)
{
    int iTimeTotal = iTimeInRun;
    for (int jRun=0; jRun<iRun; jRun++)
        iTimeTotal += _dtBinsSec[jRun].size();
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
    if ( resetWeights ) // set weights using a specific time model
        setWeightsInRun(iRun);

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

void PETRTM::setWeightsInRun(int iRun)
{
//    qDebug() << "setWeightsInRun enter" << _PETWeightingModel << isFRTMFitk4();
    double trace=0.;
    int nTimeInRun = _dtBinsSec[iRun].size();
    if ( _weights[iRun].size() != nTimeInRun )
        _weights[iRun].resize(nTimeInRun);

    for (int jt=0; jt<nTimeInRun; jt++)
    {
        if ( isFRTMFitk4() )  // isFRTMNew && _fitk4UsingFixedBPnd
        {
            _weights[iRun][jt] = SQR(_frtmConv_dCtdtE[iRun][jt]);
//            qDebug() << "use weight[" << jt << "] =" << _weights[iRun][jt];
        }
        else if ( _PETWeightingModel == Weights_Uniform )
            _weights[iRun][jt] = 1.;
        else if ( _PETWeightingModel == Weights_noUptake )
            _weights[iRun][jt] = SQR(_frtmConvDeWeightUptake[iRun][jt]);
        else
        {
            double tau = 20.2334;  // 11C, minutes
            if ( _PETWeightingModel == Weights_18F || _PETWeightingModel == Weights_18F_Noiseless )
                tau = 110.;     // 18F, minutes
            tau *= 1.442695;    // convert from half life to exponential time constant
//            double time = getTimeInRun(iRun,jt);
            double time = _timeInRun[iRun][jt];
            double dt = static_cast<double>(_dtBinsSec[iRun][jt])/60.;
            _weights[iRun][jt] = dt / qExp(time/tau);
            if ( _PETWeightingModel == Weights_11C || _PETWeightingModel == Weights_18F ) //|| _PETWeightingModel == Weights_noUptake )
                _weights[iRun][jt] *= _tissRegion[iRun][jt];
        }
        trace += _weights[iRun][jt];
    }
//    qDebug() << "trace =" << trace;
    if ( trace == 0. ) trace = 1.;
    // Normalize weights: average value should = 1; this must be true to make F statistic valid.
    // E.g., post-hoc determination of _sigma2 assumes accurate separability of _sigma2 and weights = 1/sigma2
    for (int jt=0; jt<nTimeInRun; jt++)
        _weights[iRun][jt] *= static_cast<double>(nTimeInRun)/trace;
}

bool PETRTM::isValidID(int iRun, int iType, QChar eventID)
{
    bool valid = true;
    valid &= getEventIndex(eventID) >= 0;
    if ( iType == Type_R1 && isRTM2() )
        valid = false;
    else if ( iType == Type_dCrdt && !_dCrdtIncluded )
        valid = false;
    // make sure the reference region is defined for basis functions that need it
    if ( iType == Type_R1 || iType == Type_k2 || iType == Type_dCrdt )
    {
        valid &= _refRegionRaw.size() == _nRuns;
        int nTimeInRun = _dtBinsSec[iRun].size();
        if ( valid ) valid &= _refRegionRaw[iRun].size() == nTimeInRun;
        if ( valid )
        {
            bool nonZero = false;
            for (int jt=0; jt<nTimeInRun; jt++)
                nonZero |= _refRegionRaw[iRun][jt] != 0.;
            valid &= nonZero;
        }
//        if ( iType == Type_dCrdt )
//            valid &= _tau4[iRun] != 0.;
    }
    return valid;
};

int PETRTM::countEvents()
{
//    qDebug() << "PETRTM::countEvents enter";
    // Count events
    _basisID.resize(0);
    _basisShape.resize(0);
    _challengeIndex.resize(0);
    // This method of counting events assumes that IDs are independent across type.
    // e.g., an ID cannot be a k2 and k2a event in different runs
    for ( int jRun=0; jRun<_nRuns; jRun++ )
    {
//        qDebug() << "run" << jRun << "k2=" << _k2EventID[jRun] << "valid" << isValidID(jRun, Type_k2, _k2EventID[jRun]);
        // Count k2, k2a, and R1 events
        if ( ! _basisID.contains(_k2EventID[jRun]) && isValidID(jRun, Type_k2, _k2EventID[jRun]) )
        {
            _basisID.append(_k2EventID[jRun]);
            _basisShape.append(Type_k2);
            _challengeIndex.append(-1);
        }
        if ( ! _basisID.contains(_k2aEventID[jRun]) && isValidID(jRun, Type_k2a, _k2aEventID[jRun]) )
        {
            _basisID.append(_k2aEventID[jRun]);
            _basisShape.append(Type_k2a);
            _challengeIndex.append(-1);
        }
        if ( ! _basisID.contains(_R1EventID[jRun]) && isValidID(jRun, Type_R1, _R1EventID[jRun]) && isRTM3() )
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
//        qDebug() << "PETRTM::countEvents _k2EventCoefficient" << _k2EventID[jRun] << _k2EventCoefficient[jRun];
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
//    qDebug() << "PETRTM::countEvents exit" << nCoeff;
    return nCoeff;
}

void PETRTM::setRTMModelType(QString model)
{
//    qDebug() << "PETRTM::setRTMModelType enter";
    if ( model == "SRTM3" )
        _modelRTM = RTM_SRTM3;
    else if ( model == "SRTM2" )
        _modelRTM = RTM_SRTM2;
    else if ( model == "rFRTM3" )
        _modelRTM = RTM_rFRTM3;
    else if ( model == "rFRTM2" )
        _modelRTM = RTM_rFRTM2;
    else if ( model == "SRTM2Fit")
        _modelRTM = RTM_SRTM2Fit;
    else if ( model == "rFRTM3New" )
        _modelRTM = RTM_rFRTM3New;
    else if ( model == "rFRTM2New" )
        _modelRTM = RTM_rFRTM2New;
    else
    {
        qWarning() << "Error: model not defined: " << model;
        exit(1);
    }
    setPrepared(false);
};

void PETRTM::setRTMModelType(RTMModelTypes model)
{
    _modelRTM = model;
    setPrepared(false);
};

void PETRTM::createAllBasisFunctions()
{
    int nCoeff = countEvents();
//    qInfo() << "*************** create PET-RTM basis functions **************" << nCoeff << _nRuns;
    if ( nCoeff == 0 ) return;

    int nTimeTotal = 0;
    for (int jRun=0; jRun<_nRuns; jRun++)
        nTimeTotal += _dtBinsSec[jRun].size();
    init(nTimeTotal, nCoeff);

    // FRTM convolution may be needed for basis functions or weights. Calculate them as needed
    if ( isFRTM() || _PETWeightingModel == Weights_noUptake )  calculateFRTMConvolution();

    dVector weights;    weights.resize(nTimeTotal);
    if ( _nRuns != _weights.size() ) _weights.resize(_nRuns);
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        setWeightsInRun(jRun);
        setIgnoredPoints(jRun,false,_ignoreString[jRun]);
        for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
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
        if ( _basisShape[jCoeff] == Type_challenge )
            createChallengeBasisFunction(jCoeff, basis);
        else
            createRunBasisFunction(eventID, basis);
        addOrInsertBasisFunction(jCoeff,basis);
    }
    calculatePseudoInverse();
    defineConditions(getConditionString());
    setPrepared(true);
}

void PETRTM::createRunBasisFunction(QChar eventID, dVector &basis)
{
    dMatrix basisFunction;   basisFunction.resize(_nRuns);

    // For when not new TACs, reset the tissue vector to subtract the dCR/dt term

    // Create integrals
    dMatrix tissueIntegral = _tissRegion;
    integrateByRun(tissueIntegral);
    dMatrix frtmConvolutionIntegral;
    if ( isFRTM() ) // either old or new forms
    {
        if ( isFRTMOld() )
            frtmConvolutionIntegral = _frtmConv_dCtdtE;
        else if ( isFRTMNew() )
            frtmConvolutionIntegral = _frtmConv_CtE;
        integrateByRun(frtmConvolutionIntegral);
    }

    for ( int jRun=0; jRun<_nRuns; jRun++ )
    {
        basisFunction[jRun].fill(0.,_dtBinsSec[jRun].size());

        if ( eventID == _R1EventID[jRun] )
        { // true for 3-parameter RTM (SRTM3, rFRTM3)
            for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                basisFunction[jRun][jt] = _refRegion[jRun][jt];
        }
        else if ( eventID == _k2EventID[jRun] )
        {
            for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                basisFunction[jRun][jt] = _refRegionIntegral[jRun][jt];
            if ( isFRTM() )
            {
                if ( isFRTMNew() )
                {
                    for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                        basisFunction[jRun][jt] -= tissueIntegral[jRun][jt];
                }
                else  // isFRTMOld()
                {
                    for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                        basisFunction[jRun][jt] -= frtmConvolutionIntegral[jRun][jt];
                }
            }
            // 2-parameter RTM?
            if ( isRTM2() )
            {
                double tau2Ref = 0.;
                if ( isFRTM() )
                    tau2Ref = _tau2RefFRTMFixed[jRun];
                else if ( _modelRTM == RTM_SRTM2 )
                    tau2Ref = _tau2RefSRTMFixed[jRun];
                else if ( _modelRTM == RTM_SRTM2Fit )
                {
                    double tau2RefFit  = getTau2RefSRTMCalInRun(jRun,_BPndForIterations[jRun]);
                    double scaleFactor = _tau2RefSRTMCalOffset;
                    tau2Ref = tau2RefFit * scaleFactor;
                }
                for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                    basisFunction[jRun][jt] += tau2Ref * _refRegion[jRun][jt];
            }
        }
        else if ( eventID == _k2aEventID[jRun] )
        {
            if ( isFRTMNew() )
            {
                for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                    basisFunction[jRun][jt] = frtmConvolutionIntegral[jRun][jt];
            }
            else
            { // all other forms
                for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                    basisFunction[jRun][jt] = - tissueIntegral[jRun][jt];
                if ( isFRTMOld() )
                {
                    for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                        basisFunction[jRun][jt] += frtmConvolutionIntegral[jRun][jt];
                }
            }
        } // event types
        else if ( eventID == _dCrdtEventID[jRun] )
        { // optional function to search for vessels
            for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                basisFunction[jRun][jt] = _refRegionDeriv[jRun][jt];
        }
    } // jRun

    // convert matrix "basisFunction" to concatenated vector "basis"
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        int iStartRun = jRun * _dtBinsSec[jRun].size();
        for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
            basis[jt+iStartRun] = basisFunction[jRun][jt];
    }
}

void PETRTM::createChallengeBasisFunction(int iCoeff, dVector &basis)
{
    int indexChallenge = _challengeIndex[iCoeff];
    dMatrix basisFunction;   basisFunction.resize(_nRuns);
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        basisFunction[jRun].fill(0.,_dtBinsSec[jRun].size());
        if ( isGoodChallengeInRun(indexChallenge,jRun) )
        {
            dVector shape;
            createChallengeShape(jRun, indexChallenge, shape);
            if ( isFRTMNew() )
            {
                dVector equilibrationVector = getEquilibrationVector(jRun);
                dVector gammaCt = _tissRegionRaw[jRun];
                for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                    gammaCt[jt] *= shape[jt];
                dVector convolution = convolveEquilibration(jRun, gammaCt, equilibrationVector);
                basisFunction[jRun] = convolution;
                // include fixed k4 in the basis function
                for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                    basisFunction[jRun][jt] /= _tau4[jRun];
//                for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
//                    basisFunction[jRun][jt] *= - shape[jt];  // - {Ct - dC/dt X E(BPnd,k4)} * shape(t)
            }
            else
            {
                for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                    basisFunction[jRun][jt] = _tissRegion[jRun][jt];
                if ( isFRTM() )
                {
                    for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                        basisFunction[jRun][jt] -= _frtmConv_dCtdtE[jRun][jt];  // Ct - dC/dt X E(BPnd,k4)
                }
                for ( int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                    basisFunction[jRun][jt] *= - shape[jt];  // - {Ct - dC/dt X E(BPnd,k4)} * shape(t)
            }
        }
    } // jRun
    // Integrate the basis functions by run
    integrateByRun(basisFunction);

    // convert matrix "basisFunction" to concatenated vector "basis"
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        int iStartRun = jRun * _dtBinsSec[jRun].size();
        for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
            basis[jt+iStartRun] = basisFunction[jRun][jt];
    }
}

void PETRTM::calculateFRTMConvolution()
{
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        _frtmConv_dCtdtERaw[jRun].fill(0.,_dtBinsSec[jRun].size());
        _frtmConv_CtE[jRun].fill(0.,_dtBinsSec[jRun].size());
        _frtmConvDeWeightUptake[jRun].fill(0.,_dtBinsSec[jRun].size());
    }

    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        if ( _tau4[jRun] != 0. )
        {   // tau4 != 0
            dVector equilibrationVector = getEquilibrationVector(jRun);  // = k4 * (1+BPnd)
            if ( isFRTMNew() )
            {
//                qDebug() << "calculateFRTMConvolution check off" << _fitk4UsingFixedBPnd;
                _frtmConv_CtE[jRun] = convolveEquilibration(jRun, _tissRegionRaw[jRun],   equilibrationVector);
                for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                {
                    if ( _fitk4UsingFixedBPnd )
                        _frtmConv_CtE[jRun][jt] *= _BPndForIterations[jRun];
                    else
                        _frtmConv_CtE[jRun][jt] /= _tau4[jRun];
                }
            }
            else
            {
                _frtmConv_dCtdtERaw[jRun] = convolveEquilibration(jRun, _tissRegionDeriv[jRun], equilibrationVector);
                double maxFRTMConv=0.;
                for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                    if ( qAbs(_frtmConv_dCtdtERaw[jRun][jt]) > maxFRTMConv )
                        maxFRTMConv = qAbs(_frtmConv_dCtdtERaw[jRun][jt]);
                for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
                    _frtmConvDeWeightUptake[jRun][jt] = (1. - qAbs(_frtmConv_dCtdtERaw[jRun][jt]) / maxFRTMConv);
            }
        } // tau4 != 0
    } // jRun

    _frtmConv_dCtdtE = _frtmConv_dCtdtERaw;
    if ( _smoothingScaleMin != 0. )
        // Potentially smooth the convolution to reduce noise
        fitLoessCurve(_frtmConv_dCtdtE);  // additional smoothing xxx
}
dVector PETRTM::convolveEquilibration(int iRun, dVector tissue, dVector equilibration)
{
//    dVector tissueFine = tissue;
//    dVector equilibrationFine = equilibration;
    dVector tissueFine = interpolateTissueVector(iRun, tissue);
    dVector equilibrationFine = interpolateTissueVector(iRun,equilibration);
//    qDebug() << "fine sizes" << tissueFine.size() << equilibrationFine.size();

    int nTimeFine = tissueFine.size();
    double dt = static_cast<double>(_minBin[iRun])/60./static_cast<double>(_interpolationFactor);
    dVector convolution;  convolution.fill(0.,nTimeFine);
    for (int jt=0; jt<nTimeFine; jt++)
    {
        dVector exponential;          exponential.resize(jt+1);
        for ( int jtPrime=0; jtPrime<=jt; jtPrime++ )
        {
            double timeDiff = dt * (jt - jtPrime);
            exponential[jtPrime] = qExp(-equilibrationFine[jtPrime] * timeDiff);
            // Discrete values of the integral do not represent the average value across the bin,
            // so compute a correction factor to equilibrationFine that the exponential integral matches the continuous integral.
            double correctionFactor = qExp(-equilibrationFine[jtPrime]*dt/2.) / equilibrationFine[jtPrime] / dt;
            correctionFactor *= qExp(equilibrationFine[jtPrime] * dt) - 1.;
//            qDebug() << "exponential[" << jt << "][" << jtPrime << "] =" << exponential[jtPrime] << correctionFactor;
            exponential[jtPrime] *= correctionFactor;
            convolution[jt] += tissueFine[jtPrime] * exponential[jtPrime] * dt;
        }
    } // jt
    dVector convolutionLowRes = combineFineTissueVector(iRun,convolution);
//    dVector convolutionLowRes = convolution;
//    qDebug() << "convolutionHighRes" << convolution;
//    qDebug() << "convolutionLowRes" << convolutionLowRes;
    return convolutionLowRes;
}

void PETRTM::calculateBPndOrK4ForIteration()
{
    if ( isSRTMReg() )
    {
        for (int jRun=0; jRun<_nRuns; jRun++)
            _BPndForIterations[jRun] = getBP0InRun(jRun);
    }
    else if ( isFRTMFitk4() )
    {
        for (int jRun=0; jRun<_nRuns; jRun++)
            _tau4[jRun] = _tau4[jRun];
    }
}

void PETRTM::setSmoothingScaleMin(double smoothingScale)
{
    _smoothingScaleMin = smoothingScale;
    if ( _smoothingScaleMin != 0. )
    {
        for ( int jRun=0; jRun<_nRuns; jRun++ )
            _quadLOESS[jRun].define(_timeInRun[jRun],_smoothingScaleMin);
    }
    updateReferenceRegion();
}

void PETRTM::fitLoessCurve(dMatrix &runData)
{
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        dVector fit = _quadLOESS[jRun].fit(runData[jRun]);
        for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
            runData[jRun][jt] = _quadLOESS[jRun].getFitAtCentralPoint(jt);
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
    int nTimeInRun = _dtBinsSec[iRun].size();
    shape.fill(0,nTimeInRun);
    for ( int jStim=0; jStim<_maxStimuli; jStim++)
    {
        if ( isGoodStimulusInRun(indexChallenge,jStim,iRun) )
        {
            for ( int jt=0; jt<nTimeInRun; jt++)
            {
//                double time = getTimeInRun(iRun,jt);
                double time = _timeInRun[iRun][jt];
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
                    double tau = _challengeTau[indexChallenge];
                    // double alpha = _challengeAlpha[indexChallenge];
                    if ( time >= time0 )
                        shape[jt] += (time-time0)/tau * qExp(1.-(time-time0)/tau);
                }
                else if ( _challengeShape[indexChallenge] == Challenge_Sigmoid )
                {
                    double time0 = _challengeOn[indexChallenge][jStim];
                    double tau = _challengeTau[indexChallenge];
                    // double alpha = _challengeAlpha[indexChallenge];
                    if ( time >= time0 )
                        shape[jt] += (time-time0)/tau / qSqrt((1.+ (time-time0)/tau*(time-time0)/tau  ));
                }
            } // jt
        } // isGoodStimulusInRun
    } // jStim
}
void PETRTM::integrateByRun(dMatrix &runMatrix )  // [nRuns][nTimePerRun]
{
    dMatrix copiedMatrix = runMatrix;
    if ( _dtBinsSec.size() != _nRuns )
    { // this is just to provide a "default" option; bins sizes may differ in time
        _dtBinsSec.resize(_nRuns);
        for ( int jRun=0; jRun<_nRuns; jRun++ )
        {
            int nTimeInRun = runMatrix[jRun].size();
            _dtBinsSec[jRun].fill(60,nTimeInRun);
            setTimeBinsSec(jRun,_dtBinsSec[jRun]);
        }
    }
    // Find the integral from t'=0 to t'=t for each point. Reset at the start of each run.
    for ( int jRun=0; jRun<runMatrix.size(); jRun++ )
    {
        for ( int jt=0; jt<runMatrix[jRun].size(); jt++ )
        {
            double integral = 0.;
            for ( int jtPrime=0; jtPrime<=jt; jtPrime++ )
            {
                double dt = static_cast<double>(_dtBinsSec[jRun][jtPrime])/60.;
                if ( jtPrime == 0 )
                    integral += copiedMatrix[jRun][jtPrime] * dt;
                else
                    integral += (copiedMatrix[jRun][jtPrime] + copiedMatrix[jRun][jtPrime-1])/2. * dt;
            }
            runMatrix[jRun][jt] = integral;
        }
    }
}

void PETRTM::differentiateByRun(dMatrix &runMatrix )  // [nRuns][nTimePerRun]
{
    dMatrix copiedMatrix = runMatrix;
    for (int jRun=0; jRun<_nRuns; jRun++)
    {
        for (int jt=0; jt<_dtBinsSec[jRun].size(); jt++)
        {
            double dt = static_cast<double>(_dtBinsSec[jRun][jt])/60.;
            if ( jt != 0. )
                runMatrix[jRun][jt] = (copiedMatrix[jRun][jt] - copiedMatrix[jRun][jt-1]) / dt;
            else
                runMatrix[jRun][jt] = copiedMatrix[jRun][jt] / dt;
        }
    }
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
            int nTimeInRun = _dtBinsSec[iRun].size();
            for ( int jt=0; jt<nTimeInRun; jt++)
            {
//                double time = getTimeInRun(iRun,jt);
                double time = _timeInRun[iRun][jt];
                if ( time >= onSet )
                {
                    iStart = jt;
                    break;
                }
            } // jt
//            double time0 = getTimeInRun(iRun,iStart);
            double time0 = _timeInRun[iRun][iStart];
            int iTime = iStart - _nPreForChallAv;
            for ( int jt=0; jt<_yForChallAv.size(); jt++, iTime++)
            {
                if ( iTime >= 0 && iTime < yData[iRun].size() )
                {
//                    double time = getTimeInRun(iRun,iTime);
                    double time = _timeInRun[iRun][iTime];
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
                int nTimeInRun = _dtBinsSec[iRun].size();
                for ( int jt=0; jt<nTimeInRun; jt++)
                {
//                    double time = getTimeInRun(iRun,jt);
                    double time = _timeInRun[iRun][jt];
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

/////////////////////////////////////////////////////////////////////////
