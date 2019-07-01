#include <QtWidgets>
#include <QFrame>
#include <QDebug>
#include <QVector>
#include <QObject>
#include <QFileDialog>
#include <QString>
#include <QPen>
#include <QListWidget>
#include <QListWidgetItem>
#include <QSize>

#include "plot.h"

/*
plotData::~plotData()
{
}
*/

plotData::plotData(int plotID)
{
//    _parentPage = containingPage;

    _plotID = plotID;

    _qcplot = new QCustomPlot();
    _qcplot->setObjectName(QStringLiteral("plot"));
    _qcplot->xAxis->setLabel("time");
    _qcplot->yAxis->setLabel("y");
    _qcplot->xAxis->setLabelFont(QFont("Arial", 24));
    _qcplot->yAxis->setLabelFont(QFont("Helvetica", 24));
    _qcplot->xAxis->setTickLabelFont(QFont("Arial", 20));
    _qcplot->yAxis->setTickLabelFont(QFont("Helvetica", 20));
    _qcplot->xAxis2->setVisible(true);
    _qcplot->xAxis2->setTickLabels(false);
    _qcplot->yAxis2->setVisible(true);
    _qcplot->yAxis2->setTickLabels(false);
//    _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables | QCP::iSelectAxes);
    _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    _qcplot->setAutoAddPlottableToLegend(true);
    _qcplot->legend->setVisible(false);
    _qcplot->legend->setFont(QFont("Arial", 20));

    // add the phase tracer (red circle) which sticks to the graph data (and gets updated in bracketDataSlot by timer event):
    _positionTracer = new QCPItemTracer(_qcplot);
//    phaseTracer->setInterpolating(true);
    _positionTracer->setStyle(QCPItemTracer::tsCircle);
    _positionTracer->setPen(QPen(Qt::red));
    _positionTracer->setBrush(Qt::cyan);
    _positionTracer->setSize(16);
//    _positionTracer->setVisible(false);
    _positionTracer->setGraphKey(0);
    _positionTracer->setVisible(true);

    setSelectPoints();

    connect(_qcplot, SIGNAL(axisDoubleClick(QCPAxis*,QCPAxis::SelectablePart,QMouseEvent*)), this,
            SLOT(axisLabelDoubleClick(QCPAxis*,QCPAxis::SelectablePart)));
    connect(_qcplot, SIGNAL(legendDoubleClick(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*)), this,
            SLOT(legendDoubleClick(QCPLegend*,QCPAbstractLegendItem*)));
    connect(_qcplot, SIGNAL(mousePress(QMouseEvent*)),this,SLOT(plotMousePress(QMouseEvent*)));
    connect(_qcplot, SIGNAL(mouseMove(QMouseEvent*)), this,SLOT(plotMouseMove(QMouseEvent*)));

    connect(_qcplot->xAxis, SIGNAL(rangeChanged(QCPRange)), _qcplot->xAxis2, SLOT(setRange(QCPRange)));
    connect(_qcplot->yAxis, SIGNAL(rangeChanged(QCPRange)), this, SLOT(setAxis2Range(QCPRange)));

    // setup policy and connect slot for context menu popup:
    _qcplot->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(_qcplot, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequest(QPoint)));

    if ( _plotID == 0 )
    { // for the time-series plot only
        QAction *plusAction = new QAction(this);
        plusAction->setShortcut(Qt::Key_Equal);
        connect(plusAction, SIGNAL(triggered()), this, SLOT(keyboardPlus()));
        QAction *minusAction = new QAction(this);
        minusAction->setShortcut(Qt::Key_Minus);
        connect(minusAction, SIGNAL(triggered()), this, SLOT(keyboardMinus()));

        _qcplot->addAction(plusAction);
        _qcplot->addAction(minusAction);
    }

}

void plotData::keyboardMinus()
{
    int iFile = _iFilePosition;
    int iTime = _iTimePosition - 1;
    if ( iTime < 0 )
    { // go to the previous file
        iFile--;
        int iDataCurve = getDataCurveIndex(iFile);   // data points should be 1st curve in ifile
        if ( iDataCurve < 0 )
        {
            int iDataCurve = getDataCurveIndex(_nFiles-1);
            iTime = _listOfCurves[iDataCurve].yData.size() - 1; // last time point in last file
        }
    }
    setPositionTracer(iFile, iTime);
    emit changedPointFromGraph(_iFilePosition,_iTimePosition);
}

void plotData::keyboardPlus()
{
    int iFile = _iFilePosition;
    // Go to the next time point on the same curve, if it exists
    int iDataCurve = getDataCurveIndex(iFile);
    int iTime = _iTimePosition + 1;
    if ( iTime >= _listOfCurves[iDataCurve].yData.size() )
    { // go to the next file
        iFile++;  iTime=0;
        int iDataCurve = getDataCurveIndex(iFile);
        if ( iDataCurve >= _listOfCurves.size() )
            iFile = iTime = 0; // wrap-around to the 1st file and point
    }
    setPositionTracer(iFile, iTime);
    emit changedPointFromGraph(_iFilePosition,_iTimePosition);
}

///////////////////////////////////////////////////////////////////
///////////////////////// slots ///////////////////////////////////
///////////////////////////////////////////////////////////////////
void plotData::setXZoom()
{
    _qcplot->setCursor(Qt::OpenHandCursor);
    _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    _qcplot->axisRect(0)->setRangeDrag(Qt::Horizontal);
    _qcplot->axisRect(0)->setRangeZoom(Qt::Horizontal);
    autoScale(false);
    _pressToMoveTracer = false;
}
void plotData::setYZoom()
{
    _qcplot->setCursor(Qt::OpenHandCursor);
    _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    _qcplot->axisRect(0)->setRangeDrag(Qt::Vertical);
    _qcplot->axisRect(0)->setRangeZoom(Qt::Vertical);
    autoScale(false);
    _pressToMoveTracer = false;
}
void plotData::setSelectPoints()
{
    _qcplot->setCursor(Qt::CrossCursor);
    _qcplot->setInteractions(nullptr);
    _pressToMoveTracer = true;
}

void plotData::lastTimePoint()
{
    setPositionTracer(_iFilePositionLast, _iTimePositionLast);
    emit changedPointFromGraph(_iFilePosition,_iTimePosition);
}

void plotData::setAxis2Range( QCPRange Range1 )
{
    QCPRange Range2;
    Range2.lower = Range1.lower * _yAxis2Ratio;
    Range2.upper = Range1.upper * _yAxis2Ratio;
    _qcplot->yAxis2->setRange(Range2);
}
void plotData::axisLabelDoubleClick(QCPAxis *axis, QCPAxis::SelectablePart part)
{
  // Set an axis label by double clicking on it
  if (part == QCPAxis::spAxisLabel) // only react when the actual axis label is clicked, not tick label or axis backbone
  {
    bool ok;
    QString newLabel = QInputDialog::getText(this, "QCustomPlot example", "New axis label:", QLineEdit::Normal, axis->label(), &ok);
    if (ok)
    {
      axis->setLabel(newLabel);
      _qcplot->replot();
    }
  }
}
void plotData::legendDoubleClick(QCPLegend *legend, QCPAbstractLegendItem *item)
{
  // Rename a graph by double clicking on its legend item
  Q_UNUSED(legend)
  if (item) // only react if item was clicked (user could have clicked on border padding of legend where there is no item, then item is 0)
  {
    QCPPlottableLegendItem *plItem = qobject_cast<QCPPlottableLegendItem*>(item);
    bool ok;
    QString newName = QInputDialog::getText(this, "QCustomPlot example", "New graph name:", QLineEdit::Normal, plItem->plottable()->name(), &ok);
    if (ok)
    {
      plItem->plottable()->setName(newName);
      _qcplot->replot();
    }
  }
}
void plotData::contextMenuRequest(QPoint pos)
{
  QMenu *menu = new QMenu(this);
  menu->setAttribute(Qt::WA_DeleteOnClose);

  if (_qcplot->legend->selectTest(pos, false) >= 0 && _qcplot->legend->visible() ) // context menu on legend requested
  {
      if  ( _qcplot->legend->visible() )
          menu->addAction("Hide legend", this, SLOT(hideLegend()));
      else
          menu->addAction("Show legend", this, SLOT(showLegend()));

      menu->addAction("Move to top left", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignTop|Qt::AlignLeft)));
      menu->addAction("Move to top center", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignTop|Qt::AlignHCenter)));
      menu->addAction("Move to top right", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignTop|Qt::AlignRight)));
      menu->addAction("Move to bottom right", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignBottom|Qt::AlignRight)));
      menu->addAction("Move to bottom left", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignBottom|Qt::AlignLeft)));
  }
  else if (_qcplot->selectedGraphs().size() > 0)
  {
      menu->addAction("Write selected graph", this, SLOT(writeSelectedGraph()));
      menu->addAction("Remove selected graph", this, SLOT(removeSelectedGraph()));
      menu->addAction("LineStyle dotted", this, SLOT(makeSelectedGraphDotted()));
      menu->addAction("Bigger points", this, SLOT(makeSelectedGraphBiggerPoints()));
      menu->addAction("Smaller points", this, SLOT(makeSelectedGraphSmallerPoints()));
  }
  else
      menu->addAction("Write all curves to table file", this, SLOT(writeOneGraph()));
  /*
  if (_qcplot->selectedGraphs().size() > 0)
  {
      menu->addAction("Remove selected graph", this, SLOT(removeSelectedGraph()));
      menu->addAction("LineStyle dotted", this, SLOT(makeSelectedGraphDotted()));
      menu->addAction("Bigger points", this, SLOT(makeSelectedGraphBiggerPoints()));
      menu->addAction("Smaller points", this, SLOT(makeSelectedGraphSmallerPoints()));
  }
*/

  menu->popup(_qcplot->mapToGlobal(pos));
}
void plotData::makeSelectedGraphDotted()
{
  if (_qcplot->selectedGraphs().size() > 0)
  {
      QPen currentPen = _qcplot->selectedGraphs().first()->pen();
      QCPScatterStyle currentScatter = _qcplot->selectedGraphs().first()->scatterStyle();
      QBrush currentBrush = _qcplot->selectedGraphs().first()->brush();
      currentPen.setStyle(Qt::DotLine);
//      currentPen.setColor(Qt::green);
//      currentScatter.setBrush(Qt::green);
//      currentScatter.setPen(currentPen);
//      currentScatter.setSize(5);
//      currentScatter.setShape(QCPScatterStyle::ssCircle);
      _qcplot->selectedGraphs().first()->setScatterStyle(currentScatter);
      _qcplot->selectedGraphs().first()->setPen(currentPen);
      _qcplot->replot();
  }
}
void plotData::makeSelectedGraphBiggerPoints()
{
  if (_qcplot->selectedGraphs().size() > 0)
  {
      QCPScatterStyle currentScatter = _qcplot->selectedGraphs().first()->scatterStyle();
      currentScatter.setSize(1.5*currentScatter.size());
      _qcplot->selectedGraphs().first()->setScatterStyle(currentScatter);
      _qcplot->replot();
  }
}
void plotData::makeSelectedGraphSmallerPoints()
{
  if (_qcplot->selectedGraphs().size() > 0)
  {
      QCPScatterStyle currentScatter = _qcplot->selectedGraphs().first()->scatterStyle();
      currentScatter.setSize(0.75*currentScatter.size());
      _qcplot->selectedGraphs().first()->setScatterStyle(currentScatter);
      _qcplot->replot();
  }
}
void plotData::removeSelectedGraph()
{
  if (_qcplot->selectedGraphs().size() > 0)
  {
    _qcplot->removeGraph(_qcplot->selectedGraphs().first());
    _qcplot->replot();
  }
}

void plotData::writeOneGraph()
{
    QString fileName = "/graph.dat";
    QFileDialog fileDialog;
    QString fullFileName;
    fullFileName = fileDialog.getSaveFileName(this,
                                              "Name of file",
                                              QDir::currentPath()+fileName,
                                              tr("Text files (*.roi *.dat *.txt)"));
    if ( fullFileName.isEmpty() ) return;
    writeGraph(true,fullFileName,"");
}

void plotData::writeGraph(bool newGraph, QString fileName, QString regionName)
{
    QFile file(fileName);
    QTextStream out(&file);
    if ( newGraph )
    {
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
            return;
    }
    else
    {
        if (!file.open(QIODevice::Append | QIODevice::Text))
            return;
        QString message = "write region " + regionName + " to file " + fileName;
        _statusBar->showMessage(message);
        out << "new #" + regionName + "\n";
    }

    if ( !_concatenateRuns )
    {
        // Write the headers
        out << "index " << getLabelXAxis() << " ";
        for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
        {
            plotCurve currentCurve = _listOfCurves[jCurve];
            if ( currentCurve.exportCurve )
                out << currentCurve.legend << " ";
            if ( currentCurve.errorBars != 0 )
                out << currentCurve.legend << "_err ";
            if ( currentCurve.scaleFactor != 1. )
                out << currentCurve.legend << "_scaled ";
        }
        out << "\n";
        // Write the data
        for (int jt=0; jt<_listOfCurves[0].yData.size(); jt++)
        {
            out << jt << " " << _listOfCurves[0].xPlot[jt] << " ";
            for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
            {
                plotCurve currentCurve = _listOfCurves[jCurve];
                if ( currentCurve.exportCurve )
                {
                    out << currentCurve.yData[jt] << " ";
                    if ( currentCurve.errorBars != 0 )
                        out << _listOfCurves[jt].yError[jt] << " ";
                }
            }
            out << "\n";
        }
    }
    else
    {   // Concatenate runs
        // Write the headers
        out << "index time ";
        for (int jCurve=0; jCurve<_nCurvesPerFile; jCurve++)
        {
            plotCurve currentCurve = _listOfCurves[jCurve];
            if ( currentCurve.exportCurve )
            {
                out << currentCurve.legend << " ";
                if ( currentCurve.errorBars != 0 )
                    out << currentCurve.legend << "_err ";
                if ( currentCurve.scaleFactor != 1. )
                    out << currentCurve.legend << "_scaled ";
            }
        }
        out << "\n";
        // Write the data for each file in a concatenated series
        int iTime=0;
        for ( int jFile=0; jFile<_nFiles; jFile++ )
        {
            int iDataCurve = getDataCurveIndex(jFile);
            plotCurve currentCurve = _listOfCurves[iDataCurve];
            int nTime = currentCurve.yData.size();
            for ( int jt=0; jt<nTime; jt++, iTime++ )
            {
                out << jt << " " << currentCurve.xPlot[jt] << " ";
                for (int jCurve=0; jCurve<_nCurvesPerFile; jCurve++)
                {
                    int iDataCurve = getDataCurveIndex(jFile) + jCurve;
                    plotCurve currentCurve = _listOfCurves[iDataCurve];
                    if ( currentCurve.exportCurve )
                    {
                        out << currentCurve.yData[jt] << " ";
                        if ( currentCurve.errorBars != 0 )
                            out << currentCurve.yError[jt] << " ";
                    }
                } // jCurve
                out << "\n";
            } // jt
        } // jFile
    }

    file.close();
}

void plotData::writeSelectedGraph()
{
  if (_qcplot->selectedGraphs().size() == 1)
  {
      QString fileName;
      QFileDialog fileDialog;
      fileName = fileDialog.getSaveFileName(this,
                                            "Name of file",
                                            QDir::currentPath(),
                                            "(*.dat)");
      if ( fileName.isEmpty() ) return;
      QFile file(fileName);
      QTextStream out(&file);
      if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
          return;
      out << "x y\n";
      int nTime = _qcplot->selectedGraphs().at(0)->dataCount();
      for (int jt=0; jt<nTime; jt++)
          out << _qcplot->selectedGraphs().at(0)->dataMainKey(jt)   << " "
              << _qcplot->selectedGraphs().at(0)->dataMainValue(jt) << "\n";
      file.close();
  }
}
void plotData::moveLegend()
{
  if (QAction* contextAction = qobject_cast<QAction*>(sender())) // make sure this slot is really called by a context menu action, so it carries the data we need
  {
    bool ok;
    int dataInt = contextAction->data().toInt(&ok);
    if (ok)
    {
      _qcplot->axisRect()->insetLayout()->setInsetAlignment(0, static_cast<Qt::Alignment>(dataInt));
      _qcplot->replot();
    }
  }
}

void plotData::setPositionTracer(int iFile, int iTime)
{
    if ( iFile >= _listOfCurves.size() ) return;

    static bool concatenateLast=false;
    if ( iFile == _iFilePosition && iTime == _iTimePosition && concatenateLast == _concatenateRuns )
        return;
    concatenateLast = _concatenateRuns;
    double x = setXPosition(iFile, iTime);
    if ( _iFilePositionLast != _iFilePosition ) _iFilePositionLast = _iFilePosition;
    if ( _iTimePositionLast != _iTimePosition ) _iTimePositionLast = _iTimePosition;
    _iFilePosition = iFile;
    _iTimePosition = iTime;
    setPositionTracer(iFile,x);
}

void plotData::setPositionTracer(int iFile, double xPosition)
{ // file index = curve index; xPosition is actual x value on axis, with or without concatenation
    _positionTracer->setGraphKey(xPosition);
    _positionTracer->setVisible(true);

    for ( int jCurve=0; jCurve<_listOfCurves.size(); jCurve++ )
        _listOfCurves[jCurve].isCurrentCurve = (jCurve == getDataCurveIndex(iFile));

    plotDataAndFit(false);

    QCPItemPosition *position = _positionTracer->position;
    QString message;
    message.sprintf("(scan,pt) = (%d , %d)    (x,y) = (%g , %g)",_iFilePosition+1,_iTimePosition+1,
                    position->coords().x(),position->coords().y());
    if ( _statusBar )
        _statusBar->showMessage(message);
}

void plotData::plotMousePress(QMouseEvent *event)
{
    double x = _qcplot->xAxis->pixelToCoord(event->pos().x());
    // double y = _qcplot->yAxis->pixelToCoord(event->pos().y());
    if ( event->buttons() == Qt::LeftButton && _pressToMoveTracer )
    {
        if ( !_concatenateRuns )
        {
            // then don't identify the curve; keep the same one
            setPositionTracer(_iFilePosition,x);
            int iTimePoint = whichTracerTimePoint(_iFilePosition);
            setPositionTracer(_iFilePosition,iTimePoint);  // store last time point info
            emit changedPointFromGraph(_iFilePosition,iTimePoint);
        }
        else
        {
            // Need to indentify which curve and which time point
            int iFile = whichConcatenatedTracerFile(x);
//            qDebug() << "plotData::plotMousePress iFile" << iFile;
            if ( iFile < 0 ) return;
            setPositionTracer(iFile,x);
            int iTimePoint = whichTracerTimePoint(iFile);
            setPositionTracer(iFile,iTimePoint);   // store last time point info
            emit changedPointFromGraph(iFile,iTimePoint);
        }
        _qcplot->replot();
    }
}

double plotData::setXPosition( int iFile, int iTime)
{   // This function assumes all curves for a file have the same x values
//    qDebug() << "plotData::setXPosition" << iFile << iTime;
//    qDebug() << "plotData::setXPosition index" << getDataCurveIndex(iFile);
    double xValue = _listOfCurves[getDataCurveIndex(iFile)].xData[iTime];
//    qDebug() << "plotData::setXPosition1";
    if ( iFile >= 0 && _concatenateRuns )
    {
        for (int jFile=0; jFile<iFile; jFile++)
        {
            int iDataCurve = getDataCurveIndex(jFile);
            int iLastTimeIndexLastFile = _listOfCurves[iDataCurve].xData.size() - 1;
            double deltaTime = 1.;
            if ( iLastTimeIndexLastFile != 0 )
                deltaTime = _listOfCurves[iDataCurve].xData[iLastTimeIndexLastFile] - _listOfCurves[iDataCurve].xData[iLastTimeIndexLastFile-1];
            double startTime = _listOfCurves[iDataCurve].xData[iLastTimeIndexLastFile];
            xValue += startTime + deltaTime;  // last run + 1 time step
        }
    }
//    qDebug() << "plotData::setXPosition exit" << xValue;
    return xValue;
}

int plotData::getDataCurveIndex(int iFile)
{ // data should be 1st curve
    int iCurve=0;
    for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
    {
        if ( _listOfCurves[jCurve].iFile == iFile )
        {
            iCurve = jCurve;
            break;
        }
    }
    return iCurve;
}

int plotData::whichConcatenatedTracerFile( double x )
{
    double firstX, lastX;
    for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
    {
//        qDebug() << "for curve" << jCurve << "iFile=" << _listOfCurves[jCurve].iFile << "nTime" << _listOfCurves[jCurve].xPlot.size();
        firstX = _listOfCurves[jCurve].xPlot[0];
        lastX  = _listOfCurves[jCurve].xPlot[_listOfCurves[jCurve].xPlot.size()-1];
        if ( firstX == lastX )
        {
            firstX -= 0.5;  lastX += 0.5;
        }
//        qDebug() << "x, firstX, lastX" << x << firstX << lastX << _listOfCurves[jCurve].xPlot.size();
//        if ( _listOfCurves[jCurve].xPlot.size() == 1 || (x >= firstX && x<= lastX) )
        if ( x >= firstX && x<= lastX )
            return _listOfCurves[jCurve].iFile;
    }
    return -1;
}

int plotData::whichTracerTimePoint( int iFile )
{
    _positionTracer->updatePosition();
    double xTracer = _positionTracer->position->key();
    int iDataCurve = getDataCurveIndex(iFile);
    for (int jT=0; jT<_listOfCurves[iDataCurve].xPlot.size(); jT++)
    {
        if ( xTracer == _listOfCurves[iDataCurve].xPlot[jT] )
            return jT;
    }
    qWarning() << "Programming error: time point not found for curve " << iDataCurve << " in function whichTracerTimePoint";
    exit(1);
}

void plotData::plotMouseMove(QMouseEvent *event)
{
    interpretMousePosition(event);
    // bool leftButtonMove = (event->type() == QEvent::MouseMove) && (event->buttons() == Qt::LeftButton);

    if ( _cursorOnPlot )
    {
        double x = _qcplot->xAxis->pixelToCoord(event->pos().x());
        double y = _qcplot->yAxis->pixelToCoord(event->pos().y());
//        double y2 = _qcplot->yAxis2->pixelToCoord(event->pos().y());
        QString message = QString("%1 , %2").arg(x).arg(y);
        _qcplot->setToolTip(message);
        if ( _statusBar )
            _statusBar->showMessage(message);
    }
    else
    {
        QCPItemPosition *position = _positionTracer->position;
        QString message;
        message.sprintf("(scan,pt) = (%d , %d)    (x,y) = (%g , %g)",_iFilePosition+1,_iTimePosition+1,
                        position->coords().x(),position->coords().y());
        if ( _statusBar )
            _statusBar->showMessage(message);
    }
}

void plotData::interpretMousePosition(QMouseEvent *event)
{
    qreal x = event->localPos().x();
    qreal y = event->localPos().y();
    _cursorOnPlot = false;

    if ( x >= _qcplot->axisRect()->left() && x <=_qcplot->axisRect()->rect().right() &&
         y >= _qcplot->axisRect()->top()  && y <=_qcplot->axisRect()->rect().bottom() )
        _cursorOnPlot = true;
}

double plotData::getMaxY(plotCurve *ptrCurve)
{
    double dMax = -1.e10;
    for ( int jPoint=0; jPoint<ptrCurve->yData.size(); jPoint++ )
    {
        if ( ptrCurve->yData[jPoint] > dMax ) dMax = ptrCurve->yData[jPoint];
    }
    return dMax;
}

void plotData::init()
{
    _listOfCurves.resize(0);
    setLegendOn(false);
    setLabelXAxis("time");
    setLabelYAxis("signal");
    _yAxis2Ratio = 0.;
}
void plotData::conclude(int iCurrentFile, bool singleRun)
{
    // Find current curve as first instance of current file counter
    bool found=false;
    for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
        if ( _listOfCurves[jCurve].isCurrentCurve ) found=true;

    int iFileCurrent=0;
    if ( !found )
    {
        for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
            _listOfCurves[jCurve].isCurrentCurve = false;
        for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
        {
            if ( _listOfCurves[jCurve].iFile == iCurrentFile )
            {  // set this as the current curve with the tracer
                _listOfCurves[jCurve].isCurrentCurve = true;
                iFileCurrent = _listOfCurves[jCurve].iFile;
                break;
            }
        }
    }

    _nFiles = 0;
//    qDebug() << "plotData::conclude enter" << _listOfCurves.size() << iCurrentFile;
    for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
    {
        if ( _listOfCurves[jCurve].isCurrentCurve )
            _listOfCurves[jCurve].visible = true;
        else
            _listOfCurves[jCurve].visible = !singleRun || (_listOfCurves[jCurve].iFile == iFileCurrent);

        if ( _listOfCurves[jCurve].iFile+1 > _nFiles ) _nFiles = _listOfCurves[jCurve].iFile+1;

        // Set the color for iFile=0 only (generally data)
//        qDebug() << "set color in conclude" << _listOfCurves[jCurve].legend << _listOfCurves[jCurve].isCurrentCurve << _listOfCurves[jCurve].iFile;
    }
//    qDebug() << "plotData::conclude nFiles" << _nFiles;

    // concatenate runs
    concatenateRuns();
//    qDebug() << "plotData::conclude exit";
}

void plotData::concatenateRuns()
{
//    qDebug() << "plotData::concatenateRuns enter" << _listOfCurves.size();
    for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
    {
        int iFile = _listOfCurves[jCurve].iFile;
        if ( iFile == 0 || !_concatenateRuns )
            _listOfCurves[jCurve].xPlot = _listOfCurves[jCurve].xData;
        else
        {
            int iDataLast = getDataCurveIndex(iFile-1);
//            qDebug() << "plotData::concatenateRuns jCurve iFile iDataLast" << jCurve << iFile << iDataLast;
            int nTime = _listOfCurves[jCurve].xData.size();
            double offset;
            if ( nTime == 1 )
                offset = 0.;
            else
            {
                double delta = _listOfCurves[iDataLast].xPlot[nTime-1] - _listOfCurves[iDataLast].xPlot[nTime-2];
                int iLast = _listOfCurves[iDataLast].xPlot.size()-1;
                offset = _listOfCurves[iDataLast].xPlot[iLast] + delta;
            }
            _listOfCurves[jCurve].xPlot.resize(nTime);
            for ( int jt=0; jt<nTime; jt++ )
                _listOfCurves[jCurve].xPlot[jt] = _listOfCurves[jCurve].xData[jt] + offset;
        }
//        qDebug() << "xplot size for" << _listOfCurves[jCurve].legend << " = " << _listOfCurves[jCurve].xPlot.size();
    }
//    qDebug() << "plotData::concatenateRuns exit" << _listOfCurves.size();
}

void plotData::addCurve(int iFile, QString legend)
{
    plotCurve newCurve;
    newCurve.iFile = iFile;
    newCurve.iCurve = _listOfCurves.size();
    newCurve.lineThickness = 1;
    newCurve.pointSize = 0;
    newCurve.pointStyle = QCPScatterStyle::ssNone;
    newCurve.color = Qt::black;
    newCurve.visible = true;
    newCurve.enabled = true;
    newCurve.exportCurve = true;
    newCurve.isCurrentCurve = false;
    newCurve.histogram = false;
    _listOfCurves.append(newCurve);
    setLegend(legend);
}

void plotData::setData(dVector xData, dVector yData)
{
    bVector ignoreRescale;
    ignoreRescale.fill(false,xData.size());
    dVector yError;
    yError.fill(0.,xData.size());
    setData(xData, yData, yError, ignoreRescale);
}
void plotData::setData(dVector xData, dVector yData, dVector yError)
{
    bVector ignoreRescale;
    ignoreRescale.fill(false,xData.size());
    setData(xData, yData, yError, ignoreRescale);
}
void plotData::setData(dVector xData, dVector yData, bVector ignoreRescale)
{
    dVector yError;
    yError.fill(0.,xData.size());
    setData(xData, yData, yError, ignoreRescale);
}
void plotData::setData(dVector xData, dVector yData, dVector yError, bVector ignoreRescale)
{
    _listOfCurves.last().xData         = xData;
    _listOfCurves.last().yData         = yData;
    _listOfCurves.last().yError        = yError;
    _listOfCurves.last().ignoreRescale = ignoreRescale;
}

void plotData::plotDataAndFit(bool newData)
{
    plotCurve currentCurve;
    QPen myPen;

//    qDebug() << "plotData::plotDataAndFit enter";

    // Reset x position based upon _concatenateRuns
    /*
    for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
    {
        currentCurve = _listOfCurves[jCurve];
        for (int jTime=0; jTime<_listOfCurves[jCurve].xData.size(); jTime++)
            currentCurve.xPlot[jTime] = setXPosition(currentCurve.iFile,jTime);
    }
    */

    // clear the graph (the following can return the number cleared)
    _qcplot->clearGraphs();
    QCPRange Range1 = _qcplot->yAxis->range();
    _qcplot->yAxis2->setRange(Range1);
    _qcplot->yAxis2->setTickLabels(false);
    _qcplot->yAxis2->setLabel("");

    // plot curves in reverse order so that the first curves (e.g., jCurve=0) are not occluded by later curves
//    for (int jCurve=0; jCurve<_listOfCurves.size(); jCurve++)
    for (int jCurve=_listOfCurves.size()-1; jCurve>=0; jCurve--)
    {
        currentCurve = _listOfCurves[jCurve];
        if ( currentCurve.visible && currentCurve.enabled )
        {
            // Fill color for error envelopes
            QColor fillColor;
            int hue; int saturation; int value;
            currentCurve.color.getHsv(&hue, &saturation, &value);
            fillColor.setHsv(hue,saturation/4,value);

            // Pen for lines & points
            if ( currentCurve.errorBars == 2 )
                myPen.setColor(fillColor);
            else
                myPen.setColor(currentCurve.color);
            myPen.setWidth(currentCurve.lineThickness);
            // Points
            QCPScatterStyle myScatter;
            myScatter.setShape(currentCurve.pointStyle);
            myScatter.setPen(myPen);
            myScatter.setBrush(currentCurve.color);
            myScatter.setSize(currentCurve.pointSize);

            _qcplot->addGraph();

            // set the curve name and color
            _qcplot->graph()->setName(currentCurve.legend);
            if (currentCurve.legend.compare("none") == 0  )
            {
                int nLegendItems = _qcplot->legend->itemCount();
                _qcplot->legend->removeItem(nLegendItems-1);
            }
            // Points
            _qcplot->graph()->setPen(myPen);
            _qcplot->graph()->setScatterStyle(myScatter);
            // lines
            if ( currentCurve.lineThickness == 0 )
                _qcplot->graph()->setLineStyle(QCPGraph::lsNone);
            if ( currentCurve.histogram )
                _qcplot->graph()->setLineStyle(QCPGraph::lsStepCenter);
            // Set the data
//            qDebug() << "set data" << currentCurve.legend << currentCurve.xPlot.size() << currentCurve.yData.size();
            _qcplot->graph()->setData(currentCurve.xPlot,currentCurve.yData);
            if ( _listOfCurves[jCurve].isCurrentCurve && currentCurve.yData.size() != 0 )
                _positionTracer->setGraph(_qcplot->graph());  // set this when graphs are defined
            //                if ( currentCurve.errorBars == 2 && jCurve != 0 )
            if ( currentCurve.errorBars == 1 )
            {
                // error bars
                QCPErrorBars *errorBars = new QCPErrorBars(_qcplot->xAxis, _qcplot->yAxis);
                errorBars->removeFromLegend();
                errorBars->setAntialiased(false);
                // errorBars->setDataPlottable(_qcplot->graph());
                int iGraph = _qcplot->graphCount() - 1; // iGraph-1 is  current graph
                errorBars->setDataPlottable(_qcplot->graph(iGraph));
                errorBars->setPen(myPen);
                // Set the data
                errorBars->setData(currentCurve.yError);
            }
            else if ( currentCurve.errorBars == 2 )
            {
                if ( _listOfCurves[jCurve+1].errorBars == 2 ) // jCurve+1 due to reverse order plotting; otherwise jCurve-1
                {
                    int iGraph = _qcplot->graphCount() - 2; // iGraph-1 is  current graph
                    _qcplot->graph()->setBrush(QBrush(QBrush(fillColor)));
                    _qcplot->graph()->setChannelFillGraph(_qcplot->graph(iGraph));
                }
            } // error bars
            if ( currentCurve.scaleFactor != 1. )
            {
                _qcplot->yAxis2->setVisible(true);
                _qcplot->yAxis2->setTickLabels(true);
                _qcplot->yAxis2->setTickLabelColor(currentCurve.color);
                _qcplot->yAxis2->setTickLabelFont(QFont("Helvetica", 24));
                Range1 = _qcplot->yAxis->range();
                setAxis2Range(Range1);
                _qcplot->yAxis2->setLabel(currentCurve.legend);
                _qcplot->yAxis2->setLabelFont(QFont("Helvetica", 24));
                _qcplot->yAxis2->setLabelColor((currentCurve.color));
            }
        } // visible
    } //jcurve
    if ( _autoScale && newData )
        autoScale(true);
    else
        _qcplot->replot();
//    qDebug() << "plotData::plotDataAndFit exit";
}

void plotData::autoScale(bool state)
{
//    qDebug() << "plotData::autoScale" << state;
    _autoScale = state;
    if ( _autoScale )
    {
        QCPRange xRange, yRange;
        reScaleAxes(&xRange, &yRange);
//        qDebug() << "plotData::autoScale yRange" << yRange.lower << yRange.upper;
        _qcplot->yAxis->setRange(yRange);
        setAxis2Range(yRange);
        if ( _autoScaleXRange.lower == _autoScaleXRange.upper )
            _qcplot->xAxis->setRange(xRange);
        else
            _qcplot->xAxis->setRange(_autoScaleXRange);
        _qcplot->replot();
    }
}

void plotData::reScaleAxes(QCPRange *xRange, QCPRange *yRange)
{ // don't use the QCustomPlot function rescaleAxes (_qcplot->rescaleAxes()) so we can flag ignored points in the re-scaling
    xRange->lower=1.e10;  xRange->upper=-1.e10;
    yRange->lower=1.e10;  yRange->upper=-1.e10;
//    qDebug() << "Data::reScaleAxes enter";

    for ( int jGraph=0; jGraph<_listOfCurves.size(); jGraph++ )
    {
//        qDebug() << "Data::reScaleAxes jGraph" << jGraph;
        if ( _listOfCurves[jGraph].visible && _listOfCurves[jGraph].enabled )
        {
//            qDebug() << "Data::reScaleAxes 2";
            int nPoints = _listOfCurves[jGraph].yData.size();
            for ( int jPoint=0; jPoint<nPoints; jPoint++ )
            {
//                qDebug() << "Data::reScaleAxes 3" << _listOfCurves[jGraph].ignoreRescale[jPoint];
                bool errorBars = _listOfCurves[jGraph].errorBars == 1; // for errorBars=2 (envelop, data are passed as two line curves)
                if ( !_listOfCurves[jGraph].ignoreRescale[jPoint] )
                {  // Ignore the y value of "ignored" points
//                    qDebug() << "Data::reScaleAxes 4";
                    double yMin, yMax;
                    if ( errorBars )
                    {
                        yMin = _listOfCurves[jGraph].yData[jPoint] - _listOfCurves[jGraph].yError[jPoint];
                        yMax = _listOfCurves[jGraph].yData[jPoint] + _listOfCurves[jGraph].yError[jPoint];
                    }
                    else
                        yMin = yMax = _listOfCurves[jGraph].yData[jPoint];
                    //                qDebug() << "graph" << jGraph << "point" << jPoint << "yMin, yMax" << yMin << yMax;
                    if ( yMin < yRange->lower ) yRange->lower = yMin;
                    if ( yMax > yRange->upper ) yRange->upper = yMax;
                }
                // Do not ignore the x value of "ignored" points.
                if ( _listOfCurves[jGraph].xPlot[jPoint] < xRange->lower ) xRange->lower = _listOfCurves[jGraph].xPlot[jPoint];
                if ( _listOfCurves[jGraph].xPlot[jPoint] > xRange->upper ) xRange->upper = _listOfCurves[jGraph].xPlot[jPoint];
            }
        }
    }
    // expland the x and y ranges by 10% total
    double extra = (xRange->upper - xRange->lower) * 0.05;
    xRange->lower -= extra;
    xRange->upper += extra;
    extra = (yRange->upper - yRange->lower) * 0.05;
//    qDebug() << "rescale" << yRange->lower << yRange->upper << extra;
    yRange->lower -= extra;
    yRange->upper += extra;
//    qDebug() << "Data::reScaleAxes exit";
}
